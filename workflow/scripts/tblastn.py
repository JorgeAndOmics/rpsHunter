"""Search one genome with one query protein (tBLASTn) and write its hit table.

The genome is searched translated in six frames. Each hit becomes a row with
BLAST's statistics, the hit's translation, a Tag derived from the hit itself
(ADR-004) and a flag for a stop codon inside the translation (ADR-011). The
BLAST archive is kept beside the table.

Usage (called by the Snakefile):
    python tblastn.py --species KEY --query-accession ACC \
        --query-file QUERY.fa --output-parquet HITS.parquet --output-csv HITS.csv \
        --output-asn ARCHIVE.asn

No hit gives an empty table with every column. A failing tBLASTn stops the job
(external.ToolError): a crashed search is never mistaken for a genome without hits.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import NamedTuple

import pandas as pd

import defaults
import external
import hit_table
import log
import tags

COLUMNS = (
    "Query ID",
    "Subject ID",
    "Pct Identity",
    "Alignment Length",
    "Mismatches",
    "Gap Openings",
    "Q. Start",
    "Q. End",
    "S. Start",
    "S. End",
    "E-value",
    "Bit Score",
    "Subject Sequence",
)
"""Columns of the tabular output asked of blast_formatter (FORMAT). S. Start is
greater than S. End for a hit on the minus strand."""

HIT_COLUMNS = (
    *COLUMNS,
    "Species",
    "Species_Name",
    "Query_Accession",
    "Tag",
    "Internal_Stop",
)

FORMAT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"

_INTEGER = (
    "Alignment Length",
    "Mismatches",
    "Gap Openings",
    "Q. Start",
    "Q. End",
    "S. Start",
    "S. End",
)
_FLOAT = ("Pct Identity", "E-value", "Bit Score")
_TYPES = dict.fromkeys(_INTEGER, "int64") | dict.fromkeys(_FLOAT, "float64")


class Search(NamedTuple):
    """One tBLASTn search: the programs, the genome's database, the query, the cut-off."""

    program: str
    formatter: str
    database: Path
    query: Path
    max_evalue: float


def search_command(search: Search, archive: Path) -> list[str]:
    """Return the tBLASTn command that writes an archive of the search."""
    return [
        search.program,
        "-db", str(search.database),
        "-query", str(search.query),
        "-evalue", str(search.max_evalue),
        "-outfmt", "11",
        "-out", str(archive),
    ]  # fmt: skip


def parse_hits(text: str) -> pd.DataFrame:
    """Parse blast_formatter's tabular output into a typed table with COLUMNS.

    Raises:
        ValueError: A row does not have the 13 expected fields.
    """
    rows = [line.split("\t") for line in text.splitlines() if line]
    wrong = [row for row in rows if len(row) != len(COLUMNS)]
    if wrong:
        raise ValueError(
            f"expected {len(COLUMNS)} tab-separated fields, got {len(wrong[0])}"
        )
    return pd.DataFrame(rows, columns=list(COLUMNS)).astype(_TYPES)


def describe(
    hits: pd.DataFrame, species: str, species_name: str, accession: str
) -> pd.DataFrame:
    """Add provenance, the Tag and the stop flag; a hit reported twice is kept once."""
    described = hits.assign(
        Species=species, Species_Name=species_name, Query_Accession=accession
    )
    # Duplicates go before tagging: a Tag identifies one hit, and tag_hits
    # refuses a table in which a hit appears twice.
    described = described.drop_duplicates().reset_index(drop=True)
    described["Tag"] = tags.tag_hits(described)
    described["Internal_Stop"] = hit_table.has_internal_stop(
        described["Subject Sequence"]
    )
    return described[list(HIT_COLUMNS)]


def run(search: Search, archive: Path) -> pd.DataFrame:
    """Run the search and return its hits (COLUMNS)."""
    external.run_tool(search_command(search, archive), tool="tblastn")
    formatter = [search.formatter, "-archive", str(archive), "-outfmt", FORMAT]
    return parse_hits(external.run_tool(formatter, tool="blast_formatter"))


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--species", required=True)
    parser.add_argument("--query-accession", required=True)
    parser.add_argument("--query-file", type=Path, required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-asn", type=Path, required=True)
    log.add_scope_argument(parser)
    args = parser.parse_args(argv)
    log.job_logging("tblastn", args.scope)

    search = Search(
        program=defaults.TBLASTN_CMD,
        formatter=defaults.BLAST_FORMATTER_CMD,
        database=defaults.PATH_DICT["SPECIES_DB"] / args.species / args.species,
        query=args.query_file,
        max_evalue=defaults.E_VALUE_THRESHOLD,
    )
    name = defaults.SPECIES_DICT[args.species]
    hits = describe(
        run(search, args.output_asn), args.species, name, args.query_accession
    )
    hits.to_parquet(args.output_parquet, index=False)
    hits.to_csv(args.output_csv, index=False)
    log.ok(f"{log.count(len(hits))} hits at E <= {search.max_evalue:g}")
    return 0


if __name__ == "__main__":
    sys.exit(log.run_main(main))
