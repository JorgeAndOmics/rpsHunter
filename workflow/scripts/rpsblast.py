"""Search one species' hit sequences against the families' CDD models (RPS-BLAST).

The database holds only the models listed under `domains:`, so the search is
told the size of the full CDD release (`-dbsize`). E-values are then on the
scale NCBI's CD-Search and the literature use, and do not change when the list
of families does (ADR-006). The archive (ASN.1) is kept for rpsbproc; the same
alignments are also written as a table, one row per alignment, with the model
coordinates that later steps need.

Usage (called by the Snakefile):
    python rpsblast.py --species KEY --query-accession ACC --fasta-input FASTA \
        --output-parquet TABLE --output-asn ARCHIVE --db DATABASE_PREFIX

A FASTA with no sequence (nothing passed the gates for this species) gives an
empty table and an empty archive without running the search. A search that
fails stops the job (external.ToolError).
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path
from typing import NamedTuple

import pandas as pd

import defaults
import external
import log
import tags

QUERY_ID = "Query ID"

TABLE_COLUMNS = [
    QUERY_ID,
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
    "Subject Title",
]
"""Columns of `blast_formatter -outfmt "6 std stitle"`. Q. is the hit's translated
sequence (residues); S. is the CDD model."""

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
_PROVENANCE = ("Species", "Species_Name", "Query_Accession")


class SearchSettings(NamedTuple):
    """The RPS-BLAST options that come from the config."""

    max_evalue: float
    database_residues: int
    comp_based_stats: int
    seg: str
    window_size: int


def search_command(
    program: str, database: Path, query: Path, archive: Path, settings: SearchSettings
) -> list[str]:
    """Return the RPS-BLAST command that writes an archive of the search."""
    return [
        program,
        "-db", str(database),
        "-query", str(query),
        "-evalue", str(settings.max_evalue),
        "-dbsize", str(settings.database_residues),
        "-comp_based_stats", str(settings.comp_based_stats),
        "-seg", settings.seg,
        "-window_size", str(settings.window_size),
        "-outfmt", "11",
        "-out", str(archive),
    ]  # fmt: skip


def parse_table(text: str) -> pd.DataFrame:
    """Parse tabular RPS-BLAST output into a typed table with the Tag split off.

    Args:
        text: Output of `blast_formatter -outfmt "6 std stitle"`; may be empty.

    Returns:
        TABLE_COLUMNS plus "Tag". "Query ID" no longer carries the `|tag:` suffix.

    Raises:
        ValueError: A row does not have the 13 expected fields.
    """
    rows = [line.split("\t") for line in io.StringIO(text).read().splitlines() if line]
    wrong = [row for row in rows if len(row) != len(TABLE_COLUMNS)]
    if wrong:
        raise ValueError(
            f"expected {len(TABLE_COLUMNS)} tab-separated fields, got {len(wrong[0])}: {wrong[0]}"
        )
    table = pd.DataFrame(rows, columns=TABLE_COLUMNS)
    table = table.astype(
        dict.fromkeys(_INTEGER, "int64") | dict.fromkeys(_FLOAT, "float64")
    )
    identifier = table[QUERY_ID].astype("object")
    table["Tag"] = identifier.str.extract(tags.TAG_IN_IDENTIFIER, expand=False)
    table[QUERY_ID] = identifier.str.replace(tags.TAG_IN_IDENTIFIER, "", regex=True)
    return table


def configured_settings() -> SearchSettings:
    """The search settings of the current config."""
    return SearchSettings(
        max_evalue=defaults.RPSBLAST_E_VALUE,
        database_residues=defaults.CDD_RESIDUES,
        comp_based_stats=defaults.RPSBLAST_COMP_BASED_STATS,
        seg=defaults.RPSBLAST_SEG,
        window_size=defaults.RPSBLAST_WINDOW_SIZE,
    )


def search(fasta: Path, database: Path, archive: Path) -> pd.DataFrame:
    """Run the search for one FASTA and return its alignment table."""
    if fasta.stat().st_size == 0:
        # Nothing passed the gates: a genuine empty result, known before searching.
        # (rpsblast on an empty FASTA writes an archive blast_formatter rejects.)
        archive.touch()
        return parse_table("")
    command = search_command(
        defaults.RPSBLAST_CMD, database, fasta, archive, configured_settings()
    )
    external.run_tool(command, tool="rpsblast")
    formatter = [
        defaults.BLAST_FORMATTER_CMD,
        "-archive",
        str(archive),
        "-outfmt",
        "6 std stitle",
    ]
    return parse_table(external.run_tool(formatter, tool="blast_formatter"))


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--species", required=True)
    parser.add_argument("--query-accession", required=True)
    parser.add_argument("--fasta-input", type=Path, required=True)
    parser.add_argument("--output-parquet", type=Path, required=True)
    parser.add_argument("--output-asn", type=Path, required=True)
    parser.add_argument("--db", type=Path, required=True)
    log.add_scope_argument(parser)
    args = parser.parse_args(argv)
    log.job_logging("rpsblast", args.scope)

    table = search(args.fasta_input, args.db, args.output_asn)
    table["Species"] = args.species
    table["Species_Name"] = defaults.SPECIES_DICT[args.species]
    table["Query_Accession"] = args.query_accession
    table = table[[*TABLE_COLUMNS, "Species", "Species_Name", "Tag", "Query_Accession"]]
    table.to_parquet(args.output_parquet, index=False)
    log.ok(f"{log.count(len(table))} alignments at E <= {defaults.RPSBLAST_E_VALUE:g}")
    return 0


if __name__ == "__main__":
    sys.exit(log.run_main(main))
