"""Build one species' table of domain calls from the rpsbproc report.

Each call is given its family (from the config's `domains:` block), its
evidence tier, its own position on the genome, and the part of the model it
aligns to. The logic is in domain_calls.py; this script reads the files and
writes the table.

Usage (called by the Snakefile):
    python call_domains.py --species KEY --query-accession ACC \
        --report REPORT.txt --alignments RPSBLAST.parquet --output CALLS.parquet
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

import defaults
import domain_calls
import domain_families

CALL_COLUMNS = (
    "Species",
    "Species_Name",
    "Query_Accession",
    "Tag",
    "Sequence",
    "Strand",
    "Hit_Start",
    "Hit_End",
    "Family",
    "Model",
    "Model_Name",
    "Hit_Type",
    "Evidence",
    "Evalue",
    "Bitscore",
    "Start",
    "End",
    "Residue_From",
    "Residue_To",
    "Model_From",
    "Model_To",
    "Model_Length",
    "Incomplete",
    "PSSM_ID",
    "Aligned_PSSM",
)


def cdd_tables(data_dir: Path) -> tuple[dict[int, int], dict[str, str]]:
    """Read model lengths (by PSSM id) and superfamily membership from CDD's files."""
    lengths: dict[int, int] = {}
    with (data_dir / "cddid.tbl").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            lengths[int(fields[0])] = int(fields[4])
    superfamily_of: dict[str, str] = {}
    with (data_dir / "family_superfamily_links").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            superfamily_of[fields[0]] = fields[2]
    return lengths, superfamily_of


def family_of_each_call(models: pd.Series, superfamily_of: dict[str, str]) -> pd.Series:
    """Family of each call, from the accession of the model or superfamily it reports.

    Raises:
        ValueError: A call names a model that belongs to no configured family,
            which means the search database and the config disagree.
    """
    families = defaults.FAMILIES
    by_accession = domain_families.family_by_cdd(families)
    by_accession |= domain_families.family_by_superfamily(families, superfamily_of)
    unknown = sorted(set(models) - set(by_accession))
    if unknown:
        raise ValueError(
            f"calls for models outside every family under 'domains:': {unknown}"
        )
    return models.map(by_accession)


def build(report: str, alignments: pd.DataFrame, data_dir: Path) -> pd.DataFrame:
    """Return the call table (without provenance columns) for one report."""
    lengths, superfamily_of = cdd_tables(data_dir)
    calls = domain_calls.parse_report(report)
    calls = domain_calls.add_model_coordinates(calls, alignments, lengths)
    calls = domain_calls.place_on_genome(calls)
    calls["Family"] = family_of_each_call(calls["Model"], superfamily_of)
    calls["Evidence"] = calls["Evalue"].map(domain_calls.evidence_tier).astype("object")
    return calls


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--species", required=True)
    parser.add_argument("--query-accession", required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--alignments", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    alignments = pd.read_parquet(args.alignments)
    calls = build(
        args.report.read_text(), alignments, defaults.PATH_DICT["RPSBPROC_DB"]
    )
    calls["Species"] = args.species
    calls["Species_Name"] = defaults.SPECIES_DICT[args.species]
    calls["Query_Accession"] = args.query_accession
    calls[list(CALL_COLUMNS)].to_parquet(args.output, index=False)
    print(f"{args.species}: {len(calls)} domain calls")
    return 0


if __name__ == "__main__":
    sys.exit(main())
