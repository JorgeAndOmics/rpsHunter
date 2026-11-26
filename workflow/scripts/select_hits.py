"""Select one species' tBLASTn hits by the quality gate and export them.

A hit is selected when it meets the four thresholds of the config's `blast:`
section (identity, E-value, alignment length, bit score). Nothing else decides
which hits go on to the domain searches: HMMER is run on the selected sequences
as independent evidence, not as a filter (ADR-010).

Writes three files:
  audit     every hit, with Quality_Pass and Selected
  selected  the selected hits (Parquet, and CSV beside it)
  FASTA     their translated sequences, named `<seq>:<start>-<end>|tag:<Tag>`

Usage (called by the Snakefile):
    python select_hits.py --input HITS.parquet --output-fasta FASTA \
        --output-audit AUDIT.parquet --output-selected SELECTED.parquet
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import NamedTuple

import pandas as pd

import defaults
import hit_table

SEQUENCE = "Subject Sequence"
MASKED_STOP = "X"


class QualityThresholds(NamedTuple):
    """The quality gate: a hit must meet all four."""

    min_identity: float
    max_evalue: float
    min_length: int
    min_bitscore: float


def quality_pass(hits: pd.DataFrame, thresholds: QualityThresholds) -> pd.Series:
    """Say, for each hit, whether it meets every threshold (boundaries included)."""
    return (
        (hits["Pct Identity"] >= thresholds.min_identity)
        & (hits["E-value"] <= thresholds.max_evalue)
        & (hits["Alignment Length"] >= thresholds.min_length)
        & (hits["Bit Score"] >= thresholds.min_bitscore)
    )


def searchable_sequence(aligned: str) -> str:
    """Return a hit's translation ready for a domain search.

    Alignment gaps are dropped: they stand for query residues that have no
    bases in the genome. A stop codon is masked with X instead of being removed,
    because residue n must stay 3(n-1) bases from the hit's start for domain
    positions to be right (ADR-007).
    """
    return aligned.replace(hit_table.GAP, "").replace(hit_table.STOP, MASKED_STOP)


def fasta_text(hits: pd.DataFrame) -> str:
    """Return the FASTA of the hits' searchable sequences; empty ones are left out."""
    records = []
    columns = ["Subject ID", "S. Start", "S. End", "Tag", SEQUENCE]
    for sequence_id, start, end, tag, aligned in hits[columns].itertuples(index=False):
        residues = searchable_sequence(str(aligned))
        if residues:
            records.append(f">{sequence_id}:{start}-{end}|tag:{tag}\n{residues}\n")
    return "".join(records)


def configured_thresholds() -> QualityThresholds:
    """The quality gate of the current config."""
    return QualityThresholds(
        min_identity=defaults.PERC_IDENTITY_THRESHOLD,
        max_evalue=defaults.E_VALUE_THRESHOLD,
        min_length=defaults.SEQ_LENGTH_THRESHOLD,
        min_bitscore=defaults.BITSCORE_THRESHOLD,
    )


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--output-audit", type=Path, required=True)
    parser.add_argument("--output-selected", type=Path, required=True)
    args = parser.parse_args(argv)

    hits = pd.read_parquet(args.input)
    hits["Quality_Pass"] = quality_pass(hits, configured_thresholds())
    hits["Selected"] = hits["Quality_Pass"]
    hits.to_parquet(args.output_audit, index=False)
    selected = hits[hits["Selected"]].drop(columns=["Quality_Pass", "Selected"])
    selected.to_parquet(args.output_selected, index=False)
    selected.to_csv(args.output_selected.with_suffix(".csv"), index=False)
    args.output_fasta.write_text(fasta_text(selected))
    print(f"{len(selected)} of {len(hits)} hits selected")
    return 0


if __name__ == "__main__":
    sys.exit(main())
