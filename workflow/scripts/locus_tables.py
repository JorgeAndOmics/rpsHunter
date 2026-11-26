"""Chain the run's gated hits into loci and report what each locus carries.

Writes, as Parquet and CSV:
  loci            one row per locus: position, hits, queries, families present
  locus_families  one row per locus and family: evidence tier, best E-value,
                  model coverage across the locus's exons, copies
The logic is in loci.py; this script reads the per-species tables and writes.

Usage (called by the Snakefile):
    python locus_tables.py --loci-out loci.parquet --families-out locus_families.parquet \
        --hits HITS.parquet [...] --calls CALLS.parquet [...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

import defaults
import loci

HIT_COLUMNS = {
    "Subject ID": "Sequence",
    "Pct Identity": "Identity",
    "Bit Score": "Bitscore",
}
"""Columns of the gated-hit tables, renamed to the names loci.py uses."""


def read_hits(paths: list[Path]) -> pd.DataFrame:
    """Read gated hits and express each as a start <= end interval."""
    hits = pd.concat([pd.read_parquet(path) for path in paths], ignore_index=True)
    hits = hits.rename(columns=HIT_COLUMNS)
    ends = hits[["S. Start", "S. End"]]
    return hits.assign(Start=ends.min(axis=1), End=ends.max(axis=1))


def write_table(table: pd.DataFrame, parquet: Path) -> None:
    """Write a table as Parquet and, beside it, as CSV."""
    table.to_parquet(parquet, index=False)
    table.to_csv(parquet.with_suffix(".csv"), index=False)


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--loci-out", type=Path, required=True)
    parser.add_argument("--families-out", type=Path, required=True)
    parser.add_argument("--hits", type=Path, nargs="+", required=True)
    parser.add_argument("--calls", type=Path, nargs="+", required=True)
    args = parser.parse_args(argv)

    chained = loci.chain(read_hits(args.hits), gap=defaults.LOCUS_GAP)
    calls = pd.concat([pd.read_parquet(path) for path in args.calls], ignore_index=True)
    evidence = loci.family_evidence(chained, calls)
    table = loci.locus_table(chained, evidence, family_order=list(defaults.FAMILIES))
    write_table(table, args.loci_out)
    write_table(evidence, args.families_out)
    print(
        f"{len(table)} loci from {len(chained)} gated hits (gap {defaults.LOCUS_GAP} bp)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
