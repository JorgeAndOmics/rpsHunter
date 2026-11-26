"""Combine every species' and query's domain calls into the run's domain tables.

Writes, as Parquet and CSV:
  calls      every domain call of every query and species
  domains    domain instances: overlapping calls of one model, counted once
and one GFF3 file of the instances per configured species, also for a species
with no domain (header only), so that "no domains" is never confused with "not
run".

Usage (called by the Snakefile):
    python domain_tables.py --calls-out calls.parquet --domains-out domains.parquet \
        --gff-dir DIR CALLS.parquet [CALLS.parquet ...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

import defaults
import domain_instances
import gff3


def write_table(table: pd.DataFrame, parquet: Path) -> None:
    """Write a table as Parquet and, beside it, as CSV."""
    table.to_parquet(parquet, index=False)
    table.to_csv(parquet.with_suffix(".csv"), index=False)


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--calls-out", type=Path, required=True)
    parser.add_argument("--domains-out", type=Path, required=True)
    parser.add_argument("--gff-dir", type=Path, required=True)
    parser.add_argument("inputs", type=Path, nargs="+")
    args = parser.parse_args(argv)

    calls = pd.concat(
        [pd.read_parquet(path) for path in args.inputs], ignore_index=True
    )
    instances = domain_instances.instances(calls)
    write_table(calls, args.calls_out)
    write_table(instances, args.domains_out)
    for species in defaults.SPECIES:
        of_species = instances[instances["Species"] == species]
        text = "\n".join(gff3.lines(of_species)) + "\n"
        (args.gff_dir / f"{species}.gff3").write_text(text)
    print(
        f"{len(calls)} calls, {len(instances)} domain instances, {len(defaults.SPECIES)} GFF3 files"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
