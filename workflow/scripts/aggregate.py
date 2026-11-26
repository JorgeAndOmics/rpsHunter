"""Concatenate tables of one kind (per species, or per query) into one.

Usage (called by the Snakefile):
    python aggregate.py --output ALL.parquet TABLE.parquet [TABLE.parquet ...]

Writes the Parquet file and a CSV beside it. Rows keep the order of the inputs;
nothing is filtered or changed.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

import log


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("inputs", type=Path, nargs="+")
    log.add_scope_argument(parser)
    args = parser.parse_args(argv)
    log.job_logging("aggregate", args.scope)

    combined = pd.concat(
        [pd.read_parquet(path) for path in args.inputs], ignore_index=True
    )
    combined.to_parquet(args.output, index=False)
    combined.to_csv(args.output.with_suffix(".csv"), index=False)
    log.ok(f"{log.count(len(combined))} rows from {len(args.inputs)} tables")
    return 0


if __name__ == "__main__":
    sys.exit(log.run_main(main))
