"""Combine every species' and query's domain calls into the run's domain tables.

Writes, as Parquet and CSV:
  calls      every domain call of every query and species
  domains    domain instances: overlapping calls of one model, counted once
  hmmer      HMMER's own domains (only when HMMER tables are given)
Calls and HMMER domains each carry Method_Agreement: whether the other method
reports the same family on the same hit at an overlapping position (ADR-010).
and one GFF3 file of the instances per configured species, also for a species
with no domain (header only), so that "no domains" is never confused with "not
run".

Usage (called by the Snakefile):
    python domain_tables.py --calls-out calls.parquet --domains-out domains.parquet \
        --gff-dir DIR --calls CALLS.parquet [...] \
        [--hmmer-out hmmer.parquet --hmmer HMMER.parquet ...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

import defaults
import domain_instances
import gff3
import method_agreement


def write_table(table: pd.DataFrame, parquet: Path) -> None:
    """Write a table as Parquet and, beside it, as CSV."""
    table.to_parquet(parquet, index=False)
    table.to_csv(parquet.with_suffix(".csv"), index=False)


def compared_families(hmmer_ran: bool) -> list[str]:
    """Families both methods searched for: those with CDD and Pfam models, if HMMER ran."""
    if not hmmer_ran:
        return []
    return [f.name for f in defaults.FAMILIES.values() if f.cdd and f.pfam]


def read_all(paths: list[Path]) -> pd.DataFrame:
    """Concatenate tables; no table gives an empty frame with the position columns."""
    if not paths:
        return pd.DataFrame(columns=["Tag", "Family", "Start", "End"])
    return pd.concat([pd.read_parquet(path) for path in paths], ignore_index=True)


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--calls-out", type=Path, required=True)
    parser.add_argument("--domains-out", type=Path, required=True)
    parser.add_argument("--gff-dir", type=Path, required=True)
    parser.add_argument("--hmmer-out", type=Path)
    parser.add_argument("--hmmer", type=Path, nargs="*", default=[])
    parser.add_argument("--calls", type=Path, nargs="+", required=True)
    args = parser.parse_args(argv)

    calls = read_all(args.calls)
    hmmer_domains = read_all(args.hmmer)
    compared = compared_families(hmmer_ran=bool(args.hmmer))
    instances = domain_instances.instances(calls)
    calls["Method_Agreement"] = method_agreement.agreement(
        calls, hmmer_domains, compared, "cdd_only"
    )
    if args.hmmer_out:
        hmmer_domains["Method_Agreement"] = method_agreement.agreement(
            hmmer_domains, calls, compared, "hmmer_only"
        )
        write_table(hmmer_domains, args.hmmer_out)
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
