"""Build the RPS-BLAST database from the CDD models of the configured families.

The database holds exactly the models listed under `domains:` in the config
(their `cdd` accessions), one `.smp` profile file each, taken from the unpacked
CDD release. A listed model without a profile file is an error: searching a
database that silently lacks a model would report the family as absent.

Usage (called by the Snakefile):
    python cdd_subset.py --smp-dir DIR --output-db PREFIX
Writes PREFIX.* (the database) and, beside it, smp_list.txt (the files used).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import defaults
import domain_families
import external
import log


def profile_files(accessions: list[str], smp_dir: Path) -> list[Path]:
    """Return the profile file of each accession, in order.

    Raises:
        FileNotFoundError: One or more accessions have no `<accession>.smp` in
            the folder; the message names all of them.
    """
    files = [smp_dir / f"{accession}.smp" for accession in accessions]
    missing = [path.stem for path in files if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            f"no CDD profile file in {smp_dir} for: {', '.join(missing)}. "
            "Check the accessions under 'domains:' against the CDD release in use."
        )
    return files


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--smp-dir", type=Path, required=True)
    parser.add_argument("--output-db", type=Path, required=True)
    log.add_scope_argument(parser)
    args = parser.parse_args(argv)
    log.job_logging("cdd_subset", args.scope)

    accessions = domain_families.cdd_accessions(defaults.FAMILIES)
    files = profile_files(accessions, args.smp_dir)
    listing = args.output_db.parent / "smp_list.txt"
    listing.write_text("".join(f"{path}\n" for path in files))
    command = [
        defaults.MAKEPROFILEDB_CMD,
        "-in", str(listing),
        "-out", str(args.output_db),
        "-title", "rpsHunter domain families",
    ]  # fmt: skip
    external.run_tool(command, tool="makeprofiledb")
    log.ok(f"database of {len(files)} CDD models: {args.output_db}")
    return 0


if __name__ == "__main__":
    sys.exit(log.run_main(main))
