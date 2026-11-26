"""Check a genome FASTA against the checksum pinned for its species.

Run once per genome, immediately before its BLAST database is built: every
genome the pipeline analyses passes through that step exactly once, so a file
that is not the one the config names is caught before anything is computed
from it. A species pinned by assembly accession only has nothing to compare;
its checksum is printed so that it can be added to the config.

Usage: genome_check.py --species KEY --fasta PATH
Exit status 1, with both checksums on standard error, when they differ.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import defaults
import species_config

_BLOCK_BYTES = 1 << 20  # read in 1 MiB blocks: genomes are several gigabytes


def file_sha256(path: Path) -> str:
    """Return the SHA-256 of a file as 64 hexadecimal characters."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(_BLOCK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--species", required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    args = parser.parse_args(argv)

    entry = defaults.SPECIES_ENTRIES[args.species]
    actual = file_sha256(args.fasta)
    problem = species_config.checksum_problem(args.species, entry, actual)
    if problem is not None:
        print(problem, file=sys.stderr)
        return 1
    pinned = "matches the config" if entry.sha256 else "not pinned in the config"
    print(f"{args.species}: sha256 {actual} ({pinned})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
