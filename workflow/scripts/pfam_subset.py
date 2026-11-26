"""Extract the families' Pfam profiles from the full library into one small file.

The config lists Pfam accessions without version (PF00069); the library names
its profiles (Pkinase) and versions their accessions (PF00069.30). The library's
NAME and ACC lines are read to translate, then hmmfetch copies the profiles.

Usage (called by the Snakefile):
    python pfam_subset.py --library Pfam-A.hmm --output subset.hmm

Fails if a configured accession is not in the library (DomainConfigError) or
if hmmfetch fails (external.ToolError).
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from collections.abc import Iterable, Iterator
from pathlib import Path

import defaults
import domain_families
import external
import log


def library_index(lines: Iterable[str]) -> Iterator[tuple[str, str]]:
    """Yield (profile name, versioned accession) for each profile of a library."""
    name = ""
    for line in lines:
        if line.startswith("NAME"):
            name = line.split()[1]
        elif line.startswith("ACC"):
            yield name, line.split()[1]


def build(library: Path, subset: Path) -> list[str]:
    """Write the configured families' profiles to `subset`; return their names."""
    with library.open() as handle:
        names = domain_families.pfam_profile_names(
            defaults.FAMILIES, library_index(handle)
        )
    wanted = list(names.values())
    with tempfile.TemporaryDirectory() as folder:
        keys = Path(folder) / "profiles.keys"
        keys.write_text("\n".join(wanted) + "\n")
        command = ["hmmfetch", "-f", str(library), str(keys)]
        subset.write_text(external.run_tool(command, tool="hmmfetch"))
    return wanted


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    log.add_scope_argument(parser)
    args = parser.parse_args(argv)
    log.job_logging("pfam_subset", args.scope)

    names = build(args.library, args.output)
    log.ok(f"{len(names)} Pfam profiles: {', '.join(names)}")
    return 0


if __name__ == "__main__":
    sys.exit(log.run_main(main))
