"""Check that downloaded CDD files are the release the config names.

NCBI serves only the current CDD release, so a download cannot ask for a
version. It can refuse to go on with a different one: E-values are put on the
scale of the configured release (ADR-006), and models change between releases.

Usage (called by the Snakefile after downloading `cdd.info`):
    python cdd_release.py --info cdd.info
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import defaults
import log

_VERSION = re.compile(r"^cdd version (\S+)", re.IGNORECASE)


def declared_version(info: str) -> str:
    """Return the release a `cdd.info` file declares on its first line.

    Raises:
        log.PipelineError: The text is not a cdd.info file.
    """
    found = _VERSION.match(info.strip())
    if found is None:
        raise log.PipelineError(
            "cdd.info does not state a CDD version",
            hint="the download may have failed; look at the file",
        )
    return found[1]


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--info", type=Path, required=True)
    log.add_scope_argument(parser)
    args = parser.parse_args(argv)
    log.job_logging("cdd_release", args.scope)

    served = declared_version(args.info.read_text())
    if served != defaults.CDD_VERSION:
        raise log.PipelineError(
            f"NCBI now serves CDD {served}; the config says {defaults.CDD_VERSION}",
            hint="NCBI offers only its current release. To use it, set "
            "databases.cdd.version and databases.cdd.residues to that release "
            "(docs/configuration.md); results will then be on its scale",
        )
    log.ok(f"CDD release {served}, as configured")
    return 0


if __name__ == "__main__":
    sys.exit(log.run_main(main))
