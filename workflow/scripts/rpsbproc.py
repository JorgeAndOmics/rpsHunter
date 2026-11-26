"""Post-process one RPS-BLAST archive with rpsbproc into its text report.

rpsbproc turns raw alignments into domain calls: it picks representative hits
per region, marks each as specific, non-specific or superfamily using CDD's own
thresholds, and flags calls that miss an end of the model. It is given the same
E-value cut-off as the search, so nothing the search kept is dropped unseen and
nothing above the cut-off appears (ADR-006).

Usage (called by the Snakefile):
    python rpsbproc.py --input-asn ARCHIVE --output-txt REPORT

An empty archive (a species with nothing to search) gives an empty report
without running the tool. A failing rpsbproc stops the job.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import NamedTuple

import defaults
import external
import log


class Tool(NamedTuple):
    """How to run rpsbproc: the program, CDD's annotation files, and the cut-off."""

    program: str
    data_dir: Path
    max_evalue: float


def command(tool: Tool, archive: Path, report: Path) -> list[str]:
    """Return the rpsbproc command for one archive."""
    return [
        tool.program,
        "-i", str(archive),
        "-d", str(tool.data_dir),
        "-t", "both",  # domain hits and, where CDD annotates them, site features
        "-e", str(tool.max_evalue),
        "-o", str(report),
    ]  # fmt: skip


def process(tool: Tool, archive: Path, report: Path) -> bool:
    """Write the report for one archive; return whether rpsbproc had to run."""
    if archive.stat().st_size == 0:
        report.write_text("")
        return False
    external.run_tool(command(tool, archive, report), tool="rpsbproc")
    return True


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--input-asn", type=Path, required=True)
    parser.add_argument("--output-txt", type=Path, required=True)
    log.add_scope_argument(parser)
    args = parser.parse_args(argv)
    log.job_logging("rpsbproc", args.scope)

    tool = Tool(
        program=defaults.RPSBPROC_CMD,
        data_dir=defaults.PATH_DICT["RPSBPROC_DB"],
        max_evalue=defaults.RPSBLAST_E_VALUE,
    )
    ran = process(tool, args.input_asn, args.output_txt)
    outcome = "processed" if ran else "empty archive, empty report"
    log.ok(outcome)
    return 0


if __name__ == "__main__":
    sys.exit(log.run_main(main))
