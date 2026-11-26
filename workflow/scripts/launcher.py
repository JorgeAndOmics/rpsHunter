"""Run the requested stages: one dry run for the guard, then one Snakemake call.

The command line is parsed by stages.py. This module finds the config, asks
Snakemake what would run, lets guard.py stop a run that would start a download
or a search nobody asked for, and then runs Snakemake once with every requested
target. Before any of that, preflight.py checks the config, the tools and the
input files. The exit status is Snakemake's, 2 for a usage, config or input
problem, and EXIT_GUARD when the guard stopped the run.
"""

from __future__ import annotations

import os
import sys
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

import yaml

import external
import guard
import preflight
import stages

REPO = Path(__file__).resolve().parents[2]
WORKFLOW = REPO / "workflow"
DEFAULT_CONFIG = REPO / "data" / "config" / "config.yaml"
CONFIG_VARIABLE = "RPSHUNTER_CONFIG"
EXIT_USAGE = 2
EXIT_GUARD = 3

Run = Callable[..., Any]
"""external.run_plain, or a stand-in with the same call and result shape."""
Check = Callable[[stages.Request, Path], preflight.Findings]
"""preflight.check, or a stand-in."""


class LaunchError(Exception):
    """The run cannot start; the message says why in one line."""


def config_path(flag: str | None, environment: Mapping[str, str]) -> Path:
    """Return the absolute path of the run's config: flag, else environment, else default.

    Raises:
        LaunchError: The file does not exist.
    """
    chosen = flag or environment.get(CONFIG_VARIABLE)
    path = Path(chosen).resolve() if chosen else DEFAULT_CONFIG
    if not path.is_file():
        raise LaunchError(f"config file not found: {path}")
    return path


def cores_for(request: stages.Request, config: Path) -> int:
    """Return the cores to use: --cores, else execution.num_cores, else 1."""
    if request.cores is not None:
        return request.cores
    settings = yaml.safe_load(config.read_text()) or {}
    return int((settings.get("execution") or {}).get("num_cores", 1))


def snakemake_command(request: stages.Request, cores: int) -> list[str]:
    """Return the one Snakemake call that builds every requested target."""
    return [
        sys.executable, "-m", "snakemake",
        "--snakefile", str(WORKFLOW / "Snakefile"),
        *request.targets,
        "--cores", str(cores),
        "--rerun-incomplete",
        "--show-failed-logs",  # a failed job's log is printed, not just named
        "--resources", "ncbi_api=1",  # one NCBI request at a time (HTTP 429 otherwise)
        *request.snakemake_args,
    ]  # fmt: skip


def _plan(command: list[str], environment: dict[str, str], run: Run) -> tuple[int, str]:
    """Do a dry run; return its exit status and everything it printed."""
    done = run(
        [*command, "--dry-run"],
        cwd=WORKFLOW,
        env=environment,
        capture_output=True,
        text=True,
    )
    return done.returncode, done.stdout + done.stderr


def _job_table(counts: Mapping[str, int]) -> str:
    """The jobs a run would start, one rule per line, with the total."""
    width = max(len(rule) for rule in counts)
    rows = [f"  {rule:<{width}}  {jobs:>5}" for rule, jobs in counts.items()]
    return "\n".join([*rows, f"  {'total':<{width}}  {sum(counts.values()):>5}"])


def _guard_verdict(
    request: stages.Request, plan: str, counts: Mapping[str, int]
) -> str:
    """The guard's report, or an empty string when the run may go ahead."""
    if request.allow_heavy:
        return ""
    stopped = guard.blocked(counts, request.allowed_heavy)
    return guard.report(stopped, counts, guard.first_reasons(plan)) if stopped else ""


def _checked(request: stages.Request, config: Path, check: Check) -> bool:
    """Run the checks, show what they found, and say whether the run may start."""
    findings = check(request, config)
    for note in findings.notes:
        print(f"note: {note}")
    for problem in findings.problems:
        print(f"problem: {problem}", file=sys.stderr)
    if findings.problems:
        count = len(findings.problems)
        print(f"Nothing was run: {count} problem(s) in {config}.", file=sys.stderr)
    return not findings.problems


def _launch(
    request: stages.Request, environment: dict[str, str], run: Run, check: Check
) -> int:
    """Check, plan and run; return the exit status."""
    config = config_path(request.config, environment)
    if not _checked(request, config, check):
        return EXIT_USAGE
    environment = {**environment, CONFIG_VARIABLE: str(config)}
    command = snakemake_command(request, cores_for(request, config))
    status, plan = _plan(command, environment, run)
    if status != 0:
        print(plan, file=sys.stderr)
        return int(status)
    counts = guard.job_counts(plan)
    if not counts:
        print("Nothing to do: every requested output is up to date.")
        return 0
    print(f"Jobs ({config.name}):\n{_job_table(counts)}")
    if verdict := _guard_verdict(request, plan, counts):
        print(verdict, file=sys.stderr)
        return EXIT_GUARD
    if request.dry_run:
        return 0
    return int(run(command, cwd=WORKFLOW, env=environment).returncode)


def main(
    argv: Sequence[str] | None = None,
    run: Run = external.run_plain,
    check: Check = preflight.check,
) -> int:
    """Command-line entry point; returns the process exit status."""
    request = stages.parse(sys.argv[1:] if argv is None else argv)
    try:
        return _launch(request, dict(os.environ), run, check)
    except LaunchError as problem:
        print(f"rpsHunter: {problem}", file=sys.stderr)
        return EXIT_USAGE


if __name__ == "__main__":
    sys.exit(main())
