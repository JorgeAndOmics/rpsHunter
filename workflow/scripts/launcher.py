"""Run the requested stages: one dry run for the guard, then one Snakemake call.

The command line is parsed by stages.py. This module finds the config, asks
Snakemake what would run, lets guard.py stop a run that would start a download
or a search nobody asked for, and then runs Snakemake once with every requested
target, showing the run through console.py, keeping all of its output in
`<logs>/runs/<time>.log` and writing the run's manifest (manifest.py). Before any of that, preflight.py checks the config, the tools and the
input files. The exit status is Snakemake's, 2 for a usage, config or input
problem, and EXIT_GUARD when the guard stopped the run.
"""

from __future__ import annotations

import os
import sys
import time
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any, NamedTuple

import yaml

import console
import external
import guard
import manifest
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
Stream = Callable[..., int]
"""external.run_streaming, or a stand-in."""
Check = Callable[[stages.Request, Path], preflight.Findings]
"""preflight.check, or a stand-in."""


class Tools(NamedTuple):
    """What the launcher calls outside itself; replaced in tests."""

    run: Run = external.run_plain
    stream: Stream = external.run_streaming
    check: Check = preflight.check
    record: Callable[..., Path] = manifest.record


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


def _settings(config: Path) -> dict[str, Any]:
    """The config as a mapping (the checks have already found it sound)."""
    return yaml.safe_load(config.read_text()) or {}


def version() -> str:
    """The version this copy of rpsHunter declares in CITATION.cff."""
    declared = yaml.safe_load((REPO / "CITATION.cff").read_text())
    return str(declared.get("version", "unknown"))


def _folder(settings: dict[str, Any], root: str, name: str) -> Path:
    """Return a folder under one of the config's roots.

    These folders are made when the Snakefile is read (defaults.ensure_dirs),
    which the dry run has already done.
    """
    return Path(settings["root"][root]).resolve() / name


def _screen(request: stages.Request, settings: dict[str, Any]) -> console.Console:
    """The console for this run, at the verbosity asked for."""
    chosen = request.verbosity or (settings.get("display") or {}).get(
        "verbosity", "normal"
    )
    encoding = sys.stdout.encoding or "ascii"
    style = console.style_for(os.environ, sys.stdout.isatty(), encoding)

    def write(text: str) -> None:
        sys.stdout.write(text)
        sys.stdout.flush()

    return console.Console(chosen, style, write)


class _Session(NamedTuple):
    """One run about to start: what to call, where, and what it was asked to do."""

    request: stages.Request
    config: Path
    command: list[str]
    environment: dict[str, str]
    jobs: dict[str, int]


def _execute(session: _Session, screen: console.Console, tools: Tools) -> int:
    """Run Snakemake: show the run, keep every line, write the manifest."""
    settings = _settings(session.config)
    stamp = time.strftime("%Y-%m-%d_%H%M%S")
    kept = _folder(settings, "logs_root_folder", "runs") / f"{stamp}.log"
    started, clock = time.strftime("%Y-%m-%dT%H:%M:%S"), time.monotonic()
    with kept.open("w") as record:

        def take(line: str) -> None:
            record.write(line + "\n")
            screen.feed(line)

        options = {"cwd": WORKFLOW, "env": session.environment}
        status = int(tools.stream(session.command, take, **options))
    seconds = time.monotonic() - clock
    screen.summary(status, seconds)
    if screen.warnings:
        kept.with_suffix(".warnings.txt").write_text("\n".join(screen.warnings) + "\n")
    run = manifest.Run(
        command=[Path(sys.argv[0]).name, *sys.argv[1:]],
        stages=[stage.flag for stage in session.request.stages],
        started=started,
        seconds=seconds,
        exit_status=status,
        jobs=session.jobs,
    )
    folder = _folder(settings, "results_root_folder", "manifests")
    print(
        f"Run log: {kept}\nManifest: {tools.record(run, session.config, folder, stamp)}"
    )
    return status


def _launch(request: stages.Request, environment: dict[str, str], tools: Tools) -> int:
    """Check, plan and run; return the exit status."""
    config = config_path(request.config, environment)
    if not _checked(request, config, tools.check):
        return EXIT_USAGE
    environment = {**environment, CONFIG_VARIABLE: str(config)}
    cores = cores_for(request, config)
    command = snakemake_command(request, cores)
    status, plan = _plan(command, environment, tools.run)
    if status != 0:
        print(plan, file=sys.stderr)
        return int(status)
    counts = guard.job_counts(plan)
    if not counts:
        print("Nothing to do: every requested output is up to date.")
        return 0
    print(f"Jobs ({config.name}):\n{_job_table(counts)}", flush=True)
    if verdict := _guard_verdict(request, plan, counts):
        print(verdict, file=sys.stderr)
        return EXIT_GUARD
    if request.dry_run:
        return 0
    screen = _screen(request, _settings(config))
    screen.banner(
        version(), str(config), cores, [stage.flag for stage in request.stages]
    )
    return _execute(
        _Session(request, config, command, environment, counts), screen, tools
    )


DEFAULT_TOOLS = Tools()


def main(argv: Sequence[str] | None = None, tools: Tools = DEFAULT_TOOLS) -> int:
    """Command-line entry point; returns the process exit status."""
    request = stages.parse(sys.argv[1:] if argv is None else argv)
    try:
        return _launch(request, dict(os.environ), tools)
    except LaunchError as problem:
        print(f"rpsHunter: {problem}", file=sys.stderr)
        return EXIT_USAGE


if __name__ == "__main__":
    sys.exit(main())
