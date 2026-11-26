"""One line shape for every message a job writes.

    14:02:40 WARN select_hits query_1/Homo_sapiens | 3 hits without sequence
    time     level step        scope                message

The step is the rule's name without `_species`; the scope is `<query>/<species>`,
or `all` for a job over the whole run. With many jobs writing at once, step and
scope are what make interleaved lines readable. Scripts write plain lines; only
the launcher adds colour and symbols (console.py), and it finds our lines among
Snakemake's and the tools' output with `parse`.

Levels: ERROR, the job cannot produce a correct output; WARN, the output exists
but a person should look; OK, the one headline a job prints when it finishes;
INFO, detail.
"""

from __future__ import annotations

import logging
import re
import sys
import time
import traceback
from collections.abc import Callable
from typing import NamedTuple, TextIO

OK = 25  # between INFO (20) and WARNING (30)
LEVEL_WORDS = {
    logging.DEBUG: "DEBUG",
    logging.INFO: "INFO",
    OK: "OK",
    logging.WARNING: "WARN",
    logging.ERROR: "ERROR",
    logging.CRITICAL: "ERROR",
}
EVERYTHING = "all"

_LINE = re.compile(
    r"^(?P<time>\d\d:\d\d:\d\d) (?P<level>DEBUG|INFO|OK|WARN|ERROR) "
    r"(?P<step>\S+) (?P<scope>\S+) \| (?P<message>.*)$"
)


_STREAM: list[TextIO] = [sys.stdout]
"""Where the current job's lines go; set by job_logging, read by run_main."""


class Line(NamedTuple):
    """One message, taken apart."""

    time: str
    level: str
    step: str
    scope: str
    message: str


class PipelineError(Exception):
    """A failure we know how to explain: what went wrong, and what to do about it."""

    def __init__(self, message: str, hint: str = "") -> None:
        super().__init__(message)
        self.hint = hint

    def __str__(self) -> str:
        text = super().__str__()
        return f"{text}. Fix: {self.hint}" if self.hint else text


class _LineFormatter(logging.Formatter):
    """Formats a record as one line of the shape above."""

    def __init__(self, step: str, scope: str, clock: Callable[[], str]) -> None:
        super().__init__()
        self._step, self._scope, self._clock = step, scope, clock

    def format(self, record: logging.LogRecord) -> str:
        word = LEVEL_WORDS.get(record.levelno, "INFO")
        message = " ".join(record.getMessage().split())
        return f"{self._clock()} {word} {self._step} {self._scope} | {message}"


def _now() -> str:
    """The local time as HH:MM:SS."""
    return time.strftime("%H:%M:%S")


def job_logging(
    step: str,
    scope: str = EVERYTHING,
    stream: TextIO | None = None,
    clock: Callable[[], str] = _now,
) -> None:
    """Make the root logger write this job's lines; call once, first thing in main.

    Args:
        step: The rule's name without `_species`.
        scope: `<query>/<species>`, or EVERYTHING for a job over the whole run.
        stream: Where lines go (default: standard output, which the rule sends
            to the job's log and to the launcher).
        clock: Returns the time shown; replaced in tests.
    """
    target = stream or sys.stdout
    _STREAM[0] = target
    handler = logging.StreamHandler(target)
    handler.setFormatter(_LineFormatter(step, scope, clock))
    root = logging.getLogger()
    root.handlers[:] = [handler]
    root.setLevel(logging.INFO)


def ok(message: str) -> None:
    """Write the job's headline: what it produced, with its key number."""
    logging.getLogger().log(OK, message)


def count(number: int) -> str:
    """A count with thousands separators, as messages show them."""
    return f"{number:,}"


def parse(text: str) -> Line | None:
    """Take one of our lines apart; None for anything else."""
    found = _LINE.match(text)
    return Line(**found.groupdict()) if found else None


def run_main(main: Callable[[], int]) -> int:
    """Run a script's main and turn a failure into one ERROR line and status 1.

    A PipelineError is reported as its message and fix. Any other exception is
    reported in one line and its traceback follows, for the job's log.
    """
    try:
        return main()
    except PipelineError as problem:
        logging.getLogger().error(str(problem))
    except Exception as problem:  # the job log must say why a job died, whatever it was
        logging.getLogger().error(f"{type(problem).__name__}: {problem}")
        traceback.print_exc(file=_STREAM[0])
    return 1
