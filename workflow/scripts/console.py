"""What the launcher shows of a run.

Jobs write plain lines (log.py); Snakemake adds its own. The launcher feeds
every line of that stream to a Console, which decides what reaches the screen:

  quiet    errors, warnings, the progress line, the summary
  normal   also one line per finished job
  verbose  every line, Snakemake's included

Colour and symbols carry the level and nothing else, and are used only on a
terminal that can show them (NO_COLOR and FORCE_COLOR are honoured). Nothing
here writes files: the launcher keeps the full stream in the run log.

This module is pure: it turns lines into lines.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Mapping, Sequence
from typing import NamedTuple

import log
import stages

VERBOSITIES = stages.VERBOSITIES

_PROGRESS = re.compile(r"^(\d+) of (\d+) steps \((\d+)%\) done")
_FAILED_RULE = re.compile(r"^Error in rule (\w+):")
_FAILED_LOG = re.compile(r"^\s+log: (\S+)")
_CLEAR = "\r\x1b[K"
_RESET = "\x1b[0m"
_MARKS = {  # level: (symbol, word, colour)
    "OK": ("✓", "OK", "\x1b[32m"),
    "WARN": ("⚠", "WARN", "\x1b[33m"),
    "ERROR": ("✗", "ERROR", "\x1b[1;31m"),
}
_SHOWN_FROM = {"ERROR": "quiet", "WARN": "quiet", "OK": "normal"}
"""The least verbosity at which a level is shown; other levels need `verbose`."""


class Style(NamedTuple):
    """What the terminal can show."""

    colour: bool
    symbols: bool
    live: bool
    """A terminal: a progress line can be redrawn in place."""


class Failure(NamedTuple):
    """A job Snakemake reported as failed."""

    rule: str
    log: str


def style_for(
    environment: Mapping[str, str], is_terminal: bool, encoding: str
) -> Style:
    """Choose the style for an output stream.

    Colour on a terminal, unless NO_COLOR is set; FORCE_COLOR turns it on
    anywhere. Symbols when the encoding has them.
    """
    forced = bool(environment.get("FORCE_COLOR"))
    colour = forced or (is_terminal and not environment.get("NO_COLOR"))
    return Style(colour, "utf" in encoding.lower(), is_terminal)


def duration(seconds: float) -> str:
    """A length of time as people say it."""
    whole = round(seconds)
    if whole < 60:
        return f"{whole} s"
    minutes, rest = divmod(whole, 60)
    if minutes < 60:
        return f"{minutes} min {rest} s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours} h {minutes} min"


class Console:
    """Follows one run through the lines it prints."""

    def __init__(
        self, verbosity: str, style: Style, write: Callable[[str], None]
    ) -> None:
        """Start a screen.

        Args:
            verbosity: One of VERBOSITIES.
            style: What the terminal can show.
            write: Receives each piece of text to put on the screen, line ends
                and cursor moves included.
        """
        self._level = VERBOSITIES.index(verbosity)
        self._style = style
        self._write = write
        self._errors: set[str] = set()
        self._pending_rule = ""
        self._status = ""
        self.done = 0
        self.total = 0
        self.warnings: list[str] = []
        self.failures: list[Failure] = []

    # --- what goes on the screen ---------------------------------------------------

    def _show(self, text: str) -> None:
        """Write one line, keeping the progress line (if any) below it."""
        if self._status:
            self._write(_CLEAR + text + "\n")
            self._write("\r" + self._status)
        else:
            self._write(text + "\n")

    def _mark(self, level: str) -> str:
        """The symbol or word for a level, coloured if the terminal allows."""
        symbol, word, colour = _MARKS[level]
        text = symbol if self._style.symbols else f"{word:<5}"
        return f"{colour}{text}{_RESET}" if self._style.colour else text

    def _wanted(self, level: str) -> bool:
        """Whether lines of this level are shown at the current verbosity."""
        least = _SHOWN_FROM.get(level, "verbose")
        return self._level >= VERBOSITIES.index(least)

    def banner(
        self, version: str, config: str, cores: int, flags: Sequence[str]
    ) -> None:
        """Say what is about to run."""
        self._show(
            f"rpsHunter {version} | config {config} | {cores} cores | {' '.join(flags)}"
        )

    # --- reading the run -----------------------------------------------------------

    def feed(self, line: str) -> None:
        """Take one line of the run's output."""
        ours = log.parse(line)
        if ours is not None:
            self._job_line(ours, line)
            return
        self._snakemake_line(line)
        if self._level == VERBOSITIES.index("verbose"):
            self._show(line)

    def _job_line(self, line: log.Line, raw: str) -> None:
        """One of our own lines: remember what matters and show it if wanted."""
        text = f"{line.step} {line.scope}: {line.message}"
        if line.level == "ERROR":
            if text in self._errors:
                return  # Snakemake prints a failed job's log again
            self._errors.add(text)
        if line.level == "WARN" and text not in self.warnings:
            self.warnings.append(text)
        if self._level == VERBOSITIES.index("verbose"):
            self._show(raw)
        elif self._wanted(line.level):
            self._show(
                f"{self._mark(line.level)} {line.step}  {line.scope}  {line.message}"
            )

    def _snakemake_line(self, line: str) -> None:
        """One of Snakemake's lines: follow progress and failures."""
        if found := _PROGRESS.match(line):
            self.done, self.total = int(found[1]), int(found[2])
            if self._style.live:
                self._status = f"[{self.done}/{self.total}] {found[3]}%"
                self._write("\r" + self._status)
        elif found := _FAILED_RULE.match(line):
            self._pending_rule = found[1]
        elif self._pending_rule and (found := _FAILED_LOG.match(line)):
            self.failures.append(Failure(self._pending_rule, found[1]))
            self._pending_rule = ""

    # --- the end ---------------------------------------------------------------------

    def _warning_count(self) -> str:
        """The number of warnings, in words."""
        count = len(self.warnings)
        if count == 0:
            return "no warnings"
        return "1 warning" if count == 1 else f"{count} warnings"

    def summary(self, status: int, seconds: float) -> None:
        """Say how the run ended: jobs, time, warnings, and what failed where."""
        if self._status:
            self._write(_CLEAR)
            self._status = ""
        jobs = f"{self.done} of {self.total} jobs"
        if status == 0:
            self._show(f"Done in {duration(seconds)}: {jobs}, {self._warning_count()}.")
        else:
            self._show(f"Failed after {duration(seconds)}: {jobs} finished.")
            failed = [f"  {failure.rule}: {failure.log}" for failure in self.failures]
            for line in failed or ["  See the messages above."]:
                self._show(line)
        for warning in self.warnings:
            self._show(f"  {warning}")
