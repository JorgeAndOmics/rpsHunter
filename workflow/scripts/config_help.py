"""Explain config keys from the reference documentation (`--config-help`).

The answer is read from `docs/configuration.md`, so the help cannot say
something else than the documentation: there is one text to keep right.

    ./rpsHunter --config-help              the sections
    ./rpsHunter --config-help rpsblast     one section
    ./rpsHunter --config-help rpsblast.e_value

A section is a `### `name` -- Title` heading under "Parameter Reference"; its
parameters are the rows of its tables; each table's header names the cells.

This module is pure: the launcher reads the file and prints the answer.
"""

from __future__ import annotations

import re
import textwrap
from typing import NamedTuple

_HEADING = re.compile(r"^### `(?P<name>[\w.]+)` -- (?P<title>.+)$")
_REFERENCE_START = "## Parameter Reference"
_ROW = re.compile(r"^\| `(?P<name>[\w.]+)` \|")
_UNESCAPED_BAR = re.compile(r"(?<!\\)\|")
_LINK = re.compile(r"\[([^\]]+)\]\([^)]+\)")
_WIDTH = 88


class UnknownKeyError(KeyError):
    """No such section or parameter; the message lists the ones that exist."""

    def __str__(self) -> str:
        return str(self.args[0])


class Section(NamedTuple):
    """One section of the reference."""

    title: str
    body: str


def _reference_part(reference: str) -> list[str]:
    """The lines under "Parameter Reference", up to the next second-level heading."""
    lines = reference.splitlines()
    starts = [i for i, line in enumerate(lines) if line.strip() == _REFERENCE_START]
    if not starts:
        return []
    rest = lines[starts[0] + 1 :]
    ends = [i for i, line in enumerate(rest) if line.startswith("## ")]
    return rest[: ends[0]] if ends else rest


def sections(reference: str) -> dict[str, Section]:
    """Return the config sections of the reference, by name, in document order."""
    found: dict[str, tuple[str, list[str]]] = {}
    current: list[str] = []
    for line in _reference_part(reference):
        if heading := _HEADING.match(line):
            current = []
            found[heading["name"]] = (heading["title"], current)
        else:
            current.append(line)
    return {
        name: Section(title, "\n".join(body).strip())
        for name, (title, body) in found.items()
    }


def _plain(cell: str) -> str:
    """A table cell as plain text: no backticks, links as their text, bars unescaped."""
    return _LINK.sub(r"\1", cell).replace("`", "").replace("\\|", "|").strip()


def _cells(line: str) -> list[str]:
    """The plain cells of one table line."""
    return [_plain(cell) for cell in _UNESCAPED_BAR.split(line.strip().strip("|"))]


def parameters(body: str) -> dict[str, dict[str, str]]:
    """Return the table rows of a section: parameter name to {column header: cell}.

    The header of the table each row belongs to names its cells, so tables
    with different columns (Default or Required, with or without Validation)
    are all read the same way.
    """
    rows: dict[str, dict[str, str]] = {}
    headers: list[str] = []
    for line in body.splitlines():
        if line.startswith("| Parameter |"):
            headers = _cells(line)
        elif found := _ROW.match(line):
            rows[found["name"]] = dict(zip(headers, _cells(line), strict=False))
    return rows


def _listing(known: dict[str, Section]) -> str:
    """The sections, one per line."""
    width = max(len(name) for name in known)
    rows = [f"  {name:<{width}}  {section.title}" for name, section in known.items()]
    return "\n".join(["Config sections (--config-help <section> for more):", *rows])


def _parameter(section: str, name: str, cells: dict[str, str]) -> str:
    """One parameter: each column of its row, then its description."""
    facts = {
        header: cell
        for header, cell in cells.items()
        if header not in ("Parameter", "Description")
    }
    width = max((len(header) for header in facts), default=0) + 1
    lines = [
        f"  {header.lower() + ':':<{width}}  {cell}" for header, cell in facts.items()
    ]
    meaning = textwrap.fill(
        cells.get("Description", ""),
        _WIDTH,
        initial_indent="  ",
        subsequent_indent="  ",
    )
    return "\n".join([f"{section}.{name}", *lines, meaning])


def explain(reference: str, key: str) -> str:
    """Return the help text for a key: "" (all sections), a section, or section.parameter.

    Raises:
        UnknownKeyError: The section or the parameter is not in the reference.
    """
    known = sections(reference)
    if not key:
        return _listing(known)
    name, _, parameter = key.partition(".")
    if name not in known:
        raise UnknownKeyError(
            f"no config section `{name}`; sections: {', '.join(known)}"
        )
    section = known[name]
    if not parameter:
        return f"{name}: {section.title}\n\n{section.body}"
    rows = parameters(section.body)
    if parameter not in rows:
        raise UnknownKeyError(
            f"no parameter `{parameter}` in `{name}`; it has: {', '.join(rows) or 'none'}"
        )
    return _parameter(name, parameter, rows[parameter])
