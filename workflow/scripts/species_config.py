"""The `species:` block of the configuration.

Each species is a key (used in file names) with a display name and a pin that
says exactly which genome is meant:

    species:
      Mus_musculus:
        name: Mus musculus
        assembly: GCA_000001635.9      # a versioned NCBI assembly accession
      My_isolate:
        name: My isolate
        sha256: <64 hex characters>    # a local FASTA, pinned by its checksum

At least one pin is required. A taxon name is not accepted, because NCBI resolves
it to whichever assembly is newest on the day, so the same config would mean
different genomes at different times.

This module is pure: it validates and returns values, and reads nothing.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import NamedTuple

ASSEMBLY_ACCESSION = re.compile(r"GC[AF]_\d{9}\.\d+")
"""A GenBank (GCA) or RefSeq (GCF) assembly accession with its version."""

SHA256_DIGEST = re.compile(r"[0-9a-f]{64}")

_NEW_FORM = (
    "write each species as\n"
    "  <key>:\n"
    "    name: <display name>\n"
    "    assembly: <versioned accession, e.g. GCA_000001635.9>"
)


class SpeciesConfigError(ValueError):
    """The `species:` block is not usable; the message lists every problem."""


class SpeciesEntry(NamedTuple):
    """One species: how to show it and which genome it is."""

    name: str
    assembly: str | None
    sha256: str | None


def _entry_problem(key: str, raw: object) -> str | None:
    """Return what is wrong with one species entry, or None if it is usable."""
    if not isinstance(raw, Mapping):
        return f"{key}: the name-only form is no longer accepted; {_NEW_FORM}"
    assembly, sha256 = raw.get("assembly"), raw.get("sha256")
    checks = (
        (not raw.get("name"), "has no display name (name:)"),
        (assembly is None and sha256 is None, f"names no genome; {_NEW_FORM}"),
        (
            assembly is not None and not ASSEMBLY_ACCESSION.fullmatch(str(assembly)),
            f"assembly '{assembly}' is not a versioned accession (GCA_/GCF_ with .N)",
        ),
        (
            sha256 is not None and not SHA256_DIGEST.fullmatch(str(sha256)),
            "sha256 is not a 64-character hexadecimal digest",
        ),
    )
    return next((f"{key}: {message}" for failed, message in checks if failed), None)


def parse_species(raw: object) -> dict[str, SpeciesEntry]:
    """Validate the `species:` block and return its entries in config order.

    Args:
        raw: The block as loaded from YAML; None or empty means no species.
            Anything other than a mapping is refused.

    Returns:
        Species key to entry, in the order of the config.

    Raises:
        SpeciesConfigError: An entry lacks a name or a pin, a pin is malformed,
            or an entry uses the old name-only form. Every problem is listed.
    """
    if raw is not None and not isinstance(raw, Mapping):
        raise SpeciesConfigError(
            f"species: must be a mapping of species keys; {_NEW_FORM}"
        )
    problems: list[str] = []
    parsed: dict[str, SpeciesEntry] = {}
    for key, value in (raw or {}).items():
        problem = _entry_problem(key, value)
        if problem is not None:
            problems.append(problem)
        elif isinstance(value, Mapping):
            parsed[key] = SpeciesEntry(
                name=str(value["name"]),
                assembly=_optional_text(value.get("assembly")),
                sha256=_optional_text(value.get("sha256")),
            )
    if problems:
        raise SpeciesConfigError("species: " + "; ".join(problems))
    return parsed


def _optional_text(value: object) -> str | None:
    """Return the value as text, keeping None as None."""
    return None if value is None else str(value)


def checksum_problem(key: str, entry: SpeciesEntry, actual_sha256: str) -> str | None:
    """Return what is wrong if a genome's checksum differs from its pin, else None.

    Args:
        key: Species key, for the message.
        entry: The species' config entry.
        actual_sha256: SHA-256 of the genome FASTA that is about to be used.

    Returns:
        None if the entry has no `sha256` pin or the pin matches; otherwise a
        message naming both digests.
    """
    if entry.sha256 is None or entry.sha256 == actual_sha256:
        return None
    return (
        f"{key}: the genome is not the one the config pins: its sha256 is "
        f"{actual_sha256}, the config says {entry.sha256}"
    )
