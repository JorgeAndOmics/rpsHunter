"""The `domains:` block of the configuration: which models evidence which family.

A family is a name chosen by the user and the accessions of the models that
count as evidence for it, in the two libraries the pipeline searches:

    domains:
      zinc_finger:
        cdd: [pfam00096, smart00355, COG5048]   # CDD accessions, for RPS-BLAST
        pfam: [PF00096]                         # Pfam accessions, for HMMER
        repeat: true                            # optional: occurs in tandem

Accessions are used, not names, because the libraries rename models between
releases while accessions stay. Every table then carries the family beside the
model, and no family name needs to appear in code.

This module is pure: it validates and maps, and reads nothing.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from typing import NamedTuple

PFAM_ACCESSION = re.compile(r"PF\d{5}")
"""A Pfam accession without its version (the version changes with each release)."""

_KEYS = frozenset({"cdd", "pfam", "repeat"})


class DomainConfigError(ValueError):
    """The `domains:` block is not usable; the message lists every problem."""


class Family(NamedTuple):
    """One domain family and the models that evidence it."""

    name: str
    cdd: tuple[str, ...]
    pfam: tuple[str, ...]
    repeat: bool


def _family_problems(name: str, raw: object) -> list[str]:
    """Return what is wrong with one family entry (empty if it is usable)."""
    if not isinstance(raw, Mapping):
        return [f"{name}: must be a mapping with 'cdd' and/or 'pfam' accession lists"]
    problems = [f"{name}: unknown key '{key}'" for key in raw if key not in _KEYS]
    if not raw.get("cdd") and not raw.get("pfam"):
        problems.append(f"{name}: lists no model (give 'cdd' and/or 'pfam' accessions)")
    problems += [
        f"{name}: '{accession}' is not a Pfam accession (PF plus five digits, no version)"
        for accession in raw.get("pfam") or ()
        if not PFAM_ACCESSION.fullmatch(str(accession))
    ]
    return problems


def _shared_models(families: Mapping[str, Family]) -> list[str]:
    """Return a problem for each model that more than one family lists."""
    owners: dict[str, list[str]] = {}
    for family in families.values():
        for accession in (*family.cdd, *family.pfam):
            owners.setdefault(accession, []).append(family.name)
    return [
        f"{accession} is listed by more than one family: {', '.join(names)}"
        for accession, names in owners.items()
        if len(names) > 1
    ]


def parse_families(raw: object) -> dict[str, Family]:
    """Validate the `domains:` block and return its families in config order.

    Raises:
        DomainConfigError: The block is missing or empty, a family lists no
            model or has an unknown key, a Pfam accession is malformed, or a
            model belongs to two families. Every problem is listed.
    """
    if not isinstance(raw, Mapping) or not raw:
        raise DomainConfigError("domains: at least one family must be defined")
    problems = [p for name, entry in raw.items() for p in _family_problems(name, entry)]
    families = {
        name: Family(
            name=name,
            cdd=tuple(str(a) for a in entry.get("cdd") or ()),
            pfam=tuple(str(a) for a in entry.get("pfam") or ()),
            repeat=bool(entry.get("repeat", False)),
        )
        for name, entry in raw.items()
        if isinstance(entry, Mapping)
    }
    problems += _shared_models(families)
    if problems:
        raise DomainConfigError("domains: " + "; ".join(problems))
    return families


def cdd_accessions(families: Mapping[str, Family]) -> list[str]:
    """Every CDD accession to search, once, in config order."""
    return [accession for family in families.values() for accession in family.cdd]


def family_by_cdd(families: Mapping[str, Family]) -> dict[str, str]:
    """Map each CDD accession to the family it evidences."""
    return {a: family.name for family in families.values() for a in family.cdd}


def family_by_pfam(families: Mapping[str, Family]) -> dict[str, str]:
    """Map each Pfam accession to the family it evidences."""
    return {a: family.name for family in families.values() for a in family.pfam}


def family_by_superfamily(
    families: Mapping[str, Family], superfamily_of: Mapping[str, str]
) -> dict[str, str]:
    """Map CDD superfamily accessions to families, through their member models.

    RPS-BLAST post-processing reports some hits at superfamily level (a cluster
    accession such as cl02581). Such a hit counts for the family whose models
    belong to that cluster.

    Args:
        families: The configured families.
        superfamily_of: CDD model accession to superfamily accession, for at
            least the configured models (from CDD's family_superfamily_links).

    Raises:
        DomainConfigError: Models of two families share a superfamily, so a
            superfamily-level hit could not be assigned to either.
    """
    owners: dict[str, set[str]] = {}
    for accession, family in family_by_cdd(families).items():
        superfamily = superfamily_of.get(accession)
        if superfamily is not None:
            owners.setdefault(superfamily, set()).add(family)
    shared = {s: names for s, names in owners.items() if len(names) > 1}
    if shared:
        detail = "; ".join(f"{s}: {', '.join(sorted(n))}" for s, n in shared.items())
        raise DomainConfigError(f"domains: families share a CDD superfamily ({detail})")
    return {superfamily: next(iter(names)) for superfamily, names in owners.items()}


def pfam_profile_names(
    families: Mapping[str, Family], library: Iterable[tuple[str, str]]
) -> dict[str, str]:
    """Map each configured Pfam accession to its profile name in an HMM library.

    Args:
        families: The configured families.
        library: (profile name, versioned accession) for every profile in the
            library, e.g. ("Pkinase", "PF00069.30").

    Raises:
        DomainConfigError: A configured accession is not in the library.
    """
    wanted = family_by_pfam(families)
    names = {
        accession.split(".")[0]: name
        for name, accession in library
        if accession.split(".")[0] in wanted
    }
    missing = [accession for accession in wanted if accession not in names]
    if missing:
        raise DomainConfigError(
            f"domains: Pfam accessions not found in the HMM library: {', '.join(missing)}"
        )
    return {accession: names[accession] for accession in wanted}
