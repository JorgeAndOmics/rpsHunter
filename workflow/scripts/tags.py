"""Tags: short, stable identifiers for BLAST hits.

Every tBLASTn hit gets a Tag that travels with it through the pipeline: it is
written into the FASTA identifier of the hit's sequence, read back from the
RPS-BLAST and rpsbproc output, and used to join domains to the hit they came
from.

A Tag is derived from what the hit is (query, species, sequence, subject and
query coordinates), so the same hit has the same Tag on every run and on every
machine, and hits can be compared between runs. Query coordinates are part of
the identity because one stretch of a genome routinely aligns to several places
on the query (a repeated domain such as a zinc finger): those are separate hits
with identical subject coordinates.

This module is pure: no file, clock, random or configuration access.
"""

from __future__ import annotations

import base64
import hashlib
import re
from typing import NamedTuple

import pandas as pd

TAG_LENGTH = 12
"""Characters per Tag. Base32, so 60 bits: about one chance in 10^6 of any
collision among a million hits, and tag_hits refuses a table that has one."""

TAG_IN_IDENTIFIER = re.compile(r"\|tag:(\w+)")
"""Finds the Tag in a sequence identifier of the form `<seq>:<start>-<end>|tag:<Tag>`."""

_SEPARATOR = "\x1f"  # ASCII unit separator: cannot occur in any identity field


class HitIdentity(NamedTuple):
    """What makes a BLAST hit that hit, and no other.

    Subject positions are kept as BLAST reports them: the start is greater than
    the end for a minus-strand hit, so the strand is part of the identity.
    """

    query_accession: str
    species: str
    subject_id: str
    subject_start: int
    subject_end: int
    query_start: int
    query_end: int


IDENTITY_COLUMNS = (
    "Query_Accession",
    "Species",
    "Subject ID",
    "S. Start",
    "S. End",
    "Q. Start",
    "Q. End",
)
"""Columns of a BLAST hit table, in the order of the fields of HitIdentity."""


class DuplicateHitError(ValueError):
    """Two rows of a hit table describe the same hit."""


class TagCollisionError(ValueError):
    """Two different hits were given the same Tag."""


def hit_tag(hit: HitIdentity) -> str:
    """Return the Tag of one hit: TAG_LENGTH characters from A-Z and 2-7."""
    fields = (
        hit.query_accession,
        hit.species,
        hit.subject_id,
        str(int(hit.subject_start)),
        str(int(hit.subject_end)),
        str(int(hit.query_start)),
        str(int(hit.query_end)),
    )
    digest = hashlib.sha256(_SEPARATOR.join(fields).encode()).digest()
    return base64.b32encode(digest).decode("ascii")[:TAG_LENGTH]


def tag_hits(hits: pd.DataFrame) -> pd.Series:
    """Return the Tag of every row of a hit table, on the table's index.

    Args:
        hits: A table with the columns named in IDENTITY_COLUMNS.

    Returns:
        A string Series named "Tag", aligned with `hits`.

    Raises:
        DuplicateHitError: Two rows have the same identity. Remove duplicate
            rows before tagging.
        TagCollisionError: Two different hits hash to the same Tag. The Tag
            would no longer identify a hit, so nothing is returned.
    """
    identity = hits[list(IDENTITY_COLUMNS)]
    repeated = identity[identity.duplicated(keep=False)]
    if not repeated.empty:
        example = repeated.iloc[0].to_dict()
        raise DuplicateHitError(
            f"{len(repeated)} rows describe a hit more than once, e.g. {example}"
        )
    rows = identity.itertuples(index=False, name=None)
    values = [hit_tag(HitIdentity(*row)) for row in rows]
    result = pd.Series(values, index=hits.index, name="Tag", dtype="object")
    shared = result[result.duplicated(keep=False)]
    if not shared.empty:
        raise TagCollisionError(
            f"different hits share the Tag {shared.iloc[0]}: lengthen TAG_LENGTH"
        )
    return result
