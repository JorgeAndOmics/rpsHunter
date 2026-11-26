"""Whether RPS-BLAST (CDD) and HMMER (Pfam) report the same domain at the same place.

The two methods search the same hit sequences independently. A domain reported
by one is confirmed by the other when the other reports the same family on the
same hit at an overlapping genomic position (ADR-010). Agreement without
positional overlap would only say that both methods saw the family somewhere.

This module is pure.
"""

from __future__ import annotations

from collections.abc import Collection

import pandas as pd

BOTH = "both"
NOT_COMPARED = "not_compared"
"""The family has models in only one library, or HMMER did not run."""

_KEY = ["Tag", "Family"]


def confirmed(own: pd.DataFrame, other: pd.DataFrame) -> pd.Series:
    """Say, for each row of `own`, whether `other` has an overlapping domain.

    Both tables need Tag, Family, Start and End (genomic, start <= end).
    Overlap is of at least one base. The result is on `own`'s index.
    """
    left = own[[*_KEY, "Start", "End"]].reset_index(names="_row")
    pairs = left.merge(other[[*_KEY, "Start", "End"]], on=_KEY, suffixes=("", "_other"))
    overlap = (pairs["Start"] <= pairs["End_other"]) & (
        pairs["Start_other"] <= pairs["End"]
    )
    return pd.Series(own.index.isin(pairs.loc[overlap, "_row"]), index=own.index)


def agreement(
    own: pd.DataFrame, other: pd.DataFrame, compared: Collection[str], alone: str
) -> pd.Series:
    """Label each row of `own`: BOTH, `alone`, or NOT_COMPARED.

    Args:
        own: The domains to label.
        other: The other method's domains.
        compared: Families both methods searched for.
        alone: Label for a domain of a compared family the other method lacks
            (for example "cdd_only").
    """
    labels = confirmed(own, other).map({True: BOTH, False: alone}).astype("object")
    return labels.where(own["Family"].isin(list(compared)), NOT_COMPARED)
