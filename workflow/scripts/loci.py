"""Loci: neighbouring hits taken as one gene-sized region, and what each carries.

A tBLASTn hit covers one exon. Hits that passed the gates and lie within a
chosen gap of each other on a sequence are chained into a locus (ADR-008); the
default gap of 20 kb keeps the exons of a gene together without fusing tandem
copies. Hits on either strand may share a locus.

For each locus and domain family the module reports the strength of evidence,
how much of a model is covered when the locus's exon fragments are put together
(ADR-009), and how many separate places the family occurs at. It gives no
verdict on whether a locus is a "complete" gene: which families that takes
depends on the protein family being studied, and the table lets the user decide.

This module is pure.
"""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

import domain_calls

LOCUS_COLUMNS = (
    "Species",
    "Species_Name",
    "Locus",
    "Sequence",
    "Start",
    "End",
    "Length",
    "N_Hits",
    "Queries",
    "Best_Identity",
    "Best_Bitscore",
    "Families",
    "N_Families",
)
EVIDENCE_COLUMNS = (
    "Species",
    "Locus",
    "Family",
    "Evidence",
    "Best_Evalue",
    "Coverage",
    "Best_Model",
    "Copies",
    "N_Calls",
)

_TIER_RANK = {
    domain_calls.STANDARD: 0,
    domain_calls.BORDERLINE: 1,
    domain_calls.RESCUE: 2,
}
_SCOPE = ["Species", "Sequence"]
_LOCUS = ["Species", "Locus"]


def chain(hits: pd.DataFrame, gap: int) -> pd.DataFrame:
    """Assign each hit to a locus by single linkage along its sequence.

    Args:
        hits: Gated hits with Species, Sequence, Start, End (start <= end) and
            any further columns, which are kept.
        gap: Largest distance in bases between a hit and the farthest end of
            the hits before it for the two to share a locus.

    Returns:
        The hits, ordered by position, with a "Locus" column of the form
        `<sequence>:<start>-<end>`, unique within a species.
    """
    ordered = hits.sort_values([*_SCOPE, "Start", "End"]).reset_index(drop=True)
    number = _locus_numbers(ordered, gap)
    start = ordered.groupby(number)["Start"].transform("min").astype(str)
    end = ordered.groupby(number)["End"].transform("max").astype(str)
    ordered["Locus"] = ordered["Sequence"].astype(str) + ":" + start + "-" + end
    return ordered


def _locus_numbers(ordered: pd.DataFrame, gap: int) -> pd.Series:
    """Running locus number of position-ordered hits: a new one after each gap."""
    by_scope = [ordered[column] for column in _SCOPE]
    farthest = ordered.groupby(by_scope, sort=False)["End"].cummax()
    reach = farthest.groupby(by_scope, sort=False).shift()
    return (reach.isna() | (ordered["Start"] - reach > gap)).cumsum()


def _covered_fraction(pieces: pd.DataFrame) -> float:
    """Fraction of one model covered by the union of its aligned intervals."""
    covered, reach = 0, 0
    intervals = sorted(zip(pieces["Model_From"], pieces["Model_To"], strict=True))
    for start, end in intervals:
        covered += max(0, end - max(start - 1, reach))
        reach = max(reach, end)
    return covered / int(pieces["Model_Length"].iloc[0])


def _places(pieces: pd.DataFrame) -> int:
    """Number of separate genomic places (per strand) the calls of one model fall at."""
    ordered = pieces.sort_values(["Sequence", "Strand", "Start"])
    scope = [ordered["Sequence"], ordered["Strand"]]
    farthest = ordered.groupby(scope, sort=False)["End"].cummax()
    reach = farthest.groupby(scope, sort=False).shift()
    return int((reach.isna() | (ordered["Start"] > reach)).sum())


def _one_family(calls: pd.DataFrame) -> dict[str, object]:
    """Evidence for one family in one locus, from its calls there."""
    models = dict(iter(calls.groupby("Aligned_PSSM", sort=False)))
    coverage = {pssm: _covered_fraction(pieces) for pssm, pieces in models.items()}
    best = max(coverage, key=lambda pssm: coverage[pssm])
    tiers = sorted(set(calls["Evidence"]), key=lambda tier: _TIER_RANK[tier])
    return {
        "Evidence": tiers[0],
        "Best_Evalue": float(calls["Evalue"].min()),
        "Coverage": coverage[best],
        "Best_Model": models[best]["Model"].iloc[0],
        "Copies": max(_places(pieces) for pieces in models.values()),
        "N_Calls": len(calls),
    }


def family_evidence(chained: pd.DataFrame, calls: pd.DataFrame) -> pd.DataFrame:
    """Report, per locus and family, what the locus's domain calls add up to.

    Args:
        chained: Output of chain (hits with their Locus and Tag).
        calls: Domain calls with Tag, Family, Model, Aligned_PSSM, Model_From,
            Model_To, Model_Length, Evalue, Evidence, Sequence, Strand, Start, End.

    Returns:
        One row per (locus, family) with EVIDENCE_COLUMNS. Coverage is the
        largest fraction of any one model of the family covered by the locus's
        calls taken together; Copies is the largest number of separate places
        any one model occurs at.

    Raises:
        ValueError: A call's Tag belongs to no hit of the loci, so the calls and
            the hits come from different runs.
    """
    locus_of = chained[["Tag", "Species", "Locus"]].drop_duplicates("Tag")
    # The locus decides the species; a Species column on the calls would clash.
    own = calls.drop(columns=["Species", "Locus"], errors="ignore")
    placed = own.merge(locus_of, on="Tag", how="left", validate="many_to_one")
    orphans = placed.loc[placed["Locus"].isna(), "Tag"]
    if not orphans.empty:
        raise ValueError(
            f"{len(orphans)} calls come from hits in no locus, e.g. Tag {orphans.iloc[0]}"
        )
    rows = [
        {"Species": species, "Locus": locus, "Family": family} | _one_family(group)
        for (species, locus, family), group in placed.groupby(
            [*_LOCUS, "Family"], sort=True
        )
    ]
    return pd.DataFrame(rows, columns=list(EVIDENCE_COLUMNS))


def _listed(values: pd.Series) -> str:
    """Distinct values, sorted and comma-separated."""
    return ",".join(sorted(set(values)))


def _summarise_hits(chained: pd.DataFrame) -> pd.DataFrame:
    """One row per locus from its hits, indexed by (Species, Locus)."""
    table = chained.groupby(_LOCUS, sort=False).agg(
        Species_Name=("Species_Name", "first"),
        Sequence=("Sequence", "first"),
        Start=("Start", "min"),
        End=("End", "max"),
        N_Hits=("Tag", "size"),
        Queries=("Query_Accession", _listed),
        Best_Identity=("Identity", "max"),
        Best_Bitscore=("Bitscore", "max"),
    )
    table["Length"] = table["End"] - table["Start"] + 1
    return table


def _families_per_locus(
    evidence: pd.DataFrame | None, family_order: Sequence[str]
) -> dict[tuple[str, str], list[str]]:
    """Ordered names of the families present in each locus, by (species, locus)."""
    if evidence is None or evidence.empty:
        return {}
    rank = {name: position for position, name in enumerate(family_order)}
    present: dict[tuple[str, str], list[str]] = {}
    for species, locus, family in zip(
        evidence["Species"], evidence["Locus"], evidence["Family"], strict=True
    ):
        present.setdefault((species, locus), []).append(family)
    return {
        key: sorted(names, key=lambda name: (rank.get(name, len(rank)), name))
        for key, names in present.items()
    }


def locus_table(
    chained: pd.DataFrame,
    evidence: pd.DataFrame | None = None,
    family_order: Sequence[str] = (),
) -> pd.DataFrame:
    """Return one row per locus: where it is, what found it, which families it has.

    Args:
        chained: Output of chain; needs Species_Name, Tag, Query_Accession,
            Identity and Bitscore besides the position columns.
        evidence: Output of family_evidence, or None to leave families empty.
        family_order: Family names in the order to list them (the config's).
            Families not named are listed after them, alphabetically.

    Returns:
        A table with LOCUS_COLUMNS, ordered by species, sequence and start.
    """
    if chained.empty:
        return pd.DataFrame(columns=list(LOCUS_COLUMNS))
    table = _summarise_hits(chained).reset_index()
    families = _families_per_locus(evidence, family_order)
    keys = list(zip(table["Species"], table["Locus"], strict=True))
    table["Families"] = [",".join(families.get(key, [])) for key in keys]
    table["N_Families"] = [len(families.get(key, [])) for key in keys]
    table = table.sort_values([*_SCOPE, "Start"]).reset_index(drop=True)
    return table[list(LOCUS_COLUMNS)]
