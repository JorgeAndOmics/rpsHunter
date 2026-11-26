"""Domain instances: overlapping calls of one model, counted once.

The same stretch of genome is usually called several times: by each query
protein, and by overlapping BLAST hits of one query. An instance is a maximal
run of overlapping calls of the same model on the same strand of the same
sequence; it spans all of them and carries the evidence of the best one.

Calls of different models are not merged, even within a family: a model that
covers a whole repeat array and a model for one repeat describe different
things, and merging them would hide the repeats.

This module is pure.
"""

from __future__ import annotations

import pandas as pd

KEY = ["Species", "Sequence", "Strand", "Family", "Model"]

INSTANCE_COLUMNS = (
    "Species",
    "Species_Name",
    "Sequence",
    "Strand",
    "Start",
    "End",
    "Family",
    "Model",
    "Model_Name",
    "Hit_Type",
    "Evidence",
    "Evalue",
    "Bitscore",
    "N_Calls",
    "Queries",
    "Tags",
)

_FROM_BEST_CALL = (
    "Species_Name",
    "Model_Name",
    "Hit_Type",
    "Evidence",
    "Evalue",
    "Bitscore",
)


def _listed(values: pd.Series) -> str:
    """Distinct values, sorted and comma-separated."""
    return ",".join(sorted(set(values)))


def instances(calls: pd.DataFrame) -> pd.DataFrame:
    """Merge overlapping calls of one model into instances.

    Args:
        calls: Domain calls with the KEY columns, Start, End, Evalue, Bitscore,
            Species_Name, Model_Name, Hit_Type, Evidence, Query_Accession, Tag.

    Returns:
        One row per instance with INSTANCE_COLUMNS, ordered by species,
        sequence and start. Start and End span every merged call; the columns
        describing evidence are those of the call with the lowest E-value.
    """
    ordered = calls.sort_values([*KEY, "Start", "End"]).reset_index(drop=True)
    reach = (
        ordered.groupby(KEY, sort=False)["End"]
        .cummax()
        .groupby([ordered[k] for k in KEY], sort=False)
        .shift()
    )
    ordered["Instance"] = (reach.isna() | (ordered["Start"] > reach)).cumsum()
    grouped = ordered.groupby("Instance", sort=False)
    best = ordered.loc[
        grouped["Evalue"].idxmin(), ["Instance", *_FROM_BEST_CALL]
    ].set_index("Instance")
    merged = grouped.agg(
        Species=("Species", "first"),
        Sequence=("Sequence", "first"),
        Strand=("Strand", "first"),
        Start=("Start", "min"),
        End=("End", "max"),
        Family=("Family", "first"),
        Model=("Model", "first"),
        N_Calls=("Start", "size"),
        Queries=("Query_Accession", _listed),
        Tags=("Tag", _listed),
    ).join(best)
    merged = merged.sort_values(["Species", "Sequence", "Start", "End"]).reset_index(
        drop=True
    )
    return merged[list(INSTANCE_COLUMNS)]
