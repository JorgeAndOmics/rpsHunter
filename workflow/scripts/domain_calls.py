"""What is derived for each domain call from what RPS-BLAST and rpsbproc report.

This module is pure: values in, values out.
"""

from __future__ import annotations

import re
from collections.abc import Mapping

import pandas as pd

STANDARD_MAX_EVALUE = 0.01
"""NCBI's CD-Search default: calls at or below it need no further argument."""

BORDERLINE_MAX_EVALUE = 1.0
"""NCBI: "results with E-values in the range of 1 and above should be considered
putative false positives". Between the two thresholds a call is borderline."""

STANDARD = "standard"
BORDERLINE = "borderline"
RESCUE = "rescue"


def evidence_tier(evalue: float) -> str:
    """Return how strong a call's evidence is, from its E-value on NCBI's scale.

    Returns:
        "standard" at or below 0.01, "borderline" up to 1, and "rescue" above:
        a call that only exists because the cut-off was raised on purpose.

    Raises:
        ValueError: The E-value is negative.
    """
    if evalue < 0:
        raise ValueError(f"negative E-value: {evalue}")
    if evalue <= STANDARD_MAX_EVALUE:
        return STANDARD
    if evalue <= BORDERLINE_MAX_EVALUE:
        return BORDERLINE
    return RESCUE


# --- reading rpsbproc's report --------------------------------------------------

REPORT_COLUMNS = (
    "Tag",
    "Sequence",
    "Hit_A",
    "Hit_B",
    "Hit_Type",
    "PSSM_ID",
    "Residue_From",
    "Residue_To",
    "Evalue",
    "Bitscore",
    "Model",
    "Model_Name",
    "Incomplete",
)
"""One row per call. Hit_A and Hit_B are the hit's genomic ends as tBLASTn wrote
them (A greater than B on the minus strand); Residue_From/To are positions in the
hit's translated sequence; Model is the CDD accession of the model or superfamily."""

_DEFINITION = re.compile(r"(?P<sequence>.+):(?P<a>\d+)-(?P<b>\d+)\|tag:(?P<tag>\w+)")
_INTEGERS = ("Hit_A", "Hit_B", "PSSM_ID", "Residue_From", "Residue_To")
_FLOATS = ("Evalue", "Bitscore")


def _hit_of(definition: str) -> tuple[str, str, int, int]:
    """Return (tag, sequence, a, b) from a query definition line."""
    match = _DEFINITION.fullmatch(definition)
    if match is None:
        raise ValueError(f"not a tagged hit identifier: {definition!r}")
    return match["tag"], match["sequence"], int(match["a"]), int(match["b"])


def parse_report(text: str) -> pd.DataFrame:
    """Parse an rpsbproc report into one row per domain call.

    Args:
        text: The report (tab-delimited, `-t both`); may be empty.

    Returns:
        A table with REPORT_COLUMNS, integers and floats typed. Sequences for
        which rpsbproc reports no domain contribute no row.

    Raises:
        ValueError: A query's definition line is not `<seq>:<a>-<b>|tag:<Tag>`.
    """
    rows: list[tuple[object, ...]] = []
    hit: tuple[str, str, int, int] | None = None
    in_domains = False
    for line in text.splitlines():
        fields = line.split("\t")
        if fields[0] == "QUERY":
            hit = _hit_of(fields[4])
        elif fields[0] in ("DOMAINS", "ENDDOMAINS"):
            in_domains = fields[0] == "DOMAINS"
        elif in_domains and hit is not None:
            # session, query, hit type, PSSM, from, to, E-value, bitscore, accession,
            # short name, incomplete, superfamily PSSM
            rows.append((*hit, fields[2], *fields[3:11]))
    table = pd.DataFrame(rows, columns=list(REPORT_COLUMNS))
    return table.astype(
        dict.fromkeys(_INTEGERS, "int64") | dict.fromkeys(_FLOATS, "float64")
    )


# --- placing a call on the genome (ADR-007) --------------------------------------


def genomic_span(
    hit_a: int, hit_b: int, residue_from: int, residue_to: int
) -> tuple[int, int, str]:
    """Return (start, end, strand) of residues within a hit, on the genome.

    A residue is three bases. On the plus strand positions count up from the
    hit's first base; on the minus strand (hit_a greater than hit_b) they count
    down from it. Coordinates are 1-based and inclusive, with start <= end.
    """
    if hit_a <= hit_b:
        return hit_a + 3 * (residue_from - 1), hit_a + 3 * residue_to - 1, "+"
    return hit_a - 3 * residue_to + 1, hit_a - 3 * (residue_from - 1), "-"


def place_on_genome(calls: pd.DataFrame) -> pd.DataFrame:
    """Add each call's strand, its hit's span and its own span on the genome.

    Adds Strand, Hit_Start, Hit_End (the hit, start <= end) and Start, End (the
    domain itself), and drops Hit_A and Hit_B.
    """
    forward = calls["Hit_A"] <= calls["Hit_B"]
    first, last = calls["Residue_From"], calls["Residue_To"]
    placed = calls.assign(
        Strand=forward.map({True: "+", False: "-"}).astype("object"),
        Hit_Start=calls[["Hit_A", "Hit_B"]].min(axis=1),
        Hit_End=calls[["Hit_A", "Hit_B"]].max(axis=1),
        Start=(calls["Hit_A"] + 3 * (first - 1)).where(
            forward, calls["Hit_A"] - 3 * last + 1
        ),
        End=(calls["Hit_A"] + 3 * last - 1).where(
            forward, calls["Hit_A"] - 3 * (first - 1)
        ),
    )
    return placed.drop(columns=["Hit_A", "Hit_B"])


# --- which part of the model each call covers (for ADR-009) ---------------------


def add_model_coordinates(
    calls: pd.DataFrame, alignments: pd.DataFrame, model_lengths: Mapping[int, int]
) -> pd.DataFrame:
    """Add the model positions each call aligns to, and the model's length.

    rpsbproc reports where a call lies in the hit; RPS-BLAST's table also says
    which part of the model it aligns to. The two are matched on the hit's Tag
    and the residue range. A superfamily call has no model of its own, so it
    takes the best-scoring alignment with that residue range.

    Args:
        calls: Output of parse_report.
        alignments: RPS-BLAST table (rpsblast.py) with Tag, "Subject ID"
            (`CDD:<pssm>`), "Q. Start", "Q. End", "S. Start", "S. End", "Bit Score".
        model_lengths: PSSM id to model length, from CDD's cddid.tbl.

    Returns:
        The calls with Aligned_PSSM (the PSSM id of the model the coordinates
        refer to), Model_From, Model_To and Model_Length.

    Raises:
        ValueError: A call has no alignment with its Tag and residue range.
    """
    table = alignments.assign(
        Aligned_PSSM=alignments["Subject ID"]
        .str.rsplit(":", n=1)
        .str[-1]
        .astype("int64")
    ).sort_values("Bit Score", ascending=False)
    exact = table.drop_duplicates(["Tag", "Aligned_PSSM", "Q. Start", "Q. End"])
    best = table.drop_duplicates(["Tag", "Q. Start", "Q. End"])
    columns = ["Tag", "Q. Start", "Q. End", "Aligned_PSSM", "S. Start", "S. End"]
    own = calls.merge(
        exact[columns],
        how="left",
        left_on=["Tag", "Residue_From", "Residue_To", "PSSM_ID"],
        right_on=["Tag", "Q. Start", "Q. End", "Aligned_PSSM"],
    )
    member = calls.merge(
        best[columns],
        how="left",
        left_on=["Tag", "Residue_From", "Residue_To"],
        right_on=columns[:3],
    )
    chosen = own.where(own["Aligned_PSSM"].notna(), member)
    unmatched = chosen[chosen["Aligned_PSSM"].isna()]
    if not unmatched.empty:
        raise ValueError(
            f"{len(unmatched)} calls have no RPS-BLAST alignment, e.g. Tag "
            f"{unmatched.iloc[0]['Tag']} residues {unmatched.iloc[0]['Residue_From']}-"
            f"{unmatched.iloc[0]['Residue_To']}"
        )
    return calls.assign(
        Aligned_PSSM=chosen["Aligned_PSSM"].astype("int64").to_numpy(),
        Model_From=chosen["S. Start"].astype("int64").to_numpy(),
        Model_To=chosen["S. End"].astype("int64").to_numpy(),
        Model_Length=chosen["Aligned_PSSM"]
        .astype("int64")
        .map(model_lengths)
        .to_numpy(),
    )
