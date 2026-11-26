"""What is derived for each domain call from what RPS-BLAST and rpsbproc report.

This module is pure: values in, values out.
"""

from __future__ import annotations

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
