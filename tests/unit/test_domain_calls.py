"""Each domain call says how strong its evidence is, on NCBI's E-value scale."""

import pytest

import domain_calls


@pytest.mark.parametrize(
    ("evalue", "tier"),
    [
        (0.0, "standard"),
        (1e-30, "standard"),
        (0.01, "standard"),
        (0.011, "borderline"),
        (1.0, "borderline"),
        (1.01, "rescue"),
        (1000.0, "rescue"),
    ],
)
def test_the_tier_follows_ncbi_thresholds(evalue, tier):
    assert domain_calls.evidence_tier(evalue) == tier


def test_a_negative_e_value_is_refused():
    with pytest.raises(ValueError, match="negative"):
        domain_calls.evidence_tier(-1.0)
