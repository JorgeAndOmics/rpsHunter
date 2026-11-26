"""Two methods agree on a domain only where their positions overlap."""

import pandas as pd

import method_agreement as ma


def _domains(*rows):
    return pd.DataFrame(rows, columns=["Tag", "Family", "Start", "End"])


def test_same_family_on_the_same_hit_at_an_overlapping_place_agrees():
    calls = _domains(("A", "SET", 100, 400))
    other = _domains(("A", "SET", 400, 600))  # shares one base

    assert list(ma.agreement(calls, other, {"SET"}, "cdd_only")) == ["both"]


def test_same_family_elsewhere_on_the_hit_does_not_agree():
    calls = _domains(("A", "SET", 100, 400))
    other = _domains(("A", "SET", 401, 600))

    assert list(ma.agreement(calls, other, {"SET"}, "cdd_only")) == ["cdd_only"]


def test_another_family_or_another_hit_at_the_same_place_does_not_agree():
    calls = _domains(("A", "SET", 100, 400), ("B", "SET", 100, 400))
    other = _domains(("A", "KRAB", 100, 400), ("C", "SET", 100, 400))

    assert list(ma.agreement(calls, other, {"SET", "KRAB"}, "cdd_only")) == [
        "cdd_only",
        "cdd_only",
    ]


def test_a_family_only_one_method_searches_is_not_compared():
    calls = _domains(("A", "SSXRD", 100, 400), ("A", "SET", 500, 700))
    other = _domains(("A", "SET", 500, 700))

    assert list(ma.agreement(calls, other, {"SET"}, "cdd_only")) == [
        "not_compared",
        "both",
    ]


def test_without_the_other_method_nothing_is_compared():
    calls = _domains(("A", "SET", 100, 400))

    assert list(ma.agreement(calls, _domains(), set(), "cdd_only")) == ["not_compared"]


def test_several_matches_still_give_one_label_per_row_on_the_own_index():
    calls = _domains(("A", "SET", 100, 400), ("A", "SET", 900, 950))
    calls.index = [7, 9]
    other = _domains(("A", "SET", 90, 200), ("A", "SET", 300, 500))

    labels = ma.agreement(calls, other, {"SET"}, "cdd_only")

    assert list(labels.index) == [7, 9]
    assert list(labels) == ["both", "cdd_only"]


def test_no_rows_give_no_labels():
    assert ma.agreement(_domains(), _domains(), {"SET"}, "hmmer_only").empty


def test_rows_sharing_an_index_label_are_judged_one_by_one():
    calls = _domains(
        ("A", "SET", 100, 400), ("A", "SET", 900, 950), ("A", "SET", 2000, 2100)
    )
    calls.index = [7, 7, 9]
    other = _domains(("A", "SET", 90, 200))

    assert list(ma.confirmed(calls, other)) == [True, False, False]
