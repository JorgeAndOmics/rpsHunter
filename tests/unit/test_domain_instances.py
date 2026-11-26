"""Overlapping calls of one model are one domain instance."""

import pandas as pd

import domain_instances


def _call(start, end, evalue, **overrides):
    row = {
        "Species": "Homo_sapiens",
        "Species_Name": "Homo sapiens",
        "Query_Accession": "NP_1.1",
        "Tag": f"T{start}",
        "Sequence": "chr1",
        "Strand": "+",
        "Family": "zinc_finger",
        "Model": "pfam00096",
        "Model_Name": "zf-C2H2",
        "Hit_Type": "Specific",
        "Evidence": "standard",
        "Evalue": evalue,
        "Bitscore": 30.0,
        "Start": start,
        "End": end,
    }
    return row | overrides


def test_overlapping_calls_become_one_instance_spanning_both():
    calls = pd.DataFrame([_call(100, 168, 1e-5), _call(150, 200, 1e-9)])

    result = domain_instances.instances(calls)

    assert len(result) == 1
    assert (result.loc[0, "Start"], result.loc[0, "End"]) == (100, 200)
    assert result.loc[0, "N_Calls"] == 2


def test_the_instance_carries_the_evidence_of_its_best_call():
    weak = _call(100, 168, 0.5, Evidence="borderline", Hit_Type="Non-specific")
    strong = _call(150, 200, 1e-9, Bitscore=45.0)

    result = domain_instances.instances(pd.DataFrame([weak, strong]))

    row = result.loc[0]
    assert row["Evalue"] == 1e-9
    assert row["Evidence"] == "standard"
    assert row["Hit_Type"] == "Specific"
    assert row["Bitscore"] == 45.0


def test_calls_that_do_not_overlap_stay_separate():
    calls = pd.DataFrame([_call(100, 168, 1e-5), _call(169, 237, 1e-5)])

    assert len(domain_instances.instances(calls)) == 2


def test_a_long_call_bridges_the_short_ones_it_covers():
    calls = pd.DataFrame(
        [_call(100, 400, 1e-5), _call(120, 150, 1e-5), _call(300, 380, 1e-5)]
    )

    result = domain_instances.instances(calls)

    assert len(result) == 1
    assert result.loc[0, "N_Calls"] == 3


def test_strand_sequence_species_and_model_each_keep_calls_apart():
    base = _call(100, 168, 1e-5)
    others = [
        _call(100, 168, 1e-5, Strand="-"),
        _call(100, 168, 1e-5, Sequence="chr2"),
        _call(100, 168, 1e-5, Species="Mus_musculus"),
        _call(100, 168, 1e-5, Model="smart00355", Model_Name="ZnF_C2H2"),
    ]

    assert len(domain_instances.instances(pd.DataFrame([base, *others]))) == 5


def test_the_queries_and_tags_behind_an_instance_are_listed_once_each():
    calls = pd.DataFrame(
        [
            _call(100, 168, 1e-5, Query_Accession="NP_2.1", Tag="B"),
            _call(110, 170, 1e-6, Query_Accession="NP_1.1", Tag="A"),
            _call(120, 175, 1e-7, Query_Accession="NP_2.1", Tag="C"),
        ]
    )

    row = domain_instances.instances(calls).loc[0]

    assert row["Queries"] == "NP_1.1,NP_2.1"
    assert row["Tags"] == "A,B,C"


def test_no_calls_give_an_empty_table_with_every_column():
    empty = pd.DataFrame([_call(1, 2, 1.0)]).iloc[0:0]

    result = domain_instances.instances(empty)

    assert len(result) == 0
    assert list(result.columns) == list(domain_instances.INSTANCE_COLUMNS)


def test_instances_are_ordered_by_position_within_species_and_sequence():
    calls = pd.DataFrame([_call(900, 950, 1e-5), _call(100, 168, 1e-5)])

    result = domain_instances.instances(calls)

    assert list(result["Start"]) == [100, 900]
