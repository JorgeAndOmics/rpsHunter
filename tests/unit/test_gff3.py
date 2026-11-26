"""Domain instances are written as GFF3 features, one file per species."""

import pandas as pd

import gff3

INSTANCE = {
    "Sequence": "CM000667.2",
    "Strand": "-",
    "Family": "zinc_finger",
    "Model": "pfam00096",
    "Model_Name": "zf-C2H2",
    "Hit_Type": "Specific",
    "Evidence": "standard",
    "Evalue": 1e-9,
    "Bitscore": 40.1,
    "Start": 826,
    "End": 894,
    "N_Calls": 2,
    "Queries": "NP_1.1,NP_2.1",
}


def test_an_empty_table_still_gives_a_valid_gff3_header():
    lines = gff3.lines(pd.DataFrame([INSTANCE]).iloc[0:0])

    assert lines == ["##gff-version 3"]


def test_an_instance_becomes_one_feature_line_with_nine_columns():
    lines = gff3.lines(pd.DataFrame([INSTANCE]))

    fields = lines[1].split("\t")
    assert len(fields) == 9
    assert fields[:8] == [
        "CM000667.2",
        "rpsHunter",
        "protein_match",
        "826",
        "894",
        "1e-09",
        "-",
        ".",
    ]


def test_attributes_name_the_family_model_and_evidence():
    attributes = gff3.lines(pd.DataFrame([INSTANCE]))[1].split("\t")[8]

    pairs = dict(pair.split("=") for pair in attributes.split(";"))
    assert pairs["Name"] == "zinc_finger"
    assert pairs["model"] == "pfam00096"
    assert pairs["model_name"] == "zf-C2H2"
    assert pairs["evidence"] == "standard"
    assert pairs["hit_type"] == "Specific"
    assert pairs["calls"] == "2"


def test_reserved_characters_in_values_are_percent_encoded():
    row = INSTANCE | {"Model_Name": "a;b=c,d", "Queries": "NP_1.1,NP_2.1"}

    attributes = gff3.lines(pd.DataFrame([row]))[1].split("\t")[8]

    assert "model_name=a%3Bb%3Dc%2Cd" in attributes
    assert "queries=NP_1.1,NP_2.1" in attributes  # a list keeps its commas


def test_each_feature_has_a_unique_identifier():
    rows = [INSTANCE, INSTANCE | {"Start": 900, "End": 968}]

    lines = gff3.lines(pd.DataFrame(rows))

    identifiers = [line.split("\t")[8].split(";")[0] for line in lines[1:]]
    assert len(set(identifiers)) == 2
    assert all(identifier.startswith("ID=") for identifier in identifiers)
