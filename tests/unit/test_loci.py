"""Loci are chained from gated hits; each reports the families it carries."""

import pandas as pd
import pytest

import loci


def _hit(start, end, tag, **overrides):
    row = {
        "Species": "Homo_sapiens",
        "Species_Name": "Homo sapiens",
        "Sequence": "chr1",
        "Start": start,
        "End": end,
        "Tag": tag,
        "Query_Accession": "NP_1.1",
        "Identity": 50.0,
        "Bitscore": 100.0,
    }
    return row | overrides


def _call(tag, family, model_from, model_to, **overrides):
    row = {
        "Tag": tag,
        "Family": family,
        "Aligned_PSSM": 11,
        "Model": "pfam00001",
        "Model_From": model_from,
        "Model_To": model_to,
        "Model_Length": 100,
        "Evalue": 1e-9,
        "Evidence": "standard",
        "Sequence": "chr1",
        "Strand": "+",
        "Start": 1,
        "End": 2,
    }
    return row | overrides


def test_hits_within_the_gap_share_a_locus_and_farther_ones_do_not():
    hits = pd.DataFrame(
        [_hit(1000, 1300, "A"), _hit(21300, 21500, "B"), _hit(41501, 41700, "C")]
    )

    chained = loci.chain(hits, gap=20_000)

    assert chained.loc[chained["Tag"] == "A", "Locus"].item() == "chr1:1000-21500"
    assert chained.loc[chained["Tag"] == "B", "Locus"].item() == "chr1:1000-21500"
    assert chained.loc[chained["Tag"] == "C", "Locus"].item() == "chr1:41501-41700"


def test_the_gap_is_measured_from_the_farthest_end_so_far():
    # B lies inside A; C is 15 kb beyond A's end though 24 kb beyond B's.
    hits = pd.DataFrame(
        [_hit(1000, 10_000, "A"), _hit(1200, 1400, "B"), _hit(25_000, 25_300, "C")]
    )

    assert loci.chain(hits, gap=20_000)["Locus"].nunique() == 1


def test_sequences_and_species_never_share_a_locus():
    hits = pd.DataFrame(
        [
            _hit(1000, 1300, "A"),
            _hit(1000, 1300, "B", Sequence="chr2"),
            _hit(1000, 1300, "C", Species="Mus_musculus"),
        ]
    )

    assert len(loci.chain(hits, gap=20_000).groupby(["Species", "Locus"])) == 3


def test_opposite_strands_may_share_a_locus():
    hits = pd.DataFrame([_hit(1000, 1300, "A"), _hit(1500, 1900, "B")])

    assert loci.chain(hits, gap=20_000)["Locus"].nunique() == 1


def test_no_hits_give_no_loci():
    empty = pd.DataFrame([_hit(1, 2, "A")]).iloc[0:0]

    assert len(loci.chain(empty, gap=20_000)) == 0
    assert len(loci.locus_table(loci.chain(empty, gap=20_000))) == 0


def test_the_locus_table_summarises_its_hits():
    hits = pd.DataFrame(
        [
            _hit(1000, 1300, "A", Identity=40.0, Bitscore=80.0),
            _hit(
                5000, 5600, "B", Identity=95.0, Bitscore=300.0, Query_Accession="NP_2.1"
            ),
        ]
    )

    table = loci.locus_table(loci.chain(hits, gap=20_000))

    row = table.iloc[0]
    assert (row["Start"], row["End"], row["Length"]) == (1000, 5600, 4601)
    assert row["N_Hits"] == 2
    assert row["Queries"] == "NP_1.1,NP_2.1"
    assert row["Best_Identity"] == 95.0
    assert row["Best_Bitscore"] == 300.0


def test_exon_pieces_of_one_model_add_up_to_its_coverage():
    # Three exons of a 100-residue model: 1-40, 38-70, 71-100 -> the whole model.
    calls = pd.DataFrame(
        [
            _call("A", "SET", 1, 40),
            _call("B", "SET", 38, 70),
            _call("C", "SET", 71, 100),
        ]
    )
    chained = loci.chain(
        pd.DataFrame(
            [_hit(1000, 1120, "A"), _hit(3000, 3100, "B"), _hit(6000, 6090, "C")]
        ),
        gap=20_000,
    )

    evidence = loci.family_evidence(chained, calls)

    assert evidence.loc[0, "Coverage"] == pytest.approx(1.0)
    assert evidence.loc[0, "N_Calls"] == 3


def test_a_missing_exon_leaves_the_model_partly_covered():
    calls = pd.DataFrame([_call("A", "SET", 1, 40), _call("B", "SET", 71, 100)])
    chained = loci.chain(
        pd.DataFrame([_hit(1000, 1120, "A"), _hit(3000, 3100, "B")]), gap=20_000
    )

    evidence = loci.family_evidence(chained, calls)

    assert evidence.loc[0, "Coverage"] == pytest.approx(0.70)


def test_coverage_is_that_of_the_best_covered_model_of_the_family():
    calls = pd.DataFrame(
        [
            _call("A", "SET", 1, 20),
            _call(
                "A", "SET", 1, 180, Aligned_PSSM=22, Model="cd00022", Model_Length=200
            ),
        ]
    )
    chained = loci.chain(pd.DataFrame([_hit(1000, 1600, "A")]), gap=20_000)

    evidence = loci.family_evidence(chained, calls)

    assert evidence.loc[0, "Coverage"] == pytest.approx(0.90)
    assert evidence.loc[0, "Best_Model"] == "cd00022"


def test_pieces_from_different_loci_are_not_added_together():
    calls = pd.DataFrame([_call("A", "SET", 1, 50), _call("B", "SET", 51, 100)])
    chained = loci.chain(
        pd.DataFrame([_hit(1000, 1150, "A"), _hit(90_000, 90_150, "B")]), gap=20_000
    )

    evidence = loci.family_evidence(chained, calls)

    assert sorted(evidence["Coverage"]) == pytest.approx([0.5, 0.5])


def test_the_family_carries_its_strongest_evidence_and_best_e_value():
    calls = pd.DataFrame(
        [
            _call("A", "KRAB", 1, 40, Evalue=0.5, Evidence="borderline"),
            _call("A", "KRAB", 1, 40, Evalue=1e-6, Evidence="standard"),
        ]
    )
    chained = loci.chain(pd.DataFrame([_hit(1000, 1120, "A")]), gap=20_000)

    row = loci.family_evidence(chained, calls).iloc[0]

    assert row["Evidence"] == "standard"
    assert row["Best_Evalue"] == 1e-6


def test_copies_count_separate_places_where_a_repeat_occurs():
    # One 23-residue model found at three separate places in the locus.
    calls = pd.DataFrame(
        [
            _call("A", "zinc_finger", 1, 23, Model_Length=23, Start=1000, End=1068),
            _call("A", "zinc_finger", 1, 23, Model_Length=23, Start=1084, End=1152),
            _call("B", "zinc_finger", 1, 23, Model_Length=23, Start=1090, End=1150),
            _call("B", "zinc_finger", 1, 23, Model_Length=23, Start=1168, End=1236),
        ]
    )
    chained = loci.chain(
        pd.DataFrame([_hit(1000, 1160, "A"), _hit(1080, 1240, "B")]), gap=20_000
    )

    row = loci.family_evidence(chained, calls).iloc[0]

    assert row["Copies"] == 3
    assert row["Coverage"] == pytest.approx(1.0)


def test_the_locus_table_lists_its_families_in_the_given_order():
    calls = pd.DataFrame([_call("A", "zinc_finger", 1, 23), _call("A", "KRAB", 1, 40)])
    chained = loci.chain(
        pd.DataFrame([_hit(1000, 1300, "A"), _hit(90_000, 90_100, "Z")]), gap=20_000
    )
    evidence = loci.family_evidence(chained, calls)

    table = loci.locus_table(
        chained, evidence, family_order=["KRAB", "SET", "zinc_finger"]
    )

    assert list(table["Families"]) == ["KRAB,zinc_finger", ""]
    assert list(table["N_Families"]) == [2, 0]


def test_a_call_from_a_hit_outside_every_locus_is_an_error():
    calls = pd.DataFrame([_call("UNKNOWN", "KRAB", 1, 40)])
    chained = loci.chain(pd.DataFrame([_hit(1000, 1300, "A")]), gap=20_000)

    with pytest.raises(ValueError, match="UNKNOWN"):
        loci.family_evidence(chained, calls)


def test_calls_that_already_name_their_species_are_accepted():
    calls = pd.DataFrame([_call("A", "KRAB", 1, 40, Species="Homo_sapiens")])
    chained = loci.chain(pd.DataFrame([_hit(1000, 1120, "A")]), gap=20_000)

    evidence = loci.family_evidence(chained, calls)

    assert list(evidence["Species"]) == ["Homo_sapiens"]
