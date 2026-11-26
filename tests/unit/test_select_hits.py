"""Hits are selected by the quality gate alone, and exported without losing their numbering."""

import pandas as pd

import select_hits

THRESHOLDS = select_hits.QualityThresholds(
    min_identity=20.0, max_evalue=10.0, min_length=15, min_bitscore=20.0
)


def _hits(**columns):
    base = {
        "Pct Identity": [50.0],
        "E-value": [1e-5],
        "Alignment Length": [40],
        "Bit Score": [60.0],
    }
    return pd.DataFrame(base | columns)


def test_a_hit_meeting_every_threshold_passes():
    assert list(select_hits.quality_pass(_hits(), THRESHOLDS)) == [True]


def test_each_threshold_can_fail_a_hit_on_its_own():
    failing = {
        "Pct Identity": 19.9,
        "E-value": 10.1,
        "Alignment Length": 14,
        "Bit Score": 19.9,
    }
    for column, value in failing.items():
        hits = _hits(**{column: [value]})

        assert list(select_hits.quality_pass(hits, THRESHOLDS)) == [False], column


def test_a_value_exactly_at_a_threshold_passes():
    hits = _hits(
        **{
            "Pct Identity": [20.0],
            "E-value": [10.0],
            "Alignment Length": [15],
            "Bit Score": [20.0],
        }
    )

    assert list(select_hits.quality_pass(hits, THRESHOLDS)) == [True]


def test_alignment_gaps_are_dropped_from_the_exported_sequence():
    # A gap in the subject row stands for query residues with no genomic bases.
    assert select_hits.searchable_sequence("MK--VLA") == "MKVLA"


def test_a_stop_codon_is_masked_so_later_residues_keep_their_position():
    # Removing the stop would shift every residue after it, and with them the
    # genomic position of every domain found after it, by three bases.
    assert select_hits.searchable_sequence("MKV*LA") == "MKVXLA"
    assert len(select_hits.searchable_sequence("MK-V*L*A")) == 7


def test_fasta_records_are_named_by_position_and_tag():
    hits = pd.DataFrame(
        {
            "Subject ID": ["chr1", "chr1"],
            "S. Start": [1300, 5000],
            "S. End": [1001, 5098],
            "Tag": ["AAAAAAAAAAAA", "BBBBBBBBBBBB"],
            "Subject Sequence": ["MK-V*LA", "MKV"],
        }
    )

    text = select_hits.fasta_text(hits)

    assert (
        text
        == ">chr1:1300-1001|tag:AAAAAAAAAAAA\nMKVXLA\n>chr1:5000-5098|tag:BBBBBBBBBBBB\nMKV\n"
    )


def test_no_hits_give_an_empty_fasta():
    empty = pd.DataFrame(
        columns=["Subject ID", "S. Start", "S. End", "Tag", "Subject Sequence"]
    )

    assert select_hits.fasta_text(empty) == ""


def test_a_hit_whose_sequence_is_only_gaps_is_left_out_of_the_fasta():
    hits = pd.DataFrame(
        {
            "Subject ID": ["chr1"],
            "S. Start": [1],
            "S. End": [9],
            "Tag": ["A"],
            "Subject Sequence": ["---"],
        }
    )

    assert select_hits.fasta_text(hits) == ""


def test_main_writes_the_audit_the_selected_hits_and_their_sequences(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(select_hits, "configured_thresholds", lambda: THRESHOLDS)
    hits = pd.DataFrame(
        {
            "Subject ID": ["chr1", "chr1"],
            "S. Start": [1, 500],
            "S. End": [120, 620],
            "Tag": ["AAAAAAAAAAAA", "BBBBBBBBBBBB"],
            "Subject Sequence": ["MKV*LA", "MKV"],
            "Pct Identity": [50.0, 10.0],  # the second fails on identity
            "E-value": [1e-5, 1e-5],
            "Alignment Length": [40, 40],
            "Bit Score": [60.0, 60.0],
        }
    )
    hits.to_parquet(tmp_path / "hits.parquet")

    status = select_hits.main(
        ["--input", str(tmp_path / "hits.parquet"),
         "--output-fasta", str(tmp_path / "X.fa"),
         "--output-audit", str(tmp_path / "audit.parquet"),
         "--output-selected", str(tmp_path / "selected.parquet")]
    )  # fmt: skip

    audit = pd.read_parquet(tmp_path / "audit.parquet")
    selected = pd.read_parquet(tmp_path / "selected.parquet")
    assert status == 0
    assert list(audit["Selected"]) == [True, False]
    assert list(selected["Tag"]) == ["AAAAAAAAAAAA"]
    assert "Selected" not in selected.columns
    assert (tmp_path / "X.fa").read_text() == ">chr1:1-120|tag:AAAAAAAAAAAA\nMKVXLA\n"
    assert (tmp_path / "selected.csv").is_file()


def test_the_configured_gate_comes_from_the_blast_section(monkeypatch):
    monkeypatch.setattr(select_hits.defaults, "PERC_IDENTITY_THRESHOLD", 33.0)

    assert select_hits.configured_thresholds().min_identity == 33.0
