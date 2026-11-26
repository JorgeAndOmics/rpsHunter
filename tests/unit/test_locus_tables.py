"""The run's locus tables are built from gated hits and domain calls."""

import pandas as pd

import locus_tables


def test_loci_and_family_evidence_are_written_as_parquet_and_csv(tmp_path, monkeypatch):
    families = locus_tables.defaults.domain_families.parse_families(
        {"first": {"cdd": ["pfam00001"]}}
    )
    monkeypatch.setattr(locus_tables.defaults, "FAMILIES", families)
    monkeypatch.setattr(locus_tables.defaults, "LOCUS_GAP", 20_000)
    hits = pd.DataFrame(
        {
            "Species": ["X", "X"],
            "Species_Name": ["Species x", "Species x"],
            "Subject ID": ["chr1", "chr1"],
            "S. Start": [1300, 90_000],  # the first hit is on the minus strand
            "S. End": [1001, 90_300],
            "Tag": ["A", "B"],
            "Query_Accession": ["NP_1.1", "NP_1.1"],
            "Pct Identity": [80.0, 30.0],
            "Bit Score": [200.0, 40.0],
            "Internal_Stop": [False, True],
        }
    )
    calls = pd.DataFrame(
        {
            "Tag": ["A"],
            "Family": ["first"],
            "Model": ["pfam00001"],
            "Aligned_PSSM": [11],
            "Model_From": [1],
            "Model_To": [50],
            "Model_Length": [100],
            "Evalue": [1e-9],
            "Evidence": ["standard"],
            "Sequence": ["chr1"],
            "Strand": ["-"],
            "Start": [1100],
            "End": [1249],
        }
    )
    hits.to_parquet(tmp_path / "hits.parquet")
    calls.to_parquet(tmp_path / "calls.parquet")

    status = locus_tables.main(
        ["--loci-out", str(tmp_path / "loci.parquet"),
         "--families-out", str(tmp_path / "families.parquet"),
         "--hits", str(tmp_path / "hits.parquet"),
         "--calls", str(tmp_path / "calls.parquet")]
    )  # fmt: skip

    loci = pd.read_parquet(tmp_path / "loci.parquet")
    evidence = pd.read_parquet(tmp_path / "families.parquet")
    assert status == 0
    assert list(loci["Locus"]) == ["chr1:1001-1300", "chr1:90000-90300"]
    assert list(loci["Families"]) == ["first", ""]
    assert evidence.loc[0, "Coverage"] == 0.5
    assert (tmp_path / "loci.csv").is_file()
    assert (tmp_path / "families.csv").is_file()


def test_a_species_without_any_hit_does_not_break_the_locus_tables(
    tmp_path, monkeypatch
):
    # blast.py writes a table of untyped columns when tBLASTn finds nothing.
    families = locus_tables.defaults.domain_families.parse_families(
        {"first": {"cdd": ["pfam00001"]}}
    )
    monkeypatch.setattr(locus_tables.defaults, "FAMILIES", families)
    columns = ["Species", "Species_Name", "Subject ID", "S. Start", "S. End", "Tag",
               "Query_Accession", "Pct Identity", "Bit Score", "Internal_Stop"]  # fmt: skip
    pd.DataFrame(columns=columns).to_parquet(tmp_path / "hits.parquet")
    calls = ["Tag", "Family", "Model", "Aligned_PSSM", "Model_From", "Model_To",
             "Model_Length", "Evalue", "Evidence", "Sequence", "Strand", "Start", "End"]  # fmt: skip
    pd.DataFrame(columns=calls).to_parquet(tmp_path / "calls.parquet")

    status = locus_tables.main(
        ["--loci-out", str(tmp_path / "loci.parquet"),
         "--families-out", str(tmp_path / "families.parquet"),
         "--hits", str(tmp_path / "hits.parquet"),
         "--calls", str(tmp_path / "calls.parquet")]
    )  # fmt: skip

    assert status == 0
    assert len(pd.read_parquet(tmp_path / "loci.parquet")) == 0
    assert len(pd.read_parquet(tmp_path / "families.parquet")) == 0


def test_hits_read_from_an_untyped_empty_table_are_typed():
    columns = [
        "Subject ID",
        "S. Start",
        "S. End",
        "Pct Identity",
        "Bit Score",
        "Internal_Stop",
    ]
    empty = pd.DataFrame(columns=columns)

    hits = locus_tables.typed(empty.rename(columns=locus_tables.HIT_COLUMNS))

    assert str(hits["S. Start"].dtype) == "int64"
    assert str(hits["Identity"].dtype) == "float64"
    assert hits["Internal_Stop"].dtype == bool
