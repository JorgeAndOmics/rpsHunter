"""The run's domain tables and GFF3 tracks are written for every species."""

import pandas as pd

import domain_tables


def _calls(species, start):
    return pd.DataFrame(
        [
            {
                "Species": species,
                "Species_Name": species,
                "Query_Accession": "NP_1.1",
                "Tag": f"T{start}",
                "Sequence": "chr1",
                "Strand": "+",
                "Family": "first",
                "Model": "pfam00001",
                "Model_Name": "model_one",
                "Hit_Type": "Specific",
                "Evidence": "standard",
                "Evalue": 1e-9,
                "Bitscore": 40.0,
                "Start": start,
                "End": start + 60,
            }
        ]
    )


def test_tables_and_one_gff3_per_species_are_written(tmp_path, monkeypatch):
    monkeypatch.setattr(domain_tables.defaults, "SPECIES", ["A", "B", "Empty"])
    inputs = []
    for species, start in (("A", 100), ("A", 130), ("B", 500)):
        path = tmp_path / f"{species}_{start}.parquet"
        _calls(species, start).to_parquet(path)
        inputs.append(str(path))
    gff_dir = tmp_path / "ranges"
    gff_dir.mkdir()

    status = domain_tables.main(
        ["--calls-out", str(tmp_path / "calls.parquet"),
         "--domains-out", str(tmp_path / "domains.parquet"),
         "--gff-dir", str(gff_dir), "--calls", *inputs]
    )  # fmt: skip

    assert status == 0
    assert len(pd.read_parquet(tmp_path / "calls.parquet")) == 3
    domains = pd.read_parquet(tmp_path / "domains.parquet")
    assert len(domains) == 2  # the two overlapping calls of A are one instance
    assert len(pd.read_csv(tmp_path / "domains.csv")) == 2
    assert (tmp_path / "calls.csv").is_file()
    assert (gff_dir / "A.gff3").read_text().count("\n") == 2
    assert (gff_dir / "Empty.gff3").read_text() == "##gff-version 3\n"


def test_without_hmmer_tables_no_call_is_compared(tmp_path, monkeypatch):
    monkeypatch.setattr(domain_tables.defaults, "SPECIES", ["A"])
    path = tmp_path / "A.parquet"
    _calls("A", 100).to_parquet(path)

    domain_tables.main(
        ["--calls-out", str(tmp_path / "calls.parquet"),
         "--domains-out", str(tmp_path / "domains.parquet"),
         "--gff-dir", str(tmp_path), "--calls", str(path)]
    )  # fmt: skip

    calls = pd.read_parquet(tmp_path / "calls.parquet")
    assert list(calls["Method_Agreement"]) == ["not_compared"]
    assert not (tmp_path / "hmmer.parquet").exists()


def test_hmmer_domains_are_written_and_both_tables_say_where_the_methods_agree(
    tmp_path, monkeypatch
):
    families = domain_tables.defaults.domain_families.parse_families(
        {"first": {"cdd": ["pfam00001"], "pfam": ["PF00001"]}}
    )
    monkeypatch.setattr(domain_tables.defaults, "FAMILIES", families)
    monkeypatch.setattr(domain_tables.defaults, "SPECIES", ["A"])
    calls_path, hmmer_path = tmp_path / "A.parquet", tmp_path / "A_hmmer.parquet"
    pd.concat([_calls("A", 100), _calls("A", 900)]).to_parquet(calls_path)
    hmmer_domains = pd.DataFrame(
        {
            "Species": ["A", "A"],
            "Tag": ["T100", "T900"],
            "Family": ["first", "first"],
            "Start": [120, 2000],  # the second is elsewhere on its hit
            "End": [150, 2100],
        }
    )
    hmmer_domains.to_parquet(hmmer_path)

    domain_tables.main(
        ["--calls-out", str(tmp_path / "calls.parquet"),
         "--domains-out", str(tmp_path / "domains.parquet"),
         "--gff-dir", str(tmp_path),
         "--hmmer-out", str(tmp_path / "hmmer.parquet"), "--hmmer", str(hmmer_path),
         "--calls", str(calls_path)]
    )  # fmt: skip

    calls = pd.read_parquet(tmp_path / "calls.parquet")
    written = pd.read_parquet(tmp_path / "hmmer.parquet")
    assert list(calls["Method_Agreement"]) == ["both", "cdd_only"]
    assert list(written["Method_Agreement"]) == ["both", "hmmer_only"]
    assert (tmp_path / "hmmer.csv").is_file()
