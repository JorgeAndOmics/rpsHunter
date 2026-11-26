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
         "--gff-dir", str(gff_dir), *inputs]
    )  # fmt: skip

    assert status == 0
    assert len(pd.read_parquet(tmp_path / "calls.parquet")) == 3
    domains = pd.read_parquet(tmp_path / "domains.parquet")
    assert len(domains) == 2  # the two overlapping calls of A are one instance
    assert len(pd.read_csv(tmp_path / "domains.csv")) == 2
    assert (tmp_path / "calls.csv").is_file()
    assert (gff_dir / "A.gff3").read_text().count("\n") == 2
    assert (gff_dir / "Empty.gff3").read_text() == "##gff-version 3\n"
