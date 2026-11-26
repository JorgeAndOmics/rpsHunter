"""The search database is built from exactly the configured models, or not at all."""

import pytest

import cdd_subset


def test_profile_files_are_returned_in_the_order_of_the_accessions(tmp_path):
    for accession in ("pfam00096", "cd19193"):
        (tmp_path / f"{accession}.smp").write_text("profile")

    files = cdd_subset.profile_files(["cd19193", "pfam00096"], tmp_path)

    assert [path.name for path in files] == ["cd19193.smp", "pfam00096.smp"]


def test_every_model_without_a_profile_file_is_named(tmp_path):
    (tmp_path / "pfam00096.smp").write_text("profile")

    with pytest.raises(FileNotFoundError) as raised:
        cdd_subset.profile_files(["pfam00096", "COG9999", "cd99999"], tmp_path)

    assert "COG9999" in str(raised.value)
    assert "cd99999" in str(raised.value)
    assert "pfam00096" not in str(raised.value).split("for:")[1]


def test_main_stops_before_building_when_a_model_is_missing(tmp_path, monkeypatch):
    families = cdd_subset.domain_families.parse_families({"F": {"cdd": ["cd99999"]}})
    monkeypatch.setattr(cdd_subset.defaults, "FAMILIES", families)

    with pytest.raises(FileNotFoundError, match="cd99999"):
        cdd_subset.main(
            ["--smp-dir", str(tmp_path), "--output-db", str(tmp_path / "db")]
        )

    assert not (tmp_path / "smp_list.txt").exists()


def test_main_builds_the_database_from_the_listed_profiles(tmp_path, monkeypatch):
    (tmp_path / "cd00001.smp").write_text("profile")
    families = cdd_subset.domain_families.parse_families({"F": {"cdd": ["cd00001"]}})
    monkeypatch.setattr(cdd_subset.defaults, "FAMILIES", families)
    seen = {}
    monkeypatch.setattr(
        cdd_subset.external, "run_tool", lambda command, *, tool: seen.update(tool=tool)
    )

    status = cdd_subset.main(
        ["--smp-dir", str(tmp_path), "--output-db", str(tmp_path / "db")]
    )

    assert status == 0
    assert seen["tool"] == "makeprofiledb"
    assert (tmp_path / "smp_list.txt").read_text() == f"{tmp_path / 'cd00001.smp'}\n"
