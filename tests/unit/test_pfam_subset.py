"""The Pfam subset holds exactly the profiles the configured families name."""

from pathlib import Path

import pytest

import pfam_subset

LIBRARY = "NAME  KRAB\nACC   PF01352.33\n//\nNAME  SET\nACC   PF00856.34\n//\n"


def _families(monkeypatch, raw):
    families = pfam_subset.domain_families.parse_families(raw)
    monkeypatch.setattr(pfam_subset.defaults, "FAMILIES", families)


def test_the_library_index_pairs_names_with_versioned_accessions():
    index = list(pfam_subset.library_index(LIBRARY.splitlines()))

    assert index == [("KRAB", "PF01352.33"), ("SET", "PF00856.34")]


def test_the_configured_profiles_are_fetched_by_name(tmp_path, monkeypatch):
    _families(monkeypatch, {"first": {"pfam": ["PF00856"]}})
    library = tmp_path / "Pfam-A.hmm"
    library.write_text(LIBRARY)
    asked = {}

    def fake_run(command, *, tool):
        asked["program"] = command[0]
        asked["keys"] = Path(command[-1]).read_text()
        return "NAME  SET\n//\n"

    monkeypatch.setattr(pfam_subset.external, "run_tool", fake_run)

    status = pfam_subset.main(
        ["--library", str(library), "--output", str(tmp_path / "subset.hmm")]
    )

    assert status == 0
    assert asked == {"program": "hmmfetch", "keys": "SET\n"}
    assert (tmp_path / "subset.hmm").read_text() == "NAME  SET\n//\n"


def test_an_accession_missing_from_the_library_stops_before_fetching(
    tmp_path, monkeypatch
):
    _families(monkeypatch, {"first": {"pfam": ["PF09514"]}})
    library = tmp_path / "Pfam-A.hmm"
    library.write_text(LIBRARY)

    with pytest.raises(pfam_subset.domain_families.DomainConfigError, match="PF09514"):
        pfam_subset.build(library, tmp_path / "subset.hmm")

    assert not (tmp_path / "subset.hmm").exists()
