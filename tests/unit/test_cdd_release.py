"""The CDD files on disk must be the release the config says."""

import pytest

import cdd_release

INFO = "cdd version 3.22\n71136 models in the default search set\n"


def test_the_release_is_read_from_the_first_line_of_cdd_info():
    assert cdd_release.declared_version(INFO) == "3.22"


def test_a_file_that_does_not_say_its_version_is_refused():
    with pytest.raises(cdd_release.log.PipelineError, match=r"cdd\.info"):
        cdd_release.declared_version("<html>not found</html>")


def test_the_configured_release_passes(tmp_path, monkeypatch):
    info = tmp_path / "cdd.info"
    info.write_text(INFO)
    monkeypatch.setattr(cdd_release.defaults, "CDD_VERSION", "3.22")

    assert cdd_release.main(["--info", str(info)]) == 0


def test_another_release_is_an_error_that_says_what_to_change(tmp_path, monkeypatch):
    info = tmp_path / "cdd.info"
    info.write_text(INFO.replace("3.22", "3.23"))
    monkeypatch.setattr(cdd_release.defaults, "CDD_VERSION", "3.22")

    with pytest.raises(cdd_release.log.PipelineError) as raised:
        cdd_release.main(["--info", str(info)])

    assert "3.23" in str(raised.value)
    assert "3.22" in str(raised.value)
    assert "databases.cdd" in str(raised.value)
