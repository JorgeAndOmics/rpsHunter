"""The query download fetches exactly the accession asked for, or fails."""

import os
import stat
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).parents[2] / "workflow" / "scripts" / "query_downloader.sh"


def _fake_efetch(folder: Path, reply: str, status: int = 0) -> dict[str, str]:
    """Put a stand-in `efetch` on PATH that prints `reply`; return the environment."""
    (folder / "reply.txt").write_text(reply)
    tool = folder / "efetch"
    tool.write_text(f"#!/usr/bin/env bash\ncat {folder}/reply.txt\nexit {status}\n")
    tool.chmod(tool.stat().st_mode | stat.S_IXUSR)
    return {**os.environ, "PATH": f"{folder}{os.pathsep}{os.environ['PATH']}"}


def _download(tmp_path: Path, accession: str, reply: str, status: int = 0):
    tools = tmp_path / "bin"
    tools.mkdir()
    out = tmp_path / "out"
    out.mkdir()
    command = ["bash", str(SCRIPT), accession, str(out), str(tmp_path / "record.tsv")]
    environment = _fake_efetch(tools, reply, status)
    done = subprocess.run(
        command, env=environment, capture_output=True, text=True, check=False
    )
    return done, out


def test_the_record_asked_for_is_saved_and_logged_with_its_checksum(tmp_path):
    done, out = _download(tmp_path, "NP_000001.1", ">NP_000001.1 a protein\nMKV\n")

    assert done.returncode == 0
    assert (out / "NP_000001.1.fa").read_text() == ">NP_000001.1 a protein\nMKV\n"
    record = (tmp_path / "record.tsv").read_text().split("\t")
    assert record[0] == "NP_000001.1"
    assert len(record[1]) == 64


def test_a_reply_that_is_not_the_record_fails_and_leaves_no_file(tmp_path):
    done, out = _download(tmp_path, "NP_000001.1", "<html>error</html>\n")

    assert done.returncode == 66
    assert "NP_000001.1" in done.stderr
    assert list(out.iterdir()) == []


def test_an_empty_reply_fails(tmp_path):
    done, out = _download(tmp_path, "NP_000001.1", "")

    assert done.returncode == 66
    assert list(out.iterdir()) == []


def test_a_failing_download_fails_and_leaves_no_file(tmp_path):
    done, out = _download(tmp_path, "NP_000001.1", "", status=1)

    assert done.returncode != 0
    assert list(out.iterdir()) == []


def test_a_name_or_an_accession_without_version_is_refused_before_any_download(
    tmp_path,
):
    for bad in ("kinase", "NP_000001"):
        work = tmp_path / bad
        work.mkdir()

        done, _ = _download(work, bad, ">x\n")

        assert done.returncode == 65
        assert "versioned" in done.stderr
