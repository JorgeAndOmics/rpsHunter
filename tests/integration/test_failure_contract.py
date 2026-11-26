"""A search tool that fails must fail the job, not produce an empty result.

Before this contract, a crashed tBLASTn or RPS-BLAST was logged and then
replaced by an empty table, so a species whose search had failed looked exactly
like a species with no hits. These tests run the real scripts with the search
program swapped for one that always fails (`false`). The failure is one ERROR line
in the job's output, followed by the traceback for the log.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.integration

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "workflow" / "scripts"


def _config(tmp_path: Path, failing_program: str) -> Path:
    """Write the shipped config with tmp roots and one program replaced by `false`."""
    config = yaml.safe_load((REPO / "data" / "config" / "config.yaml").read_text())
    config["root"] = {
        f"{name}_root_folder": str(tmp_path / name)
        for name in ("db", "data", "results", "logs")
    }
    config["programs"] = {failing_program: "false"}
    config["species"] = {"Toyus_plenus": {"name": "Toyus plenus", "sha256": "0" * 64}}
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    return path


def _run_script(script: str, arguments: list[str], config: Path):
    """Run a pipeline script as Snakemake would, after creating the output tree."""
    environment = {**os.environ, "RPSHUNTER_CONFIG": str(config)}
    subprocess.run(
        [sys.executable, "-c", "import defaults; defaults.ensure_dirs()"],
        cwd=SCRIPTS,
        env=environment,
        check=True,
    )
    return subprocess.run(
        [sys.executable, str(SCRIPTS / script), *arguments],
        cwd=REPO / "workflow",
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def _error_lines(output: str) -> list[str]:
    """The ERROR lines of a job's output, without their time."""
    lines = [line.split(" ", 1)[1] for line in output.splitlines() if " " in line]
    return [line for line in lines if line.startswith("ERROR ")]


def test_a_failed_tblastn_fails_the_job_and_writes_no_table(tmp_path):
    config = _config(tmp_path, "tblastn")
    query = tmp_path / "query.fa"
    query.write_text(">q\nMKV\n")
    table = tmp_path / "hits.parquet"

    done = _run_script(
        "tblastn.py",
        [
            "--species", "Toyus_plenus",
            "--query-file", str(query),
            "--query-accession", "q",
            "--output-parquet", str(table),
            "--output-csv", str(tmp_path / "hits.csv"),
            "--output-asn", str(tmp_path / "hits.asn"),
        ],
        config,
    )  # fmt: skip

    assert done.returncode != 0
    assert "tblastn exited with status 1" in done.stdout
    assert _error_lines(done.stdout)[0].startswith("ERROR tblastn all | ToolError")
    assert not table.exists()


def test_a_failed_rpsblast_fails_the_job_and_writes_no_table(tmp_path):
    config = _config(tmp_path, "rpsblast")
    sequences = tmp_path / "Toyus_plenus.fa"
    sequences.write_text(">scaffold:1-9|tag:AAAAAAAAAAAA\nMKV\n")
    table = tmp_path / "domains.parquet"

    done = _run_script(
        "rpsblast.py",
        [
            "--species", "Toyus_plenus",
            "--query-accession", "q",
            "--fasta-input", str(sequences),
            "--output-parquet", str(table),
            "--output-asn", str(tmp_path / "Toyus_plenus.asn"),
            "--db", str(tmp_path / "no_such_db"),
        ],
        config,
    )  # fmt: skip

    assert done.returncode != 0
    assert "rpsblast exited with status 1" in done.stdout
    assert not table.exists()


def test_a_genome_that_is_not_the_pinned_one_is_refused(tmp_path):
    fasta = tmp_path / "Toyus_plenus.fa"
    fasta.write_text(">scaffold\nACGTACGT\n")
    config = yaml.safe_load((REPO / "data" / "config" / "config.yaml").read_text())
    config["root"] = {
        f"{name}_root_folder": str(tmp_path / name)
        for name in ("db", "data", "results", "logs")
    }
    config["species"] = {"Toyus_plenus": {"name": "Toyus plenus", "sha256": "0" * 64}}
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))

    done = _run_script(
        "genome_check.py",
        ["--species", "Toyus_plenus", "--fasta", str(fasta)],
        config_path,
    )

    assert done.returncode == 1
    assert "not the one the config pins" in done.stdout
    assert "0" * 64 in done.stdout
    assert len(_error_lines(done.stdout)) == 1
