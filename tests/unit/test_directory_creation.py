"""Directories are created in one place, and only when asked.

`defaults.py` owns every path of the pipeline. Importing it must have no side
effect (so that printing help or inspecting a config touches nothing), and no
other script, nor the Snakefile, may create directories of its own: a path
created somewhere else is a path `defaults.py` does not know about.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import yaml

import defaults

REPO = Path(__file__).resolve().parents[2]
WORKFLOW = REPO / "workflow"
CREATES_DIRECTORY = re.compile(r"\bmkdir\b|\bmakedirs\b|\bdir\.create\b")
OWNER = WORKFLOW / "scripts" / "defaults.py"


def _config_with_roots(tmp_path: Path) -> Path:
    """Write the shipped config with its four root folders under tmp_path."""
    config = yaml.safe_load((REPO / "data" / "config" / "config.yaml").read_text())
    config["root"] = {
        f"{name}_root_folder": str(tmp_path / name)
        for name in ("db", "data", "results", "logs")
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    return path


def _run_with_defaults(code: str, config: Path) -> None:
    """Run `code` in a fresh interpreter that can import `defaults`."""
    environment = {**os.environ, "RPSHUNTER_CONFIG": str(config)}
    subprocess.run(
        [sys.executable, "-c", code],
        cwd=WORKFLOW / "scripts",
        env=environment,
        check=True,
    )


def test_only_defaults_creates_directories():
    sources = [REPO / "rpsHunter", WORKFLOW / "Snakefile"]
    sources += sorted((WORKFLOW / "scripts").rglob("*.py"))
    sources += sorted((WORKFLOW / "scripts").rglob("*.R"))
    sources += sorted((WORKFLOW / "scripts").rglob("*.sh"))
    offenders = [
        f"{path.relative_to(REPO)}:{number}"
        for path in sources
        if path != OWNER
        for number, line in enumerate(path.read_text().splitlines(), start=1)
        if CREATES_DIRECTORY.search(line)
    ]

    assert offenders == []


def test_importing_defaults_creates_nothing(tmp_path):
    config = _config_with_roots(tmp_path)

    _run_with_defaults("import defaults", config)

    assert sorted(p.name for p in tmp_path.iterdir()) == ["config.yaml"]


def test_ensure_dirs_creates_the_tree_including_per_query_folders(tmp_path):
    config = _config_with_roots(tmp_path)

    _run_with_defaults("import defaults; defaults.ensure_dirs()", config)

    assert (tmp_path / "results" / "tables" / "domains").is_dir()
    assert (tmp_path / "results" / "blast" / "Human_PRDM9").is_dir()
    assert (
        tmp_path / "results" / "tables" / "blast" / "Cow_PRDM9" / "selected"
    ).is_dir()
    assert (tmp_path / "db" / "cdd_dbs").is_dir()
    # RPS_DB is a database name prefix (cdd_dbs/Cdd.*), not a folder.
    assert not (tmp_path / "db" / "cdd_dbs" / "Cdd").exists()


def test_ensure_dirs_creates_folders_and_only_the_parent_of_a_file_prefix(
    tmp_path, monkeypatch
):
    paths = {
        "RESULTS_DIR": tmp_path / "results",
        "RPS_DB": tmp_path / "cdd_dbs" / "Cdd",
    }
    monkeypatch.setattr(defaults, "PATH_DICT", paths)
    monkeypatch.setattr(defaults, "QUERY_LABEL_LIST", ["query_1"])
    monkeypatch.setattr(
        defaults, "_PER_QUERY_PARENTS", [tmp_path / "results" / "blast"]
    )
    monkeypatch.setattr(
        defaults, "_PER_QUERY_NESTED", [tmp_path / "tables" / "{ql}" / "selected"]
    )

    defaults.ensure_dirs()

    assert (tmp_path / "results" / "blast" / "query_1").is_dir()
    assert (tmp_path / "tables" / "query_1" / "selected").is_dir()
    assert (tmp_path / "cdd_dbs").is_dir()
    assert not (tmp_path / "cdd_dbs" / "Cdd").exists()
