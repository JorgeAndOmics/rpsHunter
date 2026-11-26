"""Guards on the legacy list in pyproject.toml and .lintr.R.

Scripts that predate the quality gate are excluded from ruff, mypy, coverage
and lintr as whole files until they are rewritten. These tests keep that
arrangement honest: the copies of the list must agree, every entry must still
exist (so a rewritten or deleted file cannot linger on the list), and the list
may never grow beyond what it was when the gate was introduced.
"""

import re
import tomllib
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = "workflow/scripts/"

# The list as it stood when the gate was introduced. Entries are removed from
# pyproject.toml and .lintr.R as files are rewritten; nothing is ever added here.
ORIGINAL_PYTHON = frozenset(
    {
        "cdd_subset.py",
        "concordance.py",
        "defaults.py",
        "hit_domain_inventory.py",
        "rpsblast.py",
        "rpsbproc.py",
        "rpsbproc_parser.py",
    }
)
ORIGINAL_R = frozenset(
    {
        "completeness_detector.R",
        "contingency_sorter.R",
        "extended_plots.R",
        "merge_queries.R",
        "orf_analyser.R",
        "tests/orf_analyser_test.R",
        "tests/ranges_sorter_test.R",
    }
)


def _pyproject() -> dict:
    return tomllib.loads((REPO / "pyproject.toml").read_text())


def _script_names(paths: list[str]) -> set[str]:
    return {p.removeprefix(SCRIPTS) for p in paths if p.startswith(SCRIPTS)}


def _ruff_legacy() -> set[str]:
    return _script_names(_pyproject()["tool"]["ruff"]["extend-exclude"])


def _coverage_legacy() -> set[str]:
    return _script_names(_pyproject()["tool"]["coverage"]["run"]["omit"])


def _mypy_legacy() -> set[str]:
    patterns = _pyproject()["tool"]["mypy"]["exclude"]
    (legacy,) = [p for p in patterns if p.startswith(SCRIPTS)]
    alternatives = re.search(r"\((.*)\)", legacy)
    assert alternatives is not None
    return {f"{stem}.py" for stem in alternatives.group(1).split("|")}


def _lintr_legacy() -> set[str]:
    text = (REPO / ".lintr.R").read_text()
    return _script_names(re.findall(r'"(workflow/scripts/[^"]+\.R)"', text))


def test_python_legacy_lists_agree():
    assert _ruff_legacy() == _mypy_legacy() == _coverage_legacy()


def test_legacy_entries_still_exist():
    for name in _ruff_legacy() | _lintr_legacy():
        assert (REPO / SCRIPTS / name).is_file(), f"stale legacy entry: {name}"


def test_legacy_lists_never_grow():
    assert _ruff_legacy() <= ORIGINAL_PYTHON
    assert _lintr_legacy() <= ORIGINAL_R
