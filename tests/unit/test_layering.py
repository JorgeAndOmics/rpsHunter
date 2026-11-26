"""Pure modules hold logic only.

A pure module may not import the configuration and path layer, the logging
setup, or anything that reaches processes, files, the clock, randomness or the
network. That keeps the logic testable without a pipeline run and keeps side
effects at the edges. Add each new pure module to PURE_MODULES.

(import-linter would be the usual tool, but it analyses packages, and these
scripts are flat modules that Snakemake runs as files.)
"""

import ast
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"

PURE_MODULES = ("tags",)

FORBIDDEN = frozenset(
    {
        "defaults",
        "colored_logging",
        "subprocess",
        "os",
        "pathlib",
        "shutil",
        "tempfile",
        "logging",
        "random",
        "time",
        "datetime",
        "socket",
        "urllib",
        "requests",
        "yaml",
        "Bio",
    }
)


def _imported_top_level_names(source: str) -> set[str]:
    """Top-level package names imported anywhere in a module's source."""
    names: set[str] = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            names.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module.split(".")[0])
    return names


@pytest.mark.parametrize("module", PURE_MODULES)
def test_pure_module_imports_nothing_that_does_io(module):
    imported = _imported_top_level_names((SCRIPTS / f"{module}.py").read_text())

    assert imported & FORBIDDEN == set()


def test_the_check_sees_both_import_forms():
    source = "import os.path\nfrom Bio import SeqIO\nimport pandas as pd\n"

    assert _imported_top_level_names(source) == {"os", "Bio", "pandas"}
