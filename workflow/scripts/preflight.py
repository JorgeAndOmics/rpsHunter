"""Checks before a run: is the config sound, are the tools and inputs there?

Every check is quick and offline, none asks a question, and all problems are
reported together so that one launch shows everything there is to fix. A
problem stops the run; a note is shown and the run goes on.

Order: keys of earlier versions, the schema, the rules the schema cannot
express (families, species). Only a sound config is then used to look for the
tools and the input files of the requested stages.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import sys
from collections.abc import Callable, Collection, Iterable, Mapping
from pathlib import Path
from typing import Any, NamedTuple

import yamale
import yaml

import domain_families
import species_config
import stages

REPO = Path(__file__).resolve().parents[2]
SCHEMA = REPO / "data" / "config" / "schema.yaml"
LOCK = REPO / "data" / "config" / "environment.lock-linux-64.txt"
CONFIG_VARIABLE = "RPSHUNTER_CONFIG"
MIN_FREE_BYTES = 5 * 1024**3
"""Below this much free space in the results folder a note is shown (a model-genome
run writes about 3 GB)."""

_FAMILIES = "replaced by the `domains:` block (docs/configuration.md)"
_GONE = "no longer used; delete it"
RETIRED_KEYS: dict[tuple[str, ...], str] = {
    ("orf",): "the ORF stage was retired (docs/adr/ADR-011); delete the section",
    ("hit_domain_inventory",): "replaced by `loci.gap`",
    ("rpsblast", "target_domains"): _FAMILIES,
    ("hmmer", "profiles"): _FAMILIES,
    ("logging",): _GONE,
    (
        "random_id_length",
    ): "Tags are now derived from the hit (docs/adr/ADR-004); delete it",
    ("use_species_list",): _GONE,
    ("display", "display_requests_warning"): _GONE,
    ("display", "display_operation_info"): _GONE,
}
"""Config keys of earlier versions, by path, each with what to do instead."""

Which = Callable[[str], str | None]


class Findings(NamedTuple):
    """What the checks found: problems stop the run, notes do not."""

    problems: list[str]
    notes: list[str]


class Layout(NamedTuple):
    """Where the requested stages expect their input files."""

    genomes: Mapping[str, Path]
    queries: Mapping[str, Path]
    cdd_models: Mapping[str, Path]
    rpsbproc_data: Path
    pfam_library: Path
    hmmer_active: bool


# --- the config -----------------------------------------------------------------


def _has(config: Mapping[str, Any], path: tuple[str, ...]) -> bool:
    """Whether the config has a value at the key path."""
    node: Any = config
    for key in path:
        if not isinstance(node, Mapping) or key not in node:
            return False
        node = node[key]
    return True


def retired_key_problems(config: Mapping[str, Any]) -> list[str]:
    """Return one problem per key of an earlier version still in the config."""
    return [
        f"config key `{'.'.join(path)}`: {advice}"
        for path, advice in RETIRED_KEYS.items()
        if _has(config, path)
    ]


def schema_problems(config_file: Path, schema_file: Path) -> list[str]:
    """Return the config's violations of the schema, one per line; empty if valid."""
    try:
        data = yamale.make_data(str(config_file))
    except yaml.YAMLError as error:
        return [f"{config_file.name} is not valid YAML: {' '.join(str(error).split())}"]
    schema = yamale.make_schema(str(schema_file))
    try:
        yamale.validate(schema, data, strict=False)
    except yamale.YamaleError as error:
        return [f"config: {line}" for result in error.results for line in result.errors]
    return []


def content_problems(config: Mapping[str, Any]) -> list[str]:
    """Return violations of the rules the schema cannot express."""
    problems = []
    try:
        domain_families.parse_families(config.get("domains"))
    except domain_families.DomainConfigError as error:
        problems.append(str(error))
    try:
        species_config.parse_species(config.get("species"))
    except species_config.SpeciesConfigError as error:
        problems.append(str(error))
    return problems


# --- tools and inputs -------------------------------------------------------------


def missing_tools(
    tools: Iterable[str], configured: Mapping[str, str], which: Which
) -> list[str]:
    """Return the tools not on PATH, each looked up under its configured name."""
    names = [configured.get(tool, tool) for tool in tools]
    return [name for name in names if which(name) is None]


def _absent(files: Mapping[str, Path], what: str, fetch: str) -> list[str]:
    """One problem per missing file, naming the stage that fetches it."""
    return [
        f"{what} missing for {key}: {path} (fetch with {fetch}, or put the file there)"
        for key, path in files.items()
        if not path.is_file()
    ]


def _needs(request: stages.Request, reader: str, fetcher: str) -> bool:
    """Whether a requested stage reads inputs that this run does not also fetch."""
    flags = {stage.flag for stage in request.stages}
    return reader in flags and fetcher not in flags


def input_problems(request: stages.Request, layout: Layout) -> list[str]:
    """Return the input files the requested stages read and that are not there."""
    databases = "--download-databases"
    problems: list[str] = []
    if _needs(request, "--blast-dbs", "--download-genomes"):
        problems += _absent(layout.genomes, "genome FASTA", "--download-genomes")
    if _needs(request, "--search", "--download-queries"):
        problems += _absent(layout.queries, "query protein", "--download-queries")
    if _needs(request, "--domains", databases):
        problems += _absent(layout.cdd_models, "CDD model", databases)
        problems += _absent(
            {"rpsbproc": layout.rpsbproc_data}, "CDD annotation", databases
        )
    if layout.hmmer_active and _needs(request, "--hmmer", databases):
        problems += _absent({"HMMER": layout.pfam_library}, "Pfam library", databases)
    return problems


# --- notes ----------------------------------------------------------------------


def environment_note(lock_text: str, installed: Collection[str]) -> str | None:
    """Return a note if a package of the lock file is not installed as locked.

    Args:
        lock_text: The explicit conda lock (one package URL per line).
        installed: Installed packages as `<name>-<version>-<build>`.
    """
    locked = {
        line.split("#")[0]
        .rsplit("/", 1)[-1]
        .removesuffix(".conda")
        .removesuffix(".tar.bz2")
        for line in lock_text.splitlines()
        if line.startswith("http")
    }
    absent = sorted(locked - set(installed))
    if not absent:
        return None
    shown = ", ".join(absent[:5]) + (", ..." if len(absent) > 5 else "")
    return (
        f"the environment differs from the lock file in {len(absent)} package(s): "
        f"{shown}. Results may differ from a locked run; `make env` rebuilds it."
    )


def disk_note(free_bytes: int, folder: Path) -> str | None:
    """Return a note when the results folder's disk is nearly full."""
    if free_bytes >= MIN_FREE_BYTES:
        return None
    return f"only {free_bytes / 1024**3:.1f} GB free in {folder}"


def _installed_packages() -> set[str]:
    """Packages of the running conda environment, as `<name>-<version>-<build>`."""
    return {path.stem for path in (Path(sys.prefix) / "conda-meta").glob("*.json")}


def _notes(results: Path) -> list[str]:
    """The notes for this machine: environment against the lock, free space."""
    found = [
        environment_note(LOCK.read_text(), _installed_packages()),
        disk_note(shutil.disk_usage(results.anchor).free, results),
    ]
    return [note for note in found if note]


# --- the whole check --------------------------------------------------------------


def _settings(config_file: Path) -> Any:
    """Evaluate defaults.py for this config, apart from any copy already imported.

    defaults.py reads the config named by RPSHUNTER_CONFIG when it is executed.
    A private copy is executed here so that the answer is this config's even if
    the module was imported earlier for another.
    """
    spec = importlib.util.spec_from_file_location(
        "_preflight_defaults", Path(__file__).with_name("defaults.py")
    )
    if spec is None or spec.loader is None:
        raise ImportError("defaults.py cannot be loaded")
    module = importlib.util.module_from_spec(spec)
    before = os.environ.get(CONFIG_VARIABLE)
    os.environ[CONFIG_VARIABLE] = str(config_file)
    try:
        spec.loader.exec_module(module)
    finally:
        if before is None:
            del os.environ[CONFIG_VARIABLE]
        else:
            os.environ[CONFIG_VARIABLE] = before
    return module


def _layout(config_file: Path) -> tuple[Layout, Path]:
    """Read the pipeline's paths for a config already known to be sound."""
    defaults = _settings(config_file)
    paths = defaults.PATH_DICT
    layout = Layout(
        genomes={key: paths["SPECIES_DB"] / f"{key}.fa" for key in defaults.SPECIES},
        queries={acc: defaults.query_file_path(acc) for acc in defaults.QUERY_DICT},
        cdd_models={
            accession: paths["SMP_DIR"] / f"{accession}.smp"
            for accession in domain_families.cdd_accessions(defaults.FAMILIES)
        },
        rpsbproc_data=paths["RPSBPROC_DB"] / "cddid.tbl",
        pfam_library=paths["HMM_DB"] / "Pfam-A.hmm",
        hmmer_active=bool(defaults.HMMER_ENABLED and defaults.PFAM_ACCESSIONS),
    )
    return layout, paths["RESULTS_DIR"]


def check(
    request: stages.Request, config_file: Path, which: Which = shutil.which
) -> Findings:
    """Run every check for one launch.

    Args:
        request: The parsed command line.
        config_file: The run's config (exists).
        which: Finds a program on PATH; replaced in tests.
    """
    problems = schema_problems(config_file, SCHEMA)
    config = {} if problems else yaml.safe_load(config_file.read_text())
    problems += retired_key_problems(_raw(config_file) if problems else config)
    problems += [] if problems else content_problems(config)
    if problems:
        return Findings(problems, [])
    layout, results = _layout(config_file)
    absent = missing_tools(request.tools, config.get("programs") or {}, which)
    problems = [f"program not found on PATH: {name}" for name in absent]
    problems += input_problems(request, layout)
    return Findings(problems, _notes(results))


def _raw(config_file: Path) -> Mapping[str, Any]:
    """The config as a mapping, or an empty one when it cannot be read as such."""
    try:
        loaded = yaml.safe_load(config_file.read_text())
    except yaml.YAMLError:
        return {}
    return loaded if isinstance(loaded, Mapping) else {}
