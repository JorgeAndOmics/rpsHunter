"""The manifest of a run: what ran, with what, well enough to repeat it.

Written by the launcher after every run that executed jobs, to
`<results>/manifests/<date>_<time>.yaml`. It records the code (version and
commit), the command, a copy of the config with its checksum, the installed
versions of the analysis tools, the database releases with checksums of the
models actually searched, the genomes and the queries.

Genome FASTA files are not hashed again here (that is tens of gigabytes per
run): the checksum pinned in the config is recorded, and `genome_check.py`
verifies it whenever a BLAST database is built.
"""

from __future__ import annotations

import hashlib
import sys
from collections.abc import Callable, Collection, Iterable
from pathlib import Path
from typing import Any, NamedTuple

import yaml

import external
import preflight

REPO = Path(__file__).resolve().parents[2]
ANALYSIS_PACKAGES = (
    "blast",
    "hmmer",
    "rpsbproc",
    "snakemake-minimal",
    "python",
    "pandas",
    "pyarrow",
)
"""Conda packages whose version can change a result."""

_CHUNK = 1 << 20


class Run(NamedTuple):
    """What the launcher knows about the run it just made."""

    command: list[str]
    stages: list[str]
    started: str
    seconds: float
    exit_status: int
    jobs: dict[str, int]


def file_sha256(path: Path) -> str | None:
    """Return a file's SHA-256, or None if it is not there."""
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def code_state(
    repo: Path, run: Callable[..., Any] = external.run_plain
) -> dict[str, Any]:
    """Return the commit of this copy of the code and whether it was modified."""
    options = {"cwd": repo, "capture_output": True, "text": True}
    head = run(["git", "rev-parse", "HEAD"], **options)
    if head.returncode != 0:
        return {"commit": None, "uncommitted_changes": None}
    status = run(["git", "status", "--porcelain", "--untracked-files=no"], **options)
    return {
        "commit": head.stdout.strip(),
        "uncommitted_changes": bool(status.stdout.strip()),
    }


def package_versions(
    installed: Collection[str], names: Iterable[str]
) -> dict[str, str]:
    """Return `<version>-<build>` of each named package among those installed.

    Args:
        installed: Packages as `<name>-<version>-<build>` (conda's own naming).
        names: Package names to look for; one not installed is left out.
    """
    by_name = {package.rsplit("-", 2)[0]: package for package in installed}
    return {
        name: by_name[name].removeprefix(f"{name}-")
        for name in names
        if name in by_name
    }


def profile_accessions(profiles_text: str) -> list[str]:
    """Return the versioned accessions of the profiles in an HMM file."""
    return [
        line.split()[1] for line in profiles_text.splitlines() if line.startswith("ACC")
    ]


def _first_line(path: Path) -> str | None:
    """The first line of a small text file, or None if it is not there."""
    return path.read_text().splitlines()[0] if path.is_file() else None


def _databases(settings: Any) -> dict[str, Any]:
    """Releases and checksums of the models the run searched."""
    paths = settings.PATH_DICT
    subset = paths["HMM_OUTPUT_DIR"] / "subset.hmm"
    families = settings.FAMILIES.values()
    return {
        "cdd": {
            "version": settings.CDD_VERSION,
            "residues": settings.CDD_RESIDUES,
            "installed_release": _first_line(paths["RPSBPROC_DB"] / "cdd.info"),
            "models_sha256": {
                accession: file_sha256(paths["SMP_DIR"] / f"{accession}.smp")
                for family in families
                for accession in family.cdd
            },
        },
        "pfam": {
            "version": settings.PFAM_VERSION,
            "profiles_sha256": file_sha256(subset),
            "profiles": profile_accessions(subset.read_text())
            if subset.is_file()
            else [],
        },
    }


def _genomes(settings: Any) -> dict[str, Any]:
    """Each genome as pinned in the config, with the size of its FASTA on disk."""
    genomes = {}
    for key, entry in settings.SPECIES_ENTRIES.items():
        fasta = settings.PATH_DICT["SPECIES_DB"] / f"{key}.fa"
        genomes[key] = {
            "name": entry.name,
            "assembly": entry.assembly,
            "sha256_pinned": entry.sha256,
            "fasta_bytes": fasta.stat().st_size if fasta.is_file() else None,
        }
    return genomes


def _queries(settings: Any) -> dict[str, Any]:
    """Each query protein with its label and the checksum of its FASTA."""
    return {
        accession: {
            "label": label,
            "sha256": file_sha256(settings.query_file_path(accession)),
        }
        for accession, label in settings.QUERY_LABELS.items()
    }


def collect(run: Run, config_file: Path) -> dict[str, Any]:
    """Gather the manifest of one run as plain data."""
    settings = preflight.load_settings(config_file)
    installed = preflight.installed_packages()
    lock_note = preflight.environment_note(preflight.LOCK.read_text(), installed)
    declared = yaml.safe_load((REPO / "CITATION.cff").read_text())
    return {
        "rpshunter": {"version": str(declared.get("version")), **code_state(REPO)},
        "run": {**run._asdict(), "seconds": round(run.seconds, 1)},
        "config": {
            "file": str(config_file),
            "sha256": file_sha256(config_file),
            "settings": yaml.safe_load(config_file.read_text()),
        },
        "environment": {
            "python": sys.version.split()[0],
            "matches_lock": lock_note is None,
            "packages": package_versions(installed, ANALYSIS_PACKAGES),
        },
        "databases": _databases(settings),
        "genomes": _genomes(settings),
        "queries": _queries(settings),
    }


def write(record: dict[str, Any], folder: Path, stamp: str) -> Path:
    """Write a manifest as YAML named by the run's time stamp; return its path."""
    path = folder / f"{stamp}.yaml"
    path.write_text(yaml.safe_dump(record, sort_keys=False))
    return path


def record(run: Run, config_file: Path, folder: Path, stamp: str) -> Path:
    """Collect and write the manifest of a run; return where it was written."""
    return write(collect(run, config_file), folder, stamp)
