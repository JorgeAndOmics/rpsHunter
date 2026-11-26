"""The run's settings and paths, read once from the config file.

Every script and the Snakefile import this module, so it has to be importable
before any argument is parsed: the config is the file named by the environment
variable RPSHUNTER_CONFIG (the launcher sets it), or `data/config/config.yaml`.
One complete file; nothing is merged.

Importing has no side effect on disk. `ensure_dirs()` is the one place in the
project that creates directories.

Raises on import (the checks before a run report these first, in plain words):
    DomainConfigError, QueryConfigError, SpeciesConfigError: the `domains:`,
        `queries:` or `species:` block is not usable.
    KeyError: a required section is missing.
"""

import os
from pathlib import Path
from typing import Any

import yaml

import domain_families
import query_config
import species_config

# Configuration file path. RPSHUNTER_CONFIG selects another complete config (for
# example a machine-local one); it is an environment variable because this module
# is imported by the Snakefile and by every script, before any argument parsing.
_DEFAULT_CONFIG_FILE: Path = (
    Path(__file__).parents[2] / "data" / "config" / "config.yaml"
)
CONFIG_FILE: Path = Path(
    os.environ.get("RPSHUNTER_CONFIG") or _DEFAULT_CONFIG_FILE
).resolve()

with CONFIG_FILE.open() as f:
    config: dict[str, Any] = yaml.safe_load(f)

# BLAST
E_VALUE_THRESHOLD: float = config["blast"].get("e_value", 0.01)
PERC_IDENTITY_THRESHOLD: int = config["blast"].get("perc_identity", 60)
SEQ_LENGTH_THRESHOLD: int = config["blast"].get("seq_length", 50)
BITSCORE_THRESHOLD: int = config["blast"].get("bitscore", 70)

# RPSBLAST - Maximum sensitivity settings for short domain detection
# RPS-BLAST. E-values are on the scale of the full CDD release (ADR-006): the search
# is told that release's size, CDD_RESIDUES, because the database it runs against
# holds only the configured families' models.
RPSBLAST_E_VALUE: float = config.get("rpsblast", {}).get("e_value", 1)
CDD_RESIDUES: int = config["databases"]["cdd"]["residues"]
RPSBLAST_COMP_BASED_STATS: int = config.get("rpsblast", {}).get("comp_based_stats", 0)
RPSBLAST_SEG: str = config.get("rpsblast", {}).get("seg", "no")
RPSBLAST_WINDOW_SIZE: int = config.get("rpsblast", {}).get("window_size", 40)

# Domain families: which CDD and Pfam models evidence which family (domain_families.py).
FAMILIES: dict[str, domain_families.Family] = domain_families.parse_families(
    config.get("domains")
)
PFAM_ACCESSIONS: list[str] = list(domain_families.family_by_pfam(FAMILIES))

# Directories - using pathlib.Path
PATH_DICT: dict[str, Path] = {
    "ROOT": Path(config["root"].get("db_root_folder")).resolve()
}

# === Root Directories ===
PATH_DICT["DATA_DIR"] = Path(
    config["root"].get("data_root_folder", PATH_DICT["ROOT"] / "data")
).resolve()
PATH_DICT["RESULTS_DIR"] = Path(
    config["root"].get("results_root_folder", PATH_DICT["ROOT"] / "results")
).resolve()
PATH_DICT["LOG_DIR"] = Path(
    config["root"].get("logs_root_folder", PATH_DICT["ROOT"] / "logs")
).resolve()
PATH_DICT["WORKFLOW_DIR"] = Path(__file__).parent
# The launcher's log of each run (everything printed) and its warnings.
PATH_DICT["RUN_LOG_DIR"] = PATH_DICT["LOG_DIR"] / "runs"

# === Database Directories ===
PATH_DICT["SPECIES_DB"] = PATH_DICT["ROOT"]
PATH_DICT["HMM_DB"] = PATH_DICT["ROOT"] / "hmm_dbs"
PATH_DICT["RPS_DB"] = PATH_DICT["ROOT"] / "cdd_dbs" / "Cdd"
PATH_DICT["RPSBPROC_DB"] = PATH_DICT["ROOT"] / "rpsbproc_dbs"
PATH_DICT["SMP_DIR"] = PATH_DICT["ROOT"] / "cdd_smp"
PATH_DICT["CDD_SUBSET_DB"] = PATH_DICT["ROOT"] / "cdd_subset"

# === Data Subdirectories ===
PATH_DICT["CONFIG_DIR"] = PATH_DICT["DATA_DIR"] / "config"
PATH_DICT["INPUT_DIR"] = PATH_DICT["DATA_DIR"] / "input"
PATH_DICT["FASTA_DIR"] = PATH_DICT["INPUT_DIR"] / "fastas"

# === Results Subdirectories ===
PATH_DICT["TABLE_OUTPUT_DIR"] = PATH_DICT["RESULTS_DIR"] / "tables"
PATH_DICT["FASTA_OUTPUT_DIR"] = PATH_DICT["RESULTS_DIR"] / "fastas"
PATH_DICT["ASN_ROOT_OUTPUT_DIR"] = PATH_DICT["RESULTS_DIR"] / "asn"
PATH_DICT["RPSBPROC_OUTPUT_DIR"] = PATH_DICT["RESULTS_DIR"] / "rpsbproc"
PATH_DICT["RANGE_OUTPUT_DIR"] = PATH_DICT["RESULTS_DIR"] / "ranges"
# === Table Stage Subdirectories ===
PATH_DICT["TABLE_BLAST_DIR"] = PATH_DICT["TABLE_OUTPUT_DIR"] / "blast"
PATH_DICT["TABLE_RPSBLAST_DIR"] = PATH_DICT["TABLE_OUTPUT_DIR"] / "rpsblast"
PATH_DICT["TABLE_DOMAINS_DIR"] = PATH_DICT["TABLE_OUTPUT_DIR"] / "domains"
PATH_DICT["TABLE_LOCI_DIR"] = PATH_DICT["TABLE_OUTPUT_DIR"] / "loci"

# === Per-Species Intermediate Directories ===
# Each step that produces per-species parquets before aggregation gets its own directory.
# Aggregate rules in Snakefile combine these into results/tables/*.parquet.
PATH_DICT["BLAST_SPECIES_DIR"] = PATH_DICT["RESULTS_DIR"] / "blast"
PATH_DICT["HMM_OUTPUT_DIR"] = PATH_DICT["RESULTS_DIR"] / "hmmer"
PATH_DICT["RPSBLAST_SPECIES_DIR"] = PATH_DICT["RESULTS_DIR"] / "rpsblast"
PATH_DICT["DOMAINS_SPECIES_DIR"] = PATH_DICT["RESULTS_DIR"] / "domains"

# === ASN Subdirectories ===
PATH_DICT["ASN_TBLASTN_DIR"] = PATH_DICT["ASN_ROOT_OUTPUT_DIR"] / "tblastn"
PATH_DICT["ASN_RPSBLAST_DIR"] = PATH_DICT["ASN_ROOT_OUTPUT_DIR"] / "rpsblast"


# =============================================================================
# Query Configuration (multi-query support)
# =============================================================================


# Queries: accession to the label that names its folders (query_config.py).
QUERY_LABELS: dict[str, str] = query_config.parse_queries(config.get("queries"))
QUERIES: list[str] = list(QUERY_LABELS)
QUERY_LABEL_LIST: list[str] = list(QUERY_LABELS.values())
LABEL_TO_ACC: dict[str, str] = {
    label: accession for accession, label in QUERY_LABELS.items()
}


# =============================================================================
# Per-Query Directory Generation
# =============================================================================
# Directories that hold per-query outputs need {query_label} subdirectories.
# These must be created after QUERY_LABEL_LIST is available.

_PER_QUERY_PARENTS = [
    PATH_DICT["BLAST_SPECIES_DIR"],  # results/blast/{ql}/
    PATH_DICT["HMM_OUTPUT_DIR"],  # results/hmmer/{ql}/
    PATH_DICT["FASTA_OUTPUT_DIR"],  # results/fastas/{ql}/
    PATH_DICT["RPSBLAST_SPECIES_DIR"],  # results/rpsblast/{ql}/
    PATH_DICT["RPSBPROC_OUTPUT_DIR"],  # results/rpsbproc/{ql}/
    PATH_DICT["DOMAINS_SPECIES_DIR"],  # results/domains/{ql}/
    PATH_DICT["ASN_TBLASTN_DIR"],  # results/asn/tblastn/{ql}/
    PATH_DICT["ASN_RPSBLAST_DIR"],  # results/asn/rpsblast/{ql}/
    PATH_DICT["TABLE_BLAST_DIR"],  # results/tables/blast/{ql}/
    PATH_DICT["TABLE_RPSBLAST_DIR"],  # results/tables/rpsblast/{ql}/
]

# Also: results/tables/blast/{ql}/selected/
_PER_QUERY_NESTED = [
    PATH_DICT["TABLE_BLAST_DIR"] / "{ql}" / "selected",
]

# PATH_DICT entries that are file-name prefixes rather than folders (cdd_dbs/Cdd.*):
# only their parent folder is created.
_FILE_PREFIX_KEYS = frozenset({"RPS_DB"})


def ensure_dirs() -> None:
    """Create every directory the pipeline writes to.

    This is the only place in the project that creates directories. It is called
    once by the launcher and once when the Snakefile is parsed, so importing this
    module stays free of side effects and scripts can rely on their output folders
    existing. A new output folder is added to PATH_DICT (or to the per-query lists
    above), never created by the script that writes into it.
    """
    for key, path in PATH_DICT.items():
        folder = path.parent if key in _FILE_PREFIX_KEYS else path
        folder.mkdir(parents=True, exist_ok=True)
    for label in QUERY_LABEL_LIST:
        for parent in _PER_QUERY_PARENTS:
            (parent / label).mkdir(parents=True, exist_ok=True)
        for template in _PER_QUERY_NESTED:
            Path(str(template).replace("{ql}", label)).mkdir(
                parents=True, exist_ok=True
            )


def query_file_path(accession: str) -> Path:
    """Return where the protein FASTA of a query accession is kept."""
    return PATH_DICT["FASTA_DIR"] / f"{accession}.fa"


# Execution and requests
NUM_CORES: int = config["execution"].get("num_cores", 1)
ENTREZ_EMAIL: str = config["execution"].get("entrez_email", "")

# === HMMER Configuration ===
HMMER_ENABLED: bool = config.get("hmmer", {}).get("enabled", False)
HMMER_USE_GA: bool = config.get("hmmer", {}).get("use_gathering_threshold", False)
HMMER_EVALUE: float = config.get("hmmer", {}).get("evalue", 1e-5)
HMMER_DOM_EVALUE: float = config.get("hmmer", {}).get("dom_evalue", 1e-3)
HMMER_MIN_SCORE = config.get("hmmer", {}).get("min_score", None)
HMMER_MIN_COVERAGE: float = config.get("hmmer", {}).get("min_coverage", 0.5)
HMMER_MIN_ALN_LEN: int = config.get("hmmer", {}).get("min_alignment_length", 20)
HMMER_MAX_SENS: bool = config.get("hmmer", {}).get("max_sensitivity", False)
HMMER_BIAS_FILTER: bool = config.get("hmmer", {}).get("bias_filter", True)
HMMER_SEED: int = config.get("hmmer", {}).get("seed", 67)

# === Loci: largest gap (bp) between gated hits chained into one locus (ADR-008) ===
LOCUS_GAP: int = config.get("loci", {}).get("gap", 20000)

# === External Program Names ===
_programs = config.get("programs", {})
TBLASTN_CMD: str = _programs.get("tblastn", "tblastn")
RPSBLAST_CMD: str = _programs.get("rpsblast", "rpsblast")
BLAST_FORMATTER_CMD: str = _programs.get("blast_formatter", "blast_formatter")
MAKEBLASTDB_CMD: str = _programs.get("makeblastdb", "makeblastdb")
RPSBPROC_CMD: str = _programs.get("rpsbproc", "rpsbproc")
MAKEPROFILEDB_CMD: str = _programs.get("makeprofiledb", "makeprofiledb")

# Genomes. Every species carries a display name and a pin (assembly accession
# and/or checksum of a local FASTA); see species_config.py for the contract.
SPECIES_ENTRIES: dict[str, species_config.SpeciesEntry] = species_config.parse_species(
    config.get("species")
)
SPECIES_DICT: dict[str, str] = {
    key: entry.name for key, entry in SPECIES_ENTRIES.items()
}
SPECIES: list[str] = list(SPECIES_ENTRIES)
