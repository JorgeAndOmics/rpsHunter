"""
Defaults Configuration Script
=============================

This script loads configuration settings from a YAML file and sets up various
constants and directory paths used throughout the project.

Modules:
    - `yaml`: For parsing the YAML configuration file.
    - `pathlib.Path`: For robust path handling.

Configuration:
    - The configuration file is expected to be located at `../data/config/config.yaml`.
    - Various constants and directory paths are initialized based on the configuration file.

Usage:
    This script is intended to be imported as a module and not run directly.
"""

from pathlib import Path
from typing import Dict, Any, List

import yaml

# Configuration file path
CONFIG_FILE: Path = Path(__file__).parents[2] / 'data' / 'config' / 'config.yaml'

with open(CONFIG_FILE, 'r') as f:
    config: Dict[str, Any] = yaml.safe_load(f)

# QUERY SEQUENCE
QUERY_ACC: str = config['query'].get('accession')

# BLAST
E_VALUE_THRESHOLD: float = config['blast'].get('e_value', 0.01)
PERC_IDENTITY_THRESHOLD: int = config['blast'].get('perc_identity', 60)
SEQ_LENGTH_THRESHOLD: int = config['blast'].get('seq_length', 50)
BITSCORE_THRESHOLD: int = config['blast'].get('bitscore', 70)
ACCESSION_ID_REGEX: str = r'[A-Z]{2,}_?[0-9]+\.[0-9]{1,2}'

# Directories - using pathlib.Path
PATH_DICT: Dict[str, Path] = {
    'ROOT': Path(config['root'].get('db_root_folder')).resolve()
}

# === Root Directories ===
PATH_DICT['DATA_DIR'] = Path(config['root'].get('data_root_folder', PATH_DICT['ROOT'] / 'data')).resolve()
PATH_DICT['RESULTS_DIR'] = Path(config['root'].get('results_root_folder', PATH_DICT['ROOT'] / 'results')).resolve()
PATH_DICT['LOG_DIR'] = Path(config['root'].get('logs_root_folder', PATH_DICT['ROOT'] / 'logs')).resolve()
PATH_DICT['WORKFLOW_DIR'] = Path(__file__).parent

# === Database Directories ===
PATH_DICT['SPECIES_DB'] = PATH_DICT['ROOT']
PATH_DICT['HMM_DB'] = PATH_DICT['ROOT'] / 'hmm_dbs'
PATH_DICT['RPS_DB'] = PATH_DICT['ROOT'] / 'cdd_dbs' / 'Cdd'
PATH_DICT['RPSBPROC_DB'] = PATH_DICT['ROOT'] / 'rpsbproc_dbs'

# === Data Subdirectories ===
PATH_DICT['CONFIG_DIR'] = PATH_DICT['DATA_DIR'] / 'config'
PATH_DICT['INPUT_DIR'] = PATH_DICT['DATA_DIR'] / 'input'
PATH_DICT['FASTA_DIR'] = PATH_DICT['INPUT_DIR'] / 'fastas'
PATH_DICT['TMP_DIR'] = PATH_DICT['DATA_DIR'] / 'tmp'

# === Results Subdirectories ===
PATH_DICT['PLOT_DIR'] = PATH_DICT['RESULTS_DIR'] / 'plots'
PATH_DICT['TABLE_OUTPUT_DIR'] = PATH_DICT['RESULTS_DIR'] / 'tables'
PATH_DICT['FASTA_OUTPUT_DIR'] = PATH_DICT['RESULTS_DIR'] / 'fastas'
PATH_DICT['ASN_ROOT_OUTPUT_DIR'] = PATH_DICT['RESULTS_DIR'] / 'asn'
PATH_DICT['RPSBPROC_OUTPUT_DIR'] = PATH_DICT['RESULTS_DIR'] / 'rpsbproc'
PATH_DICT['XML_OUTPUT_DIR'] = PATH_DICT['RESULTS_DIR'] / 'xml'
PATH_DICT['RANGE_OUTPUT_DIR'] = PATH_DICT['RESULTS_DIR'] / 'ranges'
PATH_DICT['LOCI_TABLE_OUTPUT_DIR'] = PATH_DICT['TABLE_OUTPUT_DIR'] / 'loci_tables'

# === Per-Species Intermediate Directories ===
# Each step that produces per-species parquets before aggregation gets its own directory.
# Aggregate rules in Snakefile combine these into results/tables/*.parquet.
PATH_DICT['BLAST_SPECIES_DIR']    = PATH_DICT['RESULTS_DIR'] / 'blast'
PATH_DICT['ORF_OUTPUT_DIR']       = PATH_DICT['RESULTS_DIR'] / 'orf'
PATH_DICT['HMM_OUTPUT_DIR']       = PATH_DICT['RESULTS_DIR'] / 'hmmer'
PATH_DICT['RPSBLAST_SPECIES_DIR'] = PATH_DICT['RESULTS_DIR'] / 'rpsblast'
PATH_DICT['DOMAINS_SPECIES_DIR']  = PATH_DICT['RESULTS_DIR'] / 'domains'

# === ASN Subdirectories ===
PATH_DICT['ASN_TBLASTN_DIR'] = PATH_DICT['ASN_ROOT_OUTPUT_DIR'] / 'tblastn'
PATH_DICT['ASN_RPSBLAST_DIR'] = PATH_DICT['ASN_ROOT_OUTPUT_DIR'] / 'rpsblast'
PATH_DICT['ASN_RPS_PROTEIN_DIR'] = PATH_DICT['ASN_ROOT_OUTPUT_DIR'] / 'protein'

# Directory generation
for path in PATH_DICT.values():
    if isinstance(path, Path):
        path.mkdir(parents=True, exist_ok=True)

# Query configuration
QUERY_FORMAT: str = config['query'].get('format', 'fa')
QUERY_ACC: str = config['query'].get('accession')
QUERY_FILE: Path = PATH_DICT['FASTA_DIR'] / f'{QUERY_ACC}.{QUERY_FORMAT.lower()}'

# Execution and requests
NUM_CORES: int = config['execution'].get('num_cores', 1)
RANDOM_ID_LENGTH: int = config['execution'].get('random_id_length', 6)
USE_SPECIES_DICT: bool = config['execution'].get('use_species_dict', True)
RETRIEVAL_TIME_LAG: float = config['execution'].get('retrieval_time_lag', 0.3)
MAX_RETRIEVAL_ATTEMPTS: int = config['execution'].get('max_retrieval_attempts', 3)
ENTREZ_EMAIL: str = config['execution'].get('entrez_email', '')
NCBI_API_TOKEN: str = config['execution'].get('ncbi_api_token', '')

# === ORF Configuration ===
ORF_ENABLED: bool = config.get('orf', {}).get('enabled', False)
MIN_ORF_LENGTH: int = config.get('orf', {}).get('min_orf_length', 200)

# === HMMER Configuration ===
HMMER_ENABLED: bool       = config.get('hmmer', {}).get('enabled', False)
HMMER_PROFILES: list      = config.get('hmmer', {}).get('profiles', [])
HMMER_USE_GA: bool        = config.get('hmmer', {}).get('use_gathering_threshold', False)
HMMER_EVALUE: float       = config.get('hmmer', {}).get('evalue', 1e-5)
HMMER_DOM_EVALUE: float   = config.get('hmmer', {}).get('dom_evalue', 1e-3)
HMMER_MIN_SCORE           = config.get('hmmer', {}).get('min_score', None)
HMMER_MIN_COVERAGE: float = config.get('hmmer', {}).get('min_coverage', 0.5)
HMMER_MIN_ALN_LEN: int    = config.get('hmmer', {}).get('min_alignment_length', 50)
HMMER_MAX_SENS: bool      = config.get('hmmer', {}).get('max_sensitivity', False)
HMMER_BIAS_FILTER: bool   = config.get('hmmer', {}).get('bias_filter', True)
HMMER_SEED: int           = config.get('hmmer', {}).get('seed', 67)

# === External Program Names ===
_programs = config.get('programs', {})
TBLASTN_CMD: str          = _programs.get('tblastn', 'tblastn')
RPSBLAST_CMD: str         = _programs.get('rpsblast', 'rpsblast')
BLAST_FORMATTER_CMD: str  = _programs.get('blast_formatter', 'blast_formatter')
MAKEBLASTDB_CMD: str      = _programs.get('makeblastdb', 'makeblastdb')
RPSBPROC_CMD: str         = _programs.get('rpsbproc', 'rpsbproc')

# Display
DISPLAY_SNAKEMAKE_INFO: bool = config['display'].get('display_snakemake_info', False)
DISPLAY_REQUESTS_WARNING: bool = config['display'].get('display_requests_warning', False)
DISPLAY_OPERATION_INFO: bool = config['display'].get('display_operation_info', False)

# Genomes
SPECIES_DICT: Dict[str, Any] = config.get('species', {})

if not USE_SPECIES_DICT:
    SPECIES: List[str] = [f.stem for f in PATH_DICT['SPECIES_DB'].glob('*.fa')]
else:
    SPECIES: List[str] = list(SPECIES_DICT.keys())
