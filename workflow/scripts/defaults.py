"""
Defaults Configuration Script
=============================

This script loads configuration settings from a YAML file and sets up various
constants and directory paths used throughout the project.

Modules:
    - `yaml`: For parsing the YAML configuration file.
    - `os`: For handling file and directory paths.
    - `pathlib.Path`: For robust path handling.

Configuration:
    - The configuration file is expected to be located at `../data/config/config.yaml`.
    - Various constants and directory paths are initialized based on the configuration file.

Usage:
    This script is intended to be imported as a module and not run directly.
"""

from pathlib import Path
import yaml
import os
from typing import Dict, Any, List, Optional

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

# Directories
PATH_DICT: Dict[str, Any] = {
    'ROOT': os.path.abspath(config['root'].get('db_root_folder', os.path.join('..')))
}

# === Root Directories ===
PATH_DICT['DATA_DIR'] = os.path.abspath(os.path.join(config['root'].get('data_root_folder', os.path.join(PATH_DICT['ROOT'], 'data'))))
PATH_DICT['RESULTS_DIR'] = os.path.abspath(os.path.join(config['root'].get('results_root_folder', os.path.join(PATH_DICT['ROOT'], 'results'))))
PATH_DICT['LOG_DIR'] = os.path.abspath(os.path.join(config['root'].get('logs_root_folder', os.path.join(PATH_DICT['ROOT'], 'logs'))))
PATH_DICT['WORKFLOW_DIR'] = Path(__file__).parent

# === Data Subdirectories ===
PATH_DICT['CONFIG_DIR'] = os.path.abspath(os.path.join(PATH_DICT['DATA_DIR'], 'config'))
PATH_DICT['INPUT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['DATA_DIR'], 'input'))
PATH_DICT['FASTA_DIR'] = os.path.abspath(os.path.join(PATH_DICT['INPUT_DIR'], 'fastas'))
PATH_DICT['TMP_DIR'] = os.path.abspath(os.path.join(PATH_DICT['DATA_DIR'], 'tmp'))

# === Results Subdirectories ===
PATH_DICT['PLOT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'plots'))
PATH_DICT['TABLE_OUTPUT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'tables'))
PATH_DICT['FASTA_OUTPUT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'fastas'))
PATH_DICT['ASN_ROOT_OUTPUT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'asn'))
PATH_DICT['RPSBPROC_OUTPUT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'rpsbproc'))
PATH_DICT['XML_OUTPUT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'xml'))
PATH_DICT['RANGE_OUTPUT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'ranges'))

# === ASN Subdirectories ===
PATH_DICT['ASN_TBLASTN_DIR'] = os.path.abspath(os.path.join(PATH_DICT['ASN_ROOT_OUTPUT_DIR'], 'tblastn'))
PATH_DICT['ASN_RPSBLAST_DIR'] = os.path.abspath(os.path.join(PATH_DICT['ASN_ROOT_OUTPUT_DIR'], 'rpsblast'))
PATH_DICT['ASN_RPS_PROTEIN_DIR'] = os.path.abspath(os.path.join(PATH_DICT['ASN_ROOT_OUTPUT_DIR'], 'protein'))

# === Database Directories ===
PATH_DICT['ROOT_DB'] = os.path.abspath(os.path.join(PATH_DICT['ROOT'], 'local'))
PATH_DICT['SPECIES_DB'] = os.path.abspath(os.path.join(PATH_DICT['ROOT_DB'], 'blast_dbs', 'species'))
PATH_DICT['RPS_DB'] = os.path.abspath(os.path.join(PATH_DICT['ROOT_DB'], 'cdd_dbs', 'Cdd'))
PATH_DICT['RPSBPROC_DB'] = os.path.abspath(os.path.join(PATH_DICT['ROOT_DB'], 'cdd_dbs', 'rpsbproc_annot1'))

# Directory generation
for value in PATH_DICT.values():
    os.makedirs(value, exist_ok=True)

# Query configuration
QUERY_FORMAT: str = config['query'].get('format', 'fa')
QUERY_ACC: str = config['query'].get('accession')
QUERY_FILE: str = os.path.abspath(os.path.join(PATH_DICT['FASTA_DIR'], f'{QUERY_ACC}.{QUERY_FORMAT.lower()}'))

# Execution and requests
NUM_CORES: int = config['execution'].get('num_cores', 1)
USE_SPECIES_DICT: bool = config['execution'].get('use_species_dict', True)
RETRIEVAL_TIME_LAG: float = config['execution'].get('retrieval_time_lag', 0.3)
MAX_RETRIEVAL_ATTEMPTS: int = config['execution'].get('max_retrieval_attempts', 3)
MAX_THREADPOOL_WORKERS: Optional[int] = config['execution'].get('max_threadpool_workers', None)
ENTREZ_EMAIL: str = config['execution'].get('entrez_email', '')
NCBI_API_TOKEN: str = config['execution'].get('ncbi_api_token', '')

# Display
DISPLAY_SNAKEMAKE_INFO: bool = config['display'].get('display_snakemake_info', False)
DISPLAY_REQUESTS_WARNING: bool = config['display'].get('display_requests_warning', False)
DISPLAY_OPERATION_INFO: bool = config['display'].get('display_operation_info', False)

# Genomes
SPECIES_DICT: Dict[str, Any] = config.get('species', {})

if not USE_SPECIES_DICT:
    SPECIES: List[str] = [(f.split('.fa')[0]).strip() for f in os.listdir(PATH_DICT['SPECIES_DB']) if f.endswith('.fa')]
else:
    SPECIES: List[str] = SPECIES_DICT.keys()
