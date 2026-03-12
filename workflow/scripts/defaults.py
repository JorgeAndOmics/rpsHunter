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
from typing import Dict, Any, List, Union

import yaml

# Configuration file path
CONFIG_FILE: Path = Path(__file__).parents[2] / 'data' / 'config' / 'config.yaml'

with open(CONFIG_FILE, 'r') as f:
    config: Dict[str, Any] = yaml.safe_load(f)

# BLAST
E_VALUE_THRESHOLD: float = config['blast'].get('e_value', 0.01)
PERC_IDENTITY_THRESHOLD: int = config['blast'].get('perc_identity', 60)
SEQ_LENGTH_THRESHOLD: int = config['blast'].get('seq_length', 50)
BITSCORE_THRESHOLD: int = config['blast'].get('bitscore', 70)
ACCESSION_ID_REGEX: str = r'[A-Z]{2,}_?[0-9]+\.[0-9]{1,2}'

# RPSBLAST - Maximum sensitivity settings for short domain detection
RPSBLAST_E_VALUE: float = config.get('rpsblast', {}).get('e_value', 10)
RPSBLAST_COMP_BASED_STATS: int = config.get('rpsblast', {}).get('comp_based_stats', 0)
RPSBLAST_SEG: str = config.get('rpsblast', {}).get('seg', 'no')
RPSBLAST_WINDOW_SIZE: int = config.get('rpsblast', {}).get('window_size', 40)
RPSBLAST_TARGET_DOMAINS: list = config.get('rpsblast', {}).get('target_domains', [])

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
PATH_DICT['SMP_DIR'] = PATH_DICT['ROOT'] / 'cdd_smp'
PATH_DICT['CDD_SUBSET_DB'] = PATH_DICT['ROOT'] / 'cdd_subset'

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
PATH_DICT['SELECTED_OUTPUT_DIR'] = PATH_DICT['TABLE_OUTPUT_DIR'] / 'selected'

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

# =============================================================================
# Query Configuration (multi-query support)
# =============================================================================

def _make_label_safe(label: str) -> str:
    """Convert a display label to a filesystem-safe directory name."""
    return label.strip().replace(' ', '_').replace('/', '_').replace('\\', '_')


# Normalize both config formats to QUERY_DICT: Dict[str, str]  {accession: display_label}
_queries_raw: Union[dict, None] = config.get('queries', None)
_query_raw: Union[dict, None] = config.get('query', None)

if _queries_raw and isinstance(_queries_raw, dict):
    # New multi-query format: queries: {'NP_064612.2': 'human PRDM9', ...}
    QUERY_DICT: Dict[str, str] = _queries_raw
    MULTI_QUERY: bool = len(QUERY_DICT) > 1
elif _query_raw and isinstance(_query_raw, dict):
    # Legacy single-query format: query: {format: 'fa', accession: 'NP_659058.3'}
    _acc = _query_raw.get('accession', '')
    QUERY_DICT: Dict[str, str] = {_acc: _acc}
    MULTI_QUERY: bool = False
else:
    QUERY_DICT: Dict[str, str] = {}
    MULTI_QUERY: bool = False

# Ordered lists derived from QUERY_DICT
QUERIES: List[str] = list(QUERY_DICT.keys())                           # accessions
QUERY_LABELS: Dict[str, str] = {acc: _make_label_safe(label) for acc, label in QUERY_DICT.items()}  # acc → safe label
QUERY_LABEL_LIST: List[str] = list(QUERY_LABELS.values())              # safe labels in order

# Reverse lookup: safe label → accession
LABEL_TO_ACC: Dict[str, str] = {v: k for k, v in QUERY_LABELS.items()}

# Legacy singletons (backward-compatible — used by validator.py and any code not yet multi-query)
QUERY_FORMAT: str = config.get('query', {}).get('format', 'fa')
QUERY_ACC: str = QUERIES[0] if QUERIES else ''
QUERY_FILE: Path = PATH_DICT['FASTA_DIR'] / f'{QUERY_ACC}.{QUERY_FORMAT.lower()}' if QUERY_ACC else PATH_DICT['FASTA_DIR']

# =============================================================================
# Per-Query Directory Generation
# =============================================================================
# Directories that hold per-query outputs need {query_label} subdirectories.
# These must be created after QUERY_LABEL_LIST is available.

_PER_QUERY_PARENTS = [
    PATH_DICT['BLAST_SPECIES_DIR'],       # results/blast/{ql}/
    PATH_DICT['ORF_OUTPUT_DIR'],          # results/orf/{ql}/
    PATH_DICT['HMM_OUTPUT_DIR'],          # results/hmmer/{ql}/
    PATH_DICT['FASTA_OUTPUT_DIR'],        # results/fastas/{ql}/
    PATH_DICT['RPSBLAST_SPECIES_DIR'],    # results/rpsblast/{ql}/
    PATH_DICT['RPSBPROC_OUTPUT_DIR'],     # results/rpsbproc/{ql}/
    PATH_DICT['DOMAINS_SPECIES_DIR'],     # results/domains/{ql}/
    PATH_DICT['ASN_TBLASTN_DIR'],         # results/asn/tblastn/{ql}/
    PATH_DICT['ASN_RPSBLAST_DIR'],        # results/asn/rpsblast/{ql}/
    PATH_DICT['TABLE_OUTPUT_DIR'],        # results/tables/{ql}/ (aggregates + contingency)
    PATH_DICT['PLOT_DIR'],                # results/plots/{ql}/ (tile + 3D plots)
    PATH_DICT['RANGE_OUTPUT_DIR'],        # results/ranges/{ql}/ (GFF3)
]

# Also: results/tables/{ql}/selected/
_PER_QUERY_NESTED = [
    PATH_DICT['TABLE_OUTPUT_DIR'] / '{ql}' / 'selected',
]

for _ql in QUERY_LABEL_LIST:
    for _parent in _PER_QUERY_PARENTS:
        (_parent / _ql).mkdir(parents=True, exist_ok=True)
    for _nested_tmpl in _PER_QUERY_NESTED:
        Path(str(_nested_tmpl).replace('{ql}', _ql)).mkdir(parents=True, exist_ok=True)

# =============================================================================
# Extended Visualization Suite — Output File Paths
# =============================================================================
# All 7 extended plots write into results/plots/ (already created above).
# Defining the full paths here keeps defaults.py as the single source of truth
# for all path definitions, consistent with PATH_DICT.
EXT_CONC_HEATMAP_PATH: Path = PATH_DICT['PLOT_DIR'] / 'concordance_heatmap.png'
EXT_EVID_QUALITY_PATH: Path = PATH_DICT['PLOT_DIR'] / 'evidence_quality.png'
EXT_COMP_BARS_PATH:    Path = PATH_DICT['PLOT_DIR'] / 'completeness_bars.png'
EXT_SEQ_COMPLEX_PATH:  Path = PATH_DICT['PLOT_DIR'] / 'sequence_complexity.png'
EXT_CHR_DENSITY_PATH:  Path = PATH_DICT['PLOT_DIR'] / 'chromosomal_density.png'
EXT_CROSS_QUERY_PATH:  Path = PATH_DICT['PLOT_DIR'] / 'cross_query_comparison.png'
EXT_HIT_TYPE_PATH:     Path = PATH_DICT['PLOT_DIR'] / 'hit_type_distribution.png'


def query_file_path(accession: str) -> Path:
    """Return the FASTA file path for a given query accession."""
    fmt = config.get('query', {}).get('format', 'fa').lower()
    return PATH_DICT['FASTA_DIR'] / f'{accession}.{fmt}'

# Execution and requests
NUM_CORES: int = config['execution'].get('num_cores', 1)
RANDOM_ID_LENGTH: int = config['execution'].get('random_id_length', 6)
USE_SPECIES_DICT: bool = config['execution'].get('use_species_list', True)
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
MAKEPROFILEDB_CMD: str    = _programs.get('makeprofiledb', 'makeprofiledb')

# Display
DISPLAY_SNAKEMAKE_INFO: bool = config['display'].get('display_snakemake_info', False)
DISPLAY_REQUESTS_WARNING: bool = config['display'].get('display_requests_warning', False)
DISPLAY_OPERATION_INFO: bool = config['display'].get('display_operation_info', False)

# Genomes
_species_raw: Union[dict, list, None] = config.get('species', {})

if isinstance(_species_raw, dict):
    SPECIES_DICT: Dict[str, str] = _species_raw
elif isinstance(_species_raw, list):
    SPECIES_DICT: Dict[str, str] = {name: name for name in _species_raw}
else:
    SPECIES_DICT: Dict[str, str] = {}

if not USE_SPECIES_DICT:
    SPECIES: List[str] = [f.stem for f in PATH_DICT['SPECIES_DB'].glob('*.fa')]
else:
    SPECIES: List[str] = list(SPECIES_DICT.keys())


# =============================================================================
# CDD Subset Resolution
# =============================================================================

def resolve_cdd_targets(target_domains: List[str], cddid_path: Path) -> List[str]:
    """Resolve target domain names to CDD accessions via cddid.tbl.

    Matching rules (mirror HMMER prefix convention):
      - Exact match on ShortName column
      - Prefix + underscore: 'KRAB' matches 'KRAB_A-box', 'zf-C2H2' matches 'zf-C2H2_2', etc.

    Parameters
    ----------
    target_domains : list of str
        Domain ShortNames from config (e.g. ['KRAB', 'SET', 'zf-C2H2']).
    cddid_path : Path
        Path to cddid.tbl (tab-delimited: PSSM_ID, Accession, ShortName, Description, Length).

    Returns
    -------
    list of str
        Matched CDD accessions (e.g. ['pfam01352', 'cd07765', ...]).

    Raises
    ------
    RuntimeError
        If any target domain has zero matches in cddid.tbl.
    """
    import re

    # Build regex patterns: exact match OR prefix with '_'
    patterns = []
    for name in target_domains:
        escaped = re.escape(name)
        patterns.append(f'^{escaped}$|^{escaped}_')
    combined = re.compile('|'.join(f'(?:{p})' for p in patterns))

    matched_accessions = []
    matched_names = set()

    with open(cddid_path, 'r') as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            accession, short_name = parts[1], parts[2]
            if combined.search(short_name):
                matched_accessions.append(accession)
                matched_names.add(short_name)

    # Verify every requested domain had at least one match
    unmatched = []
    for name in target_domains:
        escaped = re.escape(name)
        pat = re.compile(f'^{escaped}$|^{escaped}_')
        if not any(pat.search(n) for n in matched_names):
            unmatched.append(name)
    if unmatched:
        raise RuntimeError(
            f'CDD target domain(s) not found in {cddid_path}: {unmatched}'
        )

    return matched_accessions
