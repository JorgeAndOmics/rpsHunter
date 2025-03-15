import os

# Constants
E_VALUE_THRESHOLD = 0.01
PERC_IDENTITY_THRESHOLD = 60
SEQ_LENGTH_THRESHOLD = 0
BITSCORE_THRESHOLD = 50


# Directories
ROOT_DIR = os.path.join('..')
DATA_DIR = os.path.abspath(os.path.join(ROOT_DIR, 'data'))
RESULTS_DIR = os.path.abspath(os.path.join(ROOT_DIR, 'results'))
DATA_INPUT_DIR = os.path.join(DATA_DIR, 'input')
DATA_FASTA_DIR = os.path.join(DATA_INPUT_DIR, 'fastas')
CONFIG_DIR = os.path.join(DATA_DIR, 'config')
LOG_DIR = os.path.abspath(os.path.join(ROOT_DIR, 'logs'))
WORKFLOW_DIR = os.path.abspath(os.path.join(ROOT_DIR, 'workflow'))
TMP_DIR = os.path.join(DATA_DIR, 'tmp')
PLOT_DIR = os.path.join(RESULTS_DIR, 'plots')
TABLE_OUTPUT_DIR = os.path.join(RESULTS_DIR, 'tables')
FASTA_OUTPUT_DIR = os.path.join(RESULTS_DIR, 'fastas')
ASN_ROOT_OUTPUT_DIR = os.path.join(RESULTS_DIR, 'asn')
RPSBPROC_OUTPUT_DIR = os.path.join(RESULTS_DIR, 'rpsbproc')
XML_OUTPUT_DIR = os.path.join(RESULTS_DIR, 'xml')
ASN_TBLASTN_DIR = os.path.join(ASN_ROOT_OUTPUT_DIR, 'tblastn')
ASN_RPSBLAST_DIR = os.path.join(ASN_ROOT_OUTPUT_DIR, 'rpsblast')
ASN_RPS_PROTEIN_DIR = os.path.join(ASN_ROOT_OUTPUT_DIR, 'protein')
INPUT_FASTA_DIR = os.path.join(DATA_INPUT_DIR, 'fastas')


# Database retrieval
ROOT_CONFIG_FILE = os.path.abspath(os.path.join(CONFIG_DIR, 'root_db_folder.txt'))
with open(ROOT_CONFIG_FILE, 'r') as f:
    ROOT = os.path.abspath(f.readline().strip())


# Databases
ROOT_DB = os.path.abspath(os.path.join(ROOT, 'local'))
SPECIES_DB = os.path.join(ROOT_DB, 'blast_dbs', 'species')
RPS_DB = os.path.join(ROOT_DB, 'cdd_dbs', 'Cdd')  # Includes the name of the database
RPSBPROC_DB = os.path.join(ROOT_DB, 'cdd_dbs', 'rpsbproc_annot1')


# Directory generation
os.makedirs(ROOT_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(DATA_INPUT_DIR, exist_ok=True)
os.makedirs(DATA_FASTA_DIR, exist_ok=True)
os.makedirs(CONFIG_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(WORKFLOW_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(TABLE_OUTPUT_DIR, exist_ok=True)
os.makedirs(FASTA_OUTPUT_DIR, exist_ok=True)
os.makedirs(ASN_ROOT_OUTPUT_DIR, exist_ok=True)
os.makedirs(RPSBPROC_OUTPUT_DIR, exist_ok=True)
os.makedirs(XML_OUTPUT_DIR, exist_ok=True)
os.makedirs(ASN_TBLASTN_DIR, exist_ok=True)
os.makedirs(ASN_RPSBLAST_DIR, exist_ok=True)
os.makedirs(ASN_RPS_PROTEIN_DIR, exist_ok=True)
os.makedirs(INPUT_FASTA_DIR, exist_ok=True)
os.makedirs(ROOT_DB, exist_ok=True)
os.makedirs(SPECIES_DB, exist_ok=True)
os.makedirs(RPSBPROC_DB, exist_ok=True)


# Query
QUERY_FORMAT = 'fa'
QUERY_INPUT_FILE = os.path.abspath(os.path.join(CONFIG_DIR, 'query.txt'))
with open(QUERY_INPUT_FILE, 'r') as f:
    QUERY_ACC = f.readline().strip()
QUERY_FILE = os.path.join(INPUT_FASTA_DIR, f'{QUERY_ACC}.{QUERY_FORMAT.lower()}')


# Execution and requests
MAX_THREADPOOL_WORKERS = None
USE_SPECIES_LIST: bool = True
ENTREZ_EMAIL_FILE = os.path.abspath(os.path.join(CONFIG_DIR, 'entrez_email.txt'))
with open(ENTREZ_EMAIL_FILE, 'r') as f:
    ENTREZ_EMAIL = os.path.abspath(f.readline().strip())
NCBI_API_FILE = os.path.abspath(os.path.join(CONFIG_DIR, 'ncbi_api.txt'))
with open(NCBI_API_FILE, 'r') as f:
    NCBI_API_TOKEN = os.path.abspath(f.readline().strip())


# Genomes
SPECIES_FILE = os.path.abspath(os.path.join(CONFIG_DIR, 'species.txt'))

if not USE_SPECIES_LIST:
    SPECIES: list = [(f.split('.fa')[0]).strip() for f in os.listdir(SPECIES_DB) if f.endswith('.fa')]
else:
    SPECIES: list = [line.strip() for line in open(SPECIES_FILE, 'r')]

