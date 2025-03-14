import os

E_VALUE_THRESHOLD = 0.01
PERC_IDENTITY_THRESHOLD = 60
SEQ_LENGTH_THRESHOLD = 0
BITSCORE_THRESHOLD = 50

ROOT_DB = os.path.join('/mnt', 'v', 'databases', 'local')
SPECIES_DB = os.path.join(ROOT_DB, 'blast_dbs', 'species')
RPS_DB = os.path.join(ROOT_DB, 'cdd_dbs', 'Cdd')  # Includes the name of the database
RPSBPROC_DB = os.path.join(ROOT_DB, 'cdd_dbs', 'rpsbproc_annot1')

ROOT_DIR = os.path.join('..')
LOG_DIR = os.path.join(ROOT_DIR, 'logs')
TMP_DIR = os.path.join(ROOT_DIR, 'data', 'tmp')
TABLE_OUTPUT_DIR = os.path.join(ROOT_DIR, 'results', 'output', 'tables')
FASTA_OUTPUT_DIR = os.path.join(ROOT_DIR, 'results', 'output', 'fastas')
ASN_ROOT_OUTPUT_DIR = os.path.join(ROOT_DIR, 'results', 'output', 'asn')
ASN_TBLASTN_DIR = os.path.join(ASN_ROOT_OUTPUT_DIR, 'tblastn')
ASN_RPSBLAST_DIR = os.path.join(ASN_ROOT_OUTPUT_DIR, 'rpsblast')
ASN_RPS_PROTEIN_DIR = os.path.join(ASN_ROOT_OUTPUT_DIR, 'protein')
RPSBPROC_OUTPUT_DIR = os.path.join(ROOT_DIR, 'results', 'output', 'rpsbproc')

XML_OUTPUT_DIR = os.path.join(ROOT_DIR, 'results', 'output', 'xml')
INPUT_FASTA_DIR = os.path.join(ROOT_DIR, 'data', 'input', 'fastas')

QUERY_FORMAT = 'fasta'
QUERY_SEQ = os.path.join(INPUT_FASTA_DIR, f'ADK09900.{QUERY_FORMAT.lower()}')

# Genomes
SPECIES: list = ['Desmodus_rotundus',
                 'Miniopterus_schreibersii',
                 'Tadarida_brasiliensis',
                 'Antrozous_pallidus',
                 'Molossus_molossus',
                 'Artibeus_lituratus',
                 'Eptesicus_fuscus',
                 'Myotis_myotis',
                 'Eptesicus_nilssonii',
                 'Pipistrellus_kuhlii',
                 'Rhinolophus_ferrumequinum',
                 'Saccopteryx_bilineata',
                 'Vespertilio_murinus',
                 'Plecotus_auritus',
                 'Rhinolophus_hipposideros',
                 'Phyllostomus_discolor',
                 'Myotis_daubentonii',
                 'Myotis_mystacinus',
                 'Corynorhinus_townsendii',
                 'Hipposideros_larvatus',
                 'Rhynchonycteris_naso',
                 'Saccopteryx_leptura',
                 'Molossus_alvarezi',
                 'Glossophaga_mutica',
                 'Molossus_nigricans'
                 ]