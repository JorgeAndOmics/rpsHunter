# Output Files

This document describes every output file produced by the rpsHunter pipeline, including column schemas, file formats, and directory layout.

## Directory Structure

All outputs are written under the `results/` directory:

```
results/
├── blast/              Per-species BLAST parquets
├── orf/                Per-species ORF-enriched parquets (if orf.enabled)
├── hmmer/              Per-species HMMER-enriched parquets (if hmmer.enabled)
├── fastas/             Per-species filtered FASTAs + audit parquets
├── rpsblast/           Per-species RPSBLAST parquets
├── rpsbproc/           Per-species rpsbproc text output
├── domains/            Per-species domain parquets
├── tables/
│   ├── selected/       Per-species selected (filtered) parquets
│   ├── aggregate.*     Combined enriched BLAST data (parquet + CSV)
│   ├── domains.*       Combined domain annotations (parquet + CSV)
│   ├── rpsblast.*      Combined RPSBLAST results (parquet + CSV)
│   └── contingency_table.csv
├── plots/
│   ├── tile_plot.png
│   └── scatter_3D_plot.html
├── ranges/             Per-species GFF3 files
└── asn/
    ├── tblastn/        tBLASTn ASN binary archives
    └── rpsblast/       RPSBLAST ASN binary archives
```

---

## Per-Species Intermediate Files

### blast/{species}.parquet

Raw tBLASTn results for a single species. One row per alignment hit.

| Column | Type | Description |
|--------|------|-------------|
| Query ID | str | Query protein accession |
| Subject ID | str | Genome sequence identifier (chromosome/scaffold) |
| Pct Identity | float | Percentage sequence identity |
| Alignment Length | int | Length of the alignment in residues |
| Mismatches | str | Number of mismatches |
| Gap Openings | str | Number of gap openings |
| Q. Start | str | Query sequence start position |
| Q. End | str | Query sequence end position |
| S. Start | int | Subject sequence start position |
| S. End | int | Subject sequence end position |
| E-value | float | Expect value |
| Bit Score | float | Bit score of the alignment |
| Subject Sequence | str | Aligned subject nucleotide sequence |
| Species | str | Species name (from config) |
| Tag | str | 6-character random alphanumeric identifier for traceability |

**Produced by:** `blast_species` rule (`blast.py`)

---

### orf/{species}.parquet

BLAST parquet enriched with ORF detection results. Contains all 15 BLAST columns plus 3 ORF columns. Only produced when `orf.enabled: true` in config.

| Column | Type | Description |
|--------|------|-------------|
| *(all 15 BLAST columns)* | | See `blast/{species}.parquet` above |
| ORF_Filtered | bool | `True` if the sequence contains an ORF passing the configured thresholds |
| ORF_Length | int | Length of the detected ORF in nucleotides |
| ORF_Sequence | str | Nucleotide sequence of the detected ORF |

**Produced by:** `orf_species` rule (`orf_analyser.R`)

---

### hmmer/{species}.parquet

Previous enrichment step's parquet further enriched with HMMER domain search results. Contains all upstream columns plus 5 HMM columns. Only produced when `hmmer.enabled: true` in config.

| Column | Type | Description |
|--------|------|-------------|
| *(all upstream columns)* | | From `blast/` or `orf/` depending on config |
| HMM_Filtered | bool | `True` if the sequence matched any Pfam domain passing quality filters (scalar) |
| HMM_Domain | list[str] | Pfam domain names for all quality-passing hits, sorted by E-value ascending |
| HMM_Evalue | list[float] | Per-domain E-values, same order as HMM_Domain |
| HMM_Score | list[float] | Per-domain bit scores, same order as HMM_Domain |
| HMM_Coverage | list[float] | Fraction of HMM model covered by each hit, same order as HMM_Domain |

**Produced by:** `hmmer_species` rule (`hmmer.py`)

> **Working with list-valued HMM columns in pandas**
>
> The columns `HMM_Domain`, `HMM_Evalue`, `HMM_Score`, and `HMM_Coverage` are stored as
> Python lists (serialized as numpy arrays in Parquet). When reading the file with
> `pd.read_parquet()`, each cell contains a `numpy.ndarray`. Examples:
>
> ```python
> import pandas as pd
>
> df = pd.read_parquet("results/hmmer/Desmodus_rotundus.parquet")
>
> # Access the first domain name for a given row
> df["HMM_Domain"].iloc[0][0]        # e.g. 'zf-C2H2'
>
> # Count how many domains each sequence matched
> df["HMM_Domain"].apply(len)
>
> # Explode to one row per domain hit (useful for aggregation)
> exploded = df.explode(["HMM_Domain", "HMM_Evalue", "HMM_Score", "HMM_Coverage"])
>
> # Filter to rows with a specific domain anywhere in the list
> mask = df["HMM_Domain"].apply(lambda arr: "KRAB" in arr)
> df[mask]
> ```
>
> `HMM_Filtered` remains a scalar `bool` and can be used directly in boolean indexing.

---

### fastas/{species}.fa

Filtered FASTA file containing sequences that passed all active gates (Quality, ORF, HMMER). These are the input sequences for RPSBLAST domain detection.

**Header format:**

```
>{Subject_ID}:{S.Start}-{S.End}|tag:{Tag}
```

Example:

```
>NC_065643.1:40819974-40820076|tag:x7Km2p
MKTLPRGDTSYRGTWTGEAADLGNPRWRSRGRGLALPGACDHRFHRRGPD...
```

**Produced by:** `blast_parser_species` rule (`blast_parser.py`)

---

### fastas/{species}.parquet

Audit trail parquet containing ALL rows from the enriched BLAST input (not just those that passed filtering). Includes boolean flag columns documenting which gates each row passed or failed.

| Column | Type | Description |
|--------|------|-------------|
| *(all upstream columns)* | | From the ENRICHED_BLAST_DIR source |
| Quality_Pass | bool | `True` if the row passed BLAST quality thresholds (identity, length, bitscore, e-value) |
| ORF_Pass | bool or NaN | `True` if ORF gate passed; `NaN` if ORF enrichment was disabled |
| HMM_Pass | bool or NaN | `True` if HMMER gate passed; `NaN` if HMMER enrichment was disabled |
| Selected | bool | `True` only for rows that passed all active gates and were exported to the FASTA |

**Produced by:** `blast_parser_species` rule (`blast_parser.py`)

---

### tables/selected/{species}.parquet

Rows that passed all active filter gates, with the flag columns (Quality_Pass, ORF_Pass, HMM_Pass, Selected) dropped. Schema matches the enriched BLAST parquet from the last enabled enrichment step.

**Produced by:** `blast_parser_species` rule (`blast_parser.py`)

---

### rpsblast/{species}.parquet

RPSBLAST results from searching filtered sequences against the NCBI Conserved Domain Database (or a configured CDD subset). One row per domain hit.

| Column | Type | Description |
|--------|------|-------------|
| Query ID | str | Sequence identifier (from FASTA header) |
| Subject ID | str | CDD domain accession |
| Pct Identity | float | Percentage sequence identity to the domain model |
| Alignment Length | int | Length of the alignment |
| Mismatches | str | Number of mismatches |
| Gap Openings | str | Number of gap openings |
| Q. Start | str | Query start position |
| Q. End | str | Query end position |
| S. Start | int | Subject (domain model) start position |
| S. End | int | Subject (domain model) end position |
| E-value | float | Expect value |
| Bit Score | float | Bit score |
| Subject Title | str | CDD domain title/description |
| Species | str | Species name |
| Tag | str | 6-character tag extracted from the FASTA header |

**Produced by:** `rpsblast_species` rule (`rpsblast.py`)

---

### rpsbproc/{species}.txt

Plain text output from rpsbproc post-processing. Contains structured blocks with session metadata, query information, and domain annotations.

**Block structure:**

```
SESSION ...
QUERY ...
DOMAINS
  <tab-delimited domain rows>
ENDDOMAINS
ENDQUERY
...
```

**Produced by:** `rpsbproc_species` rule (`rpsbproc.py`)

---

### domains/{species}.parquet

Parsed domain annotations from rpsbproc output. One row per domain annotation on a query sequence.

| Column | Type | Description |
|--------|------|-------------|
| Species | str | Species name |
| Session_ordinal | str | rpsbproc session index |
| Program | str | Program used (e.g. `rpsblast`) |
| Version | str | Program version |
| Database | str | Database searched |
| Score_matrix | str | Scoring matrix used |
| Evalue_threshold | str | E-value threshold applied |
| Query_ID | str | Query sequence identifier |
| Seq_type | str | Sequence type |
| Seq_length | str | Sequence length |
| Definition | str | Full definition line from the FASTA header |
| Chromosome | str | Chromosome/scaffold extracted from Query_ID |
| Start | str | Genomic start coordinate |
| End | str | Genomic end coordinate |
| Hit_type | str | Type of domain hit (specific, superfamily, etc.) |
| PSSM_ID | str | Position-Specific Scoring Matrix identifier |
| From | str | Domain alignment start on the query |
| To | str | Domain alignment end on the query |
| Evalue | str | Domain E-value |
| Bitscore | str | Domain bit score |
| Accession | str | CDD domain accession |
| Domain | str | Domain short name (e.g. `zf-C2H2`, `COG5048`, `KRAB_A-box`) |
| Incomplete | str | Completeness flag: `N` (N-terminal truncated), `C` (C-terminal truncated), `NC` (both), or empty (complete) |
| Superfamily_PSSM_ID | str | PSSM ID of the parent superfamily |
| Tag | str | 6-character alphanumeric tag extracted from Definition via regex `\|tag:(\w+)` |

**Produced by:** `rpsbproc_parser_species` rule (`rpsbproc_parser.py`)

---

## Aggregated Tables

All aggregated tables are written in dual format (Parquet + CSV) for both programmatic and human-readable access.

### tables/aggregate.parquet / aggregate.csv

Concatenation of all per-species parquets from the ENRICHED_BLAST_DIR (the last enabled enrichment step). Schema matches the source parquets --- either `blast/`, `orf/`, or `hmmer/` depending on which enrichment steps are enabled.

**Produced by:** `aggregate_blast` rule (`aggregate.py`)

---

### tables/domains.parquet / domains.csv

Concatenation of all `domains/{species}.parquet` files. Contains the full 25-column domain schema described above.

**Produced by:** `aggregate_domains` rule (`aggregate.py`)

---

### tables/rpsblast.parquet / rpsblast.csv

Concatenation of all `rpsblast/{species}.parquet` files. Contains the 15-column RPSBLAST schema described above.

**Produced by:** `aggregate_domains` rule (`aggregate.py`)

---

### tables/contingency_table.csv

Cross-tabulation of domain annotations by species, domain name, and completeness status. Uses standard comma delimiters.

| Column | Type | Description |
|--------|------|-------------|
| Species | str | Species name |
| Domain | str | Domain short name |
| Incomplete | str | Completeness category (empty, `N`, `C`, `NC`) |
| n | int | Count of annotations matching this combination |

**Produced by:** `contingency_sorter` rule (`contingency_sorter.R`)

---

## Visualizations

### plots/tile_plot.png

Domain completeness heatmap rendered as a static PNG image. Species are arranged on the y-axis, domains on the x-axis, and tile fill color represents the completeness category (complete, N-terminal truncated, C-terminal truncated, or both).

**Produced by:** `completeness_detector` rule (`completeness_detector.R`)

---

### plots/scatter_3D_plot.html

Interactive Plotly 3D scatter plot saved as a self-contained HTML file. Axes: x = chromosome position, y = domain, z = bitscore. Points are colored by species. The file can be opened in any web browser for interactive rotation, zoom, and hover inspection.

**Produced by:** `contingency_sorter` rule (`contingency_sorter.R`)

---

## Genomic Ranges

### ranges/{species}.gff3

GFF3-format files containing genomic coordinates of annotated domains after range reduction (overlapping annotations merged using GenomicRanges). One file per species.

| Field | Description |
|-------|-------------|
| seqid | Chromosome/scaffold identifier |
| source | `rpsHunter` |
| type | Feature type (varies by domain) |
| start | Genomic start coordinate |
| end | Genomic end coordinate |
| score | Domain bitscore |
| strand | Strand orientation |
| phase | `.` (not applicable) |
| attributes | GFF3 attribute string with domain metadata |

**Produced by:** `contingency_sorter` rule (`contingency_sorter.R`)

---

## ASN Binary Archives

### asn/tblastn/{species}.asn

BLAST Archive (ASN.1) binary output from tBLASTn. Retained for potential re-formatting with `blast_formatter`.

### asn/rpsblast/{species}.asn

BLAST Archive (ASN.1) binary output from RPSBLAST. Required as input for rpsbproc post-processing.

---

## Empty-Species Handling

When a species produces zero results at any stage, the pipeline creates valid but empty output files to satisfy Snakemake's output contract:

- **Empty parquets:** Written with the correct column schema and zero rows. Downstream `pd.concat` in `aggregate.py` handles these transparently.
- **Empty FASTAs:** Zero-byte files. RPSBLAST detects the empty input and writes an empty-schema parquet plus a touched ASN file.
- **Empty rpsbproc input:** rpsbproc is skipped (the output `.txt` is touched). The parser produces a zero-row parquet.

This ensures that a single species failure does not halt the entire pipeline.
