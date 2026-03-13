# Output Files

This document describes every output file produced by the rpsHunter pipeline, including column schemas, file formats, and directory layout.

## Directory Structure

All outputs are written under the `results/` directory. When using multi-query mode, per-query files are nested under `{query_label}/` subdirectories. Merged (cross-query deduplicated) outputs are written to the top-level directories.

```
results/
├── blast/{query_label}/              Per-query per-species BLAST parquets
├── orf/{query_label}/                Per-query ORF-enriched parquets (if orf.enabled)
├── hmmer/{query_label}/              Per-query HMMER-enriched parquets (if hmmer.enabled)
├── fastas/{query_label}/             Per-query filtered FASTAs + audit parquets
├── rpsblast/{query_label}/           Per-query RPSBLAST parquets
├── rpsbproc/{query_label}/           Per-query rpsbproc text output
├── domains/{query_label}/            Per-query domain parquets
├── tables/
│   ├── blast/
│   │   ├── aggregate.*               Merged enriched BLAST data
│   │   └── {query_label}/
│   │       ├── aggregate.*           Per-query enriched BLAST aggregate
│   │       └── selected/{species}.parquet  Per-query selected parquets
│   ├── rpsblast/
│   │   ├── rpsblast.*                Merged RPSBLAST results
│   │   └── {query_label}/rpsblast.*  Per-query RPSBLAST aggregate
│   ├── domains/
│   │   ├── domains.*                 Merged domain annotations (deduplicated)
│   │   └── {query_label}/domains.*   Per-query domain aggregate
│   ├── contingency/
│   │   ├── contingency_table.*       Merged contingency table
│   │   └── {query_label}/contingency_table.*  Per-query contingency table
│   ├── concordance/
│   │   ├── merged_concordance_domains.*      Merged concordance domain labels
│   │   ├── merged_concordance_sequences.*    Merged concordance sequence summary
│   │   ├── merged_concordance_summary.*      Merged concordance family summary
│   │   └── {query_label}/
│   │       ├── concordance_domains.*         Per-query concordance domain labels
│   │       ├── concordance_sequences.*       Per-query concordance sequence summary
│   │       └── concordance_summary.*         Per-query concordance family summary
│   └── inventory/
│       ├── hit_domain_inventory.*    Merged hit–domain inventory
│       ├── locus_domain_inventory.*  Merged locus–domain inventory (cross-query)
│       └── {query_label}/
│           ├── hit_domain_inventory.*    Per-query hit–domain inventory
│           └── locus_domain_inventory.*  Per-query locus–domain inventory
├── plots/
│   ├── {query_label}/
│   │   ├── tile_plot.png             Per-query tile plot
│   │   └── scatter_3D_plot.html      Per-query 3D scatter plot
│   ├── tile_plot.png                 Merged tile plot
│   ├── scatter_3D_plot.html          Merged 3D scatter plot
│   ├── concordance_heatmap.png       Extended: concordance validation
│   ├── evidence_quality.png          Extended: evidence quality landscape
│   ├── completeness_bars.png         Extended: domain completeness bars
│   ├── sequence_complexity.png       Extended: per-sequence domain complexity
│   ├── chromosomal_density.png       Extended: chromosomal domain density
│   ├── cross_query_comparison.png    Extended: cross-query domain comparison
│   └── hit_type_distribution.png     Extended: hit-type distribution
├── ranges/
│   ├── {query_label}/                Per-query GFF3 files
│   └── {species}.gff3                Merged GFF3 files
└── asn/
    ├── tblastn/{query_label}/        tBLASTn ASN binary archives
    └── rpsblast/{query_label}/       RPSBLAST ASN binary archives
```

`{query_label}` is the filesystem-safe version of the query display label (e.g., `human_PRDM9` for `'human PRDM9'`).

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
| Species | str | Species key (file-friendly, from config key) |
| Species_Name | str | Species display name (from config value, e.g. `Desmodus rotundus`) |
| Query_Accession | str | Query protein accession used for this tBLASTn search |
| Tag | str | 6-character random alphanumeric identifier for traceability |

**Produced by:** `blast_species` rule (`blast.py`)

---

### orf/{species}.parquet

BLAST parquet enriched with ORF detection results. Contains all 16 BLAST columns plus 3 ORF columns. Only produced when `orf.enabled: true` in config.

| Column | Type | Description |
|--------|------|-------------|
| *(all 16 BLAST columns)* | | See `blast/{species}.parquet` above |
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

### tables/blast/selected/{species}.parquet

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
| Species | str | Species key (file-friendly) |
| Species_Name | str | Species display name (e.g. `Desmodus rotundus`) |
| Query_Accession | str | Query protein accession |
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
| Species | str | Species key (file-friendly) |
| Species_Name | str | Species display name (e.g. `Desmodus rotundus`) |
| Query_Accession | str | Query protein accession(s). In merged outputs, comma-separated if multiple queries contributed. |
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

### tables/blast/aggregate.parquet / aggregate.csv

Concatenation of all per-species parquets from the ENRICHED_BLAST_DIR (the last enabled enrichment step). Schema matches the source parquets --- either `blast/`, `orf/`, or `hmmer/` depending on which enrichment steps are enabled.

**Produced by:** `aggregate_blast` rule (`aggregate.py`)

---

### tables/domains/{query_label}/domains.parquet / domains.csv

Per-query concatenation of all `domains/{query_label}/{species}.parquet` files. Contains the full domain schema described above.

**Produced by:** `aggregate_domains` rule (`aggregate.py`)

### tables/domains/domains.parquet / domains.csv

Cross-query merged and deduplicated domain annotations. Overlapping annotations from different queries on the same Species x Chromosome x Domain are merged using GenomicRanges `reduce()`. The `Query_Accession` column lists all contributing accessions (comma-separated).

**Produced by:** `merge_domains` rule (`merge_queries.R`)

---

### tables/rpsblast/rpsblast.parquet / rpsblast.csv

Cross-query merged concatenation of all per-query RPSBLAST aggregate tables. Contains the 16-column RPSBLAST schema described above. Per-query aggregates are at `tables/rpsblast/{query_label}/rpsblast.parquet`.

**Produced by:** `aggregate_rpsblast` (per-query) and `merge_rpsblast` (merged) rules (`aggregate.py`)

---

### tables/contingency/contingency_table.parquet / contingency_table.csv

Cross-tabulation of domain annotations by species, domain name, and completeness status. Written in dual format (Parquet + CSV).

| Column | Type | Description |
|--------|------|-------------|
| Species | str | Species key (file-friendly, e.g. `Desmodus_rotundus`) |
| Species_Name | str | Species display name (e.g. `Desmodus rotundus`) |
| Domain | str | Domain short name |
| Incomplete | str | Completeness category (`Complete`, `N-Truncated`, `C-Truncated`, `Bitruncated`) |
| Count | int | Count of annotations matching this combination |

**Produced by:** `contingency_sorter` rule (`contingency_sorter.R`)

---

## Visualizations

### plots/tile_plot.png

Domain completeness heatmap rendered as a static PNG image. Species display names are arranged on the x-axis (bold italic), domains on the y-axis, and tile fill color represents the completeness category (complete, N-terminal truncated, C-terminal truncated, or both).

**Produced by:** `completeness_detector` rule (`completeness_detector.R`)

---

### plots/scatter_3D_plot.html

Interactive Plotly 3D scatter plot saved as a self-contained HTML file. Axes: x = chromosome position, y = domain, z = bitscore. Points are colored by species. The file can be opened in any web browser for interactive rotation, zoom, and hover inspection.

**Produced by:** `contingency_sorter` rule (`contingency_sorter.R`)

---

## Extended Visualization Suite

Seven additional static PNG plots produced by the `extended_plots` rule (`extended_plots.R`). All plots read from the merged (cross-query deduplicated) output files and require that `--rpsbproc-parser`, `--contingency-parser`, and `--concordance` have been run first.

### plots/concordance_heatmap.png

Tile heatmap showing the fraction of CDD domain annotations that are confirmed by HMMER for each HMMER-checkable domain family × species combination. Fill colour encodes Concordance_Rate (0–1, viridis scale); tiles are annotated with `Confirmed/Total` counts. Faceted by Query_Accession.

**Source:** `results/tables/concordance/merged_concordance_domains.parquet`

---

### plots/evidence_quality.png

Scatter plot of CDD bitscore (y) vs. BLAST bitscore (x, log10 scale) for all HMMER-checkable domain annotations, coloured by concordance label (`confirmed` = teal, `hmmer_unmatched` = orange). Linear regression lines are overlaid per concordance group. Faceted by domain family (free scales).

**Source:** `results/tables/concordance/merged_concordance_domains.parquet`

---

### plots/completeness_bars.png

Stacked proportional bar chart (one bar per domain per species row facet) showing the breakdown of Complete / N-Truncated / C-Truncated / Bitruncated annotations. Bar segments labelled with raw counts where the segment exceeds 5 % of the bar. Domains ordered by total count descending.

**Source:** `results/tables/contingency/contingency_table.parquet`

---

### plots/sequence_complexity.png

Violin + boxplot showing the distribution of `N_CDD_Domains` (number of CDD annotations per BLAST hit sequence) per species, stratified by concordance rate quartile (Q1–Q4). Faceted by Query_Accession.

**Source:** `results/tables/concordance/merged_concordance_sequences.parquet`

---

### plots/chromosomal_density.png

Tile heatmap of log10(domain count + 1) per Species_Name × Chromosome cell. Scaffolds and unplaced contigs are excluded. Cells with more than 50 domains are labelled with the first 4 characters of the dominant domain name. Sequential blue colour scale.

**Source:** `results/tables/domains/domains.parquet`

---

### plots/cross_query_comparison.png

Grouped bar chart showing how many domain annotations per species were detected exclusively by the first query, exclusively by the second query, or by both queries. Requires at least 2 configured queries; a placeholder image is produced for single-query runs. Faceted by domain family.

**Source:** `results/tables/domains/domains.parquet`

---

### plots/hit_type_distribution.png

Stacked bar chart of CDD hit type counts (Specific / Non-specific / Superfamily) per domain per species. Domains ordered by total Specific-hit count descending. Faceted by species (free y-scale, shared x-axis).

**Source:** `results/tables/domains/domains.parquet`

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

## Concordance Tables

### tables/concordance/{query_label}/concordance_domains.parquet

Per-domain concordance labels joining CDD domain annotations to HMMER domain calls via the Tag system.

| Column | Type | Description |
|--------|------|-------------|
| *(all domain columns)* | | From `tables/domains/{query_label}/domains.parquet` |
| Domain_Family | str | CDD domain name normalized to config target via prefix matching |
| HMMER_Checkable | bool | Whether this domain family has a corresponding HMMER profile |
| Blast_Evalue | float | Source BLAST hit E-value (joined via Tag) |
| Blast_Bitscore | float | Source BLAST hit bit score (joined via Tag) |
| Blast_Identity | float | Source BLAST hit percent identity (joined via Tag) |
| Concordance | str | `confirmed` / `hmmer_unmatched` / `hmmer_not_searched` |
| HMM_Best_Evalue | float | Best HMMER E-value for this domain family on this Tag (if confirmed) |

### tables/concordance/{query_label}/concordance_sequences.parquet

Per-sequence (Tag) concordance summary.

| Column | Type | Description |
|--------|------|-------------|
| Tag | str | 6-character sequence identifier |
| Species | str | Species key |
| Species_Name | str | Species display name |
| Query_Accession | str | Query protein accession |
| N_CDD_Domains | int | Total CDD domain annotations on this sequence |
| N_CDD_Families | int | Number of distinct CDD domain families |
| N_HMM_Families | int | Number of distinct HMMER domain families |
| N_Confirmed | int | CDD domains confirmed by HMMER |
| N_Unmatched | int | CDD domains with HMMER profile but no HMMER hit |
| N_Not_Searched | int | CDD domains with no HMMER profile |
| Concordance_Rate | float | N_Confirmed / (N_Confirmed + N_Unmatched), or null |

### tables/concordance/{query_label}/concordance_summary.parquet

Per-domain-family concordance summary with aggregate statistics.

**Produced by:** `pq_concordance` rule (`concordance.py`)

| Column | Type | Description |
|--------|------|-------------|
| Query_Accession | str | Query protein accession |
| Domain_Family | str | CDD domain name normalized to config target via prefix matching |
| Total | int | Total number of annotations for this domain family |
| Confirmed | int | CDD domains confirmed by HMMER |
| Unmatched | int | CDD domains with HMMER profile but no HMMER hit |
| Not_Searched | int | CDD domains with no HMMER profile |
| Concordance_Rate | float | Confirmed / (Confirmed + Unmatched), or null if no searchable domains |
| Median_CDD_Evalue | float | Median E-value across all CDD annotations for this family |
| Median_CDD_Bitscore | float | Median bit score across all CDD annotations for this family |
| Median_HMM_Evalue | float | Median HMMER E-value for confirmed annotations, or null |

---

## Hit–Domain Inventory

### tables/inventory/{query_label}/hit_domain_inventory.parquet / hit_domain_inventory.csv

One row per BLAST hit (Tag), enriched with all CDD domains found within that hit, HMMER domains from the enrichment step, and completeness flags against the expected domain set.

| Column | Type | Description |
|--------|------|-------------|
| Tag | str | 6-character sequence identifier |
| Species | str | Species key (file-friendly) |
| Species_Name | str | Species display name |
| Query_Accession | str | Query protein accession |
| Subject_ID | str | Chromosome/scaffold identifier |
| S_Start | int | Genomic start coordinate of BLAST hit |
| S_End | int | Genomic end coordinate of BLAST hit |
| E_value | float | BLAST E-value |
| Bit_Score | float | BLAST bit score |
| Pct_Identity | float | BLAST percent identity |
| Alignment_Length | int | BLAST alignment length |
| HMMER_Domains_Found | str | Raw HMMER domain names found on this hit (semicolon-separated) |
| HMMER_Profiles_Present | str | Which `hmmer.profiles` config entries are satisfied (semicolon-separated) |
| HMMER_Profiles_Missing | str | Which `hmmer.profiles` entries are not found (semicolon-separated) |
| Complete_HMMER | bool | `True` if all `hmmer.profiles` are covered by HMMER hits |
| CDD_Domains_Found | str | Raw CDD domain names found on this hit (semicolon-separated) |
| N_CDD_Domains | int | Total number of CDD domain annotations |
| N_CDD_Families | int | Number of distinct CDD domain names |
| CDD_Targets_Present | str | Which `rpsblast.target_domains` config entries are satisfied (semicolon-separated) |
| CDD_Targets_Missing | str | Which `rpsblast.target_domains` entries are not found (semicolon-separated) |
| CDD_Domain_Detail | str | JSON array of per-domain detail dicts (`Domain`, `Evalue`, `Bitscore`, `Incomplete`, `Hit_type`, `From`, `To`) |
| Complete_CDD | bool | `True` if all `hmmer.profiles` are covered by CDD domain hits |

**Produced by:** `pq_hit_domain_inventory` rule (`hit_domain_inventory.py`)

### tables/inventory/hit_domain_inventory.parquet / hit_domain_inventory.csv

Concatenation of all per-query hit–domain inventory tables. Same schema as the per-query version.

**Produced by:** `merge_hit_domain_inventory` rule

### tables/inventory/{query_label}/locus_domain_inventory.parquet / locus_domain_inventory.csv

Locus-level aggregation of BLAST hits. Overlapping hits on the same Species x Chromosome are clustered into genomic loci (configurable via `hit_domain_inventory.locus_gap`, default: 50 kb). One row per locus, with the union of all CDD and HMMER domains across all contributing hits.

| Column | Type | Description |
|--------|------|-------------|
| Locus_ID | int | Auto-incrementing locus identifier |
| Species | str | Species key (file-friendly) |
| Species_Name | str | Species display name |
| Subject_ID | str | Chromosome/scaffold |
| Locus_Start | int | Start coordinate (min of all contributing hits) |
| Locus_End | int | End coordinate (max of all contributing hits) |
| Locus_Length | int | Locus span in bp |
| N_Hits | int | Number of BLAST hits in this locus |
| Query_Accessions | str | Contributing query accessions (semicolon-separated) |
| N_Queries | int | Number of distinct queries contributing |
| Tags | str | All Tags in this locus (semicolon-separated) |
| Best_E_value | float | Best (lowest) E-value among contributing hits |
| Best_Bit_Score | float | Best (highest) bit score among contributing hits |
| Best_Pct_Identity | float | Percent identity of the best hit |
| HMMER_Domains_Found | str | Union of raw HMMER domains across all hits (semicolon-separated) |
| HMMER_Profiles_Present | str | Which `hmmer.profiles` entries are covered (semicolon-separated) |
| HMMER_Profiles_Missing | str | Which `hmmer.profiles` entries are not covered (semicolon-separated) |
| Complete_HMMER | bool | `True` if all `hmmer.profiles` are covered by HMMER hits in this locus |
| CDD_Domains_Found | str | Union of raw CDD domains across all hits (semicolon-separated) |
| N_CDD_Domains | int | Total CDD domain annotation count across all hits |
| N_CDD_Families | int | Number of distinct CDD domain names |
| CDD_Targets_Present | str | Which `rpsblast.target_domains` entries are covered (semicolon-separated) |
| CDD_Targets_Missing | str | Which `rpsblast.target_domains` entries are not covered (semicolon-separated) |
| Complete_CDD | bool | `True` if all `hmmer.profiles` are covered by CDD hits in this locus |

**Produced by:** `pq_hit_domain_inventory` rule (`hit_domain_inventory.py`)

### tables/inventory/locus_domain_inventory.parquet / locus_domain_inventory.csv

Cross-query locus–domain inventory. All per-query aggregate tables are concatenated before clustering, so loci are formed from BLAST hits across all queries. This is the most informative table for identifying genomic loci with complete expected domain sets.

**Produced by:** `merge_hit_domain_inventory` rule (`hit_domain_inventory.py`)

---

## Empty-Species Handling

When a species produces zero results at any stage, the pipeline creates valid but empty output files to satisfy Snakemake's output contract:

- **Empty parquets:** Written with the correct column schema and zero rows. Downstream `pd.concat` in `aggregate.py` handles these transparently.
- **Empty FASTAs:** Zero-byte files. RPSBLAST detects the empty input and writes an empty-schema parquet plus a touched ASN file.
- **Empty rpsbproc input:** rpsbproc is skipped (the output `.txt` is touched). The parser produces a zero-row parquet.

This ensures that a single species failure does not halt the entire pipeline.
