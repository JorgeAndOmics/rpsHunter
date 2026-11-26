# Configuration Reference

This document is the complete parameter reference for rpsHunter's `config.yaml` file. Every configurable option is listed with its type, default value, validation constraint, and behavioral description.

## Config File Location

The primary configuration file is located at:

```
data/config/config.yaml
```

This YAML file is read at import time by `workflow/scripts/defaults.py`, which exposes all values as typed Python constants used throughout the pipeline.

## Schema Validation

Configuration is validated against `data/config/schema.yaml` using [Yamale](https://github.com/23andMe/Yamale), a strict YAML schema validator. Validation runs automatically before pipeline execution unless bypassed with the `--skip-validation` CLI flag (see [cli-reference.md](cli-reference.md)).

The schema enforces types, value ranges, regex patterns, and required/optional status for every parameter. If validation fails, the pipeline exits with a descriptive error before any computation begins.

---

## Parameter Reference

### `blast` -- tBLASTn Filtering Thresholds

These thresholds filter tBLASTn hits during the homology search stage. Hits that fail any threshold are excluded from downstream processing. Relaxing these values captures weaker homology signals at the cost of increased noise.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `e_value` | `float` | `0.01` | `num(min=0)` | Maximum E-value for a tBLASTn hit to pass quality filtering. Lower values are more stringent. |
| `perc_identity` | `float` | `60` | `num(min=0, max=100)` | Minimum percent identity between query and subject. Values below 35 introduce significant noise. |
| `seq_length` | `int` | `50` | `int(min=1)` | Minimum alignment length in amino acids. Short domains (e.g., SSXRD at 24aa) require lowering this. |
| `bitscore` | `float` | `70` | `num(min=0)` | Minimum bit-score. Context-dependent; weak but real homologs can score as low as 45. |

### `rpsblast` -- RPSBLAST / CDD Domain Detection

Controls the RPSBLAST search against the NCBI Conserved Domain Database (CDD). The `target_domains` list enables subset database construction, searching only the specified domains rather than the full CDD.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `e_value` | `float` | `1` | `num(min=0)` | Highest E-value kept, on NCBI's scale (the search is told the size of the full CDD, see `databases`). CD-Search's default is 0.01; NCBI treats 0.01 to 1 as borderline and 1 and above as putative false positives. Each call reports its tier in the `Evidence` column: `standard` (at most 0.01), `borderline` (up to 1) or `rescue` (above 1, only present if this value is raised). The same cut-off is given to rpsbproc. See [ADR-006](adr/ADR-006-rpsblast-evalue-scale-and-cutoff.md). |
| `comp_based_stats` | `int` | `0` | `int(min=0, max=4)` | Composition-based statistics mode. `0` disables (more sensitive); `1`-`4` enable progressive correction. See BLAST+ documentation. |
| `seg` | `str` | `'no'` | `str(matches='^(yes\|no)$')` | SEG low-complexity filtering. `'no'` disables, allowing short/repetitive domains through. `'yes'` masks low-complexity regions. |
| `window_size` | `int` | `40` | `int(min=0)`, optional | Window size for multi-hit algorithm. `0` switches to single-hit mode (more sensitive, slower). Omit or set to `40` for default behavior. |

### `databases` -- Reference Database Releases

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `cdd.version` | `str` | yes | CDD release the models come from, for the record (for example `'3.22'`). |
| `cdd.residues` | `int` | yes | Total residues of that full CDD release (`15821058` for 3.22; `blastdbcmd -db Cdd -info` prints it). Passed to RPS-BLAST as `-dbsize`. The search runs against the configured families' models only, and an E-value depends on database size: without this, the same alignment would score about a thousand times better than in NCBI's CD-Search, and its E-value would change whenever a family was added or removed. |

### `domains` -- Domain Families

What to look for inside the hits. Each family is a name of your choosing and the accessions of the
models that count as evidence for it ([ADR-012](adr/ADR-012-domain-families-by-accession.md)).

```yaml
domains:
  KRAB:
    cdd: [pfam01352, smart00349, cd07765]
    pfam: [PF01352]
  zinc_finger:
    cdd: [pfam00096, smart00355, pfam13465, COG5048]
    pfam: [PF00096]
    repeat: true
```

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `cdd` | list of `str` | one of the two | CDD accessions (`pfam01352`, `cd19193`, `smart00349`, `COG5048`, ...). The RPS-BLAST database is built from exactly these models; an accession without a profile in the CDD release stops the run. |
| `pfam` | list of `str` | one of the two | Pfam accessions without version (`PF01352`). Searched with HMMER when `hmmer.enabled` is true. |
| `repeat` | `bool` | no | `true` for a domain that occurs in tandem copies. |

Accessions are used because databases rename models between releases: CDD's `COG5048` changed its
short name, and a config that listed names stopped working. To find accessions, look the protein up
in NCBI's Conserved Domain Database or in InterPro. A model may belong to one family only. Calls
reported at superfamily level are assigned to the family whose models belong to that superfamily.
Every domain table carries a `Family` column beside the model.

### `hmmer` -- HMMER as Independent Evidence

Searches the selected hit sequences with the Pfam profiles of the families under `domains:`
(`pfam:` accessions). HMMER annotates and never filters: every hit that passes the quality gate
(`blast:`) goes to RPS-BLAST whether or not HMMER matches it
([ADR-010](adr/ADR-010-scope-and-hmmer-evidence.md)). Its domains are written to
`tables/domains/hmmer_domains.*` at their own genomic position, each with `Above_GA` (whether it
clears Pfam's curated gathering threshold). CDD calls and HMMER domains both carry
`Method_Agreement`: `both` when the other method reports the same family on the same hit at an
overlapping position, `cdd_only` / `hmmer_only` when it does not, `not_compared` when the family has
models in only one library or HMMER is off. HMMER runs only if at least one family lists a Pfam
accession.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `enabled` | `bool` | `true` | `bool()` | Master toggle. When `false`, HMMER is not run, no HMMER table is written and every call is `not_compared`. |
| `use_gathering_threshold` | `bool` | `false` | `bool()` | Use Pfam's curated `--cut_ga` thresholds instead of manual E-value thresholds. When `true`, `evalue` and `dom_evalue` are ignored. |
| `evalue` | `float` | `1e-5` | `num(min=0)`, optional | Sequence-level E-value threshold for hmmsearch. Ignored when `use_gathering_threshold` is `true`. |
| `dom_evalue` | `float` | `1e-3` | `num(min=0)`, optional | Per-domain E-value threshold. Controls sensitivity for individual domain hits within multi-domain sequences. Ignored when `use_gathering_threshold` is `true`. |
| `min_score` | `float\|null` | `null` | `any(num(min=0), null())`, optional | Minimum sequence bit score to report (`hmmsearch -T`). `null` disables it. Ignored when `use_gathering_threshold` is `true`. |
| `min_coverage` | `float` | `0.5` | `num(min=0, max=1)` | Minimum fraction of the HMM profile aligned to the target sequence. Range 0.0-1.0. |
| `min_alignment_length` | `int` | `20` | `int(min=0)` | Minimum alignment length in amino acids. Filters short, unreliable matches. |
| `max_sensitivity` | `bool` | `false` | `bool()` | Enable hmmsearch `--max` flag, which disables all heuristic filters. Significantly slower but maximally sensitive. |
| `bias_filter` | `bool` | `true` | `bool()` | Enable the bias composition filter. Disabling (`false`) increases sensitivity for biased-composition sequences at the cost of more false positives. |
| `seed` | `int` | `67` | `int(min=0)` | Random seed for hmmsearch reproducibility. |

### `loci` -- Locus Definition

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `gap` | `int` | `20000` | `int(min=0)`, optional | Largest distance in bases between hits for them to be chained into one locus. Only hits that passed the quality gate are chained, on either strand. 20 kb keeps the exons of a gene together (the domain-bearing exons of human PRDM9 are up to 11 kb apart) without fusing tandem copies (40 kb apart in *Desmodus*); raise it for a gene family with longer introns. See [ADR-008](adr/ADR-008-loci.md). |

A locus reports which families it carries, each with its evidence tier, the fraction of a model
covered when the locus's exon fragments are put together, and the number of separate places the
family occurs at ([ADR-009](adr/ADR-009-completeness-across-exons.md)). It gives no verdict of
"complete": which families make a whole gene depends on the protein family studied.

### `programs` -- External Program Names (Optional)

Override BLAST+ and related program binary names. This entire section is optional; when omitted, bioconda-standard names are used. Useful when your system provides modern BLAST+ naming conventions (e.g., `rpsblast+` instead of `rpsblast`).

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `tblastn` | `str` | `'tblastn'` | `str()`, optional | Binary name for tBLASTn. |
| `rpsblast` | `str` | `'rpsblast'` | `str()`, optional | Binary name for RPSBLAST. Use `'rpsblast+'` on some systems. |
| `blast_formatter` | `str` | `'blast_formatter'` | `str()`, optional | Binary name for BLAST formatter. |
| `makeblastdb` | `str` | `'makeblastdb'` | `str()`, optional | Binary name for makeblastdb. |
| `rpsbproc` | `str` | `'rpsbproc'` | `str()`, optional | Binary name for rpsbproc post-processor. |
| `makeprofiledb` | `str` | `'makeprofiledb'` | `str()`, optional | Binary name for makeprofiledb (CDD subset construction). |

### `query` -- Query Protein (single query)

Specifies a single query protein sequence used for tBLASTn homology searching. This is the backward-compatible format. See also `queries` below for multi-query support.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `format` | `str` | `'fa'` | `str()` | File format extension for the query sequence file (e.g., `'fa'`, `'fasta'`). |
| `accession` | `str` | -- | `str(matches='[A-Z]{2,}_?[0-9]+\.[0-9]{1,2}')` | NCBI protein accession with version number. Must match the pattern `XX_123456.1` or `XXX123456.1`. Example: `'NP_064612.2'`. |

### `queries` -- Multi-Query Proteins

A YAML mapping of NCBI protein accessions to display labels. Each query runs the full pipeline independently (tBLASTn through domain detection), and results are merged with cross-query domain deduplication at the end.

```yaml
queries:
  'NP_064612.2': 'human PRDM9'
  'NP_659058.3': 'mouse PRDM9'
```

**Validation:** `map(str(), key=str(), required=False)` -- any number of string key-value pairs.

- **Key:** NCBI protein accession (used for Entrez download and `Query_Accession` column in outputs)
- **Value:** Display label (used for plot titles; filesystem-safe version used for directory names, e.g., `human_PRDM9`)

When `queries` is present, it takes precedence over `query`. Both formats are supported for backward compatibility. The pipeline normalizes both to an internal `QUERY_DICT: Dict[str, str]`.

**Output directory structure:** Per-query intermediate files are stored under `{query_label}/` subdirectories (e.g., `results/blast/human_PRDM9/`). Aggregate tables are organized by stage under `results/tables/` (e.g., `results/tables/blast/`, `results/tables/domains/`). Merged outputs (cross-query deduplicated) are written at the root of each stage subdirectory (e.g., `results/tables/domains/domains.parquet`).

**Deduplication:** After all queries complete, `merge_queries.R` merges overlapping domain annotations across queries using GenomicRanges. The `Query_Accession` column in merged outputs lists all contributing query accessions (comma-separated).

### `execution` -- Runtime Parameters

Controls parallelism and NCBI identification. An NCBI API key is never written in the config: export it as `NCBI_API_KEY` and the download tools pick it up.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `num_cores` | `int` | `16` | `int(min=1)` | Maximum number of CPU cores for Snakemake. Per-species rules run in parallel up to this limit; HMMER consumes all cores per job (serialized). |
| `entrez_email` | `str` | -- | `str()` | Email address for NCBI Entrez API identification. Required by NCBI usage policies. |

### `display` -- Output Verbosity

Controls informational messages during pipeline execution.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `display_snakemake_info` | `bool` | `true` | `bool()` | Show Snakemake execution details in console output. |

### `root` -- Directory Paths

Absolute paths to the four root directories used by the pipeline. All intermediate and output directories are derived from these paths (see [architecture.md](architecture.md) for the full directory tree).

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `db_root_folder` | `str` | -- | `str()` | Root directory for all databases: BLAST nucleotide DBs, Pfam HMMs, CDD files, and rpsbproc data. Genome FASTA files (`<key>.fa`) are stored here as well. |
| `data_root_folder` | `str` | -- | `str()` | Root for pipeline input data. Contains `config/` (this config file, schema, environment spec) and `input/fastas/` (query protein files). |
| `results_root_folder` | `str` | -- | `str()` | Root for all pipeline outputs: per-species intermediates, aggregated tables, plots, GFF3 ranges, ASN archives. |
| `logs_root_folder` | `str` | -- | `str()` | Root for execution logs. Each pipeline stage writes timestamped logs here. |

### `species` -- Target Genomes

Each species is a key (underscore-separated, used in file names and as the Snakemake wildcard)
with a display name and a pin that states exactly which genome is meant.

```yaml
species:
  Desmodus_rotundus:
    name: Desmodus rotundus
    assembly: GCA_022682495.2
  My_isolate:
    name: My isolate
    sha256: 3f1c...   # 64 hexadecimal characters
```

| Parameter | Type | Required | Validation | Description |
|-----------|------|----------|------------|-------------|
| `name` | `str` | yes | `str()` | Display name used in tables and figures. |
| `assembly` | `str` | one of the two | `GCA_`/`GCF_` + 9 digits + `.N` | Versioned NCBI assembly accession. `--download-genomes` fetches exactly this version; a taxon name or an accession without a version is refused. |
| `sha256` | `str` | one of the two | 64 hex characters | SHA-256 of the genome FASTA. Use it alone for a genome you provide yourself as `<db_root_folder>/<key>.fa`, or together with `assembly` to state exactly which file was used. |

A `sha256` is checked against the FASTA immediately before that genome's BLAST database is built,
and the run stops if they differ; for a species without one, the checksum is printed in the job's
log so that it can be added. At least one of `assembly` and `sha256` is required
([ADR-005](adr/ADR-005-genomes-pinned-by-assembly.md)). The former form, `'Key': 'Display name'`,
is refused with a message showing the new one. Each download appends a line to
`<data_root_folder>/config/download.txt`: species key, accession, SHA-256 of the FASTA, UTC time
and the `datasets` version.

The species list drives the entire pipeline DAG. Each species generates an independent parallel track through all per-species rules. Adding or removing entries here is the primary way to scale the analysis.

---

## Configuration Examples

### Relaxed BLAST Thresholds for Weak Domain Capture

When searching for domains with weak homology to the query protein (e.g., KRAB, SET, SSXRD), standard thresholds filter them out. Relaxing all four BLAST parameters captures these hits at the cost of more noise in the initial hit set. Downstream HMMER filtering compensates for the relaxed entry criteria.

```yaml
blast:
  e_value: 0.01
  perc_identity: 35
  seq_length: 24
  bitscore: 45

hmmer:
  enabled: true
  profiles: ['KRAB', 'SET', 'SSXRD', 'zf-C2H2']
  evalue: 1e-5
  dom_evalue: 1e-3
  min_coverage: 0.5
  min_alignment_length: 20

orf:
  enabled: false    # Disabled -- short domains lack ORFs
```

### Targeted CDD Subset Searching

Instead of searching the full CDD (>60,000 domains), restrict RPSBLAST to specific domains of interest. The pipeline resolves ShortNames against `cddid.tbl` using prefix matching: `'KRAB'` matches both `'KRAB'` and `'KRAB_A-box'`. This reduces runtime and eliminates irrelevant domain annotations.

```yaml
rpsblast:
  e_value: 2
  comp_based_stats: 0
  seg: 'no'
  window_size: 0
  target_domains:
    - 'KRAB'
    - 'SET'
    - 'SSXRD'
    - 'zf-C2H2'
    - 'COG5048'
    - 'zf-H2C2'
    - 'ZnF_C2H2'
```

### Combined ORF + HMMER Filtering

For maximum specificity, enable both ORF detection and HMMER domain filtering. Sequences must contain an ORF of the specified minimum length *and* match at least one of the specified Pfam profiles to reach RPSBLAST. This configuration is aggressive -- expect significant data reduction. Short domains (KRAB, SET, SSXRD) may be entirely filtered out by the ORF requirement.

```yaml
orf:
  enabled: true
  min_orf_length: 200
  start_codon: 'ATG'
  longest_orf: true

hmmer:
  enabled: true
  profiles: ['zf-C2H2']
  use_gathering_threshold: false
  evalue: 1e-5
  dom_evalue: 1e-3
  min_coverage: 0.5
  min_alignment_length: 20
  seed: 67
```

### Multi-Query for Divergent Species

When analyzing distantly related species (e.g., bats and mice, ~85-90 MY divergence), a single query protein may miss divergent domains. Using multiple query proteins from different species increases sensitivity for detecting all domain families.

```yaml
# Multi-query: each query runs the full pipeline independently
queries:
  'NP_064612.2': 'human PRDM9'
  'NP_659058.3': 'mouse PRDM9'

# query: section is ignored when queries: is present
# query:
#   format: 'fa'
#   accession: 'NP_659058.3'
```

The pipeline creates per-query subdirectories for all intermediate files, then merges domain annotations across queries with overlap-based deduplication. The merged output's `Query_Accession` column lists all contributing accessions (comma-separated) for each domain annotation.

### Minimal Configuration for a New Species

To add a new species to an existing analysis, append an entry to the `species:` map with its
display name and its versioned assembly accession (find it on the assembly's NCBI page). To use a
genome you already have, place it at `<db_root_folder>/<key>.fa` and pin it with `sha256` instead.

```yaml
species:
  Desmodus_rotundus:
    name: Desmodus rotundus
    assembly: GCA_022682495.2
  Myotis_lucifugus:               # new species
    name: Myotis lucifugus
    assembly: GCF_000147115.1
```

Then run the full pipeline starting from genome download:

```bash
conda run -n rpsHunter ./rpsHunter --download-genomes --skip-validation
conda run -n rpsHunter ./rpsHunter --blast-dbs --skip-validation
conda run -n rpsHunter ./rpsHunter --blast --skip-validation
# ... continue through remaining stages
```

Snakemake will detect that the existing three species already have outputs and only process the new species through each stage.
