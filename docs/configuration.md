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
| `e_value` | `float` | `10` | `num(min=0)` | Maximum E-value for RPSBLAST hits. Higher values capture borderline domain hits (e.g., SET at E~1.03). |
| `comp_based_stats` | `int` | `0` | `int(min=0, max=4)` | Composition-based statistics mode. `0` disables (more sensitive); `1`-`4` enable progressive correction. See BLAST+ documentation. |
| `seg` | `str` | `'no'` | `str(matches='^(yes\|no)$')` | SEG low-complexity filtering. `'no'` disables, allowing short/repetitive domains through. `'yes'` masks low-complexity regions. |
| `window_size` | `int` | `40` | `int(min=0)`, optional | Window size for multi-hit algorithm. `0` switches to single-hit mode (more sensitive, slower). Omit or set to `40` for default behavior. |
| `target_domains` | `list[str]` | `[]` | `list(str())`, optional | CDD ShortName list for subset database construction. Uses prefix matching with `_` (e.g., `'KRAB'` matches `'KRAB'` and `'KRAB_A-box'`). When non-empty, `cdd_subset.py` builds a targeted database from `cddid.tbl`. When empty or omitted, RPSBLAST searches the full CDD. |

### `orf` -- Open Reading Frame Analysis

Controls the ORF detection stage, which enriches BLAST hits with ORF information and optionally filters sequences that lack ORFs. When enabled, this step runs between BLAST and HMMER in the pipeline. See [filtering.md](filtering.md) for details on how ORF filtering interacts with HMMER gating.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `enabled` | `bool` | `false` | `bool()` | Master toggle. When `false`, the ORF stage is skipped entirely and no `ORF_Filtered` column is added to per-species parquets. |
| `min_orf_length` | `int` | `200` | `int(min=1)` | Minimum ORF length in nucleotides. Lower values (e.g., 50) capture short domains like KRAB/SET but increase false positives. |
| `start_codon` | `str` | `''` | `str(matches='^([ATGC]{3})?$')` | Required start codon (e.g., `'ATG'`). Empty string disables start codon filtering, accepting ORFs beginning at any codon. |
| `longest_orf` | `bool` | `true` | `bool()` | When `true`, retain only the longest ORF per reading frame. When `false`, all ORFs meeting the length threshold are reported. |

### `hmmer` -- HMMER Domain Filtering

Controls hmmsearch against the Pfam database. When enabled, sequences are enriched with HMM domain annotations and optionally filtered by domain match. HMMER runs after ORF analysis (if enabled) and before BLAST parsing. All quality-passing domain hits per sequence are retained as list-valued columns.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `enabled` | `bool` | `true` | `bool()` | Master toggle. When `false`, the HMMER stage is skipped and no `HMM_*` columns are added. |
| `profiles` | `list[str]` | `[]` | `list(str())`, optional | Pfam profile names to search. When non-empty, `hmmfetch` extracts a subset HMM before searching. When empty, the full Pfam database is searched. |
| `use_gathering_threshold` | `bool` | `false` | `bool()` | Use Pfam's curated `--cut_ga` thresholds instead of manual E-value thresholds. When `true`, `evalue` and `dom_evalue` are ignored. |
| `evalue` | `float` | `1e-5` | `num(min=0)`, optional | Sequence-level E-value threshold for hmmsearch. Ignored when `use_gathering_threshold` is `true`. |
| `dom_evalue` | `float` | `1e-3` | `num(min=0)`, optional | Per-domain E-value threshold. Controls sensitivity for individual domain hits within multi-domain sequences. Ignored when `use_gathering_threshold` is `true`. |
| `min_score` | `float\|null` | `null` | `any(num(min=0), null())`, optional | Minimum bit-score filter. `null` disables this filter. |
| `min_coverage` | `float` | `0.5` | `num(min=0, max=1)` | Minimum fraction of the HMM profile aligned to the target sequence. Range 0.0-1.0. |
| `min_alignment_length` | `int` | `20` | `int(min=0)` | Minimum alignment length in amino acids. Filters short, unreliable matches. |
| `max_sensitivity` | `bool` | `false` | `bool()` | Enable hmmsearch `--max` flag, which disables all heuristic filters. Significantly slower but maximally sensitive. |
| `bias_filter` | `bool` | `true` | `bool()` | Enable the bias composition filter. Disabling (`false`) increases sensitivity for biased-composition sequences at the cost of more false positives. |
| `seed` | `int` | `67` | `int(min=0)` | Random seed for hmmsearch reproducibility. |

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

### `query` -- Query Protein

Specifies the query protein sequence used for tBLASTn homology searching.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `format` | `str` | `'fa'` | `str()` | File format extension for the query sequence file (e.g., `'fa'`, `'fasta'`). |
| `accession` | `str` | -- | `str(matches='[A-Z]{2,}_?[0-9]+\.[0-9]{1,2}')` | NCBI protein accession with version number. Must match the pattern `XX_123456.1` or `XXX123456.1`. Example: `'NP_064612.2'`. |

### `execution` -- Runtime Parameters

Controls parallelism, NCBI data retrieval behavior, and internal identifiers.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `num_cores` | `int` | `16` | `int(min=1)` | Maximum number of CPU cores for Snakemake. Per-species rules run in parallel up to this limit; HMMER consumes all cores per job (serialized). |
| `use_species_list` | `bool` | `true` | `bool()` | When `true`, species are read from the `species:` section of this config. When `false`, species are inferred from `.fa` files in the database root directory. |
| `random_id_length` | `int` | `6` | `int(min=1)` | Length of the alphanumeric Tag appended to each BLAST hit for traceability through the pipeline. |
| `retrieval_time_lag` | `float` | `0.3` | `num(min=0)` | Delay in seconds between NCBI Entrez requests to avoid rate limiting. |
| `max_retrieval_attempts` | `int` | `9` | `int(min=1)` | Maximum retry attempts for failed NCBI data retrievals. |
| `entrez_email` | `str` | -- | `str()` | Email address for NCBI Entrez API identification. Required by NCBI usage policies. |

### `display` -- Output Verbosity

Controls informational messages during pipeline execution.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `display_snakemake_info` | `bool` | `true` | `bool()` | Show Snakemake execution details in console output. |
| `display_requests_warning` | `bool` | `false` | `bool()` | Show HTTP request warnings (e.g., SSL verification). |
| `display_operation_info` | `bool` | `false` | `bool()` | Show detailed per-operation progress messages. |

### `root` -- Directory Paths

Absolute paths to the four root directories used by the pipeline. All intermediate and output directories are derived from these paths (see [architecture.md](architecture.md) for the full directory tree).

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `db_root_folder` | `str` | -- | `str()` | Root directory for all databases: BLAST nucleotide DBs, Pfam HMMs, CDD files, and rpsbproc data. Species genome `.fa` files are also located here when `use_species_list` is `false`. |
| `data_root_folder` | `str` | -- | `str()` | Root for pipeline input data. Contains `config/` (this config file, schema, environment spec) and `input/fastas/` (query protein files). |
| `results_root_folder` | `str` | -- | `str()` | Root for all pipeline outputs: per-species intermediates, aggregated tables, plots, GFF3 ranges, ASN archives. |
| `logs_root_folder` | `str` | -- | `str()` | Root for execution logs. Each pipeline stage writes timestamped logs here. |

### `logging` -- Console Log Styling

Configures `coloredlogs` styling for console output. Two sub-sections control level-based and field-based styling respectively.

**`logging.level_styles`** -- per-severity colors:

| Level | Properties | Validation |
|-------|-----------|------------|
| `debug` | `color` | `str()` |
| `info` | `color`, `bold` | `str()`, `str(matches='^(yes\|no)$')` |
| `warning` | `color` | `str()` |
| `error` | `color`, `bold` | `str()`, `str(matches='^(yes\|no)$')` |
| `critical` | `color`, `bold`, `background` | `str()`, `str(matches='^(yes\|no)$')`, `str()` |

**`logging.field_styles`** -- per-field colors:

| Field | Properties | Validation |
|-------|-----------|------------|
| `asctime` | `color` | `str()` |
| `levelname` | `color` | `str()` |
| `name` | `color` | `str()` |

### `species` -- Target Genomes

A YAML mapping of species identifiers to display names. Keys use underscore-separated binomial names (used as filenames and Snakemake wildcards); values are human-readable names with spaces.

```yaml
species:
  'Desmodus_rotundus': 'Desmodus rotundus'
  'Antrozous_pallidus': 'Antrozous pallidus'
```

**Validation:** `map(str(), key=str())` -- any number of string key-value pairs.

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

### Minimal Configuration for a New Species

To add a new species to an existing analysis, append an entry to the `species:` map using the underscore-separated binomial name as the key. Ensure the genome `.fna` file exists in the database root directory (downloaded via `--download-genomes` or placed manually). No other configuration changes are required -- the new species inherits all threshold and filtering settings.

```yaml
species:
  'Desmodus_rotundus': 'Desmodus rotundus'
  'Antrozous_pallidus': 'Antrozous pallidus'
  'Molossus_molossus': 'Molossus molossus'
  'Myotis_lucifugus': 'Myotis lucifugus'        # New species added
```

Then run the full pipeline starting from genome download:

```bash
conda run -n rpsHunter ./rpsHunter --download-genomes --skip-validation
conda run -n rpsHunter ./rpsHunter --blast-dbs --skip-validation
conda run -n rpsHunter ./rpsHunter --blast --skip-validation
# ... continue through remaining stages
```

Snakemake will detect that the existing three species already have outputs and only process the new species through each stage.
