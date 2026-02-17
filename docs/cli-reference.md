# rpsHunter CLI Reference

## Invocation

```bash
conda run -n rpsHunter ./rpsHunter <flag> [--skip-validation]
```

The `rpsHunter` launcher script (project root) forwards arguments to
`workflow/scripts/rpsHunter.py`, which dispatches a `snakemake` subprocess via
`run_snakemake_rule()`.

**One flag per invocation.** Each CLI flag maps to a single Snakemake target
rule. Passing multiple pipeline flags in the same invocation causes unknown-flag
errors because argparse forwards unrecognised tokens to Snakemake, which rejects
them. Run each stage as a separate command.

```bash
# Correct
conda run -n rpsHunter ./rpsHunter --hmmer --skip-validation
conda run -n rpsHunter ./rpsHunter --blast --skip-validation

# Wrong -- the second flag is forwarded to Snakemake and rejected
conda run -n rpsHunter ./rpsHunter --hmmer --blast --skip-validation
```

---

## Flag Reference

| Flag | Short | Snakemake Target | Description |
|------|-------|------------------|-------------|
| `--setup-databases` | | `setup_databases` | Download CDD, rpsbproc annotation data, and SMP files. |
| `--download-genomes` | | `genome_downloader` | Fetch reference genomes from NCBI Datasets for all configured species. |
| `--download-query` | | `query_downloader` | Download the query protein sequence from NCBI. |
| `--blast-dbs` | | `blast_db_generator` | Build BLAST nucleotide databases (one per species genome). |
| `--blast` | | `blast_parser` | Filter BLAST hits, export per-species FASTAs, write selected tables, and produce the aggregate table. **Note:** this triggers the _parsing_ step, not the tBLASTn search itself. tBLASTn runs automatically as an upstream Snakemake dependency. |
| `--orf-analyser` | | `orf_analyser` | Run ORF detection and enrich per-species BLAST parquets. No-op when `orf.enabled: false` in config. |
| `--hmmer` | | `hmmer` | Run HMMER domain filtering on BLAST hits. No-op when `hmmer.enabled: false` in config. Serialised execution (one species at a time, all configured cores per job). |
| `--rpsblast` | | `rpsblaster` | Search filtered FASTAs against NCBI CDD (or a subset database when `rpsblast.target_domains` is set). |
| `--rpsbproc` | | `rpsbproc` | Post-process RPSBLAST ASN output with rpsbproc to produce domain annotations. |
| `--rpsbproc-parser` | | `rpsbproc_parser` | Parse rpsbproc text output into structured parquets and aggregate into `tables/domains.parquet`. |
| `--completeness-detector` | | `completeness_detector` | Generate tile plot (domain completeness heatmap). |
| `--contingency-parser` | | `contingency_sorter` | Generate 3D scatter plot, contingency table (CSV), and GFF3 genomic coordinate files. |
| `--concordance` | | `concordance` | Run multi-method concordance analysis (HMMER vs CDD agreement). Produces per-query and merged concordance tables. |
| `--skip-validation` | `-skp` | _(none)_ | Skip the pre-run validation suite. Does not dispatch any Snakemake rule. |

---

## Full Pipeline Execution Order

Run each line as a separate command. Snakemake resolves upstream dependencies
within each stage automatically.

```bash
# 1. Database setup
conda run -n rpsHunter ./rpsHunter --setup-databases

# 2. Data acquisition
conda run -n rpsHunter ./rpsHunter --download-genomes
conda run -n rpsHunter ./rpsHunter --download-query

# 3. BLAST database construction
conda run -n rpsHunter ./rpsHunter --blast-dbs

# 4. ORF analysis (only if orf.enabled: true)
conda run -n rpsHunter ./rpsHunter --orf-analyser

# 5. HMMER filtering (only if hmmer.enabled: true)
#    Automatically triggers tBLASTn as an upstream dependency.
conda run -n rpsHunter ./rpsHunter --hmmer

# 6. BLAST parsing, filtering, and aggregation
#    Reads from ENRICHED_BLAST_DIR (hmmer > orf > blast, based on enabled flags).
conda run -n rpsHunter ./rpsHunter --blast

# 7. Domain detection
conda run -n rpsHunter ./rpsHunter --rpsblast
conda run -n rpsHunter ./rpsHunter --rpsbproc
conda run -n rpsHunter ./rpsHunter --rpsbproc-parser

# 8. Visualisation and export
conda run -n rpsHunter ./rpsHunter --completeness-detector
conda run -n rpsHunter ./rpsHunter --contingency-parser

# 9. Concordance analysis (optional, requires --blast and --rpsbproc-parser)
conda run -n rpsHunter ./rpsHunter --concordance
```

If both ORF and HMMER are disabled, skip steps 4 and 5. The `--blast` step will
read directly from `results/blast/`.

Append `--skip-validation` to any command to bypass the pre-run validation suite.

---

## Dependency Resolution

Each CLI flag triggers a Snakemake _target rule_. Snakemake automatically pulls
in every upstream rule required to satisfy that target's inputs. The table below
shows the per-species rules expanded by each target and their key upstream
dependencies.

| Target Rule | Per-Species Rules | Key Upstream Dependencies |
|-------------|-------------------|---------------------------|
| `setup_databases` | _(none)_ | `cdd_database_downloader`, `rpsbproc_data_downloader`, `cdd_smp_downloader` |
| `genome_downloader` | `genome_downloader_setup` x N | _(network)_ |
| `query_downloader` | `query_downloader_setup` | _(network)_ |
| `blast_db_generator` | `blast_db_generator_setup` x N | Genome FASTA files |
| `orf_analyser` | `orf_species` x N | `blast_species` |
| `hmmer` | `hmmer_species` x N | `blast_species` (or `orf_species` if ORF enabled) |
| `blast_parser` | `blast_parser_species` x N, `aggregate_blast` | ENRICHED_BLAST_DIR parquets |
| `rpsblaster` | `rpsblast_species` x N | `blast_parser_species` FASTAs, CDD database |
| `rpsbproc` | `rpsbproc_species` x N | `rpsblast_species` ASN files, rpsbproc data |
| `rpsbproc_parser` | `rpsbproc_parser_species` x N, `aggregate_domains` | `rpsbproc_species` text files |
| `completeness_detector` | `pq_completeness_detector` x Q, `merged_completeness_detector` | Per-query + merged `aggregate_domains`, `aggregate_blast`, `aggregate_rpsblast` |
| `contingency_sorter` | `pq_contingency_sorter` x Q, `merged_contingency_sorter` | Per-query + merged `aggregate_domains`, `aggregate_blast`, `aggregate_rpsblast` |
| `concordance` | `pq_concordance` x Q, `merge_concordance` | Per-query `aggregate_blast`, `aggregate_domains` |

**Q** denotes the number of configured queries (1 for single query, N for multi-query).

**ENRICHED_PATTERN** is resolved once at Snakefile parse time:

- `results/hmmer/{query_label}/{species}.parquet` when `hmmer.enabled: true`
- `results/orf/{query_label}/{species}.parquet` when only `orf.enabled: true`
- `results/blast/{query_label}/{species}.parquet` when both are disabled

---

## Pre-Run Validation

When `--skip-validation` is _not_ provided, the pipeline runs a validation suite
before dispatching any Snakemake rule. The suite checks (in order):

1. **Config schema** -- validates `data/config/config.yaml` against
   `data/config/schema.yaml` using Yamale.
2. **Query accession** -- confirms the configured protein accession exists in
   NCBI's protein database via Entrez.
3. **Genome FASTAs** -- verifies that `.fa` files exist and contain valid FASTA
   records for each configured species (skipped when `use_species_list: true`
   and species are defined in config).
4. **External programs** -- checks that all required programs are on PATH:
   - Always checked: `tblastn`, `rpsblast`, `blast_formatter`, `makeblastdb`,
     `rpsbproc`, `datasets`
   - Checked when `rpsblast.target_domains` is set: `makeprofiledb`
   - Checked when `hmmer.enabled: true`: `hmmsearch`, `hmmfetch`, `hmmpress`
5. **NCBI API key** -- prompts for `NCBI_API_KEY` if not set in the environment.

After all checks pass, the user is prompted to confirm execution (`Proceed [Y/n]`).

To bypass the entire suite:

```bash
conda run -n rpsHunter ./rpsHunter --blast --skip-validation
```

---

## Safe Rerun Recipes

Snakemake only reruns rules whose output files are missing or older than their
inputs. To force a rerun of specific stages, delete their output files and then
execute the relevant CLI flags.

**Never use `snakemake --forceall`.** This forces rerun of all rules including
genome downloads, which deletes existing genome FASTA files and re-downloads
gigabytes of data.

### Rerun from HMMER onward

Use this after changing HMMER parameters in `config.yaml` or modifying
`hmmer.py`.

```bash
rm -rf results/hmmer/*.parquet \
       results/fastas/*.fa results/fastas/*.parquet \
       results/tables/selected/ \
       results/tables/aggregate.* \
       results/rpsblast/*.parquet results/asn/rpsblast/*.asn \
       results/rpsbproc/*.txt \
       results/domains/*.parquet \
       results/tables/domains.* results/tables/rpsblast.* \
       results/tables/contingency_table.csv \
       results/plots/tile_plot.png results/plots/scatter_3D_plot.html \
       results/ranges/*.gff3

conda run -n rpsHunter ./rpsHunter --hmmer --skip-validation
conda run -n rpsHunter ./rpsHunter --blast --skip-validation
conda run -n rpsHunter ./rpsHunter --rpsblast --skip-validation
conda run -n rpsHunter ./rpsHunter --rpsbproc --skip-validation
conda run -n rpsHunter ./rpsHunter --rpsbproc-parser --skip-validation
conda run -n rpsHunter ./rpsHunter --completeness-detector --skip-validation
conda run -n rpsHunter ./rpsHunter --contingency-parser --skip-validation
```

### Rerun from BLAST parsing onward

Use this after changing BLAST filtering thresholds (`blast.e_value`,
`blast.perc_identity`, etc.) or toggling ORF/HMMER enrichment flags.

```bash
rm -rf results/fastas/*.fa results/fastas/*.parquet \
       results/tables/selected/ \
       results/tables/aggregate.* \
       results/rpsblast/*.parquet results/asn/rpsblast/*.asn \
       results/rpsbproc/*.txt \
       results/domains/*.parquet \
       results/tables/domains.* results/tables/rpsblast.* \
       results/tables/contingency_table.csv \
       results/plots/tile_plot.png results/plots/scatter_3D_plot.html \
       results/ranges/*.gff3

conda run -n rpsHunter ./rpsHunter --blast --skip-validation
conda run -n rpsHunter ./rpsHunter --rpsblast --skip-validation
conda run -n rpsHunter ./rpsHunter --rpsbproc --skip-validation
conda run -n rpsHunter ./rpsHunter --rpsbproc-parser --skip-validation
conda run -n rpsHunter ./rpsHunter --completeness-detector --skip-validation
conda run -n rpsHunter ./rpsHunter --contingency-parser --skip-validation
```

### Rerun visualisation only

Use this after modifying R visualisation scripts or wanting to regenerate plots
from existing domain data.

```bash
rm -f results/plots/tile_plot.png \
      results/plots/scatter_3D_plot.html \
      results/tables/contingency_table.csv \
      results/ranges/*.gff3

conda run -n rpsHunter ./rpsHunter --completeness-detector --skip-validation
conda run -n rpsHunter ./rpsHunter --contingency-parser --skip-validation
```

### Rerun domain annotation only

Use this after changing RPSBLAST parameters or the CDD subset configuration.

```bash
rm -rf results/rpsblast/*.parquet results/asn/rpsblast/*.asn \
       results/rpsbproc/*.txt \
       results/domains/*.parquet \
       results/tables/domains.* results/tables/rpsblast.* \
       results/tables/contingency_table.csv \
       results/plots/tile_plot.png results/plots/scatter_3D_plot.html \
       results/ranges/*.gff3

conda run -n rpsHunter ./rpsHunter --rpsblast --skip-validation
conda run -n rpsHunter ./rpsHunter --rpsbproc --skip-validation
conda run -n rpsHunter ./rpsHunter --rpsbproc-parser --skip-validation
conda run -n rpsHunter ./rpsHunter --completeness-detector --skip-validation
conda run -n rpsHunter ./rpsHunter --contingency-parser --skip-validation
```

---

## Conditional Pipeline Steps

Two enrichment stages are gated by config flags and become no-ops when disabled:

| Config Key | Default | Effect When Disabled |
|------------|---------|----------------------|
| `orf.enabled` | `true` | `--orf-analyser` produces no output. ENRICHED_BLAST_DIR skips `results/orf/`. |
| `hmmer.enabled` | `true` | `--hmmer` produces no output. ENRICHED_BLAST_DIR skips `results/hmmer/`. |

When both are disabled, `--blast` reads directly from `results/blast/` (raw
tBLASTn parquets with no enrichment columns).

---

## Parallelism

Per-species rules (`blast_species`, `orf_species`, `blast_parser_species`,
`rpsblast_species`, `rpsbproc_species`, `rpsbproc_parser_species`) declare
`threads: 1` and run N species jobs simultaneously (where N is bounded by
`execution.num_cores`).

The `hmmer_species` rule declares `threads: workflow.cores`, so HMMER runs one
species at a time using all available cores. This prevents CPU over-subscription
since `hmmsearch` is internally parallelised.

---

## Output Locations

All per-species outputs now include a `{query_label}` subdirectory for multi-query isolation. When using a single `query:` config, there is one query label directory.

| Stage | Per-Query Per-Species Output | Per-Query Aggregate | Merged Output |
|-------|------------------------------|---------------------|---------------|
| tBLASTn | `results/blast/{ql}/{species}.parquet` | `results/tables/{ql}/aggregate.parquet` | `results/tables/aggregate.parquet` |
| ORF | `results/orf/{ql}/{species}.parquet` | _(enriches blast)_ | |
| HMMER | `results/hmmer/{ql}/{species}.parquet` | _(enriches blast)_ | |
| BLAST parser | `results/fastas/{ql}/{species}.fa` | `results/tables/{ql}/aggregate.parquet` | `results/tables/aggregate.parquet` |
| RPSBLAST | `results/rpsblast/{ql}/{species}.parquet` | `results/tables/{ql}/rpsblast.parquet` | `results/tables/rpsblast.parquet` |
| rpsbproc | `results/rpsbproc/{ql}/{species}.txt` | _(text)_ | |
| Domain parser | `results/domains/{ql}/{species}.parquet` | `results/tables/{ql}/domains.parquet` | `results/tables/domains.parquet` (deduplicated) |
| Tile plot | | `results/plots/{ql}/tile_plot.png` | `results/plots/tile_plot.png` |
| 3D plot | | `results/plots/{ql}/scatter_3D_plot.html` | `results/plots/scatter_3D_plot.html` |
| Contingency | | `results/tables/{ql}/contingency_table.csv` | `results/tables/contingency_table.csv` |
| GFF3 | | `results/ranges/{ql}/{species}.gff3` | `results/ranges/{species}.gff3` |
| Concordance | | `results/concordance/{ql}/` | `results/concordance/` |

`{ql}` = query label (e.g., `human_PRDM9`, `mouse_PRDM9`)
