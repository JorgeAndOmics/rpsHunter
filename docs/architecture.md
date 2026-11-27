# Architecture

rpsHunter is a Snakemake workflow (`workflow/Snakefile`) driven by a thin launcher (`rpsHunter`).
Each rule runs one script from `workflow/scripts/`. This page describes how the parts fit; what
each step decides scientifically is in [`filtering.md`](filtering.md), and every output file in
[`outputs.md`](outputs.md).

## Data flow

```
query proteins            genomes (FASTA)        CDD models      Pfam profiles
      |                        |                     |                 |
      |                  makeblastdb          cdd_subset_builder   pfam_subset
      |                        |                     |                 |
      +------ tBLASTn ---------+                     |                 |
                 |  blast/{query}/{species}          |                 |
            select_hits   (quality gate)             |                 |
                 |  fastas/, tables/blast/.../selected                 |
                 +----------- RPS-BLAST -------------+                 |
                 |                |  rpsblast/, asn/rpsblast/          |
                 |             rpsbproc                                |
                 |                |  rpsbproc/                         |
                 |           call_domains                              |
                 |                |  domains/{query}/{species}         |
                 +------------------------------- hmmsearch -----------+
                 |                |                    |  hmmer/{query}/{species}
                 |           domain_tables <-----------+
                 |                |  tables/domains/{calls,domains,hmmer_domains}, ranges/*.gff3
                 +---------- locus_tables
                                  |  tables/loci/{loci,locus_families}
                             figure_pdfs
                                  |  plots/{overview,domains,loci,architecture}.pdf
```

Per-query, per-species jobs run in parallel; the two table rules at the end read every species and
query at once.

## Modules

Scripts are flat modules (no package). **Pure** modules hold the logic: they take and return
tables, read nothing and write nothing, and may not import I/O libraries. **Script** modules parse
arguments, read inputs, call the pure functions and write outputs. `tests/unit/test_layering.py`
enforces the split.

| Module | Kind | What it does |
|---|---|---|
| `defaults.py` | config | Reads the one config file; exposes paths and settings; `ensure_dirs()` is the only creator of directories |
| `species_config.py`, `domain_families.py` | pure | Validate the `species:` and `domains:` blocks |
| `tags.py` | pure | A hit's Tag: a hash of its identity, stable across runs ([ADR-004](adr/ADR-004-content-derived-tags.md)) |
| `external.py` | I/O | Runs an external tool; a non-zero exit raises `ToolError` |
| `genome_check.py` | script | Verifies a genome's checksum before its BLAST database is built |
| `tblastn.py` | script | tBLASTn per query and species |
| `hit_table.py` | pure | `Internal_Stop` |
| `select_hits.py` | script | Quality gate; selected table, audit table, FASTA |
| `cdd_subset.py`, `pfam_subset.py` | script | Build the families' model libraries |
| `rpsblast.py`, `rpsbproc.py` | script | Domain search and NCBI post-processing |
| `domain_calls.py` | pure | Parse rpsbproc's report; evidence tier; genomic and model coordinates |
| `call_domains.py` | script | One call table per query and species |
| `hmmer.py` | script | HMMER domains per query and species |
| `domain_instances.py`, `gff3.py`, `method_agreement.py` | pure | Merge calls into instances; GFF3 lines; agreement labels |
| `domain_tables.py` | script | Run-level domain tables and GFF3 files |
| `loci.py` | pure | Chain hits into loci; family evidence per locus |
| `locus_tables.py` | script | Run-level locus tables |
| `aggregate.py` | script | Concatenate per-species tables |
| `log.py` | I/O | The one line shape every job writes |
| `stages.py`, `guard.py` | pure | The launcher's stage table; the heavy-rule guard |
| `preflight.py`, `launcher.py` | script | Checks before a run; the one Snakemake call |
| `console.py` | pure | What the launcher shows of a run |
| `figures.R`, `plots/` | R | The figure PDFs: `style.R` (visual system), `tables.R` (counts), `pages.R` (pages) |
| `utils/log.R` | R | The log line shape, for R |
| `manifest.py` | script | The manifest of a run |
| `query_config.py` | pure | Validates the `queries:` block |
| `cdd_release.py` | script | Refuses downloaded CDD files of another release than configured |

Every script is under the quality gate. The list of files exempt by name
([ADR-002](adr/ADR-002-quality-gate-legacy-list.md)) reached zero on 2026-11-27 and a test keeps it there.

## Contracts

- **One config per run.** `RPSHUNTER_CONFIG` names a complete config file; nothing is merged.
- **A failing tool fails the job.** An empty result (no hit, no domain) is a valid, empty table with
  its full set of columns; a tool error is never turned into one.
- **Every job writes a log and a benchmark** under `logs/jobs/<rule>/` and `logs/benchmarks/<rule>/`.
- **Config changes rerun what they affect.** Each rule lists the config sections it depends on
  (`tracked(...)` in the Snakefile), so changing a threshold reruns that step and what follows.
- **Tags join tables.** Every hit has a 12-character Tag derived from its query, species, sequence
  and coordinates. Calls, HMMER domains and loci refer to hits by Tag.
- **Coordinates** are 1-based, inclusive, with `Start <= End`, on the genome's own sequence names;
  strand is a separate column.
- **Directories** are created only by `defaults.ensure_dirs()`.

## Reproducibility

The conda environment is locked (`data/config/environment.lock-linux-64.txt`). Genomes are pinned by
assembly accession or checksum ([ADR-005](adr/ADR-005-genomes-pinned-by-assembly.md)), Pfam by
release, and the CDD release is verified. HMMER is seeded. Two runs on the same inputs give
byte-identical tables, and every run leaves a manifest of what it was made with.

## Tests

`make check` runs lint, types and tests. Unit tests cover each module; the integration test builds
a small synthetic genome with planted gene copies, fragments and random sequence, runs the whole
pipeline on it against a miniature database, and checks that what was planted is what is reported.
