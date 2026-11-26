# Command-line reference

```bash
./rpsHunter [stages] [-n] [--cores N] [--config FILE] [--allow-heavy] [-- snakemake options]
```

Run inside the `rpsHunter` conda environment. `./rpsHunter --help` lists the stages by phase, with
examples, and needs no config.

## Stages

Any number of stages can be given in one call; they run in the order below, as one Snakemake run.

| Phase | Flag | What it builds |
|---|---|---|
| Setup | `--download-databases` | CDD, the rpsbproc annotation files and Pfam |
| Setup | `--download-genomes` | The genomes pinned by assembly accession |
| Setup | `--download-queries` | The query proteins |
| Search | `--blast-dbs` | BLAST databases of the genomes, after verifying their checksums |
| Search | `--search` | tBLASTn, hit selection, hit tables |
| Domains | `--hmmer` | HMMER domains of the selected hits (nothing to do if `hmmer.enabled` is false or no family lists a Pfam accession) |
| Domains | `--domains` | RPS-BLAST, rpsbproc, domain calls and instances, GFF3 tracks |
| Tables | `--loci` | Loci and the families each carries |
| | `--all` | Every stage after Setup |

## Options

| Option | Effect |
|---|---|
| `-n`, `--dry-run` | Run the checks, list the jobs, run nothing |
| `--cores N` | Cores to use; default `execution.num_cores` |
| `--config FILE` | Config for this run. Otherwise `$RPSHUNTER_CONFIG`, otherwise `data/config/config.yaml`. One complete file; nothing is merged |
| `--verbosity quiet\|normal\|verbose` | What the screen shows; default `display.verbosity`, else `normal`. See [console_style.md](console_style.md) |
| `--allow-heavy` | Let a stage start a download, database build or search it does not own (see the guard) |

## What happens before anything runs

1. **Checks** (offline, no questions, a second or two). All problems are listed together and
   nothing runs if there is one:
   - config keys of earlier versions, each with its replacement;
   - the schema (`data/config/schema.yaml`) and the rules on families and species;
   - the programs the requested stages call, on `PATH`;
   - the input files those stages read (genomes, queries, CDD models, Pfam), unless the same call
     downloads them.

   Two notes do not stop a run: the environment differing from the lock file, and less than 5 GB
   free for the results.
2. **A dry run**, whose job list is printed.
3. **The guard.** A stage builds whatever it needs upstream, which is convenient until a moved
   results folder or an edited setting makes `--loci` download genomes or redo every search. The
   downloads, the BLAST database build and tBLASTn therefore start only when their own stage is
   among those requested (`--all` covers the database build and the search, never a download).
   Otherwise the run stops, names the rule, Snakemake's reason and the stage to add.

## Exit status

| Status | Meaning |
|---|---|
| 0 | Done, or nothing to do |
| 2 | Wrong usage, a problem found by the checks, or a missing config |
| 3 | Stopped by the guard |
| other | Snakemake's own status: a job failed. The summary names the failed rule and its log under `logs/jobs/<rule>/`; everything the run printed is in `logs/runs/` |

## Flags of earlier versions

Using one stops with its replacement.

| Old | Now |
|---|---|
| `--setup-databases` | `--download-databases` |
| `--download-query` | `--download-queries` |
| `--blast` | `--search` |
| `--rpsblast`, `--rpsbproc`, `--rpsbproc-parser` | `--domains` |
| `--hit-domain-inventory` | `--loci` |
| `--concordance` | `Method_Agreement`, written by `--domains` |
| `--contingency-parser` | GFF3 tracks are written by `--domains` |
| `--orf-analyser`, `--completeness-detector`, `--extended-plots` | removed |
| `--skip-validation`, `-skp` | removed: the checks no longer prompt or use the network |

## Passing options to Snakemake

Everything after `--` goes to Snakemake unchanged, for example
`./rpsHunter --domains -- --rerun-triggers mtime`. `--forceall` is refused: it would rebuild the
databases and download the genomes again.

## What reruns when

Rules record the config sections they depend on. After editing the config, the next call reruns
exactly the affected steps:

| Changed | Reruns from |
|---|---|
| `blast.e_value` | tBLASTn (request `--search`) |
| `blast.perc_identity`, `seq_length`, `bitscore` | hit selection |
| `domains:` | model libraries, both domain searches |
| `rpsblast:` | RPS-BLAST |
| `hmmer:` | HMMER, then the domain tables |
| `loci.gap` | locus tables |

## Parallelism

Per-species jobs use one thread each and run side by side up to the number of cores. HMMER takes
all cores for one species at a time, since `hmmsearch` is parallel itself.
