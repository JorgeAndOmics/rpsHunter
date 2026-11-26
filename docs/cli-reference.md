# Command-line reference

```bash
RPSHUNTER_CONFIG=path/to/config.yaml ./rpsHunter --<stage> [--skip-validation] [snakemake options]
```

Run inside the `rpsHunter` conda environment. `RPSHUNTER_CONFIG` names the one complete config file
of the run (default: `data/config/config.yaml`).

> The launcher is being redesigned (a single run-everything command, renamed flags, checks before
> launch). Until then it has the two limits described below.

## Stages

Each flag builds its outputs and anything upstream that is missing or out of date.

| Flag | What it builds |
|---|---|
| `--setup-databases` | Downloads CDD and the rpsbproc annotation files (Pfam is fetched by `--hmmer` when it is missing) |
| `--download-genomes` | Downloads the genomes pinned by assembly accession |
| `--download-query` | Downloads the query proteins |
| `--blast-dbs` | BLAST databases of the genomes (after verifying their checksums) |
| `--blast` | tBLASTn, hit selection, hit tables |
| `--hmmer` | HMMER domains of the selected hits (nothing to do if `hmmer.enabled` is false or no family lists a Pfam accession) |
| `--rpsblast` | RPS-BLAST of the selected hits against the families' CDD models |
| `--rpsbproc` | rpsbproc reports |
| `--rpsbproc-parser` | Domain calls, domain instances, HMMER domain table, GFF3 tracks |
| `--loci` | Loci and the families each carries |

A full run, in order:

```bash
for stage in setup-databases download-genomes download-query blast-dbs \
             blast hmmer rpsblast rpsbproc rpsbproc-parser loci; do
    ./rpsHunter --$stage --skip-validation || break
done
```

Because each stage builds what it needs, a later stage can be called directly once the databases
and genomes are in place; the loop above is the explicit form.

## Limits of the current launcher

1. **One stage per call.** Only the first stage flag is used; anything after it is passed to
   Snakemake.
2. **Use `--skip-validation`.** The built-in validation prompts for input and queries NCBI on every
   run; it is being replaced.

## Passing options to Snakemake

Anything the launcher does not recognise goes to Snakemake. The most useful:

| Option | Effect |
|---|---|
| `-n` | Dry run: list the jobs that would run. Do this first and read the list |
| `--rerun-incomplete` | Redo jobs that were interrupted |

Do not use `--forceall`: it would also rebuild databases and download genomes again.

## What reruns when

Rules record the config sections they depend on. After editing the config, the next call reruns
exactly the affected steps:

| Changed | Reruns from |
|---|---|
| `blast.e_value` | tBLASTn |
| `blast.perc_identity`, `seq_length`, `bitscore` | hit selection |
| `domains:` | model libraries, both domain searches |
| `rpsblast:` | RPS-BLAST |
| `hmmer:` | HMMER, then the domain tables |
| `loci.gap` | locus tables |

## Parallelism

Per-species jobs use one thread each and run side by side up to `execution.num_cores`. HMMER takes
all cores for one species at a time, since `hmmsearch` is parallel itself.

## Exit status and logs

The launcher exits non-zero if any job fails. A failed job's log is printed and kept under
`logs/jobs/<rule>/`; see [outputs.md](outputs.md#logs-and-benchmarks).
