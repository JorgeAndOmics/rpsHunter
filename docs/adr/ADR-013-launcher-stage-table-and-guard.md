# ADR-013: One stage table, checks that always run, and a guard on heavy rules

- **Status**: Accepted (part of the modernisation plan approved 2026-11-26)
- **Date**: 2026-11-27
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: no

## Context

The launcher was two hand-written argument parsers, one calling the other, with one Snakemake call
per flag and only the first flag honoured. Its validation prompted for input and queried NCBI on
every run, so every documented command carried `--skip-validation` and nothing was ever checked.
Because a Snakemake target builds what it needs upstream, a harmless-looking call could start a
genome download or redo every tBLASTn search when a results folder had moved or a setting changed.

## Decision

1. The stages are rows of one table (`workflow/scripts/stages.py`): flag, phase, Snakemake
   targets, help line, the heavy rules the stage owns and the tools it calls. The parser, the run
   order, the guard and the tool check are derived from it. Any number of stages run as one
   Snakemake call; `--all` is every stage after Setup.
2. Flags are renamed for what they do (`--search`, `--domains`, `--download-queries`,
   `--download-databases`). A flag of an earlier version stops with its replacement.
3. Checks run before every launch and cannot be skipped: retired config keys, the schema, the
   family and species rules, the programs and the input files of the requested stages. They are
   offline, ask nothing, and report every problem at once.
4. The launcher does a dry run first. Downloads, BLAST database builds and tBLASTn start only when
   their own stage is requested, or with `--allow-heavy`; otherwise the run stops and explains.
   `--forceall` is never passed to Snakemake.
5. The config is `--config`, else `RPSHUNTER_CONFIG`, else the shipped file: one complete file.

## Consequences

- Positive: one command runs an analysis; mistakes are caught in seconds, before any compute; a
  download or a search can no longer start as a side effect.
- Negative: scripts and notes using the old flags must be updated (each old flag says how).
- Neutral: the design follows the sibling project RetroSeek's, reimplemented from its documented
  contract, so the two tools behave alike at the command line.

## Alternatives considered

- **Keep one stage per call**: simplest to keep, but it is the reason nobody could run the whole
  analysis with one command.
- **Guard only the genome download**: the search is the next most expensive surprise (ten minutes
  on the model genomes, hours on more), so it is guarded too.
