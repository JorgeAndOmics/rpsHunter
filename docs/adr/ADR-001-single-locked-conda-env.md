# ADR-001: One locked conda environment, no container

- **Status**: Accepted
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: no

## Context

`data/config/environment.yaml` pinned a handful of Python libraries and left every analysis tool
(BLAST+, HMMER, rpsbproc, R, Bioconductor) to whatever the solver picked on the day. Two people
creating the environment a month apart could get different tool versions without knowing it, and no
record existed of what any earlier run had used. The environment also referenced Anaconda's
`pkgs/main` channels and omitted five tools the scripts call.

## Decision

rpsHunter runs in a single conda environment for all rules, built from conda-forge and bioconda
only. The environment file states intent; an explicit lock,
`data/config/environment.lock-linux-64.txt`, records the exact build of every package and is what
`make env` installs. The lock is produced with conda's own `conda list --explicit --md5`. No
container image is maintained.

## Consequences

- Positive: a clone plus `make env` reproduces the tool versions exactly on any linux-64 machine;
  the lock is plain text, diffable, and needs no extra tool.
- Negative: the lock covers linux-64 only; other platforms must solve from the environment file and
  may get different versions. Servers without conda are not served.
- Neutral: changing an analysis tool version is a deliberate act: edit the environment file,
  regenerate the lock, and record the change in an ADR because it can change results.

## Alternatives considered

- **`conda-lock`**: multi-platform locks, but a new dependency for a pipeline that targets Linux only.
- **Per-rule environments (`--use-conda`)**: isolates tools, at the cost of many environments to
  solve and keep in step; nothing here conflicts enough to need it.
- **A container image**: the strongest portability, but an image to build, test and keep in step
  with the lock. Deferred, not rejected.

## Revisit trigger

A user needs to run on a machine without conda, a second platform must be supported, or two tools
require conflicting dependencies.

## References

- `data/config/environment.yaml`, `data/config/environment.lock-linux-64.txt`, `Makefile` (`env`).
