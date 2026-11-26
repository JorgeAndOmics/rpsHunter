# ADR-002: Quality gate with a whole-file legacy list

- **Status**: Accepted
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: no

## Context

rpsHunter had no linting, type checking or tests. Introducing a strict gate (ruff with a wide rule
set, mypy strict, lintr with cyclomatic complexity 8, ShellCheck, pytest) against the existing code
reports several hundred findings, and 14 of 60 Python functions exceed the complexity limit. The
code is to be rewritten stage by stage with tests written first, which will take many changes. A
gate that is red for months teaches everyone to ignore it.

## Evidence

At commit c93989e, with the gate's rules applied to every file: 478 ruff findings in 14 Python
files, 185 lintr findings in 5 R scripts, 6 ShellCheck findings in 2 shell scripts; 30 of 60 Python
functions over 30 lines and 14 over complexity 8.

## Decision

`make check` is the single gate and must pass on every commit. Scripts written before the gate are
named in a legacy list and excluded as whole files from ruff, mypy, coverage, lintr and ShellCheck
until they are rewritten. The list may only shrink: the commit that rewrites a file removes its
entry, and a test (`tests/unit/test_legacy_list.py`) fails if an entry is added, if the copies of
the list disagree, or if an entry names a file that no longer exists. Findings are never silenced
line by line.

## Consequences

- Positive: the gate is green from the first day and means something: everything outside the list
  meets the full standard. Progress is visible as the list shortens.
- Negative: legacy files can still be edited without the gate noticing, so changes to them rely on
  review and on the regression run until they are rewritten.
- Neutral: complexity and size are still measured for every file by `make quality`, so the debt is
  tracked even while it is excluded from the gate.

## Alternatives considered

- **Fix everything first**: a dedicated campaign before any other work. Rejected because the stages
  are about to be redesigned; polishing code that is then rewritten is wasted effort.
- **Per-line suppressions (`noqa`, `nolint`)**: hundreds of markers that hide which findings are
  real and tend to outlive their reason.
- **A looser rule set now, tightened later**: every tightening would reopen files already touched.

## Revisit trigger

The legacy list is empty (remove the mechanism and its test), or it has not shrunk across two
releases (the rewrite has stalled and needs a plan).

## References

- `pyproject.toml` (legacy list), `.lintr.R`, `Makefile`, `.qlty/qlty.toml`.
