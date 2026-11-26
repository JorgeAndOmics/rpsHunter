# Architectural decision records

Short records of decisions that shape rpsHunter, one file each, numbered in order and never
renumbered. An ADR is written when a choice is hard to reverse, when it changes results, or when the
reason would otherwise be lost.

Start from [`ADR-000-template.md`](ADR-000-template.md). A decision that changes scientific output
(a threshold, a method, a tool or database version) must carry measured before-and-after numbers on
the model genomes and is applied only once its status is Accepted. Commits that implement a decision
cite its number.

## Index

| #   | Title | Status | Affects results |
|-----|-------|--------|-----------------|
| 001 | [One locked conda environment, no container](ADR-001-single-locked-conda-env.md) | Accepted | no |
| 002 | [Quality gate with a whole-file legacy list](ADR-002-quality-gate-legacy-list.md) | Accepted | no |
| 003 | [MIT licence](ADR-003-mit-licence.md) | Accepted | no |
| 004 | [Tags derived from the hit, not random](ADR-004-content-derived-tags.md) | Accepted | yes |
