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
| 005 | [Genomes pinned by assembly accession or checksum](ADR-005-genomes-pinned-by-assembly.md) | Accepted | no |
| 006 | [RPS-BLAST E-values on the full-CDD scale, with a cut-off of 1](ADR-006-rpsblast-evalue-scale-and-cutoff.md) | Proposed | yes |
| 007 | [Report each domain at its own genomic position](ADR-007-domain-coordinates.md) | Proposed | yes |
| 008 | [Loci from gate-passing hits, chained within 20 kb](ADR-008-loci.md) | Proposed | yes |
| 009 | [Completeness measured per locus across exons](ADR-009-completeness-across-exons.md) | Proposed | yes |
| 010 | [State what the tool finds; HMMER as evidence, not a gate](ADR-010-scope-and-hmmer-evidence.md) | Proposed | yes |
| 011 | [Retire the ORF gate; flag internal stop codons](ADR-011-retire-orf-gate.md) | Proposed | no |
| 012 | [Domain families defined by CDD and Pfam accession](ADR-012-domain-families-by-accession.md) | Proposed | yes |
