# ADR-009: Completeness measured per locus across exons, not from per-hit flags

- **Status**: Accepted (2026-11-27)
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: yes

## Context

`rpsbproc` flags a call incomplete (N, C or NC) when 20% or more of the model is missing at an end.
The pipeline reports these flags as domain completeness: "Complete", "Truncated", "Bitruncated" in
the tile plot, the completeness bars and the contingency table. Each call comes from one tBLASTn
hit, and a hit covers one exon.

## Evidence

Measured on the five model genomes (three PRDM9 queries; CDD 3.22, Pfam 38.2, BLAST 2.17.0) against a
truth set built from RefSeq annotation of human, mouse and *Desmodus rotundus*: Pfam gathering
thresholds on the annotated proteins, mapped through exons to genome coordinates.

Every domain below is complete in the annotated protein, so a flag can only come from how the genome
was searched. Share of matching calls flagged incomplete (subset E <= 1e-3):

| True domain | Human | Mouse | *Desmodus* |
|---|---|---|---|
| SET (always split across exons) | 100% of 27 | 100% of 24 | 100% of 69 |
| Zinc finger split across exons | 99.5% of 364 | 98.9% of 372 | 99.2% of 477 |
| Zinc finger in one exon | 51.1% of 4,399 | 50.3% of 4,611 | 57.0% of 4,089 |
| KRAB split across exons | 12.4% of 129 | 29.4% of 109 | 2.6% of 78 |
| KRAB in one exon | 0% of 4 | 0% of 8 | 0% of 3 |

Adding the PRDM9-specific SET model (cd19193, absent today because its name starts with "PR-SET")
does not help in real genomes: SET stays a superfamily-level call in all five, because no single
exon covers enough of the model. On the intronless toy gene the same model gives a specific hit.

## Decision

Stop presenting the per-hit flag as truncation. For each locus and model, combine the calls from its
exons in model coordinates (RPS-BLAST reports the aligned model positions) and report the fraction
of the model covered; "complete" means covered beyond a threshold (default 80%, matching the flag's
own definition). The per-hit flag stays in the per-hit table under its own name.

## Consequences

- Positive: a multi-exon domain can be recognised as complete; completeness figures describe
  biology instead of intron positions.
- Negative: a new aggregation step; coverage across exons can be inflated when two paralogous
  fragments fall in one locus (limited by ADR-008).
- Neutral: the three completeness figures are redrawn from the new measure (figures milestone).

## Alternatives considered

- **Spliced protein-to-genome alignment** (miniprot, exonerate) before the domain search: gives whole
  proteins and would make flags meaningful, but replaces the search stage; a larger change, noted
  for the backlog.
- **Keep the flags and warn**: leaves the main figures misleading.

## Revisit trigger

The search stage produces spliced gene models.

## References

- Workbench `validation/out/incomplete_*.txt`; NCBI rpsbproc README.
