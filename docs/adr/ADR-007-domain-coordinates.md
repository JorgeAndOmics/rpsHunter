# ADR-007: Report each domain at its own genomic position

- **Status**: Proposed
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: yes

## Context

Every domain call is written with the genomic start and end of the whole tBLASTn hit it was found
in. Downstream steps (range merging, GFF3, contingency counts) treat those as the domain's position.

## Evidence

Measured on the five model genomes (three PRDM9 queries; CDD 3.22, Pfam 38.2, BLAST 2.17.0) against a
truth set built from RefSeq annotation of human, mouse and *Desmodus rotundus*: Pfam gathering
thresholds on the annotated proteins, mapped through exons to genome coordinates.

- On the toy genome, all domains of a planted gene are reported at the gene's full span.
- For zinc-finger calls matching one single-exon annotated domain, the median total boundary error
  is 411 bp (human), 528 bp (mouse) and 198 bp (*Desmodus*) with the hit's coordinates, and 48, 45
  and 54 bp with the domain's own position, computed from the call's residues within the translated
  hit (From, To). For the single-finger models (pfam00096, smart00355) the median error is 0 bp; the
  remainder comes from models that span several fingers (COG5048: 87 bp).
- Sensitivity and on-target share are unchanged by the placement.

## Decision

Compute each call's genomic start and end from its position in the translated hit and the hit's
strand, and report those as the domain's coordinates. Keep the hit's coordinates in separate columns.
Range merging, GFF3 export and counts use the domain coordinates.

## Consequences

- Positive: positions are those of the domain; overlapping-range merges stop fusing different
  domains that merely share a hit; GFF3 tracks become usable in a genome browser.
- Negative: merged-domain counts change, because merging by hit span over-merged.
- Neutral: the two strict expected-failure tests on coordinates are removed with the change.

## Alternatives considered

- **Keep hit coordinates and document it**: leaves every positional output misleading.

## Revisit trigger

The search stage changes to spliced alignments, where a hit spans introns.

## References

- Workbench `validation/out/funnel_*.txt`; `tests/integration/test_toy_pipeline.py`.
