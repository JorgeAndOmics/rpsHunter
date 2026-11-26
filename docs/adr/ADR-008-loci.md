# ADR-008: Loci from gate-passing hits, chained within 20 kb

- **Status**: Accepted (2026-11-27)
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: yes

## Context

Hits are chained into loci when they lie within 50 kb of each other on the same sequence, using
every tBLASTn hit, including those that failed the gates. A locus is called complete when all
target families occur anywhere inside it.

## Evidence

Measured on the five model genomes (three PRDM9 queries; CDD 3.22, Pfam 38.2, BLAST 2.17.0) against a
truth set built from RefSeq annotation of human, mouse and *Desmodus rotundus*: Pfam gathering
thresholds on the annotated proteins, mapped through exons to genome coordinates.

- A genome of random DNA gets a locus (12 chance hits), and two separate fragments 17 kb apart are
  pooled into one locus called complete (toy genome).
- Mouse: one of the two complete loci is a 414 kb cluster of 15 annotated zinc-finger genes with no
  PRDM9-like gene.
- *Desmodus*: four tandem PRDM9-like genes about 40 kb apart are pooled into one 163 kb locus; another
  locus spans 1 Mb and 41 annotated genes. Loci reach 1.1 Mb.
- Chaining gap against annotation (complete loci; of which holding exactly one annotated
  PRDM9-family gene; of which pooling several), gate-passing hits, calls at subset E <= 0.01:

  | Gap | Human | Mouse | *Desmodus* |
  |---|---|---|---|
  | 10 kb | 0 | 1; 1; 0 | 15; 14; 0 |
  | 20 kb | 2; 2; 0 | 1; 1; 0 | 15; 15; 0 |
  | 30 kb | 2; 2; 0 | 1; 1; 0 | 14; 13; 1 |
  | 50 kb | 2; 2; 0 | 1; 1; 0 | 12; 9; 3 |

  At 10 kb human PRDM9 and PRDM7 are split (their exons are up to 11 kb apart); at 30 kb and above
  tandem copies start to pool.
- Using all hits or only gate-passing hits gives the same complete loci; only the latter stops
  chance hits forming loci.

## Decision

Build loci from hits that passed the gates, with a default chaining gap of 20 kb
(`hit_domain_inventory.locus_gap`), and report per locus which families are present, with the
evidence tier of ADR-006 and the coverage of ADR-009. A locus carries no verdict of its own: the
single "complete" flag is dropped, because what makes a gene complete depends on the family being
studied (the PRDM9 literature uses KRAB plus SET), and the listed families let the user decide.

## Consequences

- Positive: on the three annotated genomes every complete locus corresponds to exactly one annotated
  PRDM9-family gene; no locus from noise.
- Negative: a gene with an intron longer than 20 kb between domain-bearing exons is split; the gap
  stays configurable for that reason.
- Neutral: locus counts rise (human 708 to 896) because clusters are no longer fused.

## Alternatives considered

- **Strand-aware chaining**: no gain measured, and one *Desmodus* gene is lost at 20 kb.
- **Keep 50 kb**: three pooled loci in *Desmodus*.

## Revisit trigger

A target gene family with introns longer than 20 kb, or a move to gene-model-based loci.

## References

- Workbench `validation/out/locus_gap_*.txt`, `loci_*.txt`.
