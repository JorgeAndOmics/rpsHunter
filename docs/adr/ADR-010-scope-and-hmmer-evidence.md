# ADR-010: State what the tool finds, and make HMMER evidence instead of a gate

- **Status**: Accepted (2026-11-27)
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: yes (second part)

## Context

rpsHunter is described as detecting protein domains in genomes. In fact it searches with query
proteins and annotates domains inside the hits. Separately, hits without a HMMER match are removed
before RPS-BLAST, and agreement between HMMER and CDD is then reported as "concordance".

## Evidence

Measured on the five model genomes (three PRDM9 queries; CDD 3.22, Pfam 38.2, BLAST 2.17.0) against a
truth set built from RefSeq annotation of human, mouse and *Desmodus rotundus*: Pfam gathering
thresholds on the annotated proteins, mapped through exons to genome coordinates.

Annotated domains still covered after each stage (human):

| Family | Annotated | tBLASTn hit | After gates | CDD call |
|---|---|---|---|---|
| KRAB | 381 | 95 | 95 | 90 |
| SSXRD | 10 | 2 | 2 | 2 |
| SET | 49 | 11 | 9 | 9 |
| Zinc fingers | 5,800 | 5,469 | 5,469 | 5,469 |

Mouse and *Desmodus* show the same shape (KRAB 408 to 84 and 267 to 65 at tBLASTn). Losses happen
at the search: a KRAB that does not resemble PRDM9's KRAB is never seen. The gates remove 9% of
hits and, of annotated domains, only SET exons (2 of 11 in human, 3 of 11 in mouse). Because every
sequence that reaches CDD already has a HMMER match, "concordance" cannot fall below what the gate
enforces, and it requires no positional agreement.

Tightening tBLASTn instead would cost sensitivity: at E <= 1 the search reaches 35 KRAB domains in
human instead of 95, and 8 at 1e-3. The published PRDM9 surveys search with the zinc-finger array
removed from the query and anchor on SET (Baker 2017; Cavassim 2022).

## Decision

1. Document the scope: rpsHunter finds regions similar to the query proteins and describes their
   domains; counts are of query-like loci, not of a domain family genome-wide. Keep the tBLASTn
   E-value at 10 and say why.
   The shipped PRDM9 example keeps its full-length queries; the documentation explains that a
   query containing a repeat array attracts the whole repeat family, and how to trim a query
   when only orthologs are wanted, as the published surveys do.
2. Run HMMER on every quality-passing hit and record its result as columns; do not use it to decide
   which hits reach RPS-BLAST. Report agreement between the two methods per call, requiring
   positional overlap, under the name "method agreement".

## Consequences

- Positive: the two methods become independent evidence and their agreement means something; no
  annotated SET exon is discarded.
- Negative: about 9% more sequences go to RPS-BLAST (more run time); agreement figures drop from
  near-total to an honest value.
- Neutral: the config key `hmmer.enabled` then switches an annotation, not a filter.

## Alternatives considered

- **Keep the gate and rename the metric**: cheaper, but the metric stays uninformative.
- **Drop HMMER**: loses a second line of evidence that costs little.

## Revisit trigger

A use case where the HMMER gate's reduction in run time matters more than independence.

## References

- Workbench `validation/out/funnel_*.txt`, `tblastn_cutoff_*.txt`, `literature.md`.
