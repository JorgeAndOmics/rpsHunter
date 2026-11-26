# ADR-006: RPS-BLAST E-values on the full-CDD scale, with a cut-off of 1

- **Status**: Proposed
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: yes

## Context

RPS-BLAST runs against a subset of CDD (117 profiles for the PRDM9 example) with `-evalue 1000`. An
E-value depends on database size, so against the subset every E-value is smaller than NCBI's
CD-Search would report for the same alignment, and it changes whenever the list of target domains
changes. The config comment justifies 1000 by the need "to capture mouse SSXRD".

## Evidence

Measured on the five model genomes (three PRDM9 queries; CDD 3.22, Pfam 38.2, BLAST 2.17.0) against a
truth set built from RefSeq annotation of human, mouse and *Desmodus rotundus*: Pfam gathering
thresholds on the annotated proteins, mapped through exons to genome coordinates.

- The same 12,576 alignments searched against the subset and against full CDD have identical bit
  scores and E-values that differ by a median factor of 1,051 (range 709 to 1,600): the ratio of
  the databases' sizes, 15,821,058 to 15,051 residues. The configured 1000 is about 10^6 on NCBI's
  scale.
- NCBI: CD-Search default is 0.01; "results with E-values in the range of 1 and above should be
  considered putative false positives".
- Human, calls kept at the configured cut-off: 9,863 SET calls of which 31 lie on an annotated SET
  domain (0.3%); 111 loci carry "KRAB + SET". At a subset E-value of 1e-3 (about 1 on NCBI's scale):
  263 SET calls reduce to 42 at 0.01, and loci with KRAB + SET fall to 2 (PRDM9 and PRDM7).
- Sensitivity does not change: the true domains recovered in human are 90 KRAB, 2 SSXRD, 9 SET and
  5,469 zinc fingers at every cut-off from 1000 down to a subset E-value of 0.01.
- Mouse Prdm9 keeps its SSXRD down to a subset E-value of 1e-5; the stated reason for 1000 does not
  hold with current data.
- Loci carrying all four families, by NCBI-scale cut-off (loci rebuilt as in ADR-008):

  | Species | Today | E <= 1 | E <= 0.01 |
  |---|---|---|---|
  | Human | 2 | 2 | 1 (PRDM7 lost) |
  | Mouse | 2 (one false) | 1 | 1 |
  | *Desmodus* | 16 | 12 | 4 |
  | *Antrozous* | 3 | 2 | 0 |
  | *Molossus* | 3 | 1 | 0 |

- Diverged but genuine domains sit between the two: seven annotated *Desmodus* PRDM9-like proteins
  have an N-terminal KRAB with HMMER E-values of 1e-4 to 1e-6, below Pfam's gathering threshold.

## Decision

Run RPS-BLAST with `-dbsize` set to the residue count of the full CDD release in use, so that
E-values are on NCBI's scale and do not depend on the target list. Default `rpsblast.e_value` to 1,
the upper end of what NCBI treats as borderline, and report for each call whether it is at or below
0.01 (NCBI's default) or between 0.01 and 1 (borderline). Pass the same cut-off to `rpsbproc`.

## Consequences

- Positive: E-values mean what they mean in CD-Search and in the literature; noise calls (99.7% of
  SET calls) disappear; the false complete locus in mouse disappears; no annotated domain is lost.
- Negative: fewer calls than before in every table; results are not comparable with earlier runs.
- Neutral: a user who needs a rescue search (the literature uses 100 to 1000 only for SSXRD inside an
  already confirmed ortholog) can still raise the cut-off in the config.

## Alternatives considered

- **Cut-off 0.01**: NCBI's default; loses PRDM7's completeness in human and most diverged bat copies.
- **Keep 1000**: no sensitivity gain measured, and two orders of magnitude more false calls.
- **Filter on bit score or on Specific hits only**: independent of database size, but specific hits
  are rare for domains split across exons (ADR-009).

## Revisit trigger

A CDD release changes the full database size (update `-dbsize` with the pin), or a family is found
whose genuine members score above 1.

## References

- Workbench `validation/out/`: `scores_*.tsv`, `threshold_effect.txt`, `proposal.txt`.
- https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml ; Cavassim 2022, PNAS; Baker 2017, eLife.
