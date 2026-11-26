# ADR-011: Retire the ORF gate; flag internal stop codons instead

- **Status**: Proposed
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: no (the stage is disabled by default)

## Context

An optional stage keeps only hits that contain an open reading frame found by ORFik, with
`orf.min_orf_length` documented as nucleotides and an empty `orf.start_codon` documented as "any
codon". It is disabled in the shipped config because it removes short domains.

## Evidence

- ORFik counts `minimumLength` in codons: a 7-codon ORF is found with `minimumLength = 7` and not
  with 8 (checked in the pinned ORFik 1.26.0). The configured 50 therefore means 150 bases.
- With no start codon given, ORFik uses `ATG|TTG|CTG`, not any codon.
- A tBLASTn hit covers one exon. Internal exons have no start codon, so requiring an ORF with a
  start inside the hit discards them by construction.
- Among gate-passing hits, 31.6% (human), 38.9% (mouse) and 22.1% (*Desmodus*) contain an internal
  stop codon: a useful pseudogene signal that the pipeline does not report.

## Decision

Remove the ORF stage, its config section and the ORFik dependency. Add a per-hit column stating
whether the translated hit contains an internal stop codon, and carry it into the locus table.

## Consequences

- Positive: one less stage, one less heavy dependency; a simple, correct pseudogene indicator.
- Negative: users of `orf.enabled: true` get a retired-key message instead of the old behaviour.
- Neutral: documentation of the three ORF keys is deleted.

## Alternatives considered

- **Fix the units and start-codon handling and keep the gate**: still wrong for internal exons.

## Revisit trigger

The search stage produces full coding sequences, for which an ORF check is meaningful.

## References

- `workflow/scripts/orf_analyser.R`; ORFik `findORFs` documentation.
