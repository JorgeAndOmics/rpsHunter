# ADR-004: Tags derived from the hit, not random

- **Status**: Accepted
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: yes

## Context

Each tBLASTn hit carries a Tag that links it to its sequence, its RPS-BLAST alignments and its
domains. Tags were six random characters from A-Z and 0-9, drawn afresh on every run, with nothing
to check that two hits had not drawn the same one. The hit-domain inventory and the concordance
tables join on the Tag alone, across all species and queries.

Two problems follow. No two runs can be compared, because every identifier differs. And with
36^6 (about 2.2 billion) possible Tags, n hits are expected to produce about n^2 / (2 x 36^6)
colliding pairs: six for the 166,167 hits of the model run.

## Evidence

Two runs of the unchanged pipeline (v0.9.0) on the five model genomes with three PRDM9 queries,
same environment and databases, 166,167 hits each:

| | Run 1 | Run 2 |
|---|---|---|
| Colliding Tag pairs | 5 | 10 |
| Collisions within one query and species | 0 | 0 |

Every collision paired hits from different species or queries. Comparing the two runs by content
with Tags masked: BLAST, RPS-BLAST, domain, contingency and GFF3 outputs were identical, while 37
inventory and concordance files differed. In those, a colliding hit carried the domains and scores
of the unrelated hit it shared a Tag with: domain lists and counts on a few dozen rows, BLAST
scores on 8 concordance rows, and one summary statistic (median HMMER E-value for mouse zf-C2H2,
4.95e-07 against 5.00e-07).

Of the 166,167 hits, 80,091 share their subject coordinates with another hit of the same query and
species (the same genome span aligned to different parts of the query), and none are identical once
query coordinates are included.

## Decision

A Tag is the first 12 base32 characters of a SHA-256 hash of the hit's identity: query accession,
species, sequence identifier, subject start and end as BLAST reports them (so the strand is
included), and query start and end. Duplicate rows are dropped before tagging, and tagging fails
loudly if two rows describe the same hit or if two different hits receive the same Tag. The logic
lives in `workflow/scripts/tags.py`, a pure module. The config key `execution.random_id_length` is
removed.

## Consequences

- Positive: the same hit has the same Tag on every run and machine, so runs can be compared row by
  row; joins on Tag can no longer attach one hit's data to another; a collision stops the run
  instead of corrupting it silently.
- Negative: Tags are 12 characters instead of 6, and Tags from earlier runs do not correspond to
  the new ones.
- Neutral: with 60 bits, the chance of any collision among a million hits is about one in a
  million; the check makes the remaining case an error, not a wrong result.

## Alternatives considered

- **Sequential numbers per species**: short and readable, but a hit's number changes whenever a
  threshold or database change adds or removes an earlier hit, so Tags would not be comparable
  between runs with different settings.
- **Seeded random Tags**: reproducible only while the number and order of hits stay the same, and
  still open to collisions.
- **Identity without query coordinates**: 80,091 hits would be indistinguishable from another hit.

## Revisit trigger

`TagCollisionError` is ever raised (lengthen `TAG_LENGTH`), or the definition of a hit changes, for
example if overlapping hits are merged before tagging.

## References

- `workflow/scripts/tags.py`, `tests/unit/test_tags.py`.
- ADR-002 (quality gate): `tags.py` is the first module outside the legacy list.
