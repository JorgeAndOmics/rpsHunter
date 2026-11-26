# ADR-012: Domain families defined in the config by CDD and Pfam accession

- **Status**: Proposed
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: yes

## Context

Target domains are two lists of names: `rpsblast.target_domains` (matched against CDD short names,
exactly or as a prefix followed by an underscore) and `hmmer.profiles`. Family names are also
written into the plotting and concordance code.

## Evidence

- CDD short names are not stable: COG5048 was renamed `ZMS1`, and the shipped config stopped
  resolving against CDD 3.22.
- The prefix rule selects 121 models for seven names: `SET` brings in dozens of `SET_*` models, and
  misses `PR-SET_PRDM7_9`, the model NCBI uses to annotate PRDM9's SET domain.
- The same domain is reported under many names (the baseline lists 20 zinc-finger model names),
  which the code then regroups with hard-coded patterns.
- A superfamily not asked for (`DUF4371`) appears in the output.

## Decision

Replace the two lists with one `domains:` block. Each family has a name, the CDD accessions and the
Pfam accessions that count as evidence for it, and optionally whether it is a repeat. The subset
database is built from exactly those accessions, pinned to the CDD release; every table carries the
family beside the model; no family name appears in code outside the config and the tests.

## Consequences

- Positive: a config keeps working across CDD releases; the family a call belongs to is explicit;
  the tool is usable for other protein families without code changes.
- Negative: the user lists accessions, not names; the shipped example grows.
- Neutral: the PRDM9 example includes cd19193 (see ADR-009 for why it rarely fires on real genes).

## Alternatives considered

- **Keep names and add a rename table**: fragile for the same reason.
- **Whole superfamilies by cluster id**: simple, but pulls in every member model.

## Revisit trigger

A need to search full CDD without a target list.

## References

- Baseline run tables; NCBI GenPept annotation of NP_064612.2 (PR-SET_PRDM7_9, CDD:380970).
