# ADR-005: Genomes pinned by assembly accession or checksum

- **Status**: Accepted
- **Date**: 2026-11-26
- **Deciders**: Jorge Gonzalez Garcia
- **Affects results**: no

## Context

A species was a key and a display name. To download its genome the pipeline passed the key to NCBI
as a taxon name, took the first assembly of the best available level, and always requested the
latest version. Which sequences were analysed therefore depended on the day of the download, and
nothing recorded the answer beyond a line with the accession: no version check, no checksum.

## Evidence

The five model genomes were on disk without any record of their assembly. Recovering it showed how
easily a pin goes wrong: looking up the first chromosome of *Antrozous pallidus* returned assembly
GCA_027563665.2, but the file is version .1. The two versions share their chromosomes and differ in
65 unplaced scaffolds (93 sequences against 88). Only comparing every sequence accession in the
file with NCBI's sequence report for each version identified it. The human file holds the 194
sequences of the primary assembly, while its assembly, GCA_000001405.29, lists 709 including
alternate loci and patches.

## Decision

Every species entry carries a display name and a pin: a versioned NCBI assembly accession
(`assembly`), the SHA-256 of its FASTA (`sha256`), or both. At least one is required. The downloader
takes the accession from the config and fetches exactly that version; it accepts no taxon name and
no accession without a version, checks that the archive holds the assembly requested, and records
the accession, the FASTA's checksum, the time and the tool version. A species pinned by checksum
only is never downloaded: its FASTA is provided by the user. A checksum pin is compared with the
file just before that genome's BLAST database is built (`workflow/scripts/genome_check.py`), the
one step every genome passes through once, and a mismatch stops the run. Validation of the config
block lives in `workflow/scripts/species_config.py`, a pure module.

## Consequences

- Positive: a config determines the genomes; two people running it months apart analyse the same
  sequences, or find out that they do not.
- Negative: adding a species means looking up its assembly accession; the convenience of typing a
  species name is gone. Configs written in the old form must be edited (the error says how).
- Neutral: an accession pins the whole assembly. A file that holds only part of one, such as a
  primary assembly without alternate loci, is fully described only by its checksum.

## Alternatives considered

- **Accession optional, resolved one recorded**: friendlier, but the same config would still mean
  different genomes at different times.
- **A separate manifest file listing the accessions**: keeps the config short, at the cost of two
  files that must agree.
- **Pin by taxon and date**: NCBI offers no way to ask for "the best assembly as of a date".

## Revisit trigger

A need to analyse a defined subset of an assembly (for example chromosomes only) reproducibly from
a fresh download, which a whole-assembly accession cannot express.

## References

- `workflow/scripts/species_config.py`, `workflow/scripts/genome_downloader.sh`,
  `docs/configuration.md` (species).
