# Changelog

All notable changes to rpsHunter are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `make check`: one quality gate running ruff, mypy, lintr, ShellCheck, pytest and testthat
  ([ADR-002](docs/adr/ADR-002-quality-gate-legacy-list.md)).
- Explicit conda lock, `data/config/environment.lock-linux-64.txt`, and `make env` to build from it
  ([ADR-001](docs/adr/ADR-001-single-locked-conda-env.md)).
- `RPSHUNTER_CONFIG` environment variable to run with a config other than
  `data/config/config.yaml`.
- Every job writes its output to `<logs>/jobs/<rule>/<query>/<species>.log` and a Snakemake
  benchmark (wall time, memory) to `<logs>/benchmarks/...`; the log of a failed job is shown.
- An end-to-end test on synthetic genomes with planted domains and a miniature database
  (`tests/integration/test_toy_pipeline.py`), run in CI on every push.
- `CITATION.cff`, this changelog, and architectural decision records under `docs/adr/`.

### Fixed

- **Stop codons no longer shift domain positions.** A stop codon inside a hit's translation was
  deleted before the domain searches, moving every later residue, and the genomic position of
  every domain after it, by three bases per stop. It is now masked with `X`
  ([ADR-007](docs/adr/ADR-007-domain-coordinates.md), implementation note).
- `hmmer.min_score` passed hmmsearch a flag that does not exist (`--score`); it now uses `-T`.
- A failed `tblastn`, `rpsblast` or `blast_formatter` now fails its job. Previously the error
  was logged and an empty table written, so a crashed search looked like a genome without hits.
- Changing a threshold in the config now reruns the rules that depend on it. The scripts read
  their settings directly, which Snakemake could not see, so nothing reran.
- A species with no sequence left after the gates no longer triggers a hidden `blast_formatter`
  failure in the RPS-BLAST step.

### Changed

- **Tags are derived from the hit instead of being random**
  ([ADR-004](docs/adr/ADR-004-content-derived-tags.md)): 12 characters, identical on every run.
  Random 6-character Tags collided (5 to 10 pairs per 166,167 hits), and each collision attached
  another hit's domains and scores to a hit in the inventory and concordance tables.
  The config key `execution.random_id_length` is removed.
- License is MIT ([ADR-003](docs/adr/ADR-003-mit-licence.md)). The `LICENSE` file previously held
  the GPL v2 text while the README said MIT.
- **HMMER is independent evidence, not a gate**
  ([ADR-010](docs/adr/ADR-010-scope-and-hmmer-evidence.md)). Hits are selected by the quality gate
  alone and all of them reach RPS-BLAST; HMMER searches the same sequences and writes its own
  domains to `tables/domains/hmmer_domains.*` with genomic positions and `Above_GA`. Calls and
  HMMER domains carry `Method_Agreement` (`both`, `cdd_only` / `hmmer_only`, `not_compared`),
  which requires the same family on the same hit at an overlapping position. The hit tables no
  longer have `HMM_*` columns. A consequence: a locus made of chance hits is now listed, and is
  told apart by carrying no family (`N_Families` 0).
- **Loci are chained from gate-passing hits within 20 kb and report what they carry**
  ([ADR-008](docs/adr/ADR-008-loci.md), [ADR-009](docs/adr/ADR-009-completeness-across-exons.md)).
  `tables/loci/loci.*` lists each locus with the families present; `tables/loci/locus_families.*`
  gives, per locus and family, the evidence tier, the fraction of a model covered when the locus's
  exon fragments are put together, and the number of copies. There is no "complete" flag: loci used
  to pool neighbouring genes over 50 kb (one false complete locus in mouse, four tandem genes
  counted as one in *Desmodus*), and completeness was read from per-hit flags that follow exon
  boundaries. Config key `hit_domain_inventory.locus_gap` is replaced by `loci.gap` (default 20000);
  launcher flag `--loci`.
- **Each domain is reported at its own genomic position**
  ([ADR-007](docs/adr/ADR-007-domain-coordinates.md)), computed from where it lies in the translated
  hit; the hit's span is kept in `Hit_Start` and `Hit_End`. Positions used to be those of the whole
  BLAST hit (median error for zinc fingers 411 bp in human; 48 bp with the domain's own position).
- **New domain tables.** `tables/domains/calls.*` holds every call with its family, evidence tier,
  coordinates and the part of the model it covers; `tables/domains/domains.*` holds domain
  instances (overlapping calls of one model on one strand, counted once). Numbers are stored as
  numbers. `ranges/<species>.gff3` is written for every species, also one without domains. The
  per-query domain tables and the R cross-query merge are gone.
- **RPS-BLAST E-values are on NCBI's scale, with a cut-off of 1**
  ([ADR-006](docs/adr/ADR-006-rpsblast-evalue-scale-and-cutoff.md)). The search is told the size of
  the full CDD (`databases.cdd.residues`), `rpsblast.e_value` defaults to 1 instead of 1000, the same
  cut-off goes to rpsbproc, and each call carries an `Evidence` tier (`standard`, `borderline`,
  `rescue`). The previous 1000 against the subset database was about a million on NCBI's scale: in
  human, 9,863 SET calls of which 31 lay on an annotated SET domain. No annotated domain is lost.
- **Domain families are defined by accession** ([ADR-012](docs/adr/ADR-012-domain-families-by-accession.md)):
  one `domains:` block lists, per family, the CDD and Pfam accessions that evidence it, replacing
  `rpsblast.target_domains` and `hmmer.profiles` (lists of names, one of which CDD had renamed).
  The RPS-BLAST database is built from exactly those models; domain tables gain a `Family` column.
  The example now searches 12 CDD models instead of the 121 the name rule selected, and includes
  the PRDM9-specific SET model that the rule missed.
- **Species entries name their genome** ([ADR-005](docs/adr/ADR-005-genomes-pinned-by-assembly.md)):
  each species is `name:` plus a versioned `assembly:` accession and/or a `sha256:` of its FASTA.
  The genome downloader fetches exactly that assembly; taxon names and "latest" are gone. The
  former `'Key': 'Display name'` form is refused with a message showing the new one.
- An NCBI API key is read from the `NCBI_API_KEY` environment variable only, and is no longer
  placed on a job's command line.
- All text files are stored with LF line endings.
- The environment no longer references Anaconda's `pkgs/main` and `pkgs/r` channels and now lists
  the tools the scripts call: `entrez-direct`, `jq`, `pandoc`, `unzip`, `wget`.
- Supported platform stated as Linux (including Windows through WSL); the macOS claim is withdrawn.

### Removed

- **The ORF stage** ([ADR-011](docs/adr/ADR-011-retire-orf-gate.md)): `orf_analyser.R`, the `orf`
  config section, the `--orf-analyser` flag and the ORFik dependency. It demanded an open reading
  frame with a start codon inside every hit, and a hit covers one exon; its minimum length was
  counted in codons while documented as nucleotides. In its place every hit has an
  `Internal_Stop` column (a stop codon inside its translation) and every locus counts such hits
  in `N_Hits_With_Stop`.
- The stages downstream of domain parsing, pending their rebuild on ADR-006 to ADR-012: the
  contingency table and GFF3 export, the hit-domain and locus inventories, the HMMER-versus-CDD
  concordance tables, and the tile, 3D and extended plots, with their launcher flags
  (`--completeness-detector`, `--contingency-parser`, `--hit-domain-inventory`, `--concordance`,
  `--extended-plots`).
- Config keys that nothing read: the `logging` section, `execution.retrieval_time_lag`,
  `execution.max_retrieval_attempts`, `display.display_requests_warning`,
  `display.display_operation_info`; and `execution.use_species_list`, whose "discover genomes from
  the folder" mode cannot be pinned.
- IDE settings (`.idea/`) and two scratch R scripts under `workflow/scripts/tests/` that held
  machine-specific paths and no assertions.

## [0.9.0] - 2026-11-26

State of the pipeline as developed up to March 2026, tagged as the baseline before the
modernisation work. Multi-query tBLASTn, optional ORF and HMMER enrichment, RPS-BLAST against a CDD
subset, `rpsbproc` parsing, cross-query merging, concordance and hit-domain inventory tables, and
ten figures. Not reproducible as tagged: reference databases and most tools are unpinned.

[Unreleased]: https://github.com/JorgeAndOmics/rpsHunter/compare/v0.9.0...HEAD
[0.9.0]: https://github.com/JorgeAndOmics/rpsHunter/releases/tag/v0.9.0
