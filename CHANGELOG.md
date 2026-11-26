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
- `CITATION.cff`, this changelog, and architectural decision records under `docs/adr/`.

### Changed

- **Tags are derived from the hit instead of being random**
  ([ADR-004](docs/adr/ADR-004-content-derived-tags.md)): 12 characters, identical on every run.
  Random 6-character Tags collided (5 to 10 pairs per 166,167 hits), and each collision attached
  another hit's domains and scores to a hit in the inventory and concordance tables.
  The config key `execution.random_id_length` is removed.
- License is MIT ([ADR-003](docs/adr/ADR-003-mit-licence.md)). The `LICENSE` file previously held
  the GPL v2 text while the README said MIT.
- All text files are stored with LF line endings.
- The environment no longer references Anaconda's `pkgs/main` and `pkgs/r` channels and now lists
  the tools the scripts call: `entrez-direct`, `jq`, `pandoc`, `unzip`, `wget`.
- Supported platform stated as Linux (including Windows through WSL); the macOS claim is withdrawn.

### Removed

- IDE settings (`.idea/`) and two scratch R scripts under `workflow/scripts/tests/` that held
  machine-specific paths and no assertions.

## [0.9.0] - 2026-11-26

State of the pipeline as developed up to March 2026, tagged as the baseline before the
modernisation work. Multi-query tBLASTn, optional ORF and HMMER enrichment, RPS-BLAST against a CDD
subset, `rpsbproc` parsing, cross-query merging, concordance and hit-domain inventory tables, and
ten figures. Not reproducible as tagged: reference databases and most tools are unpinned.

[Unreleased]: https://github.com/JorgeAndOmics/rpsHunter/compare/v0.9.0...HEAD
[0.9.0]: https://github.com/JorgeAndOmics/rpsHunter/releases/tag/v0.9.0
