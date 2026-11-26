"""The launcher's stages in one table, and the command-line parser built from it.

Every stage flag of `./rpsHunter` is one row of STAGES: its phase, the Snakemake
targets it asks for, its help line, the heavy rules it is allowed to run and
the tools it needs. The parser, the run order, the heavy-rule guard and the
tool check are all derived from the table, so they cannot drift apart.

This module is pure: it imports nothing from the pipeline, so `./rpsHunter -h`
works without a config.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class Stage:
    """One launcher stage.

    Attributes:
        flag: The command-line flag, e.g. "--search".
        phase: The help section it is listed under (a key of PHASES).
        targets: Snakemake target rules the stage asks for.
        help: One line for `./rpsHunter -h`.
        heavy: Heavy rules this stage exists to run. A run that would start one
            of them without its stage being requested is stopped by the guard.
        tools: Executables the stage calls, checked before a run starts.
    """

    flag: str
    phase: str
    targets: tuple[str, ...]
    help: str
    heavy: tuple[str, ...] = ()
    tools: tuple[str, ...] = ()

    @property
    def dest(self) -> str:
        """The attribute argparse stores the flag under ("--blast-dbs" -> "blast_dbs")."""
        return self.flag.lstrip("-").replace("-", "_")


PHASES: dict[str, str] = {
    "Setup": "Fetch inputs. Run once; never started by another stage.",
    "Search": "Find the regions of each genome that resemble the queries.",
    "Domains": "Describe the protein domains in those regions.",
    "Tables": "Put exons back together into loci.",
}

STAGES: tuple[Stage, ...] = (
    Stage(
        "--download-databases",
        "Setup",
        ("setup_databases", "pfam_hmm_downloader"),
        "Download CDD, the rpsbproc annotation files and Pfam (about 8 GB, 45 GB unpacked).",
        heavy=(
            "cdd_database_downloader",
            "rpsbproc_data_downloader",
            "cdd_smp_downloader",
            "pfam_hmm_downloader",
        ),
        tools=("wget",),
    ),
    Stage(
        "--download-genomes",
        "Setup",
        ("genome_downloader",),
        "Download the genomes pinned by assembly accession under `species:`.",
        heavy=("genome_downloader_setup",),
        tools=("datasets", "jq", "unzip"),
    ),
    Stage(
        "--download-queries",
        "Setup",
        ("query_downloader",),
        "Download the query proteins listed under `queries:` from NCBI.",
        heavy=("query_downloader_setup",),
        tools=("efetch",),
    ),
    Stage(
        "--blast-dbs",
        "Search",
        ("blast_db_generator",),
        "Verify each genome's checksum and build its BLAST database.",
        heavy=("blast_db_generator_setup",),
        tools=("makeblastdb",),
    ),
    Stage(
        "--search",
        "Search",
        ("blast_parser",),
        "Run tBLASTn with every query on every genome and select hits by quality.",
        heavy=("blast_species",),
        tools=("tblastn", "blast_formatter"),
    ),
    Stage(
        "--hmmer",
        "Domains",
        ("hmmer",),
        "Search the selected hits with the families' Pfam profiles (independent evidence).",
        tools=("hmmsearch", "hmmfetch"),
    ),
    Stage(
        "--domains",
        "Domains",
        ("rpsbproc_parser",),
        "Call domains with RPS-BLAST and rpsbproc; write domain tables and GFF3 tracks.",
        tools=("rpsblast", "blast_formatter", "makeprofiledb", "rpsbproc"),
    ),
    Stage(
        "--loci",
        "Tables",
        ("loci",),
        "Chain selected hits into loci and report the families each carries.",
    ),
)

HEAVY_RULES: frozenset[str] = frozenset(
    rule for stage in STAGES for rule in stage.heavy
)
"""Rules that download, or cost minutes per genome. Never started as a side effect."""

_INTO_DOMAINS = "merged into --domains"
_NO_CHECKS_TO_SKIP = (
    "removed: the checks before a run no longer prompt or use the network"
)

RETIRED: dict[str, str] = {
    "--setup-databases": "renamed: use --download-databases",
    "--download-query": "renamed: use --download-queries",
    "--blast": "renamed: use --search",
    "--rpsblast": _INTO_DOMAINS,
    "--rpsbproc": _INTO_DOMAINS,
    "--rpsbproc-parser": _INTO_DOMAINS,
    "--orf-analyser": "removed: the ORF stage was retired (docs/adr/ADR-011)",
    "--completeness-detector": "removed: see tables/loci (docs/outputs.md)",
    "--contingency-parser": "removed: GFF3 tracks are written by --domains",
    "--hit-domain-inventory": "removed: use --loci",
    "--concordance": "removed: Method_Agreement is written by --domains",
    "--extended-plots": "removed: figures are being rebuilt",
    "--skip-validation": _NO_CHECKS_TO_SKIP,
    "-skp": _NO_CHECKS_TO_SKIP,
}
"""Flags of earlier versions, each with what to do instead."""

VERBOSITIES = ("quiet", "normal", "verbose")

FORBIDDEN_SNAKEMAKE = ("--forceall", "-F")
"""Would rebuild databases and download genomes again; never passed on."""

EXAMPLES = """\
Examples:
  ./rpsHunter --all -n                 list what a full analysis would run
  ./rpsHunter --all                    databases, search, domains and loci
  ./rpsHunter --domains --loci         redo from the domain calls after a config change
  ./rpsHunter --all --config my.yaml   use another config file
  ./rpsHunter --config-help rpsblast   explain a config section

Setup stages (downloads) run only when asked for by name.
Config: --config, else $RPSHUNTER_CONFIG, else data/config/config.yaml."""


@dataclass(frozen=True)
class Request:
    """What one launcher call asks for."""

    stages: tuple[Stage, ...]
    dry_run: bool
    cores: int | None
    config: str | None
    allow_heavy: bool
    verbosity: str | None
    snakemake_args: tuple[str, ...]

    @property
    def targets(self) -> tuple[str, ...]:
        """Snakemake targets of the requested stages, once each, in stage order."""
        return tuple(dict.fromkeys(t for stage in self.stages for t in stage.targets))

    @property
    def allowed_heavy(self) -> frozenset[str]:
        """Heavy rules the requested stages own."""
        return frozenset(rule for stage in self.stages for rule in stage.heavy)

    @property
    def tools(self) -> tuple[str, ...]:
        """Executables the requested stages call, once each."""
        return tuple(dict.fromkeys(t for stage in self.stages for t in stage.tools))


def build_parser() -> argparse.ArgumentParser:
    """Return the launcher's parser: one flag per stage, grouped by phase."""
    parser = argparse.ArgumentParser(
        prog="rpsHunter",
        description="Find regions of genomes that resemble query proteins and "
        "describe the protein domains inside them.",
        epilog=EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
    )
    for phase, summary in PHASES.items():
        group = parser.add_argument_group(phase, summary)
        for stage in (s for s in STAGES if s.phase == phase):
            group.add_argument(stage.flag, action="store_true", help=stage.help)
    run = parser.add_argument_group("Run")
    run.add_argument(
        "--all",
        action="store_true",
        help="Every stage after Setup, in order.",
    )
    run.add_argument(
        "-n", "--dry-run", action="store_true", help="List the jobs; run nothing."
    )
    run.add_argument(
        "--cores", type=int, help="Cores to use (default: execution.num_cores)."
    )
    run.add_argument("--config", help="Config file for this run.")
    run.add_argument(
        "--verbosity",
        choices=VERBOSITIES,
        help="What the screen shows (default: display.verbosity, else normal).",
    )
    run.add_argument(
        "--allow-heavy",
        action="store_true",
        help="Let a stage start a search or database build it does not own.",
    )
    for flag in RETIRED:
        parser.add_argument(flag, action="store_true", help=argparse.SUPPRESS)
    return parser


def _refuse_retired(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Stop with the replacement of every retired flag that was used."""
    used = [
        f"{flag}: {advice}"
        for flag, advice in RETIRED.items()
        if getattr(args, flag.lstrip("-").replace("-", "_"))
    ]
    if used:
        parser.error("flags of earlier versions:\n  " + "\n  ".join(used))


def _split(argv: Sequence[str]) -> tuple[list[str], tuple[str, ...]]:
    """Separate the launcher's arguments from those meant for Snakemake (after `--`)."""
    words = list(argv)
    if "--" not in words:
        return words, ()
    cut = words.index("--")
    return words[:cut], tuple(words[cut + 1 :])


def parse(argv: Sequence[str]) -> Request:
    """Turn a command line into a Request.

    Exits with status 2 and a message when no stage is given, a flag is unknown
    or retired, or a forbidden option is passed on to Snakemake.
    """
    own, passed_on = _split(argv)
    parser = build_parser()
    args = parser.parse_args(own)
    _refuse_retired(parser, args)
    forbidden = [word for word in passed_on if word in FORBIDDEN_SNAKEMAKE]
    if forbidden:
        parser.error(
            f"{forbidden[0]} is never passed to Snakemake: it would rebuild the "
            "databases and download the genomes again"
        )
    chosen = tuple(
        stage
        for stage in STAGES
        if getattr(args, stage.dest) or (args.all and stage.phase != "Setup")
    )
    if not chosen:
        parser.error(
            "no stage given; try --all -n to see what a full analysis would run"
        )
    return Request(
        stages=chosen,
        dry_run=args.dry_run,
        cores=args.cores,
        config=args.config,
        allow_heavy=args.allow_heavy,
        verbosity=args.verbosity,
        snakemake_args=passed_on,
    )
