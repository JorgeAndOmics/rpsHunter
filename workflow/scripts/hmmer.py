"""Search one species' selected hit sequences with the families' Pfam profiles.

HMMER is a second, independent look at the same sequences RPS-BLAST searches.
It decides nothing: every domain it reports is written to a table with its own
position on the genome, and whether the two methods agree is worked out later
from positions (ADR-010). Each domain says whether it clears Pfam's curated
gathering threshold (Above_GA), so permissive search settings stay visible.

Usage (called by the Snakefile):
    python hmmer.py --species KEY --query-accession ACC --fasta FASTA \
        --profiles SUBSET.hmm --output TABLE.parquet --threads N

A FASTA with no sequence gives an empty table without running the search. A
failing hmmsearch stops the job (external.ToolError).
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import NamedTuple

import pandas as pd

import defaults
import domain_calls
import domain_families
import external
import log

RAW_COLUMNS = (
    "Tag",
    "Sequence",
    "Hit_A",
    "Hit_B",
    "Profile",
    "Pfam",
    "Residue_From",
    "Residue_To",
    "Model_From",
    "Model_To",
    "Model_Length",
    "Evalue",
    "Score",
)
"""One row per HMMER domain. Hit_A and Hit_B are the hit's genomic ends as
tBLASTn wrote them; residues count along the hit's translation."""

DOMAIN_COLUMNS = (
    "Species",
    "Query_Accession",
    "Tag",
    "Family",
    "Profile",
    "Pfam",
    "Sequence",
    "Strand",
    "Start",
    "End",
    "Hit_Start",
    "Hit_End",
    "Residue_From",
    "Residue_To",
    "Model_From",
    "Model_To",
    "Model_Length",
    "Evalue",
    "Score",
    "Above_GA",
)

_INTEGERS = (
    "Hit_A",
    "Hit_B",
    "Residue_From",
    "Residue_To",
    "Model_From",
    "Model_To",
    "Model_Length",
)
_FLOATS = ("Evalue", "Score")
# Fields of hmmsearch --domtblout, by position (HMMER 3.4 user guide, "tabular output").
_TARGET, _PROFILE, _ACCESSION, _MODEL_LENGTH = 0, 3, 4, 5
_I_EVALUE, _SCORE, _HMM_FROM, _HMM_TO, _ALI_FROM, _ALI_TO = 12, 13, 15, 16, 17, 18
_FIELDS = 22


class SearchSettings(NamedTuple):
    """The hmmsearch options that come from the config."""

    use_gathering_threshold: bool
    evalue: float
    dom_evalue: float
    min_score: float | None
    max_sensitivity: bool
    bias_filter: bool
    seed: int
    threads: int


def _reporting_flags(settings: SearchSettings) -> list[str]:
    """What hmmsearch reports: Pfam's gathering thresholds, or explicit cut-offs."""
    if settings.use_gathering_threshold:
        return ["--cut_ga"]
    flags = ["-E", str(settings.evalue), "--domE", str(settings.dom_evalue)]
    if settings.min_score is not None:
        flags += ["-T", str(settings.min_score)]
    return flags


def search_command(
    profiles: Path, fasta: Path, domain_table: Path, settings: SearchSettings
) -> list[str]:
    """Return the hmmsearch command that writes a per-domain table."""
    command = ["hmmsearch", *_reporting_flags(settings)]
    if settings.max_sensitivity:
        command.append("--max")
    if not settings.bias_filter:
        command.append("--nobias")
    return [
        *command,
        "--cpu", str(settings.threads),
        "--seed", str(settings.seed),
        "--domtblout", str(domain_table),
        str(profiles),
        str(fasta),
    ]  # fmt: skip


def gathering_thresholds(profiles_text: str) -> dict[str, float]:
    """Return each profile's per-domain gathering threshold (bit score) by name.

    A profile's `GA` line holds the sequence and the domain threshold, in that
    order; a domain is compared with the second.
    """
    thresholds: dict[str, float] = {}
    name = ""
    for line in profiles_text.splitlines():
        if line.startswith("NAME"):
            name = line.split()[1]
        elif line.startswith("GA"):
            thresholds[name] = float(line.split()[2].rstrip(";"))
    return thresholds


def _row(fields: list[str]) -> dict[str, object]:
    """One domain from the fields of a domain-table line."""
    target = fields[_TARGET]
    hit = domain_calls.DEFINITION.fullmatch(target)
    if hit is None:
        raise ValueError(f"target name is not `<seq>:<a>-<b>|tag:<Tag>`: {target}")
    return {
        "Tag": hit["tag"],
        "Sequence": hit["sequence"],
        "Hit_A": hit["a"],
        "Hit_B": hit["b"],
        "Profile": fields[_PROFILE],
        "Pfam": fields[_ACCESSION].split(".")[0],
        "Residue_From": fields[_ALI_FROM],
        "Residue_To": fields[_ALI_TO],
        "Model_From": fields[_HMM_FROM],
        "Model_To": fields[_HMM_TO],
        "Model_Length": fields[_MODEL_LENGTH],
        "Evalue": fields[_I_EVALUE],
        "Score": fields[_SCORE],
    }


def parse_domains(text: str) -> pd.DataFrame:
    """Parse hmmsearch's per-domain table into a typed table with RAW_COLUMNS.

    The E-value and score are the domain's own (independent E-value), not the
    whole sequence's.

    Raises:
        ValueError: A line has too few fields, or its target is not named
            `<seq>:<a>-<b>|tag:<Tag>`.
    """
    lines = [line.split() for line in text.splitlines() if line and line[0] != "#"]
    short = [fields for fields in lines if len(fields) < _FIELDS]
    if short:
        raise ValueError(f"expected at least {_FIELDS} fields, got {len(short[0])}")
    table = pd.DataFrame([_row(fields) for fields in lines], columns=list(RAW_COLUMNS))
    return table.astype(
        dict.fromkeys(_INTEGERS, "int64") | dict.fromkeys(_FLOATS, "float64")
    )


def keep_substantial(
    domains: pd.DataFrame, min_coverage: float, min_length: int
) -> pd.DataFrame:
    """Keep domains covering enough of their model over enough residues."""
    covered = domains["Model_To"] - domains["Model_From"] + 1
    length = domains["Residue_To"] - domains["Residue_From"] + 1
    keep = (covered / domains["Model_Length"] >= min_coverage) & (length >= min_length)
    return domains[keep].reset_index(drop=True)


def describe(
    domains: pd.DataFrame,
    family_of: Mapping[str, str],
    thresholds: Mapping[str, float],
) -> pd.DataFrame:
    """Add each domain's family, genomic position and gathering-threshold verdict.

    Raises:
        ValueError: A domain's profile belongs to no configured family, so the
            profile library and the config do not match.
    """
    unknown = sorted(set(domains["Pfam"]) - set(family_of))
    if unknown:
        raise ValueError(f"profiles of no configured family: {', '.join(unknown)}")
    placed = domain_calls.place_on_genome(domains)
    return placed.assign(
        Family=placed["Pfam"].map(family_of),
        Above_GA=placed["Score"] >= placed["Profile"].map(thresholds),
    )


def configured_settings(threads: int) -> SearchSettings:
    """The search settings of the current config."""
    return SearchSettings(
        use_gathering_threshold=defaults.HMMER_USE_GA,
        evalue=defaults.HMMER_EVALUE,
        dom_evalue=defaults.HMMER_DOM_EVALUE,
        min_score=defaults.HMMER_MIN_SCORE,
        max_sensitivity=defaults.HMMER_MAX_SENS,
        bias_filter=defaults.HMMER_BIAS_FILTER,
        seed=defaults.HMMER_SEED,
        threads=threads,
    )


def search(fasta: Path, profiles: Path, settings: SearchSettings) -> pd.DataFrame:
    """Run the search for one FASTA and return its domains (RAW_COLUMNS)."""
    if fasta.stat().st_size == 0:
        return parse_domains("")
    with tempfile.TemporaryDirectory() as folder:
        domain_table = Path(folder) / "domains.tbl"
        command = search_command(profiles, fasta, domain_table, settings)
        external.run_tool(command, tool="hmmsearch")
        return parse_domains(domain_table.read_text())


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--species", required=True)
    parser.add_argument("--query-accession", required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    log.add_scope_argument(parser)
    args = parser.parse_args(argv)
    log.job_logging("hmmer", args.scope)

    found = search(args.fasta, args.profiles, configured_settings(args.threads))
    kept = keep_substantial(
        found, defaults.HMMER_MIN_COVERAGE, defaults.HMMER_MIN_ALN_LEN
    )
    family_of = domain_families.family_by_pfam(defaults.FAMILIES)
    table = describe(kept, family_of, gathering_thresholds(args.profiles.read_text()))
    table = table.assign(Species=args.species, Query_Accession=args.query_accession)
    table[list(DOMAIN_COLUMNS)].to_parquet(args.output, index=False)
    log.ok(
        f"{log.count(len(table))} HMMER domains ({log.count(len(found))} before filters)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(log.run_main(main))
