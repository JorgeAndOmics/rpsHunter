"""End-to-end run of the pipeline on the toy genomes, checked against planted truth.

The toy genomes (tests/fixtures/toy/build_toy_genome.py) carry copies of the
human PRDM9 protein at known positions, and the miniature databases under
tests/fixtures/toy/db hold just the profiles that recognise its domains. One
module-scoped run feeds every test here, because the run takes about a minute.

Needs the external tools (blast, hmmer, rpsbproc, snakemake, R) on PATH, so the
tests carry the `integration` marker.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pytest
import yaml

pytestmark = pytest.mark.integration

REPO = Path(__file__).resolve().parents[2]
TOY = REPO / "tests" / "fixtures" / "toy"
PLANTED_PROTEIN = "NP_064612.2"

QUERIES = {"NP_064612.2": "Human PRDM9", "NP_659058.3": "Mouse PRDM9"}
QUERY_LABELS = ("Human_PRDM9", "Mouse_PRDM9")
SPECIES = {
    "Toyus_plenus": "Toyus plenus",
    "Toyus_fractus": "Toyus fractus",
    "Toyus_vacuus": "Toyus vacuus",
}
TARGET_DOMAINS = ["KRAB", "SSXRD", "PR-SET", "SET", "zf-C2H2", "zf-H2C2", "ZMS1"]

# Domain names (CDD short names) that count as evidence for each family.
FAMILY_MODELS = {
    "KRAB": {"KRAB", "KRAB_A-box"},
    "SSXRD": {"SSXRD"},
    "SET": {"PR-SET_PRDM7_9", "SET"},
    "zinc_fingers": {"zf-C2H2", "zf-H2C2_2", "ZMS1"},
}

# Launcher flags in pipeline order (the v0.9.0 launcher runs one stage per call).
STAGES = (
    "blast-dbs",
    "hmmer",
    "blast",
    "rpsblast",
    "rpsbproc",
    "rpsbproc-parser",
    "contingency-parser",
    "hit-domain-inventory",
    "concordance",
)


@dataclass(frozen=True)
class ToyRun:
    """Where a finished toy run left its inputs and outputs."""

    results: Path
    truth: pd.DataFrame
    query_domains: pd.DataFrame


def _load_generator():
    """Import the toy genome generator, which lives beside the fixtures."""
    path = TOY / "build_toy_genome.py"
    spec = importlib.util.spec_from_file_location("build_toy_genome", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # dataclasses resolve types through sys.modules
    spec.loader.exec_module(module)
    return module


def _write_config(work: Path) -> Path:
    """Write a complete config for the toy run, derived from the shipped one."""
    config = yaml.safe_load((REPO / "data" / "config" / "config.yaml").read_text())
    config["queries"] = QUERIES
    config["species"] = SPECIES
    config["rpsblast"]["target_domains"] = TARGET_DOMAINS
    config["execution"]["num_cores"] = 2
    config["root"] = {
        "db_root_folder": str(work / "db"),
        "data_root_folder": str(work / "data"),
        "results_root_folder": str(work / "results"),
        "logs_root_folder": str(work / "logs"),
    }
    path = work / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return path


def _prepare_workspace(work: Path) -> Path:
    """Lay out databases, genomes and queries as the pipeline expects; return truth."""
    db_root = work / "db"
    shutil.copytree(TOY / "db", db_root)
    # The Pfam download rule declares this file as an output next to Pfam-A.hmm;
    # without it Snakemake would schedule the download.
    (db_root / "hmm_dbs" / "md5sum.txt").write_text("fixture: no download\n")
    protein = TOY / "queries" / f"{PLANTED_PROTEIN}.fa"
    truth = _load_generator().build(protein, db_root)
    query_dir = work / "data" / "input" / "fastas"
    query_dir.mkdir(parents=True)
    for accession in QUERIES:
        shutil.copy(TOY / "queries" / f"{accession}.fa", query_dir)
    return truth


def _run_stage(stage: str, config: Path) -> None:
    """Run one launcher stage; fail the test with the output tail if it fails."""
    command = [
        sys.executable,
        str(REPO / "rpsHunter"),
        f"--{stage}",
        "--skip-validation",
    ]
    environment = {**os.environ, "RPSHUNTER_CONFIG": str(config)}
    done = subprocess.run(
        command, env=environment, capture_output=True, text=True, check=False
    )
    if done.returncode != 0:
        tail = "\n".join((done.stdout + done.stderr).splitlines()[-40:])
        pytest.fail(f"stage --{stage} exited {done.returncode}:\n{tail}")


@pytest.fixture(scope="module")
def toy_run(tmp_path_factory: pytest.TempPathFactory) -> ToyRun:
    work = tmp_path_factory.mktemp("toy")
    truth_path = _prepare_workspace(work)
    config = _write_config(work)
    for stage in STAGES:
        _run_stage(stage, config)
    return ToyRun(
        results=work / "results",
        truth=pd.read_csv(truth_path, sep="\t"),
        query_domains=pd.read_csv(
            TOY / "queries" / "query_domains.tsv", sep="\t", comment="#"
        ),
    )


def _domains(run: ToyRun) -> pd.DataFrame:
    """The merged domain table, with its coordinates as numbers.

    v0.9.0 writes Start and End as text; they are converted here so the tests
    about content do not all fail for that one reason. The text columns have a
    test of their own below.
    """
    table = pd.read_parquet(run.results / "tables" / "domains" / "domains.parquet")
    for column in ("Start", "End"):
        table[column] = pd.to_numeric(table[column])
    return table


def _loci(run: ToyRun) -> pd.DataFrame:
    path = run.results / "tables" / "inventory" / "locus_domain_inventory.parquet"
    return pd.read_parquet(path)


def _overlapping(table: pd.DataFrame, plant: pd.Series, columns: tuple) -> pd.DataFrame:
    """Rows of `table` on the plant's scaffold whose interval overlaps the plant."""
    scaffold, start, end = columns
    low = table[[start, end]].min(axis=1)
    high = table[[start, end]].max(axis=1)
    on_scaffold = table[scaffold] == plant["scaffold"]
    return table[on_scaffold & (high >= plant["start"]) & (low <= plant["end"])]


def _domains_over(run: ToyRun, plant: pd.Series) -> pd.DataFrame:
    """Domain rows of the plant's species that overlap the planted copy."""
    domains = _domains(run)
    same_species = domains[domains["Species"] == plant["species"]]
    return _overlapping(same_species, plant, ("Chromosome", "Start", "End"))


def _annotated(run: ToyRun, plant: pd.Series) -> pd.DataFrame:
    """The planted protein's annotated domains, with residues planted per domain."""
    domains = run.query_domains[run.query_domains["accession"] == plant["protein"]]
    first = domains["aa_start"].clip(lower=plant["aa_start"])
    last = domains["aa_end"].clip(upper=plant["aa_end"])
    planted = (last - first + 1).clip(lower=0)
    return domains.assign(planted=planted)


def _families_mostly_planted(run: ToyRun, plant: pd.Series) -> set[str]:
    """Families with at least half of their annotated residues in the copy."""
    annotated = _annotated(run, plant)
    length = annotated["aa_end"] - annotated["aa_start"] + 1
    return set(annotated.loc[annotated["planted"] * 2 >= length, "family"])


def _families_touched(run: ToyRun, plant: pd.Series) -> set[str]:
    """Families with at least one annotated residue in the copy."""
    annotated = _annotated(run, plant)
    return set(annotated.loc[annotated["planted"] > 0, "family"])


def _found_families(domains: pd.DataFrame) -> set[str]:
    """Families for which at least one of their models appears in `domains`."""
    names = set(domains["Domain"])
    return {family for family, models in FAMILY_MODELS.items() if names & models}


def _describe(plant: pd.Series) -> str:
    return f"{plant['species']} {plant['scaffold']}:{plant['start']}-{plant['end']}"


# --- behaviour that holds today -------------------------------------------------


def test_every_query_and_species_has_a_blast_table(toy_run: ToyRun):
    for label in QUERY_LABELS:
        for species in SPECIES:
            assert (toy_run.results / "blast" / label / f"{species}.parquet").is_file()


def test_planted_copies_are_found_by_tblastn_on_the_right_strand(toy_run: ToyRun):
    columns = ("Subject ID", "S. Start", "S. End")
    for _, plant in toy_run.truth.iterrows():
        path = toy_run.results / "blast" / "Human_PRDM9" / f"{plant['species']}.parquet"
        hits = _overlapping(pd.read_parquet(path), plant, columns)

        assert not hits.empty, f"no hit over {_describe(plant)}"
        forward = hits["S. Start"] < hits["S. End"]
        assert forward.all() if plant["strand"] == "+" else (~forward).all()


def test_chance_hits_in_random_dna_do_not_pass_the_gates(toy_run: ToyRun):
    # tBLASTn at the configured E-value reports chance matches even in random
    # DNA; none may survive into the selected table or the domain table.
    for label in QUERY_LABELS:
        selected = toy_run.results / "tables" / "blast" / label / "selected"

        assert len(pd.read_parquet(selected / "Toyus_vacuus.parquet")) == 0
    assert "Toyus_vacuus" not in set(_domains(toy_run)["Species"])


def test_each_planted_copy_shows_every_family_it_contains(toy_run: ToyRun):
    for _, plant in toy_run.truth.iterrows():
        found = _found_families(_domains_over(toy_run, plant))

        assert _families_mostly_planted(toy_run, plant) <= found, _describe(plant)


def test_no_domain_is_reported_away_from_a_planted_copy(toy_run: ToyRun):
    domains = _domains(toy_run)
    near_a_plant = pd.Series(False, index=domains.index)
    for _, plant in toy_run.truth.iterrows():
        near_a_plant.loc[_domains_over(toy_run, plant).index] = True

    assert near_a_plant.all(), domains.loc[~near_a_plant].to_string()


# --- behaviour that should hold and does not yet ------------------------------
# Each test states the correct behaviour and is expected to fail on the current
# code. `strict=True` turns an unexpected pass into a failure, so the marker has
# to be removed in the change that fixes the behaviour. Reasons name the decision
# record that will settle each one.

AWAITING_SCIENCE_AUDIT = "result-affecting; to be settled by an ADR after the audit"


@pytest.mark.xfail(
    strict=True,
    reason=f"chance matches are reported as domains at the loose RPS-BLAST "
    f"E-value: {AWAITING_SCIENCE_AUDIT}",
)
def test_no_family_is_reported_that_the_copy_lacks(toy_run: ToyRun):
    for _, plant in toy_run.truth.iterrows():
        found = _found_families(_domains_over(toy_run, plant))

        assert found <= _families_touched(toy_run, plant), _describe(plant)


@pytest.mark.xfail(
    strict=True,
    reason=f"domains carry the coordinates of the whole BLAST hit, not their "
    f"own: {AWAITING_SCIENCE_AUDIT}",
)
def test_a_domain_is_reported_at_its_own_position(toy_run: ToyRun):
    plant = toy_run.truth.iloc[0]  # the whole protein on the forward strand
    annotated = toy_run.query_domains
    krab = annotated[
        (annotated["accession"] == plant["protein"]) & (annotated["family"] == "KRAB")
    ].iloc[0]
    expected_start = plant["start"] + 3 * (krab["aa_start"] - 1)
    expected_end = plant["start"] + 3 * krab["aa_end"] - 1
    slack = 90  # bp; profile matches rarely end on the annotated residue
    over = _domains_over(toy_run, plant)
    reported = over[over["Domain"].isin(FAMILY_MODELS["KRAB"])]

    assert not reported.empty
    assert (reported["Start"] >= expected_start - slack).all()
    assert (reported["End"] <= expected_end + slack).all()


@pytest.mark.xfail(
    strict=True,
    reason="the merged domain table stores numeric columns as text (backlog)",
)
def test_numeric_columns_of_the_domain_table_are_numbers(toy_run: ToyRun):
    table = pd.read_parquet(toy_run.results / "tables" / "domains" / "domains.parquet")
    for column in ("Start", "End", "From", "To", "Evalue", "Bitscore"):
        assert pd.api.types.is_numeric_dtype(table[column]), column


@pytest.mark.xfail(
    strict=True,
    reason=f"loci are built from hits that failed the gates: {AWAITING_SCIENCE_AUDIT}",
)
def test_random_dna_yields_no_locus(toy_run: ToyRun):
    assert "Toyus_vacuus" not in set(_loci(toy_run)["Species"])


@pytest.mark.xfail(
    strict=True,
    reason=f"separate fragments are chained into one locus and their domains "
    f"pooled: {AWAITING_SCIENCE_AUDIT}",
)
def test_two_partial_copies_are_not_one_complete_locus(toy_run: ToyRun):
    loci = _loci(toy_run)
    fractus = loci[loci["Species"] == "Toyus_fractus"]

    assert not fractus["Complete_CDD"].any()


@pytest.mark.xfail(
    strict=True,
    reason="no GFF3 is written for a species without domains (backlog)",
)
def test_every_species_has_a_range_file(toy_run: ToyRun):
    for species in SPECIES:
        assert (toy_run.results / "ranges" / f"{species}.gff3").is_file()
