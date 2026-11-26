"""End-to-end run of the pipeline on the toy genomes, checked against planted truth.

The toy genomes (tests/fixtures/toy/build_toy_genome.py) carry copies of the
human PRDM9 protein at known positions, and the miniature databases under
tests/fixtures/toy/db hold just the profiles that recognise its domains. One
module-scoped run feeds every test here, because the run takes about a minute.

Needs the external tools (blast, hmmer, rpsbproc, snakemake, R) on PATH, so the
tests carry the `integration` marker.
"""

from __future__ import annotations

import hashlib
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
# The families the toy run looks for, by the accessions of the mini database.
DOMAINS = {
    "KRAB": {"cdd": ["pfam01352"], "pfam": ["PF01352"]},
    "SSXRD": {"cdd": ["pfam09514"], "pfam": ["PF09514"]},
    "SET": {"cdd": ["cd19193", "cd08161"], "pfam": ["PF00856"]},
    "zinc_fingers": {
        "cdd": ["pfam00096", "pfam13465", "COG5048"],
        "pfam": ["PF00096"],
        "repeat": True,
    },
}

# Launcher flags in pipeline order (the v0.9.0 launcher runs one stage per call).
STAGES = (
    "blast-dbs",
    "hmmer",
    "blast",
    "rpsblast",
    "rpsbproc",
    "rpsbproc-parser",
    "loci",
)


@dataclass(frozen=True)
class ToyRun:
    """Where a finished toy run left its inputs and outputs."""

    results: Path
    config: Path
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


def _species_block(db_root: Path) -> dict:
    """The toy species for the config, each pinned by the checksum of its FASTA."""
    return {
        key: {
            "name": name,
            "sha256": hashlib.sha256((db_root / f"{key}.fa").read_bytes()).hexdigest(),
        }
        for key, name in SPECIES.items()
    }


def _write_config(work: Path) -> Path:
    """Write a complete config for the toy run, derived from the shipped one."""
    config = yaml.safe_load((REPO / "data" / "config" / "config.yaml").read_text())
    config["queries"] = QUERIES
    config["species"] = _species_block(work / "db")
    config["domains"] = DOMAINS
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
        config=config,
        truth=pd.read_csv(truth_path, sep="\t"),
        query_domains=pd.read_csv(
            TOY / "queries" / "query_domains.tsv", sep="\t", comment="#"
        ),
    )


def _domains(run: ToyRun) -> pd.DataFrame:
    """The table of domain instances of the whole run."""
    return pd.read_parquet(run.results / "tables" / "domains" / "domains.parquet")


def _loci(run: ToyRun) -> pd.DataFrame:
    return pd.read_parquet(run.results / "tables" / "loci" / "loci.parquet")


def _locus_families(run: ToyRun, species: str, family: str) -> pd.DataFrame:
    """Rows of the per-locus family evidence for one species and family."""
    table = pd.read_parquet(run.results / "tables" / "loci" / "locus_families.parquet")
    return table[(table["Species"] == species) & (table["Family"] == family)]


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
    return _overlapping(same_species, plant, ("Sequence", "Start", "End"))


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
    """Families reported among the given domain rows."""
    return set(domains["Family"])


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


def test_planted_copies_have_no_internal_stop_codon(toy_run: ToyRun):
    # The planted copies are exact back-translations, so their hits are clean;
    # the flag exists for pseudogenes and frameshifts.
    columns = ("Subject ID", "S. Start", "S. End")
    plant = toy_run.truth.iloc[0]
    path = toy_run.results / "blast" / "Human_PRDM9" / f"{plant['species']}.parquet"
    hits = _overlapping(pd.read_parquet(path), plant, columns)

    assert hits["Internal_Stop"].dtype == bool
    assert not hits.loc[hits["Pct Identity"] == 100, "Internal_Stop"].any()


def test_chance_hits_in_random_dna_yield_no_domain(toy_run: ToyRun):
    # tBLASTn at the configured E-value reports chance matches even in random
    # DNA, and some pass the quality gate. They are searched like any other hit
    # (ADR-010: nothing but the quality gate selects) and must yield no domain.
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


def _jobs_planned_after_changing(run: ToyRun, section: str, key: str, value) -> str:
    """Dry-run the search-and-gate stage with one config value changed.

    Returns Snakemake's plan. The changed config is written beside the original
    and nothing is executed, so the finished run is left as it was.
    """
    config = yaml.safe_load(run.config.read_text())
    config[section][key] = value
    changed = run.config.with_name(f"config_{section}_{key}.yaml")
    changed.write_text(yaml.safe_dump(config, sort_keys=False))
    command = [sys.executable, str(REPO / "rpsHunter"), "--blast", "--skip-validation"]
    done = subprocess.run(
        [*command, "-n"],
        env={**os.environ, "RPSHUNTER_CONFIG": str(changed)},
        capture_output=True,
        text=True,
        check=True,
    )
    return done.stdout + done.stderr


def _scheduled(plan: str, rule: str) -> bool:
    """Whether the dry-run plan's job table lists the rule."""
    return any(line.split()[:1] == [rule] for line in plan.splitlines())


def test_every_job_leaves_a_log_and_a_benchmark(toy_run: ToyRun):
    logs = toy_run.results.parent / "logs"
    for label in QUERY_LABELS:
        for species in SPECIES:
            name = Path("blast_species") / label / species
            assert (logs / "jobs" / name.with_suffix(".log")).is_file()
            assert (logs / "benchmarks" / name.with_suffix(".tsv")).is_file()


def test_changing_a_gate_threshold_reruns_the_gate_and_not_the_search(toy_run: ToyRun):
    plan = _jobs_planned_after_changing(toy_run, "blast", "bitscore", 55)

    assert _scheduled(plan, "select_hits_species")
    assert not _scheduled(plan, "blast_species")


def test_changing_the_search_e_value_reruns_the_search(toy_run: ToyRun):
    plan = _jobs_planned_after_changing(toy_run, "blast", "e_value", 5)

    assert _scheduled(plan, "blast_species")


def test_an_unchanged_config_schedules_no_search_or_gate_job(toy_run: ToyRun):
    unchanged = toy_run.config.read_text()
    current = yaml.safe_load(unchanged)["blast"]["bitscore"]
    plan = _jobs_planned_after_changing(toy_run, "blast", "bitscore", current)

    assert not _scheduled(plan, "blast_species")
    assert not _scheduled(plan, "select_hits_species")


def test_every_call_carries_the_family_of_its_model(toy_run: ToyRun):
    path = toy_run.results / "domains" / "Human_PRDM9" / "Toyus_plenus.parquet"
    calls = pd.read_parquet(path)

    assert not calls.empty
    assert set(calls["Family"]) <= set(DOMAINS)
    by_model = calls.drop_duplicates("Model").set_index("Model")["Family"]
    assert by_model["pfam09514"] == "SSXRD"
    assert by_model["cd19193"] == "SET"


# Fixed by ADR-006 (E-values on NCBI's scale, cut-off 1): chance matches used to be
# reported as domains.
def test_no_family_is_reported_that_the_copy_lacks(toy_run: ToyRun):
    for _, plant in toy_run.truth.iterrows():
        found = _found_families(_domains_over(toy_run, plant))

        assert found <= _families_touched(toy_run, plant), _describe(plant)


# --- behaviour that should hold and does not yet ------------------------------
# Each test states the correct behaviour and is expected to fail on the current
# code. `strict=True` turns an unexpected pass into a failure, so the marker has
# to be removed in the change that fixes the behaviour. Reasons name the decision
# record that will settle each one.

AWAITING_SCIENCE_AUDIT = "result-affecting; to be settled by an ADR after the audit"


# Fixed by ADR-007: domains used to carry the coordinates of the whole BLAST hit.
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
    reported = over[over["Family"] == "KRAB"]

    assert not reported.empty
    assert (reported["Start"] >= expected_start - slack).all()
    assert (reported["End"] <= expected_end + slack).all()


# Fixed with ADR-007: the domain table used to store its numbers as text.
def test_numeric_columns_of_the_domain_table_are_numbers(toy_run: ToyRun):
    table = pd.read_parquet(toy_run.results / "tables" / "domains" / "domains.parquet")
    for column in ("Start", "End", "Evalue", "Bitscore", "N_Calls"):
        assert pd.api.types.is_numeric_dtype(table[column]), column


# ADR-008 and ADR-010: a locus is a query-like region; one made of chance hits
# is listed, and is told apart by carrying no family.
def test_random_dna_yields_no_locus_with_a_family(toy_run: ToyRun):
    loci = _loci(toy_run)
    vacuus = loci[loci["Species"] == "Toyus_vacuus"]

    assert (vacuus["N_Families"] == 0).all()
    assert set(vacuus["Families"]) <= {""}


def test_each_whole_copy_is_one_locus_carrying_every_family(toy_run: ToyRun):
    loci = _loci(toy_run)
    plenus = loci[(loci["Species"] == "Toyus_plenus") & (loci["N_Families"] > 0)]

    assert len(plenus) == 2
    assert set(plenus["Families"]) == {"KRAB,SSXRD,SET,zinc_fingers"}


def test_a_whole_copy_has_its_set_domain_fully_covered(toy_run: ToyRun):
    coverage = _locus_families(toy_run, "Toyus_plenus", "SET")["Coverage"]

    assert len(coverage) == 2
    assert (coverage >= 0.8).all()


# Fixed by ADR-008 and ADR-009: two fragments used to be pooled into a locus that
# was called complete. They may still share a locus (they are 16 kb apart), but
# the locus now says how much of each model it holds instead of giving a verdict.
def test_fragments_of_a_copy_do_not_add_up_to_a_whole_set_domain(toy_run: ToyRun):
    coverage = _locus_families(toy_run, "Toyus_fractus", "SET")["Coverage"]

    assert not coverage.empty
    assert (coverage < 0.5).all()


# ADR-010: HMMER is independent evidence. It reports its own domains and filters nothing.
def _hmmer_domains(run: ToyRun) -> pd.DataFrame:
    return pd.read_parquet(run.results / "tables" / "domains" / "hmmer_domains.parquet")


def test_hmmer_reports_the_planted_families_at_their_own_position(toy_run: ToyRun):
    found = _hmmer_domains(toy_run)
    plenus = found[found["Species"] == "Toyus_plenus"]

    assert {"KRAB", "SET", "zinc_fingers"} <= set(plenus["Family"])
    assert (plenus["Start"] >= plenus["Hit_Start"]).all()
    assert (plenus["End"] <= plenus["Hit_End"]).all()
    assert plenus["Above_GA"].any()


def test_the_two_methods_agree_on_the_planted_set_domains(toy_run: ToyRun):
    calls = pd.read_parquet(toy_run.results / "tables" / "domains" / "calls.parquet")
    planted = calls[(calls["Species"] == "Toyus_plenus") & (calls["Family"] == "SET")]

    assert not planted.empty
    assert "both" in set(planted["Method_Agreement"])
    assert set(calls["Method_Agreement"]) <= {"both", "cdd_only", "not_compared"}


def test_every_selected_hit_is_searched_whether_or_not_hmmer_matched_it(
    toy_run: ToyRun,
):
    selected = toy_run.results / "tables" / "blast" / "Human_PRDM9" / "selected"
    hits = pd.read_parquet(selected / "Toyus_plenus.parquet")
    fasta = toy_run.results / "fastas" / "Human_PRDM9" / "Toyus_plenus.fa"

    assert "HMM_Filtered" not in hits.columns
    assert fasta.read_text().count(">") == len(hits)


# Fixed with ADR-007: a species without domains used to get no GFF3 file.
def test_every_species_has_a_range_file(toy_run: ToyRun):
    for species in SPECIES:
        assert (toy_run.results / "ranges" / f"{species}.gff3").is_file()
