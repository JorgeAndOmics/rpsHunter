"""The guard reads a Snakemake dry run and stops heavy rules nobody asked for."""

from pathlib import Path

import guard

DRY_RUNS = Path(__file__).parents[1] / "fixtures" / "dry_runs"
WORK = (DRY_RUNS / "work_to_do.txt").read_text()  # a real Snakemake 9.27 dry run
NOTHING = (DRY_RUNS / "nothing_to_do.txt").read_text()


def test_jobs_are_counted_per_rule_from_a_real_dry_run():
    counts = guard.job_counts(WORK)

    assert counts["blast_species"] == 15
    assert counts["locus_tables"] == 1
    assert "total" not in counts
    assert sum(counts.values()) == 77  # the table is printed twice; it is read once


def test_a_dry_run_with_nothing_to_do_has_no_jobs():
    assert guard.job_counts(NOTHING) == {}


def test_text_that_is_not_a_dry_run_has_no_jobs():
    assert guard.job_counts("") == {}
    assert (
        guard.job_counts("Rules with missing metadata: genome_downloader_setup\n") == {}
    )


def test_the_first_reason_of_each_rule_is_kept():
    reasons = guard.first_reasons(WORK)

    assert reasons["blast_species"].startswith("Missing output files")


def test_a_heavy_rule_outside_the_requested_stages_is_blocked():
    counts = guard.job_counts(WORK)

    assert guard.blocked(counts, allowed=frozenset()) == ["blast_species"]
    assert guard.blocked(counts, allowed=frozenset({"blast_species"})) == []


def test_ordinary_rules_are_never_blocked():
    assert guard.blocked({"locus_tables": 1, "rpsblast_species": 15}, frozenset()) == []


def test_the_report_names_the_rule_its_jobs_its_reason_and_the_stage_that_owns_it():
    counts = {"blast_species": 15, "genome_downloader_setup": 1}
    reasons = {"blast_species": "Missing output files: /work/x.parquet"}

    text = guard.report(["blast_species", "genome_downloader_setup"], counts, reasons)

    assert "blast_species (15 jobs)" in text
    assert "Missing output files" in text
    assert "--search" in text
    assert "genome_downloader_setup (1 job)" in text
    assert "--download-genomes" in text
    assert "Nothing was run" in text


def test_a_changed_setting_and_a_newer_input_get_different_advice():
    changed = guard.advice("Params have changed since last execution")
    newer = guard.advice("Updated input files: /work/Pfam-A.hmm")
    missing = guard.advice("Missing output files: /work/x")

    assert len({changed, newer, missing}) == 3
    assert "config" in missing
