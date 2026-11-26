"""A run's manifest says what ran, with what, well enough to repeat it."""

import hashlib
from pathlib import Path
from types import SimpleNamespace

import yaml

import manifest

REPO = Path(__file__).parents[2]


def test_a_file_is_identified_by_its_sha256(tmp_path):
    path = tmp_path / "x.fa"
    path.write_bytes(b">q\nMKV\n")

    assert manifest.file_sha256(path) == hashlib.sha256(b">q\nMKV\n").hexdigest()
    assert manifest.file_sha256(tmp_path / "absent") is None


def test_the_code_is_identified_by_commit_and_whether_it_was_modified():
    def git(command, **_options):
        answer = "abc123\n" if "rev-parse" in command else " M workflow/Snakefile\n"
        return SimpleNamespace(returncode=0, stdout=answer)

    assert manifest.code_state(REPO, git) == {
        "commit": "abc123",
        "uncommitted_changes": True,
    }


def test_a_copy_that_is_not_a_git_checkout_says_so():
    def git(command, **_options):
        return SimpleNamespace(returncode=128, stdout="")

    assert manifest.code_state(REPO, git) == {
        "commit": None,
        "uncommitted_changes": None,
    }


def test_tool_versions_come_from_the_installed_packages():
    installed = {
        "blast-2.17.0-h66d330f_0",
        "hmmer-3.4-hb6cb901_4",
        "python-dateutil-2.9-x",
    }

    versions = manifest.package_versions(
        installed, ("blast", "hmmer", "rpsbproc", "python")
    )

    assert versions == {"blast": "2.17.0-h66d330f_0", "hmmer": "3.4-hb6cb901_4"}


def test_profile_accessions_are_read_with_their_versions():
    text = "NAME  KRAB\nACC   PF01352.33\n//\nNAME  SET\nACC   PF00856.34\n//\n"

    assert manifest.profile_accessions(text) == ["PF01352.33", "PF00856.34"]


def test_the_manifest_is_written_as_yaml_named_by_the_run(tmp_path):
    record = {"run": {"exit_status": 0}, "genomes": {"X": {"assembly": "GCA_1.1"}}}

    path = manifest.write(record, tmp_path, "2026-11-27_120000")

    assert path == tmp_path / "2026-11-27_120000.yaml"
    assert yaml.safe_load(path.read_text()) == record


def test_a_whole_manifest_describes_the_run_the_config_and_the_inputs(tmp_path):
    config = yaml.safe_load((REPO / "data" / "config" / "config.yaml").read_text())
    config["root"] = {key: str(tmp_path / key) for key in config["root"]}
    config_file = tmp_path / "config.yaml"
    config_file.write_text(yaml.safe_dump(config, sort_keys=False))
    run = manifest.Run(
        command=["./rpsHunter", "--all"],
        stages=["--search", "--loci"],
        started="2026-11-27T12:00:00",
        seconds=75.4,
        exit_status=0,
        jobs={"blast_species": 15},
    )

    record = manifest.collect(run, config_file)

    assert record["run"]["stages"] == ["--search", "--loci"]
    assert record["run"]["jobs"] == {"blast_species": 15}
    assert record["config"]["sha256"] == manifest.file_sha256(config_file)
    assert record["config"]["settings"]["hmmer"]["seed"] == config["hmmer"]["seed"]
    assert record["databases"]["cdd"]["version"] == "3.22"
    assert record["databases"]["pfam"]["version"] == "38.2"
    assert set(record["genomes"]) == set(config["species"])
    first = next(iter(config["species"]))
    assert record["genomes"][first]["assembly"] == config["species"][first]["assembly"]
    assert record["genomes"][first]["fasta_bytes"] is None  # not on this disk
    assert set(record["queries"]) == set(config["queries"])
    assert "commit" in record["rpshunter"]
    assert yaml.safe_load(yaml.safe_dump(record)) == record  # plain data only
