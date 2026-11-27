"""Checks before a run: quick, offline, and every problem said in one go."""

from pathlib import Path

import yaml

import preflight
import stages

REPO = Path(__file__).parents[2]
SHIPPED = REPO / "data" / "config" / "config.yaml"
SCHEMA = REPO / "data" / "config" / "schema.yaml"


def _config(**changes):
    config = yaml.safe_load(SHIPPED.read_text())
    config.update(changes)
    return config


def _written(tmp_path, config):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return path


# --- keys of earlier versions ---------------------------------------------------


def test_the_shipped_config_has_no_retired_key():
    assert preflight.retired_key_problems(_config()) == []


def test_a_retired_section_is_named_with_what_replaces_it():
    problems = preflight.retired_key_problems(_config(orf={"enabled": True}))

    assert len(problems) == 1
    assert "orf" in problems[0]
    assert "ADR-011" in problems[0]


def test_keys_that_lived_under_execution_are_found_there():
    config = _config()
    config["execution"]["random_id_length"] = 6
    config["execution"]["max_retrieval_attempts"] = 3

    problems = preflight.retired_key_problems(config)

    assert any("execution.random_id_length" in p and "ADR-004" in p for p in problems)
    assert any("execution.max_retrieval_attempts" in p for p in problems)


def test_a_retired_key_inside_a_section_is_found_by_its_path():
    config = _config()
    config["rpsblast"]["target_domains"] = ["KRAB"]
    config["hit_domain_inventory"] = {"locus_gap": 50000}

    problems = preflight.retired_key_problems(config)

    assert any("rpsblast.target_domains" in p and "domains" in p for p in problems)
    assert any("hit_domain_inventory" in p and "loci.gap" in p for p in problems)


# --- schema and content ---------------------------------------------------------


def test_the_shipped_config_passes_the_schema(tmp_path):
    assert preflight.schema_problems(SHIPPED, SCHEMA) == []


def test_a_wrong_type_is_reported_with_its_key(tmp_path):
    config = _config()
    config["blast"]["perc_identity"] = "high"

    problems = preflight.schema_problems(_written(tmp_path, config), SCHEMA)

    assert any("perc_identity" in problem for problem in problems)


def test_a_file_that_is_not_yaml_is_one_problem_not_a_traceback(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text("blast: [unclosed\n")

    problems = preflight.schema_problems(path, SCHEMA)

    assert len(problems) == 1
    assert "config.yaml" in problems[0]


def test_family_and_species_rules_beyond_the_schema_are_checked():
    config = _config()
    config["domains"]["second"] = dict(
        config["domains"]["SET"]
    )  # same models, two families
    first = next(iter(config["species"]))
    config["species"][first] = {"name": "No pin"}

    problems = preflight.content_problems(config)

    assert any("more than one family" in problem for problem in problems)
    assert any(first in problem for problem in problems)


def test_the_shipped_config_has_no_content_problem():
    assert preflight.content_problems(_config()) == []


# --- tools ----------------------------------------------------------------------


def test_missing_tools_are_listed_and_present_ones_are_not():
    installed = {"tblastn": "/bin/tblastn"}

    missing = preflight.missing_tools(("tblastn", "rpsbproc"), {}, installed.get)

    assert missing == ["rpsbproc"]


def test_a_program_renamed_in_the_config_is_looked_up_by_its_configured_name():
    installed = {"tblastn-2.17": "/opt/tblastn-2.17"}

    missing = preflight.missing_tools(
        ("tblastn",), {"tblastn": "tblastn-2.17"}, installed.get
    )

    assert missing == []


# --- input files ----------------------------------------------------------------


def _layout(tmp_path, present=()):
    files = {
        "genomes": {"Homo_sapiens": tmp_path / "Homo_sapiens.fa"},
        "queries": {"NP_1.1": tmp_path / "NP_1.1.fa"},
        "cdd_models": {"pfam00096": tmp_path / "pfam00096.smp"},
        "rpsbproc_data": tmp_path / "cddid.tbl",
        "pfam_library": tmp_path / "Pfam-A.hmm",
    }
    for name in present:
        (tmp_path / name).write_text("x")
    return preflight.Layout(hmmer_active=True, **files)


def _problems(tmp_path, argv, present=()):
    return preflight.input_problems(stages.parse(argv), _layout(tmp_path, present))


def test_each_stage_needs_the_inputs_it_reads(tmp_path):
    assert "Homo_sapiens.fa" in _problems(tmp_path, ["--blast-dbs"])[0]
    assert "NP_1.1.fa" in _problems(tmp_path, ["--search"])[0]
    assert "Pfam-A.hmm" in _problems(tmp_path, ["--hmmer"])[0]
    # a CDD model, rpsbproc's data, and Pfam: --domains also writes HMMER's table
    assert len(_problems(tmp_path, ["--domains"])) == 3
    assert _problems(tmp_path, ["--loci"]) == []


def test_inputs_that_are_there_raise_nothing(tmp_path):
    present = (
        "Homo_sapiens.fa",
        "NP_1.1.fa",
        "pfam00096.smp",
        "cddid.tbl",
        "Pfam-A.hmm",
    )

    assert _problems(tmp_path, ["--all"], present) == []


def test_a_missing_input_names_the_stage_that_fetches_it(tmp_path):
    assert "--download-genomes" in _problems(tmp_path, ["--blast-dbs"])[0]
    assert "--download-queries" in _problems(tmp_path, ["--search"])[0]
    assert "--download-databases" in _problems(tmp_path, ["--domains"])[0]


def test_an_input_the_same_run_downloads_is_not_missing(tmp_path):
    assert _problems(tmp_path, ["--download-genomes", "--blast-dbs"]) == []
    assert _problems(tmp_path, ["--download-databases", "--hmmer", "--domains"]) == []


def test_a_configured_species_tree_must_exist_when_figures_are_drawn(tmp_path):
    layout = _layout(tmp_path)._replace(species_tree=tmp_path / "species.nwk")

    assert (
        "species.nwk"
        in preflight.input_problems(stages.parse(["--figures"]), layout)[0]
    )
    assert preflight.input_problems(stages.parse(["--loci"]), layout) == []
    (tmp_path / "species.nwk").write_text("(A,B);")
    assert preflight.input_problems(stages.parse(["--figures"]), layout) == []


def test_pfam_is_not_needed_when_hmmer_is_off(tmp_path):
    layout = _layout(tmp_path)._replace(hmmer_active=False)

    assert preflight.input_problems(stages.parse(["--hmmer"]), layout) == []


# --- notes ----------------------------------------------------------------------

LOCK = """\
# This file may be used to create an environment using:
@EXPLICIT
https://conda.anaconda.org/bioconda/linux-64/blast-2.17.0-h66d330f_0.conda#abc
https://conda.anaconda.org/bioconda/linux-64/hmmer-3.4-hb6cb901_4.tar.bz2#def
"""


def test_an_environment_matching_the_lock_gives_no_note():
    installed = {"blast-2.17.0-h66d330f_0", "hmmer-3.4-hb6cb901_4"}

    assert preflight.environment_note(LOCK, installed) is None


def test_a_package_that_differs_from_the_lock_is_named():
    installed = {"blast-2.16.0-h66d330f_0", "hmmer-3.4-hb6cb901_4", "extra-1.0-0"}

    note = preflight.environment_note(LOCK, installed)

    assert note is not None
    assert "blast-2.17.0" in note
    assert "make env" in note


def test_little_free_space_is_a_note_with_the_figure():
    assert preflight.disk_note(2 * 1024**3, Path("/results")) is not None
    assert "2.0 GB" in str(preflight.disk_note(2 * 1024**3, Path("/results")))
    assert preflight.disk_note(200 * 1024**3, Path("/results")) is None


# --- the whole check --------------------------------------------------------------


def test_a_broken_config_stops_at_the_config_and_lists_every_problem(tmp_path):
    config = _config(orf={"enabled": False})
    config["blast"]["perc_identity"] = "high"

    findings = preflight.check(
        stages.parse(["--loci"]), _written(tmp_path, config), which=lambda _: "x"
    )

    assert len(findings.problems) >= 2
    assert any("orf" in problem for problem in findings.problems)
    assert any("perc_identity" in problem for problem in findings.problems)


def test_a_sound_config_with_its_inputs_missing_reports_the_inputs(tmp_path):
    config = _config()
    config["root"] = {key: str(tmp_path / key) for key in config["root"]}

    findings = preflight.check(
        stages.parse(["--blast-dbs"]), _written(tmp_path, config), which=lambda _: "x"
    )

    assert findings.problems
    assert all("--download-genomes" in problem for problem in findings.problems)
    assert all(
        str(tmp_path) in problem for problem in findings.problems
    )  # this config's paths


def test_a_missing_tool_is_a_problem(tmp_path):
    config = _config()
    config["root"] = {key: str(tmp_path / key) for key in config["root"]}

    findings = preflight.check(
        stages.parse(["--loci", "--hmmer"]),
        _written(tmp_path, config),
        which=lambda _: None,
    )

    assert any("hmmsearch" in problem for problem in findings.problems)
