"""The launcher's stages come from one table; the parser is built from it."""

import pytest

import stages


def test_every_stage_has_a_distinct_flag_a_known_phase_targets_and_help():
    flags = [stage.flag for stage in stages.STAGES]

    assert len(set(flags)) == len(flags)
    for stage in stages.STAGES:
        assert stage.flag.startswith("--")
        assert stage.phase in stages.PHASES
        assert stage.targets
        assert stage.help.endswith(".")


def test_one_stage_is_requested_by_its_flag():
    request = stages.parse(["--loci"])

    assert [stage.flag for stage in request.stages] == ["--loci"]
    assert request.targets == ("loci",)


def test_several_stages_run_in_table_order_whatever_the_order_typed():
    request = stages.parse(["--loci", "--search"])

    assert [stage.flag for stage in request.stages] == ["--search", "--loci"]


def test_all_asks_for_every_analysis_stage_and_no_download():
    request = stages.parse(["--all"])

    flags = [stage.flag for stage in request.stages]
    assert flags == ["--blast-dbs", "--search", "--hmmer", "--domains", "--loci"]
    assert not any(stage.phase == "Setup" for stage in request.stages)


def test_targets_are_listed_once_in_stage_order():
    request = stages.parse(["--all"])

    assert len(set(request.targets)) == len(request.targets)
    assert request.targets[-1] == "loci"


def test_a_stage_may_run_its_own_heavy_rules_and_no_others():
    request = stages.parse(["--loci"])

    assert "blast_species" not in request.allowed_heavy
    assert "blast_species" in stages.parse(["--search"]).allowed_heavy
    assert "genome_downloader_setup" not in stages.parse(["--all"]).allowed_heavy


def test_every_heavy_rule_belongs_to_a_stage():
    owned = {rule for stage in stages.STAGES for rule in stage.heavy}

    assert owned == stages.HEAVY_RULES
    assert "genome_downloader_setup" in stages.HEAVY_RULES


def test_no_stage_is_an_error_that_says_how_to_start(capsys):
    with pytest.raises(SystemExit) as stopped:
        stages.parse([])

    assert stopped.value.code == 2
    assert "--all" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("old", "new"),
    [
        ("--blast", "--search"),
        ("--rpsblast", "--domains"),
        ("--rpsbproc-parser", "--domains"),
        ("--download-query", "--download-queries"),
        ("--setup-databases", "--download-databases"),
    ],
)
def test_a_renamed_flag_names_its_replacement(old, new, capsys):
    with pytest.raises(SystemExit) as stopped:
        stages.parse([old])

    assert stopped.value.code == 2
    message = capsys.readouterr().err
    assert old in message
    assert new in message


@pytest.mark.parametrize(
    "gone", ["--orf-analyser", "--concordance", "--skip-validation", "-skp"]
)
def test_a_removed_flag_says_why_it_is_gone(gone, capsys):
    with pytest.raises(SystemExit):
        stages.parse([gone, "--loci"])

    assert gone in capsys.readouterr().err


def test_retired_flags_are_not_advertised_in_the_help(capsys):
    with pytest.raises(SystemExit) as stopped:
        stages.parse(["--help"])

    text = capsys.readouterr().out
    assert stopped.value.code == 0
    assert "--rpsbproc-parser" not in text
    assert "--search" in text
    assert "Examples" in text


def test_help_lists_stages_under_their_phase(capsys):
    with pytest.raises(SystemExit):
        stages.parse(["-h"])

    text = capsys.readouterr().out
    for phase in stages.PHASES:
        assert phase in text
    assert text.index("Setup") < text.index("Search") < text.index("Domains")


def test_run_options_are_carried_on_the_request():
    request = stages.parse(
        ["--loci", "-n", "--cores", "4", "--config", "my.yaml", "--allow-heavy"]
    )

    assert request.dry_run is True
    assert request.cores == 4
    assert request.config == "my.yaml"
    assert request.allow_heavy is True


def test_defaults_leave_cores_and_config_to_the_config_file():
    request = stages.parse(["--loci"])

    assert (request.dry_run, request.cores, request.config, request.allow_heavy) == (
        False,
        None,
        None,
        False,
    )


def test_words_after_a_double_dash_go_to_snakemake():
    request = stages.parse(["--loci", "--", "--rerun-triggers", "mtime"])

    assert request.snakemake_args == ("--rerun-triggers", "mtime")


def test_forceall_is_refused_even_after_the_double_dash(capsys):
    with pytest.raises(SystemExit):
        stages.parse(["--loci", "--", "--forceall"])

    assert "--forceall" in capsys.readouterr().err


def test_an_unknown_flag_is_an_error():
    with pytest.raises(SystemExit) as stopped:
        stages.parse(["--loci", "--nonsense"])

    assert stopped.value.code == 2
