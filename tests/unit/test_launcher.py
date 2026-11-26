"""The launcher makes one Snakemake call, after a dry run the guard has read."""

from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

import launcher
import stages

DRY_RUNS = Path(__file__).parents[1] / "fixtures" / "dry_runs"
WORK = (DRY_RUNS / "work_to_do.txt").read_text()
NOTHING = (DRY_RUNS / "nothing_to_do.txt").read_text()
RUN_OUTPUT = (DRY_RUNS.parent / "runs" / "clean.txt").read_text().splitlines()


class Runner:
    """Stands in for the programs the launcher starts.

    Records calls, answers dry runs with a text, and plays a real run back line
    by line.
    """

    def __init__(self, dry_run: str, dry_status: int = 0, run_status: int = 0) -> None:
        self.calls: list[SimpleNamespace] = []
        self.real: list[SimpleNamespace] = []
        self.manifests: list = []
        self._dry_run, self._dry_status, self._run_status = (
            dry_run,
            dry_status,
            run_status,
        )

    def __call__(self, command, **options):
        self.calls.append(SimpleNamespace(command=list(command), options=options))
        return SimpleNamespace(
            returncode=self._dry_status, stdout=self._dry_run, stderr=""
        )

    def stream(self, command, on_line, **options):
        self.real.append(SimpleNamespace(command=list(command), options=options))
        for line in RUN_OUTPUT:
            on_line(line)
        return self._run_status

    def tools(self, check=None):
        return launcher.Tools(
            run=self, stream=self.stream, check=check or _sound, record=self.record
        )

    def record(self, run, config_file, folder, stamp):
        self.manifests.append(run)
        return folder / f"{stamp}.yaml"


def _sound(_request, _config):
    """Checks that find nothing."""
    return launcher.preflight.Findings(problems=[], notes=[])


@pytest.fixture
def config(tmp_path):
    logs = tmp_path / "logs"
    (logs / "runs").mkdir(parents=True)  # the Snakefile makes it in a real run
    path = tmp_path / "config.yaml"
    roots = {
        "logs_root_folder": str(logs),
        "results_root_folder": str(tmp_path / "results"),
    }
    settings = {"execution": {"num_cores": 6}, "root": roots}
    path.write_text(yaml.safe_dump(settings))
    return path


def test_the_config_comes_from_the_flag_then_the_environment_then_the_default(
    config, tmp_path
):
    other = tmp_path / "other.yaml"
    other.write_text("execution: {}\n")
    environment = {"RPSHUNTER_CONFIG": str(other)}

    assert launcher.config_path(str(config), environment) == config
    assert launcher.config_path(None, environment) == other
    assert launcher.config_path(None, {}) == launcher.DEFAULT_CONFIG


def test_a_relative_config_path_is_made_absolute(config, monkeypatch):
    monkeypatch.chdir(config.parent)

    assert launcher.config_path("config.yaml", {}) == config


def test_a_missing_config_is_an_error_naming_the_file(tmp_path):
    with pytest.raises(launcher.LaunchError, match=r"nowhere\.yaml"):
        launcher.config_path(str(tmp_path / "nowhere.yaml"), {})


def test_cores_come_from_the_flag_or_else_the_config(config):
    assert launcher.cores_for(stages.parse(["--loci"]), config) == 6
    assert launcher.cores_for(stages.parse(["--loci", "--cores", "2"]), config) == 2


def test_the_command_names_every_target_once_and_the_fixed_options(config):
    request = stages.parse(["--domains", "--loci", "--", "--rerun-triggers", "mtime"])

    command = launcher.snakemake_command(request, cores=6)

    assert command[1:3] == ["-m", "snakemake"]
    assert command[command.index("--cores") + 1] == "6"
    assert "rpsbproc_parser" in command
    assert "loci" in command
    assert "--rerun-incomplete" in command
    assert command[-2:] == ["--rerun-triggers", "mtime"]


def test_a_run_is_one_dry_run_then_one_real_call_with_the_config_in_the_environment(
    config,
):
    runner = Runner(WORK)

    status = launcher.main(
        ["--search", "--loci", "--config", str(config)], tools=runner.tools()
    )

    assert status == 0
    assert len(runner.calls) == 1  # the dry run
    assert len(runner.real) == 1
    assert runner.real[0].options["env"]["RPSHUNTER_CONFIG"] == str(config)
    assert runner.real[0].options["cwd"] == launcher.WORKFLOW


def test_a_heavy_rule_nobody_asked_for_stops_the_run_before_it_starts(config, capsys):
    runner = Runner(WORK)

    status = launcher.main(["--loci", "--config", str(config)], tools=runner.tools())

    assert status == launcher.EXIT_GUARD
    assert runner.real == []
    assert "blast_species" in capsys.readouterr().err


def test_allow_heavy_lets_the_run_through(config):
    runner = Runner(WORK)

    status = launcher.main(
        ["--loci", "--allow-heavy", "--config", str(config)], tools=runner.tools()
    )

    assert status == 0
    assert len(runner.real) == 1


def test_a_dry_run_lists_the_jobs_and_runs_nothing(config, capsys):
    runner = Runner(WORK)

    status = launcher.main(
        ["--search", "--loci", "-n", "--config", str(config)], tools=runner.tools()
    )

    shown = capsys.readouterr().out
    assert status == 0
    assert runner.real == []
    assert "blast_species" in shown
    assert "77" in shown


def test_a_dry_run_says_when_the_guard_would_stop_the_real_run(config, capsys):
    status = launcher.main(
        ["--loci", "-n", "--config", str(config)], tools=Runner(WORK).tools()
    )

    assert status == launcher.EXIT_GUARD
    assert "--search" in capsys.readouterr().err


def test_nothing_to_do_is_said_and_is_a_success(config, capsys):
    runner = Runner(NOTHING)

    status = launcher.main(["--loci", "--config", str(config)], tools=runner.tools())

    assert status == 0
    assert runner.real == []
    assert "up to date" in capsys.readouterr().out


def test_a_failing_dry_run_shows_snakemake_s_output_and_its_status(config, capsys):
    runner = Runner("SyntaxError in Snakefile", dry_status=1)

    status = launcher.main(["--loci", "--config", str(config)], tools=runner.tools())

    assert status == 1
    assert "SyntaxError" in capsys.readouterr().err


def test_the_exit_status_of_snakemake_is_passed_on(config):
    runner = Runner(WORK, run_status=5)

    assert (
        launcher.main(
            ["--search", "--loci", "--config", str(config)], tools=runner.tools()
        )
        == 5
    )


def test_a_missing_config_is_reported_in_one_line_with_status_two(tmp_path, capsys):
    status = launcher.main(
        ["--loci", "--config", str(tmp_path / "nowhere.yaml")],
        tools=Runner(WORK).tools(),
    )

    assert status == 2
    assert "nowhere.yaml" in capsys.readouterr().err


def test_a_problem_found_by_the_checks_stops_everything_before_snakemake(
    config, capsys
):
    runner = Runner(WORK)

    def failing(_request, _config):
        problems = ["program not found on PATH: tblastn"]
        return launcher.preflight.Findings(problems=problems, notes=[])

    status = launcher.main(
        ["--search", "--config", str(config)], tools=runner.tools(failing)
    )

    assert status == 2
    assert runner.calls == []
    assert "tblastn" in capsys.readouterr().err


def test_a_note_is_shown_and_the_run_still_goes_on(config, capsys):
    runner = Runner(WORK)

    def noting(_request, _config):
        notes = ["only 2.0 GB free in /results"]
        return launcher.preflight.Findings(problems=[], notes=notes)

    arguments = ["--search", "--loci", "--config", str(config)]
    status = launcher.main(arguments, tools=runner.tools(noting))

    assert status == 0
    assert len(runner.real) == 1
    assert "2.0 GB" in capsys.readouterr().out


def test_a_run_is_shown_job_by_job_and_summed_up(config, capsys):
    runner = Runner(WORK)

    launcher.main(["--search", "--loci", "--config", str(config)], tools=runner.tools())

    shown = capsys.readouterr().out
    assert "rpsHunter " in shown  # the banner
    assert "call_domains" in shown
    assert "localrule" not in shown  # Snakemake's own lines stay in the run log
    assert "Done in" in shown


def test_everything_the_run_printed_is_kept_in_the_run_log(config):
    runner = Runner(WORK)

    launcher.main(["--search", "--loci", "--config", str(config)], tools=runner.tools())

    logs = list((config.parent / "logs" / "runs").glob("*.log"))
    assert len(logs) == 1
    assert logs[0].read_text().splitlines() == RUN_OUTPUT


def test_quiet_shows_no_job_line(config, capsys):
    runner = Runner(WORK)
    arguments = ["--search", "--loci", "--verbosity", "quiet", "--config", str(config)]

    launcher.main(arguments, tools=runner.tools())

    shown = capsys.readouterr().out
    assert "domain calls" not in shown  # the jobs' headlines
    assert "Done in" in shown


def test_the_version_comes_from_the_citation_file():
    assert launcher.version()[0].isdigit()


def test_a_run_that_executed_jobs_gets_a_manifest_with_its_stages_jobs_and_status(
    config, capsys
):
    runner = Runner(WORK, run_status=5)

    launcher.main(["--search", "--loci", "--config", str(config)], tools=runner.tools())

    (run,) = runner.manifests
    assert run.stages == ["--search", "--loci"]
    assert run.jobs["blast_species"] == 15
    assert run.exit_status == 5  # a failed run is recorded too
    assert "Manifest: " in capsys.readouterr().out


def test_a_dry_run_and_a_run_with_nothing_to_do_write_no_manifest(config):
    dry, idle = Runner(WORK), Runner(NOTHING)

    launcher.main(
        ["--search", "--loci", "-n", "--config", str(config)], tools=dry.tools()
    )
    launcher.main(["--loci", "--config", str(config)], tools=idle.tools())

    assert dry.manifests == []
    assert idle.manifests == []
