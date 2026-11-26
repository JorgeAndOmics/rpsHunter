"""The launcher's screen: what it shows of a run, and what it keeps."""

from pathlib import Path

import pytest

import console

RUNS = Path(__file__).parents[1] / "fixtures" / "runs"
CLEAN = (RUNS / "clean.txt").read_text().splitlines()  # real Snakemake 9.27 output
FAILED = (RUNS / "failed.txt").read_text().splitlines()
PLAIN = console.Style(colour=False, symbols=False, live=False)


def _screen(verbosity="normal", style=PLAIN):
    shown: list[str] = []
    return console.Console(
        verbosity, style, lambda text: shown.append(text.rstrip("\n"))
    ), shown


def _fed(lines, verbosity="normal", style=PLAIN):
    screen, shown = _screen(verbosity, style)
    for line in lines:
        screen.feed(line)
    return screen, shown


def test_a_normal_run_shows_one_line_per_finished_job_and_nothing_of_snakemake():
    _, shown = _fed(CLEAN)

    assert len(shown) == 8  # six call tables, the domain tables, the locus tables
    assert all(line.startswith("OK ") for line in shown)
    assert not any("jobid" in line or "localrule" in line for line in shown)


def test_a_job_line_reads_step_scope_and_message():
    _, shown = _fed(["00:44:33 OK call_domains query_1/Toyus_plenus | 14 domain calls"])

    assert shown == ["OK    call_domains  query_1/Toyus_plenus  14 domain calls"]


def test_quiet_shows_no_job_line_and_verbose_shows_everything():
    _, quiet = _fed(CLEAN, "quiet")
    _, verbose = _fed(CLEAN, "verbose")

    assert quiet == []
    assert len(verbose) == len(CLEAN)


def test_a_warning_is_shown_at_every_verbosity_and_kept():
    line = "09:00:00 WARN select_hits query_1/Homo_sapiens | 3 hits without sequence"
    for verbosity in ("quiet", "normal", "verbose"):
        screen, shown = _fed([line], verbosity)

        assert len(shown) == 1, verbosity
        assert screen.warnings == [
            "select_hits query_1/Homo_sapiens: 3 hits without sequence"
        ]


def test_progress_is_followed_from_snakemake_s_step_counter():
    screen, _ = _fed(CLEAN)

    assert (screen.done, screen.total) == (10, 10)


def test_a_failed_job_is_shown_once_with_its_rule_and_log():
    screen, shown = _fed(FAILED, "quiet")

    errors = [line for line in shown if line.startswith("ERROR")]
    assert len(errors) == 1  # Snakemake repeats the job's log; the line is shown once
    assert "rpsbproc exited with status 1" in errors[0]
    assert screen.failures == [
        console.Failure(
            "rpsbproc_species",
            "/work/logs/jobs/rpsbproc_species/Human_PRDM9/Toyus_vacuus.log",
        )
    ]


def test_the_summary_of_a_clean_run_counts_jobs_and_time():
    screen, shown = _fed(CLEAN)
    del shown[:]

    screen.summary(status=0, seconds=75.4)

    assert shown == ["Done in 1 min 15 s: 10 of 10 jobs, no warnings."]


def test_the_summary_lists_warnings_once_each():
    line = "09:00:00 WARN select_hits query_1/Homo_sapiens | 3 hits without sequence"
    screen, shown = _fed([line, line])
    del shown[:]

    screen.summary(status=0, seconds=2.0)

    assert shown[0].startswith("Done in 2 s")
    assert "1 warning" in shown[0]
    assert shown[1:] == ["  select_hits query_1/Homo_sapiens: 3 hits without sequence"]


def test_the_summary_of_a_failed_run_names_each_failed_rule_and_its_log():
    screen, shown = _fed(FAILED)
    del shown[:]

    screen.summary(status=1, seconds=3.0)

    assert shown[0].startswith("Failed after 3 s: 14 of 39 jobs")
    assert shown[1] == (
        "  rpsbproc_species: /work/logs/jobs/rpsbproc_species/Human_PRDM9/Toyus_vacuus.log"
    )


def test_a_failure_without_a_recognised_job_still_says_it_failed():
    screen, shown = _screen()

    screen.summary(status=1, seconds=1.0)

    assert shown[0].startswith("Failed after 1 s")
    assert "above" in shown[1]


def test_symbols_and_colour_are_added_only_when_the_style_asks():
    fancy = console.Style(colour=True, symbols=True, live=False)

    _, shown = _fed(
        ["00:44:33 OK call_domains query_1/Toyus_plenus | 14 domain calls"], style=fancy
    )

    assert shown[0].startswith("\x1b[32m✓\x1b[0m ")
    assert "\x1b[" not in shown[0].split(" ", 1)[1]  # only the symbol is coloured


def test_a_live_screen_keeps_a_progress_line_below_the_messages():
    live = console.Style(colour=False, symbols=False, live=True)
    screen, shown = _screen(style=live)

    screen.feed("3 of 10 steps (30%) done")
    screen.feed("00:44:33 OK call_domains query_1/Toyus_plenus | 14 domain calls")

    assert shown[0] == "\r[3/10] 30%"
    assert shown[1].startswith("\r\x1b[K")  # the progress line is cleared first
    assert shown[-1] == "\r[3/10] 30%"  # and drawn again below the message


@pytest.mark.parametrize(
    ("environment", "terminal", "colour"),
    [
        ({}, True, True),
        ({}, False, False),
        ({"NO_COLOR": "1"}, True, False),
        ({"FORCE_COLOR": "1"}, False, True),
    ],
)
def test_colour_follows_the_terminal_and_the_two_standard_switches(
    environment, terminal, colour
):
    style = console.style_for(environment, is_terminal=terminal, encoding="utf-8")

    assert style.colour is colour
    assert style.live is terminal


def test_symbols_need_an_encoding_that_has_them():
    assert console.style_for({}, True, "utf-8").symbols is True
    assert console.style_for({}, True, "ascii").symbols is False


def test_the_banner_says_what_runs_with_what():
    screen, shown = _screen()

    screen.banner(
        version="0.9.0",
        config="/work/config.yaml",
        cores=16,
        flags=["--search", "--loci"],
    )

    assert shown == [
        "rpsHunter 0.9.0 | config /work/config.yaml | 16 cores | --search --loci"
    ]


def test_durations_read_naturally():
    assert console.duration(0.4) == "0 s"
    assert console.duration(59.6) == "1 min 0 s"
    assert console.duration(3725) == "1 h 2 min"
