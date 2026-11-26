"""Every message of every job is one line of one shape."""

import io
import logging

import pytest

import log


def _capture(step="select_hits", scope="Human_PRDM9/Homo_sapiens"):
    stream = io.StringIO()
    log.job_logging(step, scope, stream=stream, clock=lambda: "14:02:40")
    return stream


def test_a_line_has_time_level_step_scope_and_message():
    stream = _capture()

    logging.getLogger().info("33,721 of 33,863 hits selected")

    assert stream.getvalue() == (
        "14:02:40 INFO select_hits Human_PRDM9/Homo_sapiens | 33,721 of 33,863 hits selected\n"
    )


def test_levels_are_written_as_fixed_words():
    stream = _capture()

    logging.getLogger().warning("look at this")
    logging.getLogger().error("cannot go on")
    log.ok("done")

    words = [line.split()[1] for line in stream.getvalue().splitlines()]
    assert words == ["WARN", "ERROR", "OK"]


def test_a_job_over_everything_has_the_scope_all():
    stream = io.StringIO()
    log.job_logging("locus_tables", stream=stream, clock=lambda: "09:00:00")

    log.ok("952 loci")

    assert stream.getvalue() == "09:00:00 OK locus_tables all | 952 loci\n"


def test_a_message_with_line_breaks_stays_one_line():
    stream = _capture()

    logging.getLogger().warning("first\nsecond")

    assert stream.getvalue().count("\n") == 1
    assert "first second" in stream.getvalue()


def test_setting_up_twice_does_not_write_every_line_twice():
    _capture()
    stream = _capture()

    log.ok("once")

    assert stream.getvalue().count("once") == 1


def test_a_line_can_be_read_back_into_its_parts():
    line = log.parse("14:02:40 WARN hmmer Cow_PRDM9/Mus_musculus | 3 profiles | odd")

    assert line == log.Line(
        "14:02:40", "WARN", "hmmer", "Cow_PRDM9/Mus_musculus", "3 profiles | odd"
    )


@pytest.mark.parametrize(
    "text",
    [
        "",
        "Building DAG of jobs...",
        "14:02:40 LOUD step all | x",
        "1 of 47 steps (2%) done",
    ],
)
def test_text_that_is_not_one_of_our_lines_is_not_parsed(text):
    assert log.parse(text) is None


def test_thousands_are_separated():
    assert log.count(33721) == "33,721"
    assert log.count(7) == "7"


def test_a_known_failure_is_one_error_line_with_its_fix_and_status_one():
    stream = _capture()

    def failing() -> int:
        raise log.PipelineError("no profile for PF00096", hint="check the Pfam release")

    status = log.run_main(failing)

    assert status == 1
    assert stream.getvalue() == (
        "14:02:40 ERROR select_hits Human_PRDM9/Homo_sapiens | "
        "no profile for PF00096. Fix: check the Pfam release\n"
    )


def test_an_unexpected_failure_is_one_error_line_and_the_traceback_follows():
    stream = _capture()

    def failing() -> int:
        raise KeyError("Tag")

    status = log.run_main(failing)

    lines = stream.getvalue().splitlines()
    assert status == 1
    assert lines[0].startswith(
        "14:02:40 ERROR select_hits Human_PRDM9/Homo_sapiens | KeyError"
    )
    assert any("Traceback" in line for line in lines[1:])


def test_a_main_that_succeeds_passes_its_status_on():
    _capture()

    assert log.run_main(lambda: 0) == 0
