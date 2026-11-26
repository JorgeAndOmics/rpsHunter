"""External tools either succeed or stop the job: a failure is never an empty result."""

import sys

import pytest

import external

PYTHON = sys.executable


def test_returns_what_the_tool_printed():
    output = external.run_tool([PYTHON, "-c", "print('twelve columns')"])

    assert output == "twelve columns\n"


def test_a_failing_tool_raises_with_its_name_status_and_last_words():
    script = (
        "import sys; sys.stderr.write('BLAST Database error: no alias\\n'); sys.exit(3)"
    )

    with pytest.raises(external.ToolError) as raised:
        external.run_tool([PYTHON, "-c", script], tool="tblastn")

    assert raised.value.tool == "tblastn"
    assert raised.value.returncode == 3
    assert "tblastn exited with status 3" in str(raised.value)
    assert "BLAST Database error: no alias" in str(raised.value)


def test_only_the_end_of_a_long_error_is_kept():
    script = (
        "import sys\nfor i in range(500): sys.stderr.write(f'line {i}\\n')\nsys.exit(1)"
    )

    with pytest.raises(external.ToolError) as raised:
        external.run_tool([PYTHON, "-c", script])

    message = str(raised.value)
    assert "line 499" in message
    assert "line 0\n" not in message


def test_a_tool_that_is_not_installed_raises_the_same_error():
    with pytest.raises(external.ToolError) as raised:
        external.run_tool(["no-such-tool-rpshunter", "-version"])

    assert raised.value.tool == "no-such-tool-rpshunter"
    assert raised.value.returncode is None
    assert "not found" in str(raised.value)


def test_the_tool_name_defaults_to_the_program():
    with pytest.raises(external.ToolError) as raised:
        external.run_tool([PYTHON, "-c", "raise SystemExit(2)"])

    assert raised.value.tool == PYTHON


def test_an_empty_command_is_a_programming_error():
    with pytest.raises(ValueError, match="empty command"):
        external.run_tool([])
