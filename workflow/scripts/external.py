"""Run external programs (BLAST+, HMMER, rpsbproc) with one failure rule.

Every external tool the pipeline calls goes through run_tool. A tool that
exits with a non-zero status, or that is not installed, raises ToolError; the
calling script does not catch it, so the Snakemake job fails and names the
species it failed for. Nothing may turn a failed search into an empty result:
"the tool crashed" and "the genome has no hits" are different answers.
"""

from __future__ import annotations

import subprocess
from collections.abc import Sequence

STDERR_TAIL_LINES = 20
"""Lines of the tool's error output kept in the message: the end is where
BLAST+ and HMMER say what went wrong, and a full log can run to megabytes."""


class ToolError(RuntimeError):
    """An external tool failed or could not be started.

    Attributes:
        tool: Name of the tool, for the reader of the log.
        returncode: The tool's exit status, or None if it never started.
    """

    def __init__(self, tool: str, returncode: int | None, detail: str) -> None:
        self.tool = tool
        self.returncode = returncode
        if returncode is None:
            headline = f"{tool} could not be started"
        else:
            headline = f"{tool} exited with status {returncode}"
        super().__init__(f"{headline}: {detail}" if detail else headline)


def run_tool(command: Sequence[str], *, tool: str | None = None) -> str:
    """Run an external program to completion and return its standard output.

    Args:
        command: The program and its arguments. No shell is involved.
        tool: Name to use in error messages; defaults to the program.

    Returns:
        Everything the program wrote to standard output, as text.

    Raises:
        ToolError: The program exited with a non-zero status, or was not found.
        ValueError: The command is empty.
    """
    if not command:
        raise ValueError("empty command")
    name = tool or command[0]
    try:
        done = subprocess.run(command, capture_output=True, text=True, check=False)
    except FileNotFoundError as error:
        raise ToolError(name, None, f"'{command[0]}' not found on PATH") from error
    if done.returncode != 0:
        tail = "\n".join(done.stderr.strip().splitlines()[-STDERR_TAIL_LINES:])
        raise ToolError(name, done.returncode, tail)
    return done.stdout
