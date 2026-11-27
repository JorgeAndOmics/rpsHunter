"""The heavy-rule guard: stop a run that would start a download or a search unasked.

A stage builds whatever it needs upstream. That is convenient until a moved
results folder or an edited setting makes `--loci` quietly download genomes or
rerun every tBLASTn search. Before a run the launcher does a dry run; this
module reads it and says which heavy rules (stages.HEAVY_RULES) would start
although the stage that owns them was not requested.

It reads two things from Snakemake's text: the first "Job stats" table and each
job's "reason:" line. If a later Snakemake changes them the guard finds no jobs
and lets the run through; the tests, built on a real dry run, would flag that.

This module is pure.
"""

from __future__ import annotations

import re
from collections.abc import Collection, Mapping

import stages

_RULE = re.compile(r"^(?:local)?rule (\w+):\s*$")
_REASON = re.compile(r"^\s+reason: (.+)$")
_ROW = re.compile(r"^(\w+)\s+(\d+)\s*$")
_STATS = "Job stats:"
_OWNER = {rule: stage.flag for stage in stages.STAGES for rule in stage.heavy}


def job_counts(dry_run: str) -> dict[str, int]:
    """Return jobs per rule from the first "Job stats" table of a dry run.

    Snakemake prints the table again at the end, so reading stops at the
    table's "total" row. Empty when nothing would run.
    """
    counts: dict[str, int] = {}
    in_table = False
    for line in dry_run.splitlines():
        if line.startswith(_STATS):
            in_table = True
        elif in_table and (row := _ROW.match(line)):
            if row[1] == "total":
                break
            counts[row[1]] = int(row[2])
    return counts


def first_reasons(dry_run: str) -> dict[str, str]:
    """Return the first reason Snakemake gives for running each rule."""
    reasons: dict[str, str] = {}
    rule = ""
    for line in dry_run.splitlines():
        if found := _RULE.match(line):
            rule = found[1]
        elif (found := _REASON.match(line)) and rule:
            reasons.setdefault(rule, found[1].strip())
    return reasons


def blocked(counts: Mapping[str, int], allowed: Collection[str]) -> list[str]:
    """Return the heavy rules that would run although no requested stage owns them."""
    return [
        rule for rule in counts if rule in stages.HEAVY_RULES and rule not in allowed
    ]


def advice(reason: str) -> str:
    """Say what to do about one reason Snakemake gave for running a heavy rule."""
    lowered = reason.lower()
    if "missing output" in lowered:
        return (
            "Its outputs are not where the config says they should be. Check the "
            "folders under `root:` in the config before running anything."
        )
    if "code has changed" in lowered:
        return (
            "The rule's own text changed (a new version of rpsHunter), not its "
            "inputs or settings. If its outputs are still right, keep them by "
            "running once with: -- --rerun-triggers mtime params input"
        )
    if "have changed" in lowered:
        return (
            "A setting changed since its outputs were made. If that is intended, "
            "request its stage."
        )
    if "updated" in lowered:
        return (
            "One of its inputs is newer than its outputs, or an earlier rule would "
            "run first. Look at that input or rule; this one only follows."
        )
    return "Read the reason above; if unsure, do not run."


def report(
    rules: list[str], counts: Mapping[str, int], reasons: Mapping[str, str]
) -> str:
    """Return the message shown when the guard stops a run."""
    lines = ["These rules would run although their stage was not requested:"]
    for rule in rules:
        jobs = counts.get(rule, 0)
        reason = reasons.get(rule, "no reason given")
        lines += [
            f"  {rule} ({jobs} job{'' if jobs == 1 else 's'}), stage {_OWNER[rule]}",
            f"    reason: {reason[:300]}",
            f"    {advice(reason)}",
        ]
    lines.append(
        "Nothing was run. To go ahead, add the stage named above, or --allow-heavy."
    )
    return "\n".join(lines)
