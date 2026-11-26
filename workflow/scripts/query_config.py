"""The `queries:` block of the configuration: which proteins to search with.

    queries:
      NP_000000.1: human kinase      # protein accession: label of your choice

The label names the query's folders and appears in every job's log lines, so
it is made safe for a path (spaces and slashes become underscores) and two
queries may not end up with the same one.

This module is pure.
"""

from __future__ import annotations

from collections.abc import Mapping


class QueryConfigError(ValueError):
    """The `queries:` block is not usable; the message says why."""


def safe_label(label: str) -> str:
    """Return a label usable as a folder name."""
    return label.strip().replace(" ", "_").replace("/", "_").replace("\\", "_")


def parse_queries(raw: object) -> dict[str, str]:
    """Validate the `queries:` block; return accession to safe label, in config order.

    Raises:
        QueryConfigError: The block is missing or empty, a label is empty, or
            two queries have the same label once made safe.
    """
    if not isinstance(raw, Mapping) or not raw:
        raise QueryConfigError("queries: at least one query must be listed")
    labels = {
        str(accession): safe_label(str(label)) for accession, label in raw.items()
    }
    empty = [accession for accession, label in labels.items() if not label]
    if empty:
        raise QueryConfigError(f"queries: no label for {', '.join(empty)}")
    owners: dict[str, list[str]] = {}
    for accession, label in labels.items():
        owners.setdefault(label, []).append(accession)
    shared = [f"{label} ({', '.join(a)})" for label, a in owners.items() if len(a) > 1]
    if shared:
        raise QueryConfigError(
            f"queries: the same label for two queries: {'; '.join(shared)}"
        )
    return labels
