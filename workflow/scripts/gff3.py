"""Domain instances as GFF3 features.

One feature per instance: the sequence, the domain's own span and strand, the
E-value of its best call as score, and attributes naming the family, the model
and the strength of evidence. Type `protein_match` (SO:0000349) is the Sequence
Ontology term for a match to a protein or protein model.

This module is pure: it returns the lines; writing them is the caller's job.
"""

from __future__ import annotations

import pandas as pd

HEADER = "##gff-version 3"
SOURCE = "rpsHunter"
FEATURE_TYPE = "protein_match"

# GFF3 reserves these in attribute values; they and control characters are escaped.
_RESERVED = frozenset(";=&,%")


def _value(text: object) -> str:
    """Percent-encode an attribute value as GFF3 requires."""
    return "".join(
        f"%{ord(char):02X}" if char in _RESERVED or ord(char) < 32 else char
        for char in str(text)
    )


def _attributes(row: pd.Series, number: int) -> str:
    """The ninth column for one instance."""
    queries = ",".join(_value(query) for query in str(row["Queries"]).split(","))
    pairs = (
        ("ID", f"domain_{number}"),
        ("Name", _value(row["Family"])),
        ("model", _value(row["Model"])),
        ("model_name", _value(row["Model_Name"])),
        ("hit_type", _value(row["Hit_Type"])),
        ("evidence", _value(row["Evidence"])),
        ("calls", str(row["N_Calls"])),
        ("queries", queries),
    )
    return ";".join(f"{key}={value}" for key, value in pairs)


def lines(instances: pd.DataFrame) -> list[str]:
    """Return the GFF3 lines for the instances of one species, header first.

    An empty table gives the header alone, which is a valid GFF3 file: a species
    without domains still gets its track.
    """
    result = [HEADER]
    for number, (_, row) in enumerate(instances.iterrows(), start=1):
        columns = (
            row["Sequence"],
            SOURCE,
            FEATURE_TYPE,
            str(row["Start"]),
            str(row["End"]),
            f"{row['Evalue']:g}",
            row["Strand"],
            ".",
            _attributes(row, number),
        )
        result.append("\t".join(columns))
    return result
