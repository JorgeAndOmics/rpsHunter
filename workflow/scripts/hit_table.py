"""What is derived for each tBLASTn hit from its own row.

This module is pure.
"""

from __future__ import annotations

import pandas as pd

STOP = "*"
GAP = "-"


def has_internal_stop(sequences: pd.Series) -> pd.Series:
    """Say, for each hit, whether its translation has a stop codon inside it.

    tBLASTn reports the subject as translated protein, with `*` for a stop codon
    and `-` for an alignment gap. A stop before the last residue means the
    region does not encode one uninterrupted peptide: a pseudogene, a frameshift,
    or an alignment that ran past the end of an exon. A stop as the very last
    residue is the natural end of a coding sequence and does not count.

    Args:
        sequences: Subject sequences as tBLASTn writes them.

    Returns:
        A boolean Series on the same index.
    """
    residues = sequences.fillna("").str.replace(GAP, "", regex=False)
    return residues.str[:-1].str.contains(STOP, regex=False).astype(bool)
