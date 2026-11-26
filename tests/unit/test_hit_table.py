"""A stop codon inside a hit's translation marks decayed coding sequence."""

import pandas as pd

import hit_table


def test_a_stop_before_the_last_residue_is_internal():
    sequences = pd.Series(["MKV*LLA", "MKVLLA", "MKVLLA*", "*MKV", ""])

    flags = hit_table.has_internal_stop(sequences)

    assert list(flags) == [True, False, False, True, False]


def test_alignment_gaps_do_not_hide_or_create_a_stop():
    sequences = pd.Series(["MK-V*L-LA", "MK--VLLA", "MKV-*"])

    assert list(hit_table.has_internal_stop(sequences)) == [True, False, False]


def test_the_result_is_boolean_and_keeps_the_index():
    sequences = pd.Series(["MK*V", "MKV"], index=[7, 9])

    flags = hit_table.has_internal_stop(sequences)

    assert flags.dtype == bool
    assert list(flags.index) == [7, 9]
