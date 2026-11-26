"""Each domain call says how strong its evidence is, on NCBI's E-value scale."""

import pandas as pd
import pytest

import domain_calls


@pytest.mark.parametrize(
    ("evalue", "tier"),
    [
        (0.0, "standard"),
        (1e-30, "standard"),
        (0.01, "standard"),
        (0.011, "borderline"),
        (1.0, "borderline"),
        (1.01, "rescue"),
        (1000.0, "rescue"),
    ],
)
def test_the_tier_follows_ncbi_thresholds(evalue, tier):
    assert domain_calls.evidence_tier(evalue) == tier


def test_a_negative_e_value_is_refused():
    with pytest.raises(ValueError, match="negative"):
        domain_calls.evidence_tier(-1.0)


REPORT = """#Post-RPSBLAST Processing Utility
#DATA
SESSION\t1\tRPSBLAST 2.17.0+\tfamilies\tBLOSUM62\t1
QUERY\tQuery_1\tPeptide\t513\tCM000667.2:23526232-23527770|tag:YUJR7CCRYKBA
DOMAINS
1\tQuery_1\tNon-specific\t512735\t305\t506\t2.85669e-70\t231.372\tCOG5048\tZMS1\tN\t-
1\tQuery_1\tSuperfamily\t0\t100\t119\t0.5\t18.3394\tcl40432\tSET\tC\t-
ENDDOMAINS
SITES
ENDSITES
ENDQUERY
QUERY\tQuery_2\tPeptide\t40\tCM000667.2:900-781|tag:AAAAAAAAAAAA
DOMAINS
1\tQuery_2\tSpecific\t513144\t3\t25\t1e-9\t40.1\tpfam00096\tzf-C2H2\t-\t-
ENDDOMAINS
ENDQUERY
QUERY\tQuery_3\tPeptide\t33\tCM000667.2:5000-5098|tag:BBBBBBBBBBBB
DOMAINS
ENDDOMAINS
ENDQUERY
"""


def test_a_report_gives_one_row_per_call_with_numbers_as_numbers():
    calls = domain_calls.parse_report(REPORT)

    assert len(calls) == 3
    first = calls.iloc[0]
    assert first["Tag"] == "YUJR7CCRYKBA"
    assert first["Sequence"] == "CM000667.2"
    assert first["Model"] == "COG5048"
    assert first["Hit_Type"] == "Non-specific"
    assert first["Residue_From"] == 305
    assert first["Evalue"] == pytest.approx(2.85669e-70)
    assert calls["PSSM_ID"].dtype.kind == "i"
    assert calls["Bitscore"].dtype.kind == "f"


def test_a_sequence_without_calls_gives_no_row():
    calls = domain_calls.parse_report(REPORT)

    assert "BBBBBBBBBBBB" not in set(calls["Tag"])


def test_an_empty_report_gives_an_empty_table_with_every_column():
    calls = domain_calls.parse_report("")

    assert len(calls) == 0
    assert list(calls.columns) == list(domain_calls.REPORT_COLUMNS)


def test_a_definition_line_that_is_not_a_tagged_hit_is_an_error():
    report = (
        "QUERY\tQuery_1\tPeptide\t10\tsomething else\nDOMAINS\nENDDOMAINS\nENDQUERY\n"
    )

    with pytest.raises(ValueError, match="something else"):
        domain_calls.parse_report(report)


def test_a_forward_hit_places_its_domain_three_bases_per_residue_from_the_start():
    # Hit 1001-1300 (100 residues); residues 11-20 are bases 1031-1060.
    assert domain_calls.genomic_span(1001, 1300, 11, 20) == (1031, 1060, "+")


def test_a_reverse_hit_counts_down_from_its_first_base():
    # Hit written 1300-1001: residue 1 is bases 1300-1298; residues 11-20 are 1241-1270.
    assert domain_calls.genomic_span(1300, 1001, 11, 20) == (1241, 1270, "-")


def test_the_whole_hit_maps_onto_itself():
    assert domain_calls.genomic_span(1001, 1300, 1, 100) == (1001, 1300, "+")
    assert domain_calls.genomic_span(1300, 1001, 1, 100) == (1001, 1300, "-")


def test_calls_are_placed_on_the_genome_with_hit_and_domain_coordinates():
    placed = domain_calls.place_on_genome(domain_calls.parse_report(REPORT))

    forward = placed.iloc[0]
    assert (forward["Hit_Start"], forward["Hit_End"], forward["Strand"]) == (
        23526232,
        23527770,
        "+",
    )
    assert forward["Start"] == 23526232 + 3 * 304
    assert forward["End"] == 23526232 + 3 * 506 - 1
    reverse = placed.iloc[2]
    assert (reverse["Hit_Start"], reverse["Hit_End"], reverse["Strand"]) == (
        781,
        900,
        "-",
    )
    assert (reverse["Start"], reverse["End"]) == (900 - 3 * 25 + 1, 900 - 3 * 2)


ALIGNMENTS = pd.DataFrame(
    {
        "Tag": ["YUJR7CCRYKBA", "YUJR7CCRYKBA", "YUJR7CCRYKBA", "AAAAAAAAAAAA"],
        "Subject ID": ["CDD:512735", "CDD:380970", "CDD:380914", "CDD:513144"],
        "Q. Start": [305, 100, 100, 3],
        "Q. End": [506, 119, 119, 25],
        "S. Start": [580, 60, 10, 1],
        "S. End": [785, 79, 29, 23],
        "Bit Score": [231.0, 18.3, 12.0, 40.1],
    }
)
LENGTHS = {512735: 792, 380970: 129, 380914: 72, 513144: 23}


def test_model_coordinates_come_from_the_matching_alignment():
    calls = domain_calls.parse_report(REPORT)

    with_model = domain_calls.add_model_coordinates(calls, ALIGNMENTS, LENGTHS)

    first = with_model.iloc[0]
    assert (first["Model_From"], first["Model_To"], first["Model_Length"]) == (
        580,
        785,
        792,
    )


def test_a_superfamily_call_takes_the_coordinates_of_its_best_member_alignment():
    calls = domain_calls.parse_report(REPORT)

    with_model = domain_calls.add_model_coordinates(calls, ALIGNMENTS, LENGTHS)

    superfamily = with_model.iloc[1]
    assert superfamily["Hit_Type"] == "Superfamily"
    assert (superfamily["Model_From"], superfamily["Model_To"]) == (60, 79)
    assert superfamily["Model_Length"] == 129


def test_a_call_without_a_matching_alignment_is_an_error():
    calls = domain_calls.parse_report(REPORT)

    with pytest.raises(ValueError, match="AAAAAAAAAAAA"):
        domain_calls.add_model_coordinates(calls, ALIGNMENTS.iloc[:3], LENGTHS)
