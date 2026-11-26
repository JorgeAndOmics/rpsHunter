"""One species' call table: family, evidence tier, own coordinates, model coverage."""

import pandas as pd
import pytest

import call_domains

REPORT = (
    "QUERY\tQuery_1\tPeptide\t100\tchr1:1001-1300|tag:AAAAAAAAAAAA\n"
    "DOMAINS\n"
    "1\tQuery_1\tSpecific\t11\t11\t20\t1e-9\t40.1\tpfam00001\tmodel_one\t-\t-\n"
    "1\tQuery_1\tSuperfamily\t0\t50\t70\t0.5\t18.0\tcl00009\tcluster\tC\t-\n"
    "ENDDOMAINS\nENDQUERY\n"
)
ALIGNMENTS = pd.DataFrame(
    {
        "Tag": ["AAAAAAAAAAAA", "AAAAAAAAAAAA"],
        "Subject ID": ["CDD:11", "CDD:22"],
        "Q. Start": [11, 50],
        "Q. End": [20, 70],
        "S. Start": [1, 5],
        "S. End": [10, 25],
        "Bit Score": [40.1, 18.0],
    }
)
FAMILIES = {"first": {"cdd": ["pfam00001"]}, "second": {"cdd": ["cd00022"]}}


@pytest.fixture
def cdd_dir(tmp_path, monkeypatch):
    """A CDD annotation folder with two models, and a config that lists them."""
    (tmp_path / "cddid.tbl").write_text(
        "11\tpfam00001\tmodel_one\tdescription\t10\n22\tcd00022\tmodel_two\tdescription\t40\n"
    )
    (tmp_path / "family_superfamily_links").write_text(
        "pfam00001\t11\tcl00001\t901\ncd00022\t22\tcl00009\t909\n"
    )
    families = call_domains.domain_families.parse_families(FAMILIES)
    monkeypatch.setattr(call_domains.defaults, "FAMILIES", families)
    monkeypatch.setattr(call_domains.defaults, "SPECIES_DICT", {"X": "Species x"})
    monkeypatch.setitem(call_domains.defaults.PATH_DICT, "RPSBPROC_DB", tmp_path)
    return tmp_path


def test_each_call_gets_its_family_tier_position_and_model_coverage(cdd_dir):
    calls = call_domains.build(REPORT, ALIGNMENTS, cdd_dir)

    specific = calls.iloc[0]
    assert specific["Family"] == "first"
    assert specific["Evidence"] == "standard"
    assert (specific["Start"], specific["End"], specific["Strand"]) == (1031, 1060, "+")
    assert (specific["Model_From"], specific["Model_To"], specific["Model_Length"]) == (
        1,
        10,
        10,
    )


def test_a_superfamily_call_belongs_to_the_family_of_its_cluster(cdd_dir):
    calls = call_domains.build(REPORT, ALIGNMENTS, cdd_dir)

    superfamily = calls.iloc[1]
    assert superfamily["Family"] == "second"
    assert superfamily["Evidence"] == "borderline"
    assert superfamily["Model_Length"] == 40


def test_a_call_for_a_model_outside_every_family_is_an_error(cdd_dir):
    report = REPORT.replace("pfam00001", "pfam99999")

    with pytest.raises(ValueError, match="pfam99999"):
        call_domains.build(report, ALIGNMENTS, cdd_dir)


def test_main_writes_the_table_with_provenance_and_fixed_columns(cdd_dir, tmp_path):
    report = tmp_path / "X.txt"
    report.write_text(REPORT)
    alignments = tmp_path / "X.rps.parquet"
    ALIGNMENTS.to_parquet(alignments)
    output = tmp_path / "X.calls.parquet"

    status = call_domains.main(
        ["--species", "X", "--query-accession", "NP_1.1", "--report", str(report),
         "--alignments", str(alignments), "--output", str(output)]
    )  # fmt: skip

    written = pd.read_parquet(output)
    assert status == 0
    assert list(written.columns) == list(call_domains.CALL_COLUMNS)
    assert set(written["Species_Name"]) == {"Species x"}
    assert set(written["Query_Accession"]) == {"NP_1.1"}


def test_an_empty_report_gives_an_empty_table_with_the_same_columns(cdd_dir, tmp_path):
    report = tmp_path / "X.txt"
    report.write_text("")
    alignments = tmp_path / "X.rps.parquet"
    ALIGNMENTS.iloc[0:0].to_parquet(alignments)
    output = tmp_path / "X.calls.parquet"

    call_domains.main(
        ["--species", "X", "--query-accession", "NP_1.1", "--report", str(report),
         "--alignments", str(alignments), "--output", str(output)]
    )  # fmt: skip

    written = pd.read_parquet(output)
    assert len(written) == 0
    assert list(written.columns) == list(call_domains.CALL_COLUMNS)
