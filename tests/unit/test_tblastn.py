"""tBLASTn hits of one query on one genome become one typed, tagged table."""

from pathlib import Path

import pandas as pd
import pytest

import tblastn

ROW = "NP_1.1\tchr1\t55.5\t40\t10\t1\t5\t44\t1300\t1181\t1e-9\t61.2\tMKV*L-A"
OTHER = "NP_1.1\tchr2\t30.0\t20\t14\t0\t50\t69\t100\t159\t2.5\t25.0\tMKVLA"


def test_the_search_writes_an_archive_at_the_configured_e_value():
    search = tblastn.Search(
        "tblastn", "blast_formatter", Path("db/X/X"), Path("q.fa"), 10.0
    )

    command = tblastn.search_command(search, Path("X.asn"))

    assert command[0] == "tblastn"
    assert command[command.index("-db") + 1] == "db/X/X"
    assert command[command.index("-evalue") + 1] == "10.0"
    assert command[command.index("-outfmt") + 1] == "11"
    assert command[-1] == "X.asn"


def test_every_numeric_column_is_typed():
    hits = tblastn.parse_hits(ROW + "\n")

    integers = (
        "Alignment Length",
        "Mismatches",
        "Gap Openings",
        "Q. Start",
        "Q. End",
        "S. Start",
        "S. End",
    )
    for column in integers:
        assert str(hits[column].dtype) == "int64", column
    for column in ("Pct Identity", "E-value", "Bit Score"):
        assert str(hits[column].dtype) == "float64", column
    assert hits.loc[0, "S. Start"] == 1300
    assert hits.loc[0, "Subject Sequence"] == "MKV*L-A"


def test_no_output_is_an_empty_table_with_every_column_typed():
    hits = tblastn.parse_hits("")

    assert len(hits) == 0
    assert list(hits.columns) == list(tblastn.COLUMNS)
    assert str(hits["S. Start"].dtype) == "int64"


def test_a_row_with_the_wrong_number_of_fields_is_refused():
    with pytest.raises(ValueError, match="13"):
        tblastn.parse_hits("NP_1.1\tchr1\t55.5\n")


def test_hits_get_their_provenance_a_tag_and_the_stop_flag():
    hits = tblastn.describe(
        tblastn.parse_hits(f"{ROW}\n{OTHER}\n"), "X", "Species x", "NP_1.1"
    )

    assert list(hits.columns) == list(tblastn.HIT_COLUMNS)
    assert set(hits["Species"]) == {"X"}
    assert set(hits["Species_Name"]) == {"Species x"}
    assert hits["Tag"].str.len().eq(12).all()
    assert hits["Tag"].is_unique
    assert list(hits["Internal_Stop"]) == [True, False]


def test_a_hit_reported_twice_is_kept_once():
    hits = tblastn.describe(
        tblastn.parse_hits(f"{ROW}\n{ROW}\n"), "X", "Species x", "NP_1.1"
    )

    assert len(hits) == 1


def test_an_empty_table_stays_typed_through_describe():
    hits = tblastn.describe(tblastn.parse_hits(""), "X", "Species x", "NP_1.1")

    assert list(hits.columns) == list(tblastn.HIT_COLUMNS)
    assert str(hits["S. End"].dtype) == "int64"
    assert hits["Internal_Stop"].dtype == bool


def test_the_same_hit_gets_the_same_tag_in_every_run():
    first = tblastn.describe(tblastn.parse_hits(ROW + "\n"), "X", "Species x", "NP_1.1")
    again = tblastn.describe(tblastn.parse_hits(ROW + "\n"), "X", "Species x", "NP_1.1")

    assert list(first["Tag"]) == list(again["Tag"])


def test_main_runs_both_tools_and_writes_parquet_and_csv(tmp_path, monkeypatch):
    calls = []

    def fake_run(command, *, tool):
        calls.append(tool)
        return f"{ROW}\n" if tool == "blast_formatter" else ""

    monkeypatch.setattr(tblastn.external, "run_tool", fake_run)
    monkeypatch.setattr(tblastn.defaults, "SPECIES_DICT", {"X": "Species x"})
    output = tmp_path / "X.parquet"

    status = tblastn.main(
        ["--species", "X", "--query-accession", "NP_1.1",
         "--query-file", "q.fa", "--output-parquet", str(output),
         "--output-csv", str(tmp_path / "X.csv"), "--output-asn", str(tmp_path / "X.asn")]
    )  # fmt: skip

    assert status == 0
    assert calls == ["tblastn", "blast_formatter"]
    assert len(pd.read_parquet(output)) == 1
    assert (tmp_path / "X.csv").read_text().count("\n") == 2


def test_a_failing_search_fails_the_job_and_writes_no_table(tmp_path, monkeypatch):
    def failing(command, *, tool):
        raise tblastn.external.ToolError(tool, 3, "BLAST Database error")

    monkeypatch.setattr(tblastn.external, "run_tool", failing)
    monkeypatch.setattr(tblastn.defaults, "SPECIES_DICT", {"X": "Species x"})
    output = tmp_path / "X.parquet"

    with pytest.raises(tblastn.external.ToolError):
        tblastn.main(
            ["--species", "X", "--query-accession", "NP_1.1",
             "--query-file", "q.fa", "--output-parquet", str(output),
             "--output-csv", str(tmp_path / "X.csv"), "--output-asn", str(tmp_path / "X.asn")]
        )  # fmt: skip

    assert not output.exists()
