"""RPS-BLAST runs on NCBI's E-value scale, and its table is typed and tagged."""

from pathlib import Path

import pandas as pd
import pytest

import rpsblast

SETTINGS = rpsblast.SearchSettings(
    max_evalue=1.0,
    database_residues=15_821_058,
    comp_based_stats=0,
    seg="no",
    window_size=0,
)
ROW = (
    "CM000667.2:23526232-23527770|tag:YUJR7CCRYKBA\tCDD:512735\t45.146\t206\t109\t1"
    "\t305\t506\t580\t785\t2.86e-70\t231\tCOG5048, ZMS1, Zn-finger protein\n"
)


def _command(settings=SETTINGS):
    return rpsblast.search_command(
        "rpsblast", Path("/db/families_db"), Path("q.fa"), Path("out.asn"), settings
    )


def test_the_search_states_the_size_of_the_full_database():
    command = _command()

    assert command[command.index("-dbsize") + 1] == "15821058"


def test_the_search_uses_the_configured_cut_off():
    command = _command()

    assert command[command.index("-evalue") + 1] == "1.0"


def test_the_search_writes_an_archive_for_rpsbproc():
    command = _command()

    assert command[command.index("-outfmt") + 1] == "11"
    assert command[command.index("-out") + 1] == "out.asn"


def test_a_window_size_of_zero_is_passed_on_as_the_one_hit_algorithm():
    command = _command()

    assert command[command.index("-window_size") + 1] == "0"


def test_the_table_has_numbers_as_numbers_and_the_tag_in_its_own_column():
    table = rpsblast.parse_table(ROW)

    row = table.iloc[0]
    assert row["Query ID"] == "CM000667.2:23526232-23527770"
    assert row["Tag"] == "YUJR7CCRYKBA"
    assert row["Q. Start"] == 305
    assert row["S. End"] == 785
    assert row["E-value"] == 2.86e-70
    assert table["Q. Start"].dtype.kind == "i"
    assert table["Bit Score"].dtype.kind == "f"


def test_no_alignments_gives_an_empty_table_with_every_column():
    table = rpsblast.parse_table("")

    assert len(table) == 0
    assert list(table.columns) == [*rpsblast.TABLE_COLUMNS, "Tag"]


def test_a_row_with_the_wrong_number_of_fields_is_an_error():

    with pytest.raises(ValueError, match="13"):
        rpsblast.parse_table("a\tb\tc\n")


def test_an_empty_fasta_is_not_searched(tmp_path, monkeypatch):
    fasta = tmp_path / "X.fa"
    fasta.touch()
    archive = tmp_path / "X.asn"
    monkeypatch.setattr(rpsblast.external, "run_tool", _refuse_to_run)

    table = rpsblast.search(fasta, tmp_path / "db", archive)

    assert len(table) == 0
    assert archive.read_bytes() == b""


def test_a_search_runs_rpsblast_then_formats_its_archive(tmp_path, monkeypatch):
    fasta = tmp_path / "X.fa"
    fasta.write_text(">s|tag:YUJR7CCRYKBA\nMKV\n")
    calls = []

    def fake_run(command, *, tool):
        calls.append(tool)
        return ROW if tool == "blast_formatter" else ""

    monkeypatch.setattr(rpsblast.external, "run_tool", fake_run)

    table = rpsblast.search(fasta, tmp_path / "db", tmp_path / "X.asn")

    assert calls == ["rpsblast", "blast_formatter"]
    assert list(table["Tag"]) == ["YUJR7CCRYKBA"]


def test_main_writes_the_table_with_its_provenance_columns(tmp_path, monkeypatch):
    monkeypatch.setattr(rpsblast, "search", lambda *_: rpsblast.parse_table(ROW))
    monkeypatch.setattr(rpsblast.defaults, "SPECIES_DICT", {"X": "Species x"})
    output = tmp_path / "X.parquet"
    arguments = [
        "--species", "X", "--query-accession", "NP_1.1", "--fasta-input", "X.fa",
        "--output-parquet", str(output), "--output-asn", "X.asn", "--db", "db",
    ]  # fmt: skip

    status = rpsblast.main(arguments)

    written = pd.read_parquet(output)
    assert status == 0
    assert list(written.columns[-4:]) == [
        "Species",
        "Species_Name",
        "Tag",
        "Query_Accession",
    ]
    assert written.loc[0, "Species_Name"] == "Species x"
    assert written.loc[0, "Query_Accession"] == "NP_1.1"


def _refuse_to_run(*_args, **_kwargs):
    raise AssertionError("no external tool should run")
