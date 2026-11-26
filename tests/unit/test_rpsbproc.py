"""rpsbproc is told the same cut-off as the search it post-processes."""

from pathlib import Path

import rpsbproc

TOOL = rpsbproc.Tool(
    program="rpsbproc", data_dir=Path("/db/rpsbproc_dbs"), max_evalue=1.0
)


def test_the_command_carries_the_cut_off_and_both_report_types():
    command = rpsbproc.command(TOOL, Path("in.asn"), Path("out.txt"))

    assert command[command.index("-e") + 1] == "1.0"
    assert command[command.index("-t") + 1] == "both"
    assert command[command.index("-i") + 1] == "in.asn"
    assert command[command.index("-d") + 1] == "/db/rpsbproc_dbs"
    assert command[command.index("-o") + 1] == "out.txt"


def test_an_empty_archive_gives_an_empty_report_without_running_the_tool(tmp_path):
    archive = tmp_path / "empty.asn"
    archive.touch()
    report = tmp_path / "report.txt"
    missing_tool = TOOL._replace(program="no-such-tool")

    ran = rpsbproc.process(missing_tool, archive, report)

    assert ran is False
    assert report.read_text() == ""


def test_a_non_empty_archive_is_handed_to_rpsbproc(tmp_path, monkeypatch):
    archive = tmp_path / "X.asn"
    archive.write_text("Blast4-archive")
    seen = {}

    def fake_run(command, *, tool):
        seen["command"], seen["tool"] = command, tool
        return ""

    monkeypatch.setattr(rpsbproc.external, "run_tool", fake_run)

    ran = rpsbproc.process(TOOL, archive, tmp_path / "X.txt")

    assert ran is True
    assert seen["tool"] == "rpsbproc"
    assert str(archive) in seen["command"]


def test_main_reads_its_settings_from_the_config(tmp_path, monkeypatch, capsys):
    archive = tmp_path / "X.asn"
    archive.touch()
    report = tmp_path / "X.txt"

    status = rpsbproc.main(["--input-asn", str(archive), "--output-txt", str(report)])

    assert status == 0
    assert report.exists()
    assert "empty archive" in capsys.readouterr().out
