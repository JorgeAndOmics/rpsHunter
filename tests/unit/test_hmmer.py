"""HMMER describes the selected sequences on its own; it filters nothing."""

from pathlib import Path

import pandas as pd
import pytest

import hmmer

SETTINGS = hmmer.SearchSettings(
    use_gathering_threshold=False,
    evalue=10.0,
    dom_evalue=10.0,
    min_score=None,
    max_sensitivity=False,
    bias_filter=True,
    seed=67,
    threads=2,
)
#                 target                           tlen  query  accession   qlen  E    score bias  #  of  cE   iE      score bias hmm     ali      env      acc
DOMTBL = """# comment line
chr1:1300-1001|tag:AAAAAAAAAAAA -  100 KRAB PF01352.33 42 1e-9 40.0 0.1 1 1 2e-12 1e-9 39.5 0.1 1 40 12 51 10 53 0.95 description text
chr1:1300-1001|tag:AAAAAAAAAAAA -  100 SET  PF00856.34 105 0.5 9.0 0.0 1 1 0.02 0.5 8.0 0.0 60 80 70 90 68 92 0.80 another one
"""
HMM = """HMMER3/f
NAME  KRAB
ACC   PF01352.33
LENG  42
GA    20.40 20.40;
//
NAME  SET
ACC   PF00856.34
GA    22.40 22.40;
//
"""


def _command(settings=SETTINGS):
    return hmmer.search_command(
        Path("subset.hmm"), Path("x.fa"), Path("x.domtbl"), settings
    )


def test_the_search_is_seeded_and_writes_a_domain_table():
    command = _command()

    assert command[0] == "hmmsearch"
    assert command[command.index("--seed") + 1] == "67"
    assert command[command.index("--domtblout") + 1] == "x.domtbl"
    assert command[-2:] == ["subset.hmm", "x.fa"]


def test_explicit_e_values_are_used_unless_gathering_thresholds_are_asked_for():
    explicit, gathering = (
        _command(),
        _command(SETTINGS._replace(use_gathering_threshold=True)),
    )

    assert explicit[explicit.index("-E") + 1] == "10.0"
    assert "--cut_ga" not in explicit
    assert "--cut_ga" in gathering
    assert "-E" not in gathering


def test_optional_flags_follow_the_settings():
    command = _command(
        SETTINGS._replace(max_sensitivity=True, bias_filter=False, min_score=15.0)
    )

    assert "--max" in command
    assert "--nobias" in command
    assert command[command.index("-T") + 1] == "15.0"


def test_gathering_thresholds_are_read_per_profile_for_domains():
    assert hmmer.gathering_thresholds(HMM) == {"KRAB": 20.4, "SET": 22.4}


def test_the_domain_table_is_parsed_with_the_tag_and_aligned_residues():
    domains = hmmer.parse_domains(DOMTBL)

    first = domains.iloc[0]
    assert first["Tag"] == "AAAAAAAAAAAA"
    assert first["Profile"] == "KRAB"
    assert first["Pfam"] == "PF01352"
    assert (first["Residue_From"], first["Residue_To"]) == (12, 51)
    assert (first["Model_From"], first["Model_To"], first["Model_Length"]) == (
        1,
        40,
        42,
    )
    assert first["Evalue"] == pytest.approx(1e-9)
    assert first["Score"] == pytest.approx(39.5)


def test_an_empty_domain_table_gives_an_empty_frame_with_every_column():
    domains = hmmer.parse_domains("# only comments\n")

    assert len(domains) == 0
    assert list(domains.columns) == list(hmmer.RAW_COLUMNS)


def test_short_or_barely_covering_alignments_are_dropped():
    domains = hmmer.parse_domains(DOMTBL)

    kept = hmmer.keep_substantial(domains, min_coverage=0.5, min_length=10)

    assert list(kept["Profile"]) == ["KRAB"]  # SET covers 21 of 105 model positions


def test_each_domain_gets_its_family_and_whether_it_clears_the_gathering_threshold():
    domains = hmmer.parse_domains(DOMTBL)

    described = hmmer.describe(
        domains,
        {"PF01352": "KRAB_family", "PF00856": "SET_family"},
        {"KRAB": 20.4, "SET": 22.4},
    )

    assert list(described["Family"]) == ["KRAB_family", "SET_family"]
    assert list(described["Above_GA"]) == [True, False]


def test_a_domain_of_a_profile_outside_every_family_is_an_error():
    domains = hmmer.parse_domains(DOMTBL)

    with pytest.raises(ValueError, match="PF00856"):
        hmmer.describe(domains, {"PF01352": "KRAB_family"}, {"KRAB": 20.4, "SET": 22.4})


def test_an_empty_fasta_is_not_searched(tmp_path, monkeypatch):
    fasta = tmp_path / "X.fa"
    fasta.touch()
    monkeypatch.setattr(hmmer.external, "run_tool", _refuse_to_run)

    domains = hmmer.search(fasta, tmp_path / "subset.hmm", SETTINGS)

    assert len(domains) == 0


def test_a_search_reads_back_the_domain_table_hmmsearch_wrote(tmp_path, monkeypatch):
    fasta = tmp_path / "X.fa"
    fasta.write_text(">chr1:1300-1001|tag:AAAAAAAAAAAA\nMKV\n")

    def fake_run(command, *, tool):
        Path(command[command.index("--domtblout") + 1]).write_text(DOMTBL)
        return ""

    monkeypatch.setattr(hmmer.external, "run_tool", fake_run)

    domains = hmmer.search(fasta, tmp_path / "subset.hmm", SETTINGS)

    assert len(domains) == 2


def test_main_writes_the_described_domains_with_provenance(tmp_path, monkeypatch):
    profiles = tmp_path / "subset.hmm"
    profiles.write_text(HMM)
    monkeypatch.setattr(hmmer, "search", lambda *_: hmmer.parse_domains(DOMTBL))
    families = hmmer.domain_families.parse_families(
        {"KRAB_family": {"pfam": ["PF01352"]}, "SET_family": {"pfam": ["PF00856"]}}
    )
    monkeypatch.setattr(hmmer.defaults, "FAMILIES", families)
    monkeypatch.setattr(hmmer.defaults, "HMMER_MIN_COVERAGE", 0.1)
    monkeypatch.setattr(hmmer.defaults, "HMMER_MIN_ALN_LEN", 10)
    output = tmp_path / "X.parquet"

    status = hmmer.main(
        ["--species", "X", "--query-accession", "NP_1.1", "--fasta", "X.fa",
         "--profiles", str(profiles), "--output", str(output), "--threads", "1"]
    )  # fmt: skip

    written = pd.read_parquet(output)
    assert status == 0
    assert list(written.columns) == list(hmmer.DOMAIN_COLUMNS)
    assert set(written["Species"]) == {"X"}
    assert list(written["Above_GA"]) == [True, False]


def _refuse_to_run(*_args, **_kwargs):
    raise AssertionError("no external tool should run")


def test_the_configured_settings_come_from_the_hmmer_section(monkeypatch):
    monkeypatch.setattr(hmmer.defaults, "HMMER_SEED", 5)

    settings = hmmer.configured_settings(threads=3)

    assert (settings.seed, settings.threads) == (5, 3)


def test_a_malformed_domain_table_is_refused():
    with pytest.raises(ValueError, match="fields"):
        hmmer.parse_domains("too few fields\n")
    with pytest.raises(ValueError, match="target name"):
        hmmer.parse_domains(DOMTBL.replace("|tag:AAAAAAAAAAAA", ""))


def test_a_species_without_sequences_still_gets_a_table_with_every_column(
    tmp_path, monkeypatch
):
    profiles, fasta, output = (
        tmp_path / "subset.hmm",
        tmp_path / "X.fa",
        tmp_path / "X.parquet",
    )
    profiles.write_text(HMM)
    fasta.touch()
    monkeypatch.setattr(hmmer.external, "run_tool", _refuse_to_run)

    hmmer.main(
        ["--species", "X", "--query-accession", "NP_1.1", "--fasta", str(fasta),
         "--profiles", str(profiles), "--output", str(output)]
    )  # fmt: skip

    written = pd.read_parquet(output)
    assert len(written) == 0
    assert list(written.columns) == list(hmmer.DOMAIN_COLUMNS)
