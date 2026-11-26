"""The checksum of a genome file is computed without loading it whole."""

import hashlib

import genome_check
import species_config


def test_file_sha256_matches_hashlib_on_a_file_larger_than_one_block(tmp_path):
    path = tmp_path / "data.bin"
    path.write_bytes(b"rpsHunter" * 200_000)  # 1.8 MB: spans two read blocks

    expected = hashlib.sha256(path.read_bytes()).hexdigest()
    assert genome_check.file_sha256(path) == expected


def test_file_sha256_of_an_empty_file(tmp_path):
    path = tmp_path / "empty.fa"
    path.touch()

    assert genome_check.file_sha256(path) == hashlib.sha256(b"").hexdigest()


def _pin(monkeypatch, sha256, assembly=None):
    """Make the config hold one species, X, with the given pins."""
    entry = species_config.SpeciesEntry("X", assembly, sha256)
    monkeypatch.setattr(genome_check.defaults, "SPECIES_ENTRIES", {"X": entry})


def _genome(tmp_path):
    path = tmp_path / "X.fa"
    path.write_text(">s\nACGT\n")
    return path


def test_a_genome_matching_its_pin_passes(tmp_path, monkeypatch, capsys):
    fasta = _genome(tmp_path)
    _pin(monkeypatch, genome_check.file_sha256(fasta))

    status = genome_check.main(["--species", "X", "--fasta", str(fasta)])

    assert status == 0
    assert "matches the config" in capsys.readouterr().out


def test_a_genome_differing_from_its_pin_fails_and_says_why(
    tmp_path, monkeypatch, capsys
):
    fasta = _genome(tmp_path)
    _pin(monkeypatch, "0" * 64)

    status = genome_check.main(["--species", "X", "--fasta", str(fasta)])

    assert status == 1
    assert "not the one the config pins" in capsys.readouterr().err


def test_a_species_pinned_by_accession_only_gets_its_checksum_printed(
    tmp_path, monkeypatch, capsys
):
    fasta = _genome(tmp_path)
    _pin(monkeypatch, None, assembly="GCA_000001635.9")

    status = genome_check.main(["--species", "X", "--fasta", str(fasta)])

    assert status == 0
    output = capsys.readouterr().out
    assert genome_check.file_sha256(fasta) in output
    assert "not pinned in the config" in output
