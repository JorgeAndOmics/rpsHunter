"""Each species names the genome it means, so a config alone determines the data."""

import pytest

import species_config

SHA = "c29768b9b8c97fe1d5169f5b9c5fd426d5336bcf5d621245678c30886f256ddc"


def test_a_species_with_an_assembly_accession_is_accepted():
    parsed = species_config.parse_species(
        {"Mus_musculus": {"name": "Mus musculus", "assembly": "GCA_000001635.9"}}
    )

    entry = parsed["Mus_musculus"]
    assert entry.name == "Mus musculus"
    assert entry.assembly == "GCA_000001635.9"
    assert entry.sha256 is None


def test_a_local_genome_is_pinned_by_its_checksum_instead():
    parsed = species_config.parse_species(
        {"Local_genome": {"name": "Local genome", "sha256": SHA}}
    )

    assert parsed["Local_genome"].assembly is None
    assert parsed["Local_genome"].sha256 == SHA


def test_both_pins_may_be_given_together():
    entry = {"name": "Homo sapiens", "assembly": "GCF_000001405.40", "sha256": SHA}

    parsed = species_config.parse_species({"Homo_sapiens": entry})

    assert parsed["Homo_sapiens"].assembly == "GCF_000001405.40"
    assert parsed["Homo_sapiens"].sha256 == SHA


def test_species_keep_the_order_of_the_config():
    entries = {
        key: {"name": key, "assembly": "GCA_000001635.9"} for key in ("b", "a", "c")
    }

    assert list(species_config.parse_species(entries)) == ["b", "a", "c"]


def test_no_species_is_an_empty_mapping():
    assert species_config.parse_species({}) == {}
    assert species_config.parse_species(None) == {}


def test_a_species_without_any_pin_is_refused_by_name():
    with pytest.raises(species_config.SpeciesConfigError, match="Desmodus_rotundus"):
        species_config.parse_species(
            {"Desmodus_rotundus": {"name": "Desmodus rotundus"}}
        )


@pytest.mark.parametrize(
    "accession",
    ["GCA_000001635", "GCA_000001635.latest", "Mus musculus", "NC_000067.7", "10090"],
)
def test_an_assembly_must_be_a_versioned_accession(accession):
    with pytest.raises(species_config.SpeciesConfigError, match="versioned"):
        species_config.parse_species({"X": {"name": "X", "assembly": accession}})


def test_a_checksum_must_be_a_sha256_digest():
    with pytest.raises(species_config.SpeciesConfigError, match="sha256"):
        species_config.parse_species({"X": {"name": "X", "sha256": "abc123"}})


def test_the_old_name_only_form_is_refused_with_the_new_form_spelled_out():
    with pytest.raises(species_config.SpeciesConfigError) as raised:
        species_config.parse_species({"Homo_sapiens": "Homo sapiens"})

    message = str(raised.value)
    assert "Homo_sapiens" in message
    assert "name:" in message
    assert "assembly:" in message


def test_an_entry_without_a_display_name_is_refused():
    with pytest.raises(species_config.SpeciesConfigError, match="name"):
        species_config.parse_species({"X": {"assembly": "GCA_000001635.9"}})


def test_all_problems_are_reported_together():
    with pytest.raises(species_config.SpeciesConfigError) as raised:
        species_config.parse_species(
            {"A": {"name": "A"}, "B": {"name": "B", "assembly": "latest"}}
        )

    assert "A" in str(raised.value)
    assert "B" in str(raised.value)


def test_a_species_block_that_is_a_list_is_refused_with_the_new_form():
    with pytest.raises(species_config.SpeciesConfigError, match="assembly:"):
        species_config.parse_species(["Homo_sapiens", "Mus_musculus"])


def test_a_genome_matching_its_checksum_passes():
    entry = species_config.SpeciesEntry("X", None, "a" * 64)

    assert species_config.checksum_problem("X", entry, "a" * 64) is None


def test_a_genome_that_differs_from_its_checksum_is_reported_with_both_values():
    entry = species_config.SpeciesEntry("X", None, "a" * 64)

    problem = species_config.checksum_problem("X", entry, "b" * 64)

    assert problem is not None
    assert "a" * 64 in problem
    assert "b" * 64 in problem


def test_a_species_without_a_checksum_has_nothing_to_compare():
    entry = species_config.SpeciesEntry("X", "GCA_000001635.9", None)

    assert species_config.checksum_problem("X", entry, "b" * 64) is None
