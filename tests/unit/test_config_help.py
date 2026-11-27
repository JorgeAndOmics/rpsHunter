"""--config-help answers from the configuration reference itself."""

from pathlib import Path

import pytest
import yaml

import config_help

REPO = Path(__file__).parents[2]
REFERENCE = (REPO / "docs" / "configuration.md").read_text()
SAMPLE = """\
## Parameter Reference

### `alpha` -- First Things

About alpha.

| Parameter | Type | Default | Validation | Description |
|-----------|------|---------|------------|-------------|
| `size` | `int` | `3` | `int(min=0)` | How big, with a `code` word and a [link](other.md). |
| `mode` | `str` | `'no'` | `str(matches='^(yes\\|no)$')` | Which way. |

### `beta` -- Second Things (Optional)

About beta.

## Configuration Examples

### Not a section
"""


def test_sections_are_listed_with_their_titles_in_document_order():
    listing = config_help.explain(SAMPLE, "")

    assert listing.splitlines()[1:3] == [
        "  alpha  First Things",
        "  beta   Second Things (Optional)",
    ]
    assert "Not a section" not in listing


def test_a_section_is_shown_whole():
    text = config_help.explain(SAMPLE, "alpha")

    assert text.startswith("alpha: First Things")
    assert "About alpha." in text
    assert "size" in text
    assert "mode" in text
    assert "About beta" not in text


def test_one_parameter_is_shown_with_type_default_and_meaning():
    text = config_help.explain(SAMPLE, "alpha.size")

    assert text.splitlines() == [
        "alpha.size",
        "  type:        int",
        "  default:     3",
        "  validation:  int(min=0)",
        "  How big, with a code word and a link.",
    ]


def test_an_escaped_bar_inside_a_cell_does_not_split_it():
    text = config_help.explain(SAMPLE, "alpha.mode")

    assert "  validation:  str(matches='^(yes|no)$')" in text
    assert text.splitlines()[-1] == "  Which way."


def test_an_unknown_section_or_parameter_names_what_exists():
    with pytest.raises(config_help.UnknownKeyError, match="alpha, beta"):
        config_help.explain(SAMPLE, "gamma")
    with pytest.raises(config_help.UnknownKeyError, match="size, mode"):
        config_help.explain(SAMPLE, "alpha.colour")


def test_every_section_of_the_shipped_config_is_in_the_reference():
    config = yaml.safe_load((REPO / "data" / "config" / "config.yaml").read_text())

    documented = set(config_help.sections(REFERENCE))

    assert set(config) <= documented


def test_every_documented_parameter_of_the_real_reference_can_be_explained():
    for name, section in config_help.sections(REFERENCE).items():
        for parameter in config_help.parameters(section.body):
            text = config_help.explain(REFERENCE, f"{name}.{parameter}")

            assert text.startswith(f"{name}.{parameter}")
            assert "|" not in text.splitlines()[-1][:3]


def test_a_table_with_other_columns_is_read_by_its_own_header():
    reference = SAMPLE.replace(
        "### `beta` -- Second Things (Optional)\n\nAbout beta.\n",
        "### `beta` -- Second Things (Optional)\n\n"
        "| Parameter | Type | Required | Description |\n|---|---|---|---|\n"
        "| `name` | `str` | yes | What it is called. |\n",
    )

    text = config_help.explain(reference, "beta.name")

    assert text.splitlines() == [
        "beta.name",
        "  type:      str",
        "  required:  yes",
        "  What it is called.",
    ]


def test_no_parameter_of_the_real_reference_shows_a_description_as_a_fact():
    for name, section in config_help.sections(REFERENCE).items():
        for parameter, cells in config_help.parameters(section.body).items():
            assert "Description" in cells, f"{name}.{parameter}"
            assert all(
                len(cell) < 60
                for header, cell in cells.items()
                if header != "Description"
            )
