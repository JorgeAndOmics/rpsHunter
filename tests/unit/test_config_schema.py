"""The shipped example config must satisfy the shipped schema.

A config that the repository itself cannot validate would fail every new user
at their first launch.
"""

from pathlib import Path

import yamale

CONFIG_DIR = Path(__file__).resolve().parents[2] / "data" / "config"


def test_shipped_config_validates_against_schema():
    schema = yamale.make_schema(CONFIG_DIR / "schema.yaml")
    data = yamale.make_data(CONFIG_DIR / "config.yaml")

    yamale.validate(schema, data, strict=True)  # raises YamaleError on a mismatch
