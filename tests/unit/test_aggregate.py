"""Per-species tables are put together without changing a row."""

import pandas as pd

import aggregate


def _table(species, rows):
    return pd.DataFrame({"Species": [species] * rows, "Value": list(range(rows))})


def test_tables_are_concatenated_in_the_order_given(tmp_path):
    paths = []
    for species, rows in (("B", 2), ("A", 1)):
        path = tmp_path / f"{species}.parquet"
        _table(species, rows).to_parquet(path)
        paths.append(str(path))
    output = tmp_path / "all.parquet"

    status = aggregate.main(["--output", str(output), *paths])

    combined = pd.read_parquet(output)
    assert status == 0
    assert list(combined["Species"]) == ["B", "B", "A"]
    assert list(combined.index) == [0, 1, 2]
    assert (tmp_path / "all.csv").read_text().splitlines()[0] == "Species,Value"


def test_empty_tables_keep_the_columns_and_types(tmp_path):
    path = tmp_path / "empty.parquet"
    _table("A", 0).astype({"Value": "int64"}).to_parquet(path)
    output = tmp_path / "all.parquet"

    aggregate.main(["--output", str(output), str(path), str(path)])

    combined = pd.read_parquet(output)
    assert len(combined) == 0
    assert list(combined.columns) == ["Species", "Value"]
    assert str(combined["Value"].dtype) == "int64"
