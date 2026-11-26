"""The `queries:` block maps protein accessions to labels that name folders."""

import pytest

import query_config


def test_queries_keep_their_order_and_get_folder_safe_labels():
    labels = query_config.parse_queries(
        {"NP_1.1": "human kinase", "NP_2.1": "mouse/rat kinase"}
    )

    assert labels == {"NP_1.1": "human_kinase", "NP_2.1": "mouse_rat_kinase"}
    assert list(labels) == ["NP_1.1", "NP_2.1"]


@pytest.mark.parametrize("raw", [None, {}, [], "NP_1.1"])
def test_a_config_without_queries_is_refused(raw):
    with pytest.raises(query_config.QueryConfigError, match="at least one"):
        query_config.parse_queries(raw)


def test_two_queries_whose_labels_name_the_same_folder_are_refused():
    with pytest.raises(query_config.QueryConfigError, match="same label") as raised:
        query_config.parse_queries({"NP_1.1": "my query", "NP_2.1": "my_query"})

    assert "NP_1.1" in str(raised.value)
    assert "NP_2.1" in str(raised.value)


def test_an_empty_label_is_refused():
    with pytest.raises(query_config.QueryConfigError, match="NP_1"):
        query_config.parse_queries({"NP_1.1": "  "})
