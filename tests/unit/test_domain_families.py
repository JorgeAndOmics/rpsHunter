"""Domain families are defined by the accessions of the models that evidence them."""

import pytest

import domain_families as df

RAW = {
    "KRAB": {"cdd": ["pfam01352", "cd07765"], "pfam": ["PF01352"]},
    "zinc_finger": {
        "cdd": ["pfam00096", "COG5048"],
        "pfam": ["PF00096"],
        "repeat": True,
    },
    "SSXRD": {"cdd": ["pfam09514"]},
}


def test_families_keep_their_models_and_the_order_of_the_config():
    families = df.parse_families(RAW)

    assert list(families) == ["KRAB", "zinc_finger", "SSXRD"]
    assert families["KRAB"].cdd == ("pfam01352", "cd07765")
    assert families["KRAB"].pfam == ("PF01352",)
    assert families["KRAB"].repeat is False
    assert families["zinc_finger"].repeat is True
    assert families["SSXRD"].pfam == ()


def test_each_model_maps_back_to_its_family():
    families = df.parse_families(RAW)

    assert df.family_by_cdd(families) == {
        "pfam01352": "KRAB",
        "cd07765": "KRAB",
        "pfam00096": "zinc_finger",
        "COG5048": "zinc_finger",
        "pfam09514": "SSXRD",
    }
    assert df.family_by_pfam(families) == {"PF01352": "KRAB", "PF00096": "zinc_finger"}


def test_cdd_accessions_are_listed_once_in_config_order():
    families = df.parse_families(RAW)

    assert df.cdd_accessions(families) == [
        "pfam01352",
        "cd07765",
        "pfam00096",
        "COG5048",
        "pfam09514",
    ]


@pytest.mark.parametrize("raw", [None, {}, [], "KRAB"])
def test_a_config_without_families_is_refused(raw):
    with pytest.raises(df.DomainConfigError, match="at least one family"):
        df.parse_families(raw)


def test_a_family_without_any_model_is_refused_by_name():
    with pytest.raises(df.DomainConfigError, match="SET"):
        df.parse_families({"SET": {"repeat": False}})


def test_a_model_claimed_by_two_families_is_refused():
    raw = {"A": {"cdd": ["pfam00096"]}, "B": {"cdd": ["pfam00096"]}}

    with pytest.raises(df.DomainConfigError, match="pfam00096") as raised:
        df.parse_families(raw)

    assert "A" in str(raised.value)
    assert "B" in str(raised.value)


def test_a_pfam_accession_carries_no_version_and_no_name():
    for bad in ("PF01352.33", "KRAB", "pfam01352"):
        with pytest.raises(df.DomainConfigError, match="PF"):
            df.parse_families({"KRAB": {"pfam": [bad]}})


def test_an_unknown_key_in_a_family_is_refused():
    with pytest.raises(df.DomainConfigError, match="profiles"):
        df.parse_families({"KRAB": {"cdd": ["pfam01352"], "profiles": ["KRAB"]}})


def test_all_problems_are_reported_together():
    raw = {"A": {}, "B": {"pfam": ["KRAB"]}}

    with pytest.raises(df.DomainConfigError) as raised:
        df.parse_families(raw)

    assert "A" in str(raised.value)
    assert "B" in str(raised.value)


def test_superfamilies_inherit_the_family_of_their_member_models():
    families = df.parse_families(RAW)
    links = {"pfam01352": "cl02581", "cd07765": "cl02581", "pfam00096": "cl23491"}

    assert df.family_by_superfamily(families, links) == {
        "cl02581": "KRAB",
        "cl23491": "zinc_finger",
    }


def test_a_superfamily_shared_by_two_families_is_refused():
    families = df.parse_families(RAW)
    links = {"pfam01352": "cl00001", "pfam00096": "cl00001"}

    with pytest.raises(df.DomainConfigError, match="cl00001"):
        df.family_by_superfamily(families, links)


def test_pfam_accessions_resolve_to_profile_names_whatever_the_version():
    families = df.parse_families(RAW)
    index = [("KRAB", "PF01352.33"), ("zf-C2H2", "PF00096.32"), ("SET", "PF00856.34")]

    assert df.pfam_profile_names(families, index) == {
        "PF01352": "KRAB",
        "PF00096": "zf-C2H2",
    }


def test_a_pfam_accession_missing_from_the_library_is_named():
    families = df.parse_families(RAW)

    with pytest.raises(df.DomainConfigError, match="PF00096"):
        df.pfam_profile_names(families, [("KRAB", "PF01352.33")])
