"""Tags identify a BLAST hit by what it is, so they are the same on every run."""

import re

import pandas as pd
import pytest

import tags

HIT = tags.HitIdentity(
    query_accession="NP_064612.2",
    species="Homo_sapiens",
    subject_id="CM000667.2",
    subject_start=23526232,
    subject_end=23527239,
    query_start=524,
    query_end=859,
)


def _frame(hits: list[tags.HitIdentity]) -> pd.DataFrame:
    """Build a BLAST-like hit table from identities."""
    return pd.DataFrame(hits, columns=list(tags.IDENTITY_COLUMNS))


def test_the_same_hit_always_gets_the_same_tag():
    assert tags.hit_tag(HIT) == tags.hit_tag(HIT)


def test_a_tag_is_twelve_characters_that_survive_a_fasta_header():
    tag = tags.hit_tag(HIT)

    assert re.fullmatch(r"[A-Z2-7]{12}", tag)
    found = tags.TAG_IN_IDENTIFIER.search(f"CM000667.2:1-9|tag:{tag}")
    assert found is not None
    assert found.group(1) == tag


@pytest.mark.parametrize(
    ("field", "other"),
    [
        ("query_accession", "NP_659058.3"),
        ("species", "Mus_musculus"),
        ("subject_id", "CM000668.2"),
        ("subject_start", 23526233),
        ("subject_end", 23527240),
        ("query_start", 525),
        ("query_end", 860),
    ],
)
def test_changing_any_part_of_the_hit_changes_the_tag(field, other):
    assert tags.hit_tag(HIT._replace(**{field: other})) != tags.hit_tag(HIT)


def test_strand_is_part_of_the_identity():
    # A minus-strand hit is written with start > end; swapping them is another hit.
    flipped = HIT._replace(subject_start=HIT.subject_end, subject_end=HIT.subject_start)

    assert tags.hit_tag(flipped) != tags.hit_tag(HIT)


def test_fields_cannot_run_into_each_other():
    # "ab" + "c" and "a" + "bc" must not hash alike.
    one = tags.hit_tag(HIT._replace(species="ab", subject_id="c"))
    two = tags.hit_tag(HIT._replace(species="a", subject_id="bc"))

    assert one != two


def test_tag_hits_tags_every_row_from_its_own_columns():
    second = HIT._replace(query_start=552, query_end=887)
    frame = _frame([HIT, second])

    result = tags.tag_hits(frame)

    assert list(result) == [tags.hit_tag(HIT), tags.hit_tag(second)]
    assert result.index.equals(frame.index)


def test_an_empty_hit_table_gets_no_tags():
    result = tags.tag_hits(_frame([HIT]).iloc[0:0])

    assert len(result) == 0


def test_two_rows_describing_the_same_hit_are_refused():
    with pytest.raises(tags.DuplicateHitError, match=r"CM000667\.2"):
        tags.tag_hits(_frame([HIT, HIT]))


def test_different_hits_with_the_same_tag_are_refused(monkeypatch):
    monkeypatch.setattr(tags, "hit_tag", lambda _hit: "AAAAAAAAAAAA")
    second = HIT._replace(query_start=552)

    with pytest.raises(tags.TagCollisionError, match="AAAAAAAAAAAA"):
        tags.tag_hits(_frame([HIT, second]))
