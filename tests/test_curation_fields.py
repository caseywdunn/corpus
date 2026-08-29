"""`doclang` and `pagemap`: recorded, round-tripped, and read by nothing (#214).

Curating a scanned library means writing down two things the pipeline has no
way to hold. `ocrlang` (#176) is an *instruction* to Tesseract — a pack list,
meaningful only when OCR runs. What was missing is the *fact*: this paper is
Russian; this one is 19th-c. German set in Fraktur; this one is Ancient Greek.
That is what a human or an annotation pass actually determines by looking, it
is worth recording for born-digital papers too, and it is in a different
vocabulary from Tesseract's. `pagemap` is the same idea for structure: a
`keeppages` range is unreviewable on its own, because a reader cannot tell
whether pages 1-2 were a scanner wrapper, a blank verso, or a mistake.

**Inertness is the feature, not a limitation.** Neither field is read at run
time, so neither can change output, so neither may appear in a fingerprint —
which is what makes them safe to correct. Fixing a typo in a note must not
reprocess a document. The instinct when adding a bib field is to copy the
`ocrlang` template wholesale, and that template fingerprints; the tests at the
foot of this file exist to make that copy fail loudly.
"""
from __future__ import annotations

import sqlite3

import pytest

from bib.authority import create_schema
from bib.export import export_bibtex
from bib.parser import bib_entry_to_metadata
from pipeline.stages import _expected_fingerprints_for_run


ENTRY = {
    "_key": "Stepanjants1970",
    "_type": "article",
    "title": "{Siphonophores of the Kurile-Kamchatka Trench}",
    "author": "Stepanjants, S. D.",
    "year": "1970",
    "file": "Stepanjants1970.pdf",
    "keeppages": "3--20",
    "doclang": "ru",
    "pagemap": "1--2 BHL wrapper; 3--20 Russian original; "
               "21--40 English translation; 41--43 blank",
    "ocrlang": "rus+eng",
}


def test_the_metadata_dict_carries_both_fields():
    meta = bib_entry_to_metadata(ENTRY, "Stepanjants1970.pdf")
    assert meta["doclang"] == "ru"
    assert meta["pagemap"].startswith("1--2 BHL wrapper")
    # And the instruction it sits beside is untouched.
    assert meta["ocrlang"] == "rus+eng"


def test_absent_fields_are_none_not_empty_string():
    """`None` so a NULL column round-trips as an omitted bib field."""
    meta = bib_entry_to_metadata({"_key": "X"}, "X.pdf")
    assert meta["doclang"] is None
    assert meta["pagemap"] is None


@pytest.mark.parametrize("tag", ["de-Latf", "zh-Hant", "grc", "sr-Latn"])
def test_doclang_holds_tags_no_iso_639_code_can_express(tag):
    """The cases that decide the vocabulary.

    `de-Latf` is the one that settles it: without `deu_latf` a scanned
    19th-c. German paper OCRs to whitespace, and langdetect can only ever
    emit bare `de`. "German, and set in Fraktur" is a judgment a visual pass
    can make and detection cannot express at all. The field stores the tag
    verbatim — resolving it to packs is `bcp47_to_tesseract`'s job (#215),
    and happens before a run, not during one.
    """
    meta = bib_entry_to_metadata({"_key": "X", "doclang": tag}, "X.pdf")
    assert meta["doclang"] == tag


def test_pagemap_is_free_text_and_is_not_parsed():
    """No schema, no validation. It is a note to the next reader."""
    scrawl = "front matter?? maybe 1-4; plates unpaginated at end"
    meta = bib_entry_to_metadata({"_key": "X", "pagemap": scrawl}, "X.pdf")
    assert meta["pagemap"] == scrawl


# ---------------------------------------------------------------------------
# Round-trip through the authority DB
# ---------------------------------------------------------------------------

def _seeded_db(path, **cols):
    conn = sqlite3.connect(path)
    create_schema(conn)
    conn.execute(
        "INSERT INTO works (work_id, guid_type, title, year, in_corpus, "
        "source, created_at, updated_at) "
        "VALUES ('w1', 'corpus_key', 'A paper', 1970, 1, 'corpus_paper', 0, 0)")
    for name, value in cols.items():
        conn.execute(f"UPDATE works SET {name} = ? WHERE work_id = 'w1'", (value,))
    conn.commit()
    conn.close()
    return path


def test_export_round_trips_both_fields(tmp_path):
    db = _seeded_db(tmp_path / "a.sqlite", doclang="de-Latf",
                    pagemap="1--2 JSTOR cover; 3-- article")
    text = export_bibtex(db)
    assert "doclang = {de-Latf}" in text
    assert "1--2 JSTOR cover" in text


def test_migration_re_adds_the_columns_on_an_older_db(tmp_path):
    """An existing biblio_authority.sqlite predates the columns.

    `CREATE TABLE IF NOT EXISTS` is a no-op on it, so the ALTER TABLE pass
    is the only thing that adds them — and it has two halves, the column
    list and the tuple it is unpacked into. Missing either breaks new DBs or
    old ones respectively, and only one of those shows up in a fresh test.
    """
    db = tmp_path / "old.sqlite"
    conn = sqlite3.connect(db)
    create_schema(conn)
    conn.execute("ALTER TABLE works DROP COLUMN doclang")
    conn.execute("ALTER TABLE works DROP COLUMN pagemap")
    conn.commit()
    have = {row[1] for row in conn.execute("PRAGMA table_info(works)")}
    assert not {"doclang", "pagemap"} & have, "precondition: columns removed"

    create_schema(conn)                  # re-running must migrate them back

    have = {row[1] for row in conn.execute("PRAGMA table_info(works)")}
    conn.close()
    assert {"doclang", "pagemap"} <= have


# ---------------------------------------------------------------------------
# The property that makes them usable
# ---------------------------------------------------------------------------

def test_no_stage_fingerprint_mentions_the_curation_fields():
    """Editing a note must not reprocess a document.

    `ocrlang` rewrites `processed.pdf` and so fingerprints every
    OCR-dependent stage. These two steer nothing, so if they ever reach a
    fingerprint the only thing that can follow is spurious re-runs — and,
    worse, the appearance that a correction did something.
    """
    fps = _expected_fingerprints_for_run(
        ocrlang="rus+eng",
        keeppages=None,
        taxonomy_fingerprint={"sha256": "abc"},
        lexicon_fingerprints={"anatomy": {"sha256": "def"}},
    )
    flat = repr(fps)
    assert "doclang" not in flat
    assert "pagemap" not in flat
    for stage, fp in fps.items():
        assert "doclang" not in fp, stage
        assert "pagemap" not in fp, stage


def test_the_fingerprint_builder_takes_no_curation_argument():
    """The structural guarantee behind the test above.

    Asserted directly because a future caller could otherwise pass one in
    and the value-level check would still pass on today's call sites.
    """
    import inspect
    params = inspect.signature(_expected_fingerprints_for_run).parameters
    assert "doclang" not in params
    assert "pagemap" not in params


def test_no_accessor_exists_for_either_field():
    """`ocrlang` has `entry_ocrlang()` and `ocrlang_for_pdf()` because the
    scan stage reads it. Nothing reads these, so nothing should be able to.
    """
    import bib.parser as parser
    for name in ("entry_doclang", "doclang_for_pdf",
                 "entry_pagemap", "pagemap_for_pdf"):
        assert not hasattr(parser, name), f"{name} would give the field an effect"
