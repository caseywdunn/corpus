"""Diacritic-folding author lookup (#122).

`work_authors.surname_normalized` (and the in-memory `author_to_papers`
index) are built with `bib.authority.normalize_for_key`, which strips
diacritics. The serve-side author tools used to query with a plain
`.strip().lower()`, so `get_works_by_author("Müller")` returned 0 while
`"Muller"` returned the real set. These tests pin that the query side
now uses the same normalizer:

- `Müller` ≡ `Muller` (both → "muller")
- `Mueller` stays a distinct surname (normalize_for_key is diacritic-
  strip, NOT oe-expansion)
- ASCII-only queries are unchanged (regression guard)
"""
from __future__ import annotations

import sqlite3
import time
import types
from pathlib import Path

import pytest

from bib.authority import create_schema, normalize_for_key
from mcpsrv import app as mcp_app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import get_works_by_author
from mcpsrv.tools.taxonomy import get_papers_by_author


def _insert_work(conn, work_id, surname, *, forename=None):
    now = time.time()
    conn.execute(
        """INSERT INTO works
           (work_id, guid_type, title, year, journal, doi, in_corpus,
            corpus_hash, source, confidence, bib_imported_at,
            created_at, updated_at)
           VALUES (?, 'corpus_key', ?, 2010, 'J', NULL, 0, NULL,
                   'cited_reference', 1.0, NULL, ?, ?)""",
        (work_id, f"work {work_id}", now, now),
    )
    # surname_normalized is populated exactly as the importer does it.
    conn.execute(
        """INSERT INTO work_authors
           (work_id, position, surname, surname_normalized, forename)
           VALUES (?, 0, ?, ?, ?)""",
        (work_id, surname, normalize_for_key(surname), forename),
    )


@pytest.fixture
def index(tmp_path: Path):
    db_path = tmp_path / "biblio.sqlite"
    conn = sqlite3.connect(db_path)
    create_schema(conn)
    # Two diacritic surnames, one oe-spelled, one ASCII control.
    _insert_work(conn, "w:muller1", "Müller")
    _insert_work(conn, "w:muller2", "Müller")
    _insert_work(conn, "w:mueller1", "Mueller")   # legitimately distinct
    _insert_work(conn, "w:smith1", "Smith")
    conn.commit()
    conn.close()

    fake = types.SimpleNamespace(biblio_db=BiblioAuthority(db_path))
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    yield fake
    mcp_app.set_index(original)


def test_diacritic_and_ascii_spelling_match(index):
    umlaut = {w["work_id"] for w in get_works_by_author("Müller")}
    ascii_ = {w["work_id"] for w in get_works_by_author("Muller")}
    assert umlaut == ascii_ == {"w:muller1", "w:muller2"}


def test_oe_spelling_stays_distinct(index):
    """normalize_for_key strips diacritics; it does NOT expand ü→ue, so
    'Mueller' must not collapse into the 'Müller' set."""
    mueller = {w["work_id"] for w in get_works_by_author("Mueller")}
    assert mueller == {"w:mueller1"}
    assert "w:muller1" not in mueller


def test_ascii_surname_unchanged(index):
    smith = {w["work_id"] for w in get_works_by_author("Smith")}
    assert smith == {"w:smith1"}


def test_get_papers_by_author_folds_diacritics():
    """The Grobid-metadata path (author_to_papers) is keyed with the same
    normalizer, so a diacritic query and its ASCII spelling hit the same
    key."""
    fake = types.SimpleNamespace(
        author_to_papers={normalize_for_key("Müller"): ["hash00000001"]},
        papers={"hash00000001": {"title": "T", "year": 2010, "authors": []}},
    )
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    try:
        by_umlaut = {p["hash"] for p in get_papers_by_author("Müller")}
        by_ascii = {p["hash"] for p in get_papers_by_author("Muller")}
        assert by_umlaut == by_ascii == {"hash00000001"}
    finally:
        mcp_app.set_index(original)
