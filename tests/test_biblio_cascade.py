"""Unit tests for build_biblio_authority.py reference-resolution cascade.

Covers the behavior introduced to fix GH issue #2: the cascade must not
route a titled reference onto the only (surname, year) candidate when
the titles are dissimilar, and must not cache low-confidence matches
via alias (one bad decision must not amplify into many).
"""

from __future__ import annotations

import importlib.util
import sqlite3
import sys
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_PATH = REPO_ROOT / "build_biblio_authority.py"

_spec = importlib.util.spec_from_file_location("build_biblio_authority", SCRIPT_PATH)
biblio = importlib.util.module_from_spec(_spec)
sys.modules["build_biblio_authority"] = biblio
_spec.loader.exec_module(biblio)


@pytest.fixture
def conn():
    conn = sqlite3.connect(":memory:")
    biblio.create_schema(conn)
    yield conn
    conn.close()


def _seed_work(conn, work_id, title, year, surname, position=0):
    now = time.time()
    conn.execute(
        """INSERT OR IGNORE INTO works (work_id, guid_type, title, year, journal, doi,
           corpus_hash, in_corpus, source, confidence, created_at, updated_at)
           VALUES (?, 'corpus_key', ?, ?, '', NULL, NULL, 0, 'cited_reference', 1.0, ?, ?)""",
        (work_id, title, year, now, now),
    )
    biblio.insert_authors(conn, work_id, [(surname, "")])
    # Also register the normal alias_key so alias_exact lookups can hit
    biblio.insert_alias(conn, biblio.make_alias_key(surname, year, title), work_id)
    conn.commit()


# ── Regression: Lensia / Synopsis misrouting (issue #2) ─────────────

def test_different_title_does_not_route_to_unique_candidate(conn):
    """The failure that motivated issue #2.

    Given a single (Totton, 1965) ghost for 'A synopsis of the Siphonophora',
    a reference with title 'A new species of Lensia ...' must not route
    onto the Synopsis row. ``token_set_ratio`` alone scores the pair at
    82 (misleadingly high due to shared "a of the siphonophora"
    tokens), but ``ratio`` scores it at 43 — the dual-metric guard
    rejects this pair.
    """
    _seed_work(conn, "corpus:totton|1965|a synopsis of the siphonophora",
               "A synopsis of the Siphonophora", 1965, "Totton")

    ref = {
        "title": "A new species of Lensia (Siphonophora: Diphyidae) from the coastal waters of Vancouver",
        "year": 1965,
        "authors": ["A Totton"],
        "raw": "Totton, A. K. 1965. A new species of Lensia ...",
    }
    work_id, method, score = biblio._resolve_reference(conn, ref)

    assert work_id != "corpus:totton|1965|a synopsis of the siphonophora"
    assert method == "new"


def test_rejected_match_does_not_cache_alias(conn):
    """The failure mode from issue #2 amplifies via alias caching: one
    misroute gets cached, then every subsequent ref with the same title
    hits alias_exact. After the fix, a rejected fuzzy match must leave
    the cached alias pointing at the newly-created ghost, NOT at the
    seeded candidate.
    """
    _seed_work(conn, "corpus:totton|1965|a synopsis of the siphonophora",
               "A synopsis of the Siphonophora", 1965, "Totton")

    lensia_title = "A new species of Lensia (Siphonophora: Diphyidae) from Vancouver"
    ref1 = {"title": lensia_title, "year": 1965, "authors": ["A Totton"], "raw": ""}
    ref2 = {"title": lensia_title, "year": 1965, "authors": ["A Totton"], "raw": ""}

    wid1, _, _ = biblio._resolve_reference(conn, ref1)
    wid2, method2, _ = biblio._resolve_reference(conn, ref2)
    conn.commit()

    # Both refs must resolve to the SAME new ghost (not the Synopsis row)
    assert wid1 == wid2
    assert wid1 != "corpus:totton|1965|a synopsis of the siphonophora"
    # Second ref uses alias_exact (cached during ref1's new-work creation),
    # proving the legitimate caching path still works for new ghosts.
    assert method2 == "alias_exact"


def test_no_title_still_uses_author_year_fallback(conn):
    """When the ref has no title at all we accept the unique (surname,
    year) candidate (author_year_only). The cascade's title guard only
    kicks in when a title is present.
    """
    _seed_work(conn, "corpus:only|1900|the one and only work by only",
               "The one and only work by Only", 1900, "Only")

    ref = {
        "title": "",
        "year": 1900,
        "authors": ["A Only"],
        "raw": "Only, A. 1900.",
    }
    work_id, method, _ = biblio._resolve_reference(conn, ref)
    assert work_id == "corpus:only|1900|the one and only work by only"
    assert method == "author_year_only"


def test_no_title_author_year_match_does_not_cache_alias(conn):
    """author_year_only matches must never cache an alias. A subsequent
    ref for a DIFFERENT paper by the same (surname, year) — once a
    second ghost exists — must be able to get its own ghost, not
    collide through a cached empty-title alias.
    """
    _seed_work(conn, "corpus:x|1950|first paper by x",
               "First paper by X", 1950, "X")

    ref = {"title": "", "year": 1950, "authors": ["Q X"], "raw": ""}
    biblio._resolve_reference(conn, ref)
    conn.commit()

    # No alias should be registered for "x|1950|" (empty-title key)
    empty_alias = biblio.make_alias_key("X", 1950, "")
    rows = conn.execute(
        "SELECT work_id FROM work_aliases WHERE alias_key = ?", (empty_alias,)
    ).fetchall()
    # The seed set its own alias for "first paper by x", not "", so
    # there should be no alias matching the empty-title key from the
    # author_year_only match.
    assert all(r[0] != "corpus:x|1950|first paper by x" or
               biblio.make_alias_key("X", 1950, "First paper by X") == empty_alias
               for r in rows)


def test_high_confidence_match_does_cache_alias(conn):
    """≥80 fuzzy match is the sweet spot where caching is safe: very
    likely the same work, different OCR. Subsequent refs with the same
    (now-cached) title should alias_exact to the same work.
    """
    _seed_work(conn, "corpus:author|2001|foo bar baz quux title",
               "Foo Bar Baz Quux Title", 2001, "Author")

    ref = {
        "title": "Foo Bar Baz Quux Title",   # exact match, should score 100
        "year": 2001,
        "authors": ["Author"],
        "raw": "",
    }
    work_id, method, _ = biblio._resolve_reference(conn, ref)
    assert work_id == "corpus:author|2001|foo bar baz quux title"
    assert method in {"title_fuzzy", "alias_exact"}

    # Alias for this title form must now resolve to the seed row
    alias = biblio.make_alias_key("Author", 2001, "Foo Bar Baz Quux Title")
    row = conn.execute(
        "SELECT work_id FROM work_aliases WHERE alias_key = ?", (alias,)
    ).fetchone()
    assert row is not None
    assert row[0] == "corpus:author|2001|foo bar baz quux title"


# ── Regression guard: first-author position ─────────────────────────

def test_match_restricted_to_first_author(conn):
    """A candidate where the target surname sits at a non-first author
    position must NOT be considered a match for author_year_match.
    Multi-author ghosts get resurrected wrongly if we allow any-position
    matching.
    """
    # Seed a work where Totton is the 2nd author
    now = time.time()
    conn.execute(
        """INSERT INTO works (work_id, guid_type, title, year, journal, doi,
           corpus_hash, in_corpus, source, confidence, created_at, updated_at)
           VALUES ('corpus:other|1965|coauthored paper', 'corpus_key', 'Coauthored paper', 1965, '', NULL, NULL, 0, 'cited_reference', 1.0, ?, ?)""",
        (now, now),
    )
    biblio.insert_authors(conn, "corpus:other|1965|coauthored paper",
                          [("Other", ""), ("Totton", "")])
    conn.commit()

    ref = {"title": "", "year": 1965, "authors": ["A Totton"], "raw": ""}
    result = biblio.author_year_match(conn, "Totton", 1965)
    assert result is None  # no first-author (Totton, 1965) row exists
