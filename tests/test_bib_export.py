"""Tests for bib_export (#26 export half)."""
from __future__ import annotations

import sqlite3
import time
from pathlib import Path

import pytest

from bib.export import (
    _alpha_suffix,
    _ascii_surname,
    assign_cite_keys,
    export_bibtex,
    render_entry,
)


# ---------------------------------------------------------------------------
# Cite-key helpers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("surname, expected", [
    ("Smith", "Smith"),
    ("Pugé", "Puge"),
    ("Müller", "Muller"),
    ("O'Brien", "OBrien"),
    ("van der Zee", "vanderZee"),
    ("", ""),
])
def test_ascii_surname_strips_diacritics_and_punct(surname, expected):
    assert _ascii_surname(surname) == expected


@pytest.mark.parametrize("n, expected", [
    (0, ""),
    (1, "a"),
    (2, "b"),
    (26, "z"),
    (27, "aa"),
    (28, "ab"),
    (52, "az"),
    (53, "ba"),
])
def test_alpha_suffix(n, expected):
    assert _alpha_suffix(n) == expected


def test_assign_cite_keys_unique_pairs_get_bare_keys():
    keys = assign_cite_keys([
        ("w1", 2010, "Smith"),
        ("w2", 1995, "Jones"),
        ("w3", 2020, "Lee"),
    ])
    assert keys == {"w1": "Smith2010", "w2": "Jones1995", "w3": "Lee2020"}


def test_assign_cite_keys_disambiguates_collisions():
    keys = assign_cite_keys([
        ("w1", 2010, "Smith"),
        ("w2", 2010, "Smith"),
        ("w3", 2010, "Smith"),
    ])
    assert keys["w1"] == "Smith2010"
    assert keys["w2"] == "Smith2010a"
    assert keys["w3"] == "Smith2010b"


def test_assign_cite_keys_handles_missing_year_and_surname():
    keys = assign_cite_keys([
        ("w1", None, ""),
        ("w2", None, ""),
    ])
    assert keys["w1"] == "Anonnodate"
    assert keys["w2"] == "Anonnodatea"


# ---------------------------------------------------------------------------
# render_entry
# ---------------------------------------------------------------------------


def test_render_entry_omits_empty_fields():
    text = render_entry("Smith2010", "article", {
        "author": "Smith, J",
        "title": "On siphonophores",
        "year": "2010",
        "journal": None,
        "doi": "",
    })
    assert "author = {Smith, J}" in text
    assert "title = {On siphonophores}" in text
    assert "year = {2010}" in text
    assert "journal" not in text
    assert "doi" not in text


def test_render_entry_escapes_braces():
    text = render_entry("X", "misc", {"title": "{boom}"})
    assert r"\{boom\}" in text


# ---------------------------------------------------------------------------
# export_bibtex (integration)
# ---------------------------------------------------------------------------


def _make_min_db(path: Path) -> None:
    """Minimal works/work_authors fixture matching biblio_authority.sqlite shape."""
    conn = sqlite3.connect(path)
    conn.executescript("""
        CREATE TABLE works (
            work_id TEXT PRIMARY KEY, guid_type TEXT, title TEXT, year INTEGER,
            journal TEXT, doi TEXT, bhl_item_id TEXT, bhl_part_id TEXT,
            openalex_id TEXT, corpus_hash TEXT, in_corpus INTEGER NOT NULL DEFAULT 0,
            source TEXT, confidence REAL,
            license TEXT, license_url TEXT, license_source TEXT,
            publishable INTEGER,
            serve INTEGER NOT NULL DEFAULT 1, serve_reason TEXT,
            created_at REAL, updated_at REAL
        );
        CREATE TABLE work_authors (
            work_id TEXT, position INTEGER, surname TEXT, surname_normalized TEXT,
            forename TEXT, PRIMARY KEY (work_id, position)
        );
    """)
    now = time.time()
    conn.executemany(
        """INSERT INTO works
           (work_id, guid_type, title, year, journal, doi, corpus_hash,
            in_corpus, source, confidence, created_at, updated_at)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        [
            ("10.1/aaa", "doi", "Earlier work", 1995, "Test J", "10.1/aaa",
             "aaa", 1, "corpus_paper", 1.0, now, now),
            ("10.1/bbb", "doi", "Later work", 2010, "Test J", "10.1/bbb",
             "bbb", 1, "corpus_paper", 1.0, now, now),
            # cited-only ghost — should NOT appear with corpus_only=True
            ("corpus:smith|1990|extra", "corpus_key", "Ghost", 1990,
             None, None, None, 0, "cited_reference", 1.0, now, now),
        ],
    )
    conn.executemany(
        "INSERT INTO work_authors (work_id, position, surname, surname_normalized, forename) "
        "VALUES (?, ?, ?, ?, ?)",
        [
            ("10.1/aaa", 0, "Smith", "smith", "Jane"),
            ("10.1/bbb", 0, "Jones", "jones", "K"),
            ("corpus:smith|1990|extra", 0, "Smith", "smith", "X"),
        ],
    )
    conn.commit()
    conn.close()


def test_export_corpus_only_skips_ghost(tmp_path):
    db = tmp_path / "biblio.sqlite"
    _make_min_db(db)
    bib = export_bibtex(db, corpus_only=True)
    assert "Earlier work" in bib
    assert "Later work" in bib
    assert "Ghost" not in bib


def test_export_all_works_includes_ghost(tmp_path):
    db = tmp_path / "biblio.sqlite"
    _make_min_db(db)
    bib = export_bibtex(db, corpus_only=False)
    assert "Ghost" in bib


def test_export_emits_file_field_when_documents_dir_supplied(tmp_path):
    db = tmp_path / "biblio.sqlite"
    _make_min_db(db)
    bib = export_bibtex(db, documents_dir=tmp_path / "documents")
    assert "file = {" in bib
    assert "aaa/processed.pdf" in bib


def test_export_omits_file_field_when_documents_dir_none(tmp_path):
    db = tmp_path / "biblio.sqlite"
    _make_min_db(db)
    bib = export_bibtex(db, documents_dir=None)
    assert "file = {" not in bib


def test_export_round_trip_safe_against_re_runs(tmp_path):
    """Re-export of the same DB produces byte-identical output."""
    db = tmp_path / "biblio.sqlite"
    _make_min_db(db)
    a = export_bibtex(db)
    b = export_bibtex(db)
    assert a == b


def test_export_cite_keys_use_first_author_and_year(tmp_path):
    db = tmp_path / "biblio.sqlite"
    _make_min_db(db)
    bib = export_bibtex(db)
    # Order is "year, work_id" → 1995 first, then 2010
    assert "@article{Smith1995" in bib
    assert "@article{Jones2010" in bib
