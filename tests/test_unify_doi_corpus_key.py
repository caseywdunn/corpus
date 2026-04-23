"""Unit tests for unify_doi_corpus_key.py."""

from __future__ import annotations

import importlib.util
import sqlite3
import sys
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_PATH = REPO_ROOT / "unify_doi_corpus_key.py"

_spec = importlib.util.spec_from_file_location("unify_doi_corpus_key", SCRIPT_PATH)
unify_mod = importlib.util.module_from_spec(_spec)
sys.modules["unify_doi_corpus_key"] = unify_mod
_spec.loader.exec_module(unify_mod)


def _row(**kwargs):
    """Build a dict matching the sqlite3.Row surface the helpers use."""
    defaults = {
        "work_id": "x", "title": "t", "year": 2005, "in_corpus": 0,
        "doi": None, "surname": "dunn", "cites": 0, "corpus_hash": None,
    }
    defaults.update(kwargs)
    return defaults


# ── score_title ─────────────────────────────────────────────────────

def test_score_title_identical_titles():
    set_s, ratio_s = unify_mod.score_title(
        "Molecular Phylogenetics of the Siphonophora",
        "Molecular phylogenetics of the siphonophora",
    )
    assert set_s == 100 and ratio_s >= 95


def test_score_title_empty_inputs():
    assert unify_mod.score_title("", "anything") == (0, 0)
    assert unify_mod.score_title("anything", "") == (0, 0)


# ── pick_match ──────────────────────────────────────────────────────

def test_pick_match_unambiguous():
    doi_row = _row(title="Molecular phylogenetics of the Siphonophora")
    candidates = [
        _row(work_id="corpus:dunn|2005|molecular phylogenetics of the siphonoph",
             title="Molecular Phylogenetics of the Siphonophora (Cnidaria)", in_corpus=1),
    ]
    match, status, _ = unify_mod.pick_match(
        doi_row, candidates, min_set_score=85, min_partial_ratio=75,
    )
    assert status == "matched"
    assert match["work_id"].startswith("corpus:dunn|2005|")


def test_pick_match_no_candidates():
    doi_row = _row(title="Something about jellyfish")
    match, status, _ = unify_mod.pick_match(
        doi_row, [], min_set_score=85, min_partial_ratio=75,
    )
    assert status == "no_match"
    assert match is None


def test_pick_match_no_qualifying_candidate():
    # Candidate shares surname+year but a different title
    doi_row = _row(title="Molecular phylogenetics of the Siphonophora")
    candidates = [
        _row(work_id="corpus:dunn|2005|somethingelse",
             title="Chapter 5: Marrus claudanielis, a new species"),
    ]
    _, status, _ = unify_mod.pick_match(
        doi_row, candidates, min_set_score=85, min_partial_ratio=75,
    )
    assert status == "no_match"


def test_pick_match_ambiguous():
    # Two candidates both qualify AND neither is in_corpus=1 — ambiguous, skip
    doi_row = _row(title="Molecular phylogenetics of the Siphonophora")
    candidates = [
        _row(work_id="corpus:dunn|2005|a",
             title="Molecular Phylogenetics of the Siphonophora v1"),
        _row(work_id="corpus:dunn|2005|b",
             title="Molecular Phylogenetics of the Siphonophora v2"),
    ]
    _, status, _ = unify_mod.pick_match(
        doi_row, candidates, min_set_score=85, min_partial_ratio=75,
    )
    assert status == "ambiguous"


def test_pick_match_in_corpus_wins_tie():
    """When multiple candidates qualify but exactly one is in_corpus=1,
    the in-corpus row wins — this is the Dunn 2005 case where a
    canonical corpus paper shares (surname, year) with OCR-variant
    ghosts that ghost-dedup left just below its merge threshold.
    """
    doi_row = _row(title="Molecular phylogenetics of the Siphonophora")
    candidates = [
        _row(work_id="corpus:dunn|2005|canonical",
             title="Molecular Phylogenetics of the Siphonophora", in_corpus=1),
        _row(work_id="corpus:dunn|2005|ocr_variant",
             title="Molecular phylogenetics ofthe Siphonophora"),
    ]
    match, status, _ = unify_mod.pick_match(
        doi_row, candidates, min_set_score=85, min_partial_ratio=75,
    )
    assert status == "matched"
    assert match["work_id"] == "corpus:dunn|2005|canonical"


def test_pick_match_two_in_corpus_still_ambiguous():
    """If two qualifying candidates are both in_corpus=1, there's no
    tie-break — stay ambiguous (very unlikely but don't silently pick)."""
    doi_row = _row(title="Molecular phylogenetics of the Siphonophora")
    candidates = [
        _row(work_id="corpus:dunn|2005|a",
             title="Molecular Phylogenetics of the Siphonophora", in_corpus=1),
        _row(work_id="corpus:dunn|2005|b",
             title="Molecular Phylogenetics of the Siphonophora", in_corpus=1),
    ]
    _, status, _ = unify_mod.pick_match(
        doi_row, candidates, min_set_score=85, min_partial_ratio=75,
    )
    assert status == "ambiguous"


# ── pick_survivor ───────────────────────────────────────────────────

def test_survivor_corpus_paper_wins_over_doi_ghost():
    doi_row = _row(work_id="10.1080/x", in_corpus=0)
    corpus_row = _row(work_id="corpus:dunn|2005|...", in_corpus=1)
    survivor, dup, reason = unify_mod.pick_survivor(doi_row, corpus_row)
    assert survivor == "corpus:dunn|2005|..."
    assert dup == "10.1080/x"
    assert reason == "corpus_paper_wins"


def test_survivor_doi_wins_when_both_ghosts():
    doi_row = _row(work_id="10.1080/x", in_corpus=0)
    corpus_row = _row(work_id="corpus:dunn|2005|...", in_corpus=0)
    survivor, dup, reason = unify_mod.pick_survivor(doi_row, corpus_row)
    assert survivor == "10.1080/x"
    assert dup == "corpus:dunn|2005|..."
    assert reason == "doi_wins_ghost_ghost"


def test_survivor_raises_when_both_in_corpus():
    doi_row = _row(work_id="10.1080/x", in_corpus=1)
    corpus_row = _row(work_id="corpus:dunn|2005|...", in_corpus=1)
    with pytest.raises(ValueError):
        unify_mod.pick_survivor(doi_row, corpus_row)


# ── merge_duplicate — DB integration ────────────────────────────────

def _make_db():
    conn = sqlite3.connect(":memory:")
    conn.row_factory = sqlite3.Row
    conn.executescript("""
        CREATE TABLE works (
            work_id TEXT PRIMARY KEY, guid_type TEXT NOT NULL,
            title TEXT, year INTEGER, journal TEXT, doi TEXT,
            bhl_item_id TEXT, bhl_part_id TEXT, openalex_id TEXT,
            corpus_hash TEXT, in_corpus INTEGER NOT NULL DEFAULT 0,
            source TEXT NOT NULL, confidence REAL DEFAULT 1.0,
            created_at REAL NOT NULL, updated_at REAL NOT NULL
        );
        CREATE TABLE work_authors (
            work_id TEXT NOT NULL, position INTEGER NOT NULL,
            surname TEXT NOT NULL, surname_normalized TEXT NOT NULL,
            forename TEXT,
            PRIMARY KEY (work_id, position)
        );
        CREATE TABLE citations (
            citing_work_id TEXT NOT NULL, cited_work_id TEXT NOT NULL,
            citing_corpus_hash TEXT NOT NULL,
            grobid_xml_id TEXT, raw_citation TEXT,
            match_method TEXT NOT NULL, match_score REAL DEFAULT 1.0,
            PRIMARY KEY (citing_work_id, cited_work_id, citing_corpus_hash)
        );
        CREATE TABLE work_aliases (
            alias_key TEXT NOT NULL, work_id TEXT NOT NULL,
            PRIMARY KEY (alias_key, work_id)
        );
    """)
    return conn


def _ins_work(conn, work_id, guid_type, title, in_corpus=0, doi=None, corpus_hash=None):
    now = time.time()
    conn.execute(
        """INSERT INTO works (work_id, guid_type, title, year, journal, doi,
           corpus_hash, in_corpus, source, confidence, created_at, updated_at)
           VALUES (?, ?, ?, 2005, '', ?, ?, ?,
                   'cited_reference', 1.0, ?, ?)""",
        (work_id, guid_type, title, doi, corpus_hash, in_corpus, now, now),
    )


def _ins_author(conn, work_id, position, surname):
    conn.execute(
        """INSERT INTO work_authors (work_id, position, surname, surname_normalized, forename)
           VALUES (?, ?, ?, ?, '')""",
        (work_id, position, surname, unify_mod.normalize_for_key(surname)),
    )


def test_merge_corpus_row_survives_doi_goes_into_doi_column():
    """Dunn 2005 Molecular Phylogenetics: in-corpus row survives,
    absorbs the DOI ghost's citations, and the DOI string ends up in
    its ``doi`` column."""
    conn = _make_db()
    _ins_work(conn, "corpus:dunn|2005|mp", "corpus_key",
              "Molecular Phylogenetics of the Siphonophora",
              in_corpus=1, corpus_hash="01fd6fbf6f76")
    _ins_author(conn, "corpus:dunn|2005|mp", 0, "Dunn")
    _ins_work(conn, "10.1080/10635150500354837", "doi",
              "Molecular phylogenetics of the siphonophora",
              doi="10.1080/10635150500354837")
    _ins_author(conn, "10.1080/10635150500354837", 0, "Dunn")
    # Cites split between the two rows
    for i in range(3):
        _ins_work(conn, f"p{i}", "corpus_key", "", in_corpus=1)
        conn.execute(
            f"INSERT INTO citations VALUES "
            f"('p{i}','corpus:dunn|2005|mp','h{i}','b{i}','','x',1.0)"
        )
    for i in range(3, 10):
        _ins_work(conn, f"p{i}", "corpus_key", "", in_corpus=1)
        conn.execute(
            f"INSERT INTO citations VALUES "
            f"('p{i}','10.1080/10635150500354837','h{i}','b{i}','','x',1.0)"
        )
    conn.commit()

    unify_mod.unify(conn)

    # Corpus row survives and now has the DOI
    row = conn.execute(
        "SELECT doi, corpus_hash, in_corpus FROM works WHERE work_id='corpus:dunn|2005|mp'"
    ).fetchone()
    assert row["doi"] == "10.1080/10635150500354837"
    assert row["corpus_hash"] == "01fd6fbf6f76"
    assert row["in_corpus"] == 1

    # DOI row gone
    assert conn.execute(
        "SELECT 1 FROM works WHERE work_id='10.1080/10635150500354837'"
    ).fetchone() is None

    # All 10 cites now on the corpus row
    n = conn.execute(
        "SELECT COUNT(*) FROM citations WHERE cited_work_id='corpus:dunn|2005|mp'"
    ).fetchone()[0]
    assert n == 10


def test_merge_two_ghosts_doi_row_wins():
    """When neither side is in_corpus=1, the DOI row survives
    (globally unique identifier) and absorbs the corpus_key ghost."""
    conn = _make_db()
    _ins_work(conn, "corpus:foo|2005|some title", "corpus_key",
              "Some title about widgets", in_corpus=0)
    _ins_author(conn, "corpus:foo|2005|some title", 0, "Foo")
    _ins_work(conn, "10.5555/example", "doi",
              "Some title about widgets", in_corpus=0,
              doi="10.5555/example")
    _ins_author(conn, "10.5555/example", 0, "Foo")
    conn.commit()

    unify_mod.unify(conn)

    # DOI row survives
    assert conn.execute(
        "SELECT 1 FROM works WHERE work_id='10.5555/example'"
    ).fetchone() is not None
    # Ghost corpus row gone
    assert conn.execute(
        "SELECT 1 FROM works WHERE work_id='corpus:foo|2005|some title'"
    ).fetchone() is None


def test_unify_skips_when_no_corpus_key_sibling():
    """DOI ghost with no matching corpus_key row survives untouched
    — those are real external references we just don't have a PDF
    for."""
    conn = _make_db()
    _ins_work(conn, "10.5555/lonely", "doi", "Paper about jellyfish",
              in_corpus=0, doi="10.5555/lonely")
    _ins_author(conn, "10.5555/lonely", 0, "Lonely")
    conn.commit()

    unify_mod.unify(conn)

    assert conn.execute(
        "SELECT 1 FROM works WHERE work_id='10.5555/lonely'"
    ).fetchone() is not None


def test_unify_skips_ambiguous_two_matching_candidates():
    """Two corpus_key rows both fuzzy-match the DOI — skip, don't
    guess.  Regression guard for pathological Grobid output."""
    conn = _make_db()
    _ins_work(conn, "corpus:yyz|2005|widget paper v1", "corpus_key",
              "Widget paper about jellyfish morphology", in_corpus=0)
    _ins_author(conn, "corpus:yyz|2005|widget paper v1", 0, "Yyz")
    _ins_work(conn, "corpus:yyz|2005|widget paper v2", "corpus_key",
              "Widget paper about jellyfish morphology", in_corpus=0)
    _ins_author(conn, "corpus:yyz|2005|widget paper v2", 0, "Yyz")
    _ins_work(conn, "10.5555/widget", "doi",
              "Widget paper about jellyfish morphology",
              in_corpus=0, doi="10.5555/widget")
    _ins_author(conn, "10.5555/widget", 0, "Yyz")
    conn.commit()

    counts = unify_mod.unify(conn)

    assert counts["ambiguous"] == 1
    assert counts["merged"] == 0
    # All three survive
    for wid in ("corpus:yyz|2005|widget paper v1",
                "corpus:yyz|2005|widget paper v2",
                "10.5555/widget"):
        assert conn.execute(
            "SELECT 1 FROM works WHERE work_id=?", (wid,)
        ).fetchone() is not None
