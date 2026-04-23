"""Unit tests for reconcile_corpus_to_biblio.py."""

from __future__ import annotations

import importlib.util
import sqlite3
import sys
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_PATH = REPO_ROOT / "reconcile_corpus_to_biblio.py"

# Load the reconcile module by path — the repo isn't structured as a package.
_spec = importlib.util.spec_from_file_location("reconcile_corpus_to_biblio", SCRIPT_PATH)
reconcile = importlib.util.module_from_spec(_spec)
sys.modules["reconcile_corpus_to_biblio"] = reconcile
_spec.loader.exec_module(reconcile)


# ── derive_from_filename ────────────────────────────────────────────

@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        ("Totton1965a.pdf", ("Totton", 1965)),
        ("Totton1954.pdf", ("Totton", 1954)),
        ("Totton_1954.pdf", ("Totton", 1954)),
        ("Dunn_Pugh_Haddock2005.pdf", ("Dunn", 2005)),
        ("Dunn-etal2005_Marrus.pdf", ("Dunn", 2005)),
        ("Dunn2005b_Bargmannia.pdf", ("Dunn", 2005)),
        ("Haeckel1888.pdf", ("Haeckel", 1888)),
        ("Munro_etal2019.pdf", ("Munro", 2019)),
    ],
)
def test_derive_from_filename_happy_paths(filename, expected):
    assert reconcile.derive_from_filename(filename) == expected


def test_derive_from_filename_empty():
    assert reconcile.derive_from_filename("") == ("", None)


def test_derive_from_filename_no_year():
    # Surname but no 4-digit year anywhere
    assert reconcile.derive_from_filename("Totton.pdf") == ("", None)


def test_derive_from_filename_implausible_year():
    # 4 digits but not a plausible publication year
    surname, year = reconcile.derive_from_filename("Paper0042.pdf")
    assert surname == "Paper"
    assert year is None


def test_derive_from_filename_no_leading_surname():
    # Leading non-alpha (e.g., a scan-id number) — regex requires alpha start
    assert reconcile.derive_from_filename("123456_Totton1965.pdf") == ("", None)


# ── merge behavior (integration against an in-memory DB) ────────────

def _make_db():
    conn = sqlite3.connect(":memory:")
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


def _insert_work(conn, work_id, title, year, corpus_hash=None, in_corpus=0, source="cited_reference"):
    now = time.time()
    conn.execute(
        """INSERT INTO works (work_id, guid_type, title, year, journal, doi,
           corpus_hash, in_corpus, source, confidence, created_at, updated_at)
           VALUES (?, 'corpus_key', ?, ?, '', NULL, ?, ?, ?, 1.0, ?, ?)""",
        (work_id, title, year, corpus_hash, in_corpus, source, now, now),
    )


def _insert_author(conn, work_id, position, surname):
    conn.execute(
        """INSERT INTO work_authors (work_id, position, surname, surname_normalized, forename)
           VALUES (?, ?, ?, ?, '')""",
        (work_id, position, surname, reconcile.normalize_for_key(surname)),
    )


def test_merge_redirects_citations_and_flips_ghost():
    conn = _make_db()
    # The corpus paper's garbage Phase-1 row
    _insert_work(conn, "corpus:garbage|1965|", "Table (continued)", 1965,
                 corpus_hash="abc123", in_corpus=1, source="corpus_paper")
    _insert_author(conn, "corpus:garbage|1965|", 0, "Garbage")
    # A ghost cited-reference row — correct metadata, no corpus_hash
    _insert_work(conn, "corpus:totton|1965|a synopsis", "A synopsis", 1965)
    _insert_author(conn, "corpus:totton|1965|a synopsis", 0, "Totton")
    # Two citations: one INCOMING to the ghost (other paper citing Totton 1965),
    # one OUTGOING from the Phase-1 row (Totton's own reference to someone else)
    _insert_work(conn, "other_citing", "Other", 2000, corpus_hash="hhh", in_corpus=1, source="corpus_paper")
    _insert_work(conn, "someone_else", "Something", 1900)
    conn.execute(
        "INSERT INTO citations VALUES ('other_citing', 'corpus:totton|1965|a synopsis', 'hhh', 'b1', '', 'alias_exact', 1.0)"
    )
    conn.execute(
        "INSERT INTO citations VALUES ('corpus:garbage|1965|', 'someone_else', 'abc123', 'b2', '', 'alias_exact', 1.0)"
    )
    conn.commit()

    reconcile.merge_phase1_into_ghost(
        conn, "corpus:garbage|1965|", "corpus:totton|1965|a synopsis", "abc123"
    )
    conn.commit()

    # Ghost row is now the corpus paper
    row = conn.execute(
        "SELECT in_corpus, corpus_hash, source FROM works WHERE work_id = ?",
        ("corpus:totton|1965|a synopsis",),
    ).fetchone()
    assert row == (1, "abc123", "corpus_paper")

    # Phase-1 row is gone
    assert conn.execute(
        "SELECT 1 FROM works WHERE work_id = ?", ("corpus:garbage|1965|",),
    ).fetchone() is None

    # Incoming citation still points at the ghost
    assert conn.execute(
        "SELECT 1 FROM citations WHERE cited_work_id = ? AND citing_work_id = ?",
        ("corpus:totton|1965|a synopsis", "other_citing"),
    ).fetchone() == (1,)

    # Outgoing citation redirected onto the ghost as citing_work_id
    assert conn.execute(
        "SELECT 1 FROM citations WHERE citing_work_id = ? AND cited_work_id = ?",
        ("corpus:totton|1965|a synopsis", "someone_else"),
    ).fetchone() == (1,)

    # Old citing_work_id no longer appears
    assert conn.execute(
        "SELECT 1 FROM citations WHERE citing_work_id = ?", ("corpus:garbage|1965|",),
    ).fetchone() is None


def test_merge_is_idempotent_when_ids_equal():
    conn = _make_db()
    _insert_work(conn, "x", "t", 1999, corpus_hash="h", in_corpus=1, source="corpus_paper")
    conn.commit()
    # Should be a no-op
    reconcile.merge_phase1_into_ghost(conn, "x", "x", "h")
    row = conn.execute("SELECT in_corpus, corpus_hash FROM works WHERE work_id='x'").fetchone()
    assert row == (1, "h")


# ── pick_winner ─────────────────────────────────────────────────────

def test_pick_winner_single_candidate_high_score():
    scored = [("w1", 95, "title", 10)]
    assert reconcile.pick_winner(scored, min_score=80, margin=15) == ("matched", "w1", 95)


def test_pick_winner_single_candidate_low_score():
    scored = [("w1", 55, "title", 10)]
    assert reconcile.pick_winner(scored, min_score=80, margin=15) == ("low_score", "w1", 55)


def test_pick_winner_ambiguous_when_runner_close_and_cites_similar():
    # Close title scores, similar citation counts → truly ambiguous
    scored = [("w1", 92, "t1", 10), ("w2", 85, "t2", 8)]
    status, wid, score = reconcile.pick_winner(scored, min_score=80, margin=15)
    assert status == "ambiguous"
    assert wid == "w1"
    assert score == 92


def test_pick_winner_matched_with_clear_margin():
    scored = [("w1", 95, "t1", 5), ("w2", 70, "t2", 10)]
    status, wid, _ = reconcile.pick_winner(scored, min_score=80, margin=15)
    assert status == "matched"
    assert wid == "w1"


def test_pick_winner_no_candidates():
    assert reconcile.pick_winner([], 80, 15) == ("no_candidates", None, 0)


def test_pick_winner_citation_dominance_breaks_fuzzy_tie():
    # Title scores tied at 100 but the canonical variant has 190 cites
    # vs 12 for the next — dominance > 3x → match the canonical one.
    scored = [
        ("canonical", 100, "A synopsis of the Siphonophora", 190),
        ("variant", 100, "A Synopsis of the Siphonophora. London", 12),
    ]
    status, wid, _ = reconcile.pick_winner(scored, min_score=80, margin=15)
    assert status == "matched"
    assert wid == "canonical"


def test_pick_winner_picks_canonical_even_when_shorter_title_scores_higher():
    # Real case from Totton 1954: the short-title ghost scores 100 on
    # partial_ratio because its full title fits in the source text,
    # while the long (canonical) ghost scores 98 because the full title
    # runs past the end of the search window. Citation counts tell the
    # truth: 112 vs 15. Within-margin tier → pick by cites.
    scored = [
        ("short_ghost", 100, "Siphonophora of the Indian Ocean", 15),
        ("canonical", 98, "Siphonophora of the Indian Ocean together with...", 112),
        ("other", 97, "other fragment", 1),
    ]
    status, wid, _ = reconcile.pick_winner(scored, min_score=80, margin=15)
    assert status == "matched"
    assert wid == "canonical"


def test_pick_winner_citation_dominance_respects_dominance_factor():
    # Similar counts (not 2x dominance) → stays ambiguous
    scored = [
        ("canonical", 100, "Title A", 12),
        ("variant", 100, "Title B", 11),
    ]
    status, wid, _ = reconcile.pick_winner(scored, min_score=80, margin=15)
    assert status == "ambiguous"
    assert wid == "canonical"


def test_pick_winner_lensia_case_with_moderate_dominance():
    # From the real Totton 1965 Lensia case: canonical has 14 cites,
    # closest fragment ghost has 6. 14 >= 2 × 6 → should match.
    scored = [
        ("bmnh_fragment", 100, "British Museum (Natural History)", 6),
        ("canonical", 98, "A new species of Lensia (Siphonophora Diphyidae) ...", 14),
        ("ocr_variant", 98, "A new species of LRnsia", 1),
    ]
    status, wid, _ = reconcile.pick_winner(scored, min_score=80, margin=15)
    assert status == "matched"
    assert wid == "canonical"


def test_pick_winner_citation_dominance_requires_nonzero_winner():
    # Zero-cite winner should never beat a zero-cite runner via dominance
    scored = [
        ("w1", 100, "t", 0),
        ("w2", 100, "t", 0),
    ]
    status, _, _ = reconcile.pick_winner(scored, min_score=80, margin=15)
    assert status == "ambiguous"
