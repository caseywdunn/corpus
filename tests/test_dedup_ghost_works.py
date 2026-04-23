"""Unit tests for dedup_ghost_works.py."""

from __future__ import annotations

import importlib.util
import sqlite3
import sys
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_PATH = REPO_ROOT / "dedup_ghost_works.py"

_spec = importlib.util.spec_from_file_location("dedup_ghost_works", SCRIPT_PATH)
dedup = importlib.util.module_from_spec(_spec)
sys.modules["dedup_ghost_works"] = dedup
_spec.loader.exec_module(dedup)


# ── pick_canonical + should_merge — no DB needed ────────────────────

def _row(work_id, title, cites, guid_type="corpus_key", in_corpus=0):
    """Mimic the sqlite3.Row dict-access surface the helpers expect."""
    return {
        "work_id": work_id, "title": title, "cites": cites,
        "guid_type": guid_type, "year": 1888, "surname": "haeckel",
        "in_corpus": in_corpus,
    }


def test_pick_canonical_prefers_substantive_title_high_cites():
    family = [
        _row("a", "Report on the Siphonophorae", 23),
        _row("b", "", 4),
        _row("c", "Rep. sci. Res", 2),
        _row("d", "Siphonophorae of the Challenger", 5),
    ]
    canon = dedup.pick_canonical(family, min_alpha=20)
    assert canon["work_id"] == "a"


def test_pick_canonical_none_when_all_titles_empty_or_short():
    family = [
        _row("a", "", 10),
        _row("b", "Rep. sci.", 5),
        _row("c", None, 1),
    ]
    assert dedup.pick_canonical(family, min_alpha=20) is None


def test_should_merge_empty_title():
    cand = _row("x", "", 1)
    canon = _row("c", "Report on the Siphonophorae", 23)
    merge, reason = dedup.should_merge(cand, canon)
    assert merge is True
    assert reason == "empty_title"


def test_should_merge_short_fragment():
    cand = _row("x", "Rep. sci. Res", 2)
    canon = _row("c", "Report on the Siphonophorae", 23)
    merge, reason = dedup.should_merge(cand, canon)
    assert merge is True
    assert reason == "fragment"


def test_should_merge_fuzzy_variant():
    cand = _row("x", "Report on the Siphonophorae. Report on the Scientific Results of HMS Challenger", 5)
    canon = _row("c", "Report on the Siphonophorae", 23)
    merge, reason = dedup.should_merge(cand, canon)
    assert merge is True
    assert reason.startswith("fuzzy(")


def test_should_not_merge_genuinely_different_work():
    # Haeckel 1888's Prodromus Medusarum vs. Challenger Siphonophorae
    # report — same surname+year but different underlying works.
    cand = _row("x", "Prodromus Systemas Medusarum Jena", 1)
    canon = _row("c", "Report on the Siphonophorae", 23)
    merge, reason = dedup.should_merge(cand, canon)
    assert merge is False
    assert reason.startswith("distinct(")


def test_should_not_merge_totton_synopsis_vs_lensia():
    # Regression: the pair that motivated issue #2 — Totton 1965's
    # A Synopsis vs. his Lensia paper — must stay separate.
    cand = _row("x",
                "A new species of Lensia Siphonophora Diphyidae from the coastal waters of Vancouver",
                14)
    canon = _row("c", "A synopsis of the Siphonophora", 190)
    merge, reason = dedup.should_merge(cand, canon)
    assert merge is False


# ── merge_ghost_into_canonical — DB integration ─────────────────────

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


def _insert_work(conn, work_id, title, in_corpus=0):
    now = time.time()
    conn.execute(
        """INSERT INTO works (work_id, guid_type, title, year, journal, doi,
           corpus_hash, in_corpus, source, confidence, created_at, updated_at)
           VALUES (?, 'corpus_key', ?, 1888, '', NULL, NULL, ?,
                   'cited_reference', 1.0, ?, ?)""",
        (work_id, title, in_corpus, now, now),
    )


def _insert_author(conn, work_id, position, surname):
    conn.execute(
        """INSERT INTO work_authors (work_id, position, surname, surname_normalized, forename)
           VALUES (?, ?, ?, ?, '')""",
        (work_id, position, surname, dedup.normalize_for_key(surname)),
    )


def test_merge_redirects_incoming_citations_and_drops_dup():
    conn = _make_db()
    _insert_work(conn, "dup", "Rep. sci. Res.")
    _insert_author(conn, "dup", 0, "Haeckel")
    _insert_work(conn, "canon", "Report on the Siphonophorae")
    _insert_author(conn, "canon", 0, "Haeckel")
    _insert_work(conn, "citer1", "Someone's paper", in_corpus=1)
    _insert_work(conn, "citer2", "Another paper", in_corpus=1)
    conn.execute(
        "INSERT INTO citations VALUES ('citer1','dup','h1','b1','','alias_exact',1.0)"
    )
    conn.execute(
        "INSERT INTO citations VALUES ('citer2','canon','h2','b2','','alias_exact',1.0)"
    )
    conn.commit()

    dedup.merge_ghost_into_canonical(conn, "dup", "canon")
    conn.commit()

    # Dup row gone
    assert conn.execute(
        "SELECT 1 FROM works WHERE work_id='dup'"
    ).fetchone() is None

    # Both citations now point at the canonical
    n = conn.execute(
        "SELECT COUNT(*) FROM citations WHERE cited_work_id='canon'"
    ).fetchone()[0]
    assert n == 2

    # No stragglers pointing at dup
    assert conn.execute(
        "SELECT COUNT(*) FROM citations WHERE cited_work_id='dup'"
    ).fetchone()[0] == 0


def test_merge_is_idempotent_when_ids_equal():
    conn = _make_db()
    _insert_work(conn, "x", "title")
    _insert_author(conn, "x", 0, "X")
    conn.commit()
    dedup.merge_ghost_into_canonical(conn, "x", "x")
    assert conn.execute("SELECT 1 FROM works WHERE work_id='x'").fetchone() is not None


# ── dedup() end-to-end on an in-memory DB ───────────────────────────

def test_dedup_merges_fragment_and_variant_but_preserves_distinct_work():
    conn = _make_db()
    # Canonical: Haeckel 1888 Challenger report
    _insert_work(conn, "canon", "Report on the Siphonophorae")
    _insert_author(conn, "canon", 0, "Haeckel")
    # Fragment (short title): should merge
    _insert_work(conn, "frag", "Rep. sci. Res.")
    _insert_author(conn, "frag", 0, "Haeckel")
    # OCR variant (fuzzy match): should merge
    _insert_work(conn, "variant",
                 "Report on the Siphonophorae collected by HMS Challenger during 1873")
    _insert_author(conn, "variant", 0, "Haeckel")
    # Distinct work (Prodromus Medusarum): should stay.  Needs ≥ 20
    # alpha chars and some incoming cite so load_families sees the
    # family as size ≥ 2 after counting.
    _insert_work(conn, "distinct", "Prodromus Systemas Medusarum Jena")
    _insert_author(conn, "distinct", 0, "Haeckel")
    # Add incoming cites — each from a distinct citing paper so the
    # citations PK (citing, cited, citing_corpus_hash) doesn't collapse
    # post-merge rows. In the real DB the three merged refs come from
    # three different citing papers, which is what this mirrors.
    for i, cited in enumerate(
        ["canon", "canon", "canon", "frag", "variant", "distinct"]
    ):
        _insert_work(conn, f"p{i}", "", in_corpus=1)
        conn.execute(
            f"INSERT INTO citations VALUES "
            f"('p{i}','{cited}','h{i}','b{i}','','x',1.0)"
        )
    conn.commit()

    counts = dedup.dedup(conn)

    assert counts["merged"] == 2       # frag + variant
    assert counts["kept_separate"] == 1  # distinct
    # canon now has its 3 original cites + 1 from frag + 1 from variant
    canon_cites = conn.execute(
        "SELECT COUNT(*) FROM citations WHERE cited_work_id='canon'"
    ).fetchone()[0]
    assert canon_cites == 5
    # distinct Prodromus still standalone
    assert conn.execute(
        "SELECT 1 FROM works WHERE work_id='distinct'"
    ).fetchone() is not None


def test_dedup_skips_short_surname():
    """Grobid-misparse ghosts under single-letter 'surnames' (p, a, c…)
    must not be collapsed — would destroy data across unrelated works."""
    conn = _make_db()
    _insert_work(conn, "a1", "Some title A")
    _insert_author(conn, "a1", 0, "P")
    _insert_work(conn, "a2", "Some title B")
    _insert_author(conn, "a2", 0, "P")
    conn.commit()
    counts = dedup.dedup(conn)
    assert counts["merged"] == 0
    # Both rows should remain
    assert conn.execute("SELECT COUNT(*) FROM works").fetchone()[0] == 2


def test_dedup_merges_ghost_into_in_corpus_canonical():
    """Regression: when a family contains an in_corpus=1 row (the real
    paper) alongside ghost OCR variants, the ghosts should merge INTO
    the corpus-paper row, not into a sibling ghost with a noisier
    title. Models the post-reconciliation state for Totton 1965."""
    conn = _make_db()
    # The in-corpus row (merged from Totton1965a.pdf in #1 reconcile).
    _insert_work(conn, "canon", "A synopsis of the Siphonophora", in_corpus=1)
    _insert_author(conn, "canon", 0, "Haeckel")  # using Haeckel for year=1888
    # Ghost OCR variant — noisier title, fewer direct cites
    _insert_work(conn, "variant_a",
                 "A synopsis of the Siphonophora London British Museum")
    _insert_author(conn, "variant_a", 0, "Haeckel")
    # Unrelated citer to give the canonical some cites
    for i, cited in enumerate(["canon", "canon", "canon", "variant_a"]):
        _insert_work(conn, f"p{i}", "", in_corpus=1)
        conn.execute(
            f"INSERT INTO citations VALUES "
            f"('p{i}','{cited}','h{i}','b{i}','','x',1.0)"
        )
    conn.commit()

    dedup.dedup(conn)

    # Canonical retains its in_corpus=1 status and absorbs the variant
    row = conn.execute(
        "SELECT in_corpus FROM works WHERE work_id='canon'"
    ).fetchone()
    assert row["in_corpus"] == 1
    # Variant merged away
    assert conn.execute(
        "SELECT 1 FROM works WHERE work_id='variant_a'"
    ).fetchone() is None


def test_dedup_never_merges_in_corpus_row_as_candidate():
    """In-corpus rows must never be removed. If a family has two
    in_corpus=1 rows (rare but possible after reconciliation), neither
    is merged out — only the ghosts around them get touched."""
    conn = _make_db()
    _insert_work(conn, "canon", "A synopsis of the Siphonophora", in_corpus=1)
    _insert_author(conn, "canon", 0, "Haeckel")
    # Second in_corpus=1 row — distinct work, same (surname, year)
    _insert_work(conn, "other", "A different synopsis", in_corpus=1)
    _insert_author(conn, "other", 0, "Haeckel")
    # Ghost fragment
    _insert_work(conn, "frag", "Rep. sci.")
    _insert_author(conn, "frag", 0, "Haeckel")
    conn.commit()

    dedup.dedup(conn)

    # Both corpus rows must survive
    for wid in ("canon", "other"):
        assert conn.execute(
            "SELECT 1 FROM works WHERE work_id=?", (wid,)
        ).fetchone() is not None


def test_dedup_ignores_doi_keyed_rows():
    """DOI-keyed rows are issue-#4 territory; this pass must leave
    them alone even when they share surname+year with a corpus-key
    ghost."""
    conn = _make_db()
    _insert_work(conn, "canon", "Report on the Siphonophorae")
    _insert_author(conn, "canon", 0, "Haeckel")
    _insert_work(conn, "frag", "Rep. sci.")
    _insert_author(conn, "frag", 0, "Haeckel")
    # DOI row for same underlying work
    conn.execute(
        """INSERT INTO works (work_id, guid_type, title, year, journal, doi,
           corpus_hash, in_corpus, source, confidence, created_at, updated_at)
           VALUES ('10.5962/bhl.title.3968', 'doi',
                   'Report on the Siphonophorae collected by Challenger',
                   1888, '', '10.5962/bhl.title.3968', NULL, 0,
                   'cited_reference', 1.0, ?, ?)""",
        (time.time(), time.time()),
    )
    _insert_author(conn, "10.5962/bhl.title.3968", 0, "Haeckel")
    conn.commit()

    dedup.dedup(conn)

    # DOI row must survive untouched
    assert conn.execute(
        "SELECT 1 FROM works WHERE work_id='10.5962/bhl.title.3968'"
    ).fetchone() is not None
