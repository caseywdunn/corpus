"""Bib-provenance preservation through import + reconciliation (#100).

A reference present in the user-edited .bib must stay authoritative
`bib` provenance (bib_imported_at IS NOT NULL) — including (a) when an
unchanged .bib is re-imported, (b) after a Phase-1 corpus row is
reconciled (merged) onto a pre-existing ghost cited-reference row, and
(c) for curated references matched by work_id when they carry neither
corpus_hash nor DOI.
"""
from __future__ import annotations

import sqlite3
import time

import pytest

from bib.authority import create_schema
from bib.importer import apply_entry, find_matching_work_id
from bib.reconcile import merge_phase1_into_ghost


def _insert_work(conn, work_id, *, title, year=None, journal=None, doi=None,
                 source="cited_reference", in_corpus=0, corpus_hash=None,
                 bib_imported_at=None):
    now = time.time()
    conn.execute(
        """INSERT INTO works
           (work_id, guid_type, title, year, journal, doi, in_corpus,
            corpus_hash, source, confidence, bib_imported_at,
            created_at, updated_at)
           VALUES (?, 'corpus_key', ?, ?, ?, ?, ?, ?, ?, 1.0, ?, ?, ?)""",
        (work_id, title, year, journal, doi, in_corpus, corpus_hash,
         source, bib_imported_at, now, now),
    )


def _bib_imported_at(conn, work_id):
    row = conn.execute(
        "SELECT bib_imported_at FROM works WHERE work_id = ?", (work_id,),
    ).fetchone()
    return row[0] if row else None


@pytest.fixture
def conn():
    c = sqlite3.connect(":memory:")
    create_schema(c)
    yield c
    c.close()


# (a) unchanged re-import still stamps bib provenance ----------------------


def test_unchanged_reimport_stamps_bib_imported_at(conn):
    _insert_work(conn, "corpus:dunn|2005|x", title="Marrus claudanielis",
                 year=2005, bib_imported_at=None)
    # Entry matches the row (title only, no author) → empty diff.
    changed = apply_entry(conn, "corpus:dunn|2005|x",
                          {"title": "Marrus claudanielis"})
    assert changed == 0                       # no field-level change
    assert _bib_imported_at(conn, "corpus:dunn|2005|x") is not None  # but stamped


# (c) match by work_id field when corpus_hash + DOI are absent -------------


def test_find_matching_by_work_id_field(conn):
    _insert_work(conn, "corpus:totton|1965|synopsis",
                 title="A synopsis of the Siphonophora", year=1965)
    wid, method = find_matching_work_id(
        conn, {"work_id": "corpus:totton|1965|synopsis"},
    )
    assert (wid, method) == ("corpus:totton|1965|synopsis", "work_id")


def test_find_matching_prefers_corpus_hash_then_doi_then_work_id(conn):
    _insert_work(conn, "w:byhash", title="By hash", corpus_hash="hash00000000")
    assert find_matching_work_id(conn, {"corpus_hash": "hash00000000"}) == (
        "w:byhash", "corpus_hash")
    assert find_matching_work_id(conn, {"work_id": "nope"}) == (None, "no_match")


# (b) reconciliation carries bib authority onto the surviving ghost --------


def test_reconcile_carries_bib_fields_and_stamp_forward(conn):
    # Ghost: GROBID-derived cited reference, no bib stamp.
    _insert_work(conn, "ghost:1", title="Grobid Title", year=1999,
                 journal="Wrong J", source="cited_reference")
    # Phase-1: the same paper, freshly bib-imported with curated fields.
    _insert_work(conn, "phase1:1", title="Curated Title", year=2000,
                 journal="Right J", doi="10.1/right", source="corpus_paper",
                 in_corpus=1, bib_imported_at=time.time())

    merge_phase1_into_ghost(conn, "phase1:1", "ghost:1", "corpushash01")

    row = conn.execute(
        "SELECT title, year, journal, doi, source, in_corpus, corpus_hash, "
        "bib_imported_at FROM works WHERE work_id = 'ghost:1'"
    ).fetchone()
    assert row[7] is not None                      # bib provenance preserved
    assert (row[0], row[1], row[2], row[3]) == (
        "Curated Title", 2000, "Right J", "10.1/right")   # bib fields win
    assert row[4] == "corpus_paper" and row[5] == 1 and row[6] == "corpushash01"
    # Phase-1 row removed.
    assert conn.execute(
        "SELECT 1 FROM works WHERE work_id = 'phase1:1'").fetchone() is None


def test_reconcile_without_bib_keeps_ghost_fields(conn):
    _insert_work(conn, "ghost:2", title="Ghost Title", year=1999,
                 source="cited_reference")
    _insert_work(conn, "phase1:2", title="Phase1 Title", year=2000,
                 source="corpus_paper", in_corpus=1, bib_imported_at=None)

    merge_phase1_into_ghost(conn, "phase1:2", "ghost:2", "corpushash02")

    row = conn.execute(
        "SELECT title, source, bib_imported_at FROM works "
        "WHERE work_id = 'ghost:2'"
    ).fetchone()
    # No bib authority on either side → ghost fields untouched, still no stamp.
    assert row[0] == "Ghost Title"
    assert row[1] == "corpus_paper"
    assert row[2] is None
