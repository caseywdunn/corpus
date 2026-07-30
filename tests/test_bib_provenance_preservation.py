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


# (d) the CALLER must reach the stamping branch (#142 → #152) --------------
#
# `test_unchanged_reimport_stamps_bib_imported_at` above calls apply_entry
# directly, which is why #100 looked complete: the no-diff branch worked,
# but `import_bibtex` early-`continue`d before ever calling it, so the
# branch was unreachable in production. A real round-trip re-import
# stamped only entries that happened to have field edits — 20 of 19,834 in
# the reported corpus — leaving the rest in the reconciliation-warning
# tier (#152), and an authority-DB rebuild re-broke it every time.
#
# These go through the public entry point so the wiring is covered, not
# just the branch.


def _write_bib(tmp_path, body: str):
    p = tmp_path / "lib.bib"
    p.write_text(body, encoding="utf-8")
    return p


def _seed_db(tmp_path, *, corpus_hash="aabbccddeeff", title="Marrus claudanielis",
             year=2005):
    db = tmp_path / "biblio_authority.sqlite"
    c = sqlite3.connect(db)
    create_schema(c)
    _insert_work(
        c, f"corpus:dunn|{year}|marrus", title=title, year=year,
        source="corpus_paper", in_corpus=1, corpus_hash=corpus_hash,
        bib_imported_at=None,
    )
    c.commit()
    c.close()
    return db


def test_import_bibtex_stamps_unchanged_entries(tmp_path):
    """The #142 regression test: an entry that matches the DB exactly must
    still come out of `import_bibtex` stamped into the bib tier."""
    from bib.importer import import_bibtex

    db = _seed_db(tmp_path)
    bib = _write_bib(tmp_path, """@article{dunn2005,
  corpus_hash = {aabbccddeeff},
  title = {Marrus claudanielis},
  year = {2005},
}
""")
    counters = import_bibtex(db, bib)
    assert counters["no_changes"] == 1, counters   # genuinely no field diff
    assert counters["changed"] == 0, counters

    c = sqlite3.connect(db)
    stamped = _bib_imported_at(c, "corpus:dunn|2005|marrus")
    c.close()
    assert stamped is not None, (
        "import_bibtex left bib_imported_at NULL on a matched-but-unchanged "
        "entry — #142 has regressed and #152's warnings are back"
    )


def test_unchanged_reimport_lands_in_the_silent_bib_tier(tmp_path):
    """#152's acceptance criterion, end to end: after re-importing an
    unchanged .bib, the work's provenance tier is `bib`, which is the tier
    `format_citations` emits *without* a reconciliation warning."""
    from mcpsrv.indexes import BiblioAuthority

    from bib.importer import import_bibtex

    db = _seed_db(tmp_path)

    # Before: a corpus_paper row with no bib stamp reads as reconciled, so
    # format_citations warns on it.
    before = BiblioAuthority(db).provenance("corpus:dunn|2005|marrus")
    assert before != "bib", before

    bib = _write_bib(tmp_path, """@article{dunn2005,
  corpus_hash = {aabbccddeeff},
  title = {Marrus claudanielis},
  year = {2005},
}
""")
    import_bibtex(db, bib)
    # ...and again, to prove idempotence — the reported symptom was that a
    # *re*-import of an already-clean .bib changed nothing.
    import_bibtex(db, bib)

    after = BiblioAuthority(db).provenance("corpus:dunn|2005|marrus")
    assert after == "bib", (
        f"provenance is {after!r}, so format_citations still emits "
        "'generated via reconciliation, check if correct' on a work the "
        "user curated by hand (#152)"
    )


def test_dry_run_does_not_stamp(tmp_path):
    """--dry-run must stay side-effect free, including for the new stamp."""
    from bib.importer import import_bibtex

    db = _seed_db(tmp_path)
    bib = _write_bib(tmp_path, """@article{dunn2005,
  corpus_hash = {aabbccddeeff},
  title = {Marrus claudanielis},
  year = {2005},
}
""")
    import_bibtex(db, bib, dry_run=True)
    c = sqlite3.connect(db)
    stamped = _bib_imported_at(c, "corpus:dunn|2005|marrus")
    c.close()
    assert stamped is None, "dry-run wrote bib_imported_at"
