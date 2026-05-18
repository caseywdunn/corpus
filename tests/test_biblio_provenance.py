"""Tests for BiblioAuthority.provenance() — the tier classifier
that drives the format_citation MCP tool's warning footnotes (#79).

The three tiers must be sharply distinguished — they're what stops
the LLM from emitting a tier-3 (unresolved) citation as if it were
tier-1 (human-curated bib). Drift here directly produces the
hallucination class the tool exists to prevent.
"""
from __future__ import annotations

import sqlite3
import time
from pathlib import Path

import pytest

from bib.authority import create_schema
from mcpsrv.indexes import BiblioAuthority


def _make_authority_db(path: Path) -> sqlite3.Connection:
    """Production schema via create_schema, plus a writable handle
    so each test can seed the rows it needs. BiblioAuthority opens
    the same file read-only."""
    conn = sqlite3.connect(path)
    create_schema(conn)
    return conn


def _insert_work(
    conn: sqlite3.Connection,
    work_id: str,
    *,
    source: str,
    bib_imported_at: float | None = None,
    in_corpus: int = 0,
    corpus_hash: str | None = None,
) -> None:
    now = time.time()
    conn.execute(
        """INSERT INTO works
           (work_id, guid_type, source, in_corpus, corpus_hash,
            bib_imported_at, created_at, updated_at)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
        (work_id, "corpus_key", source, in_corpus, corpus_hash,
         bib_imported_at, now, now),
    )


def _insert_citation(
    conn: sqlite3.Connection,
    citing_work_id: str,
    cited_work_id: str,
    *,
    method: str,
    score: float,
) -> None:
    conn.execute(
        """INSERT INTO citations
           (citing_work_id, cited_work_id, citing_corpus_hash,
            match_method, match_score)
           VALUES (?, ?, ?, ?, ?)""",
        (citing_work_id, cited_work_id, "h00deadbeef00", method, score),
    )


@pytest.fixture
def authority(tmp_path: Path):
    """Yields (writable_conn, BiblioAuthority(read-only))."""
    db_path = tmp_path / "biblio.sqlite"
    conn = _make_authority_db(db_path)
    yield conn, BiblioAuthority(db_path)


def test_unresolved_when_work_missing(authority):
    conn, ba = authority
    conn.commit()
    assert ba.provenance("no:such:work") == "unresolved"


def test_bib_tier_when_bib_imported_at_set(authority):
    """Any row touched by the importer is tier-1 regardless of its
    source or how it was originally resolved — the human edit is
    the trump card."""
    conn, ba = authority
    _insert_work(conn, "corpus:dunn|2005", source="cited_reference",
                 bib_imported_at=time.time())
    conn.commit()
    assert ba.provenance("corpus:dunn|2005") == "bib"


def test_corpus_paper_is_grobid_reconciled_by_default(authority):
    """A corpus_paper row without bib import is tier-2: we have the
    PDF and Grobid extracted its metadata, but a human hasn't
    blessed the fields. Default to the reconciliation warning."""
    conn, ba = authority
    _insert_work(conn, "corpus:siebert|2011", source="corpus_paper",
                 in_corpus=1, corpus_hash="abcabcabcabc")
    conn.commit()
    assert ba.provenance("corpus:siebert|2011") == "grobid_reconciled"


def test_cited_reference_with_doi_exact_is_grobid_reconciled(authority):
    conn, ba = authority
    _insert_work(conn, "10.1/foo", source="cited_reference")
    _insert_work(conn, "10.1/citer", source="corpus_paper", in_corpus=1)
    _insert_citation(conn, "10.1/citer", "10.1/foo",
                     method="doi_exact", score=1.0)
    conn.commit()
    assert ba.provenance("10.1/foo") == "grobid_reconciled"


def test_cited_reference_with_alias_exact_is_grobid_reconciled(authority):
    conn, ba = authority
    _insert_work(conn, "corpus:totton|1965", source="cited_reference")
    _insert_work(conn, "10.1/citer", source="corpus_paper", in_corpus=1)
    _insert_citation(conn, "10.1/citer", "corpus:totton|1965",
                     method="alias_exact", score=1.0)
    conn.commit()
    assert ba.provenance("corpus:totton|1965") == "grobid_reconciled"


def test_cited_reference_with_low_fuzzy_is_unresolved(authority):
    """Below the 0.9 reconciled threshold — the cascade had to fall
    back to a heuristic match. Per #2, fuzzy < 0.9 is exactly the
    class of "matched" reference where amalgamation happens."""
    conn, ba = authority
    _insert_work(conn, "corpus:ghost|1900", source="cited_reference")
    _insert_work(conn, "10.1/citer", source="corpus_paper", in_corpus=1)
    _insert_citation(conn, "10.1/citer", "corpus:ghost|1900",
                     method="title_fuzzy", score=0.7)
    conn.commit()
    assert ba.provenance("corpus:ghost|1900") == "unresolved"


def test_cited_reference_with_author_year_only_is_unresolved(authority):
    """The author_year_only fallback explicitly amplified the
    Lensia/Synopsis misrouting in #2. Must be tier-3 regardless of
    its 0.6 'confidence'."""
    conn, ba = authority
    _insert_work(conn, "corpus:lensia|1965", source="cited_reference")
    _insert_work(conn, "10.1/citer", source="corpus_paper", in_corpus=1)
    _insert_citation(conn, "10.1/citer", "corpus:lensia|1965",
                     method="author_year_only", score=0.6)
    conn.commit()
    assert ba.provenance("corpus:lensia|1965") == "unresolved"


def test_cited_reference_with_no_citations_is_unresolved(authority):
    """A cited_reference row with no citation edges pointing at it
    (somehow — a `new_work` ghost from issue #1's ghost-merge path
    before any incoming citation lands) is unresolved."""
    conn, ba = authority
    _insert_work(conn, "corpus:orphan|1900", source="cited_reference")
    conn.commit()
    assert ba.provenance("corpus:orphan|1900") == "unresolved"


def test_bib_trumps_low_confidence_citation(authority):
    """A row that was poorly resolved by the cascade BUT subsequently
    touched by bib/importer.py is tier-1 — the human edit is what
    matters now."""
    conn, ba = authority
    _insert_work(conn, "corpus:edited|1965",
                 source="cited_reference",
                 bib_imported_at=time.time())
    _insert_work(conn, "10.1/citer", source="corpus_paper", in_corpus=1)
    _insert_citation(conn, "10.1/citer", "corpus:edited|1965",
                     method="author_year_only", score=0.6)
    conn.commit()
    assert ba.provenance("corpus:edited|1965") == "bib"
