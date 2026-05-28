"""Integration tests for the format_citation MCP tool (#79).

Wires bib.format + BiblioAuthority.provenance + the warning-text
table the LLM is instructed to preserve verbatim. The contract
this test pins is the public output shape: an LLM client that
follows the default_instructions guidance pastes the ``formatted``,
``inline``, and ``warning`` strings literally into its document,
so drift in any of them is user-visible.
"""
from __future__ import annotations

import sqlite3
import time
import types
from pathlib import Path

import pytest

from bib.authority import create_schema
from mcpsrv import app as mcp_app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import _PROVENANCE_WARNING, format_citation


def _make_authority_db(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    create_schema(conn)
    return conn


def _insert_work(
    conn: sqlite3.Connection,
    work_id: str,
    *,
    title: str,
    year: int,
    journal: str | None = None,
    doi: str | None = None,
    source: str = "cited_reference",
    in_corpus: int = 0,
    corpus_hash: str | None = None,
    bib_imported_at: float | None = None,
) -> None:
    now = time.time()
    conn.execute(
        """INSERT INTO works
           (work_id, guid_type, title, year, journal, doi, in_corpus,
            corpus_hash, source, confidence, bib_imported_at,
            created_at, updated_at)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (work_id, "corpus_key", title, year, journal, doi, in_corpus,
         corpus_hash, source, 1.0, bib_imported_at, now, now),
    )


def _insert_author(
    conn: sqlite3.Connection,
    work_id: str,
    position: int,
    surname: str,
    forename: str | None = None,
) -> None:
    conn.execute(
        """INSERT INTO work_authors
           (work_id, position, surname, surname_normalized, forename)
           VALUES (?, ?, ?, ?, ?)""",
        (work_id, position, surname, surname.strip().lower(), forename),
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
def index(tmp_path: Path):
    """Build a tiny biblio_authority.sqlite, wrap with BiblioAuthority,
    inject it into the MCP module-global so the tool's
    _need_index() call succeeds.

    Three works seeded — one per provenance tier — so the tool can
    be exercised against each branch in a single test pass.
    """
    db_path = tmp_path / "biblio.sqlite"
    conn = _make_authority_db(db_path)

    # Tier 1: human-curated via bib import.
    _insert_work(conn, "corpus:dunn|2005|marrus",
                 title="Marrus claudanielis, a new species",
                 year=2005, journal="Invertebrate Biology",
                 doi="10.1111/j.1744-7410.2005.00012.x",
                 source="corpus_paper", in_corpus=1, corpus_hash="dunnhash000",
                 bib_imported_at=time.time())
    _insert_author(conn, "corpus:dunn|2005|marrus", 0, "Dunn", "Casey W.")

    # Tier 2: grobid_reconciled corpus_paper without bib import.
    _insert_work(conn, "corpus:siebert|2011|hapless",
                 title="Functional morphology of a hapless siphonophore",
                 year=2011, journal="Marine Biology",
                 source="corpus_paper", in_corpus=1, corpus_hash="siebrt00000")
    _insert_author(conn, "corpus:siebert|2011|hapless", 0, "Siebert", "Stefan")
    _insert_author(conn, "corpus:siebert|2011|hapless", 1, "Pugh", "Philip R.")

    # Tier 3: unresolved cited_reference (only a low-fuzzy citation
    # points at it).
    _insert_work(conn, "corpus:lensia|1965|unknown",
                 title="A Lensia species redescription",
                 year=1965, source="cited_reference")
    _insert_author(conn, "corpus:lensia|1965|unknown", 0, "Lensia", None)
    _insert_citation(conn, "corpus:siebert|2011|hapless",
                     "corpus:lensia|1965|unknown",
                     method="author_year_only", score=0.6)

    conn.commit()

    # Inject into the MCP global so _need_index() returns our test
    # index. We only need biblio_db wired — the tool doesn't touch
    # taxonomy_db, papers, etc.
    fake_index = types.SimpleNamespace(biblio_db=BiblioAuthority(db_path))
    original = mcp_app._INDEX
    mcp_app.set_index(fake_index)
    yield fake_index
    mcp_app.set_index(original)


def test_resolves_work_id_with_bib_tier(index):
    out = format_citation(work_id="corpus:dunn|2005|marrus")
    assert out["provenance"] == "bib"
    assert out["warning"] == ""
    assert "Dunn, C. W." in out["formatted"]
    assert "Marrus claudanielis" in out["formatted"]
    assert "*Invertebrate Biology*" in out["formatted"]
    assert "https://doi.org/10.1111/j.1744-7410.2005.00012.x" in out["formatted"]
    assert out["inline"] == "(Dunn, 2005)"


def test_resolves_work_id_with_grobid_reconciled_tier(index):
    out = format_citation(work_id="corpus:siebert|2011|hapless")
    assert out["provenance"] == "grobid_reconciled"
    assert out["warning"] == \
        "* generated via reconciliation in corpus, check if correct"
    assert "Siebert, S., & Pugh, P. R." in out["formatted"]
    assert out["inline"] == "(Siebert & Pugh, 2011)"


def test_resolves_work_id_with_unresolved_tier(index):
    out = format_citation(work_id="corpus:lensia|1965|unknown")
    assert out["provenance"] == "unresolved"
    assert out["warning"] == \
        "* reference not present in bibliography, check if correct"


def test_resolves_paper_hash(index):
    """corpus_hash is the bridge from in-corpus PDFs to authority rows."""
    out = format_citation(paper_hash="dunnhash000")
    assert out["work_id"] == "corpus:dunn|2005|marrus"
    assert out["provenance"] == "bib"


def test_resolves_query_free_text_unique(index):
    """A query that matches exactly one work returns the formatted
    citation — no disambiguation needed."""
    out = format_citation(query="Dunn 2005")
    assert out["work_id"] == "corpus:dunn|2005|marrus"
    assert out["provenance"] == "bib"
    assert "Marrus claudanielis" in out["formatted"]


def test_resolves_query_with_title_fragment_narrows(index):
    out = format_citation(query="Siebert 2011 hapless")
    assert out["work_id"] == "corpus:siebert|2011|hapless"


def test_query_not_found_returns_error_not_fabrication(index):
    """The contract: when a work isn't in the corpus, do NOT
    synthesize a citation. Return an explicit not_found so the
    model says 'not in the corpus' in its prose."""
    out = format_citation(query="NotAnAuthor 2099")
    assert out["error"] == "not_found"
    assert out["queried"] == "NotAnAuthor 2099"
    assert out["parsed_author"] == "NotAnAuthor"
    assert out["parsed_year"] == 2099
    assert "formatted" not in out


def test_ambiguous_query_returns_match_list(index):
    """Seed a second work that would also match 'Dunn 2005'."""
    conn = sqlite3.connect(index.biblio_db.conn.execute(
        "PRAGMA database_list").fetchone()[2])
    _insert_work(conn, "corpus:dunn|2005|other",
                 title="An unrelated 2005 paper",
                 year=2005, source="cited_reference")
    _insert_author(conn, "corpus:dunn|2005|other", 0, "Dunn", "Other")
    conn.commit()
    conn.close()

    out = format_citation(query="Dunn 2005")
    assert out["error"] == "ambiguous"
    assert len(out["matches"]) == 2
    assert {m["work_id"] for m in out["matches"]} == {
        "corpus:dunn|2005|marrus",
        "corpus:dunn|2005|other",
    }


def test_rejects_missing_resolver_args(index):
    out = format_citation()
    assert "exactly one of" in out["error"]


def test_rejects_multiple_resolver_args(index):
    out = format_citation(query="Dunn 2005", work_id="corpus:dunn|2005|marrus")
    assert "exactly one of" in out["error"]


def test_rejects_unknown_style(index):
    out = format_citation(work_id="corpus:dunn|2005|marrus", style="vancouver")
    assert "unknown style" in out["error"]
    assert "author-year" in out["supported_styles"]


def test_warning_table_is_complete():
    """The _PROVENANCE_WARNING table must cover every tier the
    provenance helper returns. A new tier added to BiblioAuthority.
    provenance() without a matching warning entry would crash the
    tool — this test forces them to evolve together."""
    expected_tiers = {"bib", "grobid_reconciled", "unresolved"}
    assert set(_PROVENANCE_WARNING.keys()) == expected_tiers
