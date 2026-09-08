"""`get_missing_references` withholds rows that cannot be leads (#155).

The tool ranks candidates by raw citation-string count, and its top ranks
were dominated by unreconciled variants of works already in the corpus.
v1.3 closed the resolver-safe half — DOI normalization, cross-block
author-set matching — which collapsed the known false-positive clusters
and left roughly 96 title/year-only review leads that no automated rule
can adjudicate without either a similarity threshold loose enough to
merge distinct works or a per-block LLM pass.

This is the cheap slice PLAN.md asked for, plus saying the rest out loud.
A row with **no title and no year** leaves nothing to search for: it is a
mis-parsed reference string counted as a work, not a lead. Measured on
the reference corpus, that is 477 of 6,953 rows at the default threshold
(6.9%) — and it matters because they outrank real gaps.
``corpus:|unknown|``, an empty-titled node with 30 citations, sat 11th in
the default output, above Bigelow 1906, which the audit confirmed is
genuinely missing.
"""
from __future__ import annotations

import sqlite3
import types

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import get_missing_references


@pytest.fixture
def served(tmp_path):
    """An authority DB with a real lead, a titled-but-undated one, a dated
    -but-untitled one, and the parse-debris shape."""
    from bib.authority import create_schema, insert_authors, insert_work

    db = tmp_path / "biblio_authority.sqlite"
    conn = sqlite3.connect(db)
    create_schema(conn)

    rows = [
        # (work_id, title, year, n_citations)
        ("10.5479/si.00963801.44-1946.1", "Medusae and Siphonophorae", 1906, 23),
        ("corpus:beklemishev|1969|principles", "Principles of Comparative Anatomy", 1969, 9),
        ("corpus:notitle|1888|", "", 1888, 12),        # dated, untitled
        ("corpus:nodate|unknown|x", "A Title With No Year", None, 11),  # titled, undated
        ("corpus:|unknown|", "", None, 30),            # parse debris, top-ranked
        ("corpus:c|unknown|", "", None, 15),           # parse debris
    ]
    for wid, title, year, _n in rows:
        insert_work(conn, wid, "corpus_key", title=title, year=year, journal="",
                    doi="", corpus_hash=None, in_corpus=False,
                    source="cited_reference", confidence=0.5)
        insert_authors(conn, wid, [("Someone", "")])

    # `citations` is unique on (citing, cited, citing_hash), so a count of
    # N means N distinct citing papers — which is the metric the tool
    # reports. Make that many citers.
    n_citers = max(n for _w, _t, _y, n in rows)
    for i in range(n_citers):
        insert_work(conn, f"citer{i}", "corpus_key", title=f"Citing {i}",
                    year=2020, journal="", doi="", corpus_hash=f"{i:012d}",
                    in_corpus=True, source="corpus_paper", confidence=1.0)
    for wid, _t, _y, n in rows:
        for i in range(n):
            conn.execute(
                """INSERT INTO citations
                   (citing_work_id, cited_work_id, citing_corpus_hash,
                    match_method, match_score)
                   VALUES (?, ?, ?, ?, ?)""",
                (f"citer{i}", wid, f"{i:012d}", "test", 1.0),
            )
    conn.commit()
    conn.close()

    original = mcp_app._INDEX
    mcp_app.set_index(types.SimpleNamespace(biblio_db=BiblioAuthority(db)))
    yield
    mcp_app.set_index(original)


def _ids(rows):
    return [r["work_id"] for r in rows]


# --- the slice ----------------------------------------------------------


def test_parse_debris_is_withheld_however_often_it_is_cited(served):
    """`corpus:|unknown|` has 30 citations — the highest count in the
    fixture — and is still not a lead, because there is nothing to search
    for."""
    out = get_missing_references(min_citations=2, limit=50)
    assert "corpus:|unknown|" not in _ids(out)
    assert "corpus:c|unknown|" not in _ids(out)


def test_a_real_lead_survives_and_outranks_what_is_left(served):
    out = get_missing_references(min_citations=2, limit=50)
    assert _ids(out)[0] == "10.5479/si.00963801.44-1946.1"


def test_a_title_alone_is_enough_to_keep_a_row(served):
    """Undated is not unusable — you can search for a title."""
    assert "corpus:nodate|unknown|x" in _ids(
        get_missing_references(min_citations=2, limit=50))


def test_a_year_alone_is_enough_to_keep_a_row(served):
    """Untitled but dated still narrows a search, and it is what an
    OCR-damaged entry for a known work looks like — exactly the case
    #155's residual 96 leads are made of."""
    assert "corpus:notitle|1888|" in _ids(
        get_missing_references(min_citations=2, limit=50))


def test_only_title_and_year_together_being_absent_withholds(served):
    out = _ids(get_missing_references(min_citations=1, limit=50))
    assert set(out) == {
        "10.5479/si.00963801.44-1946.1",
        "corpus:beklemishev|1969|principles",
        "corpus:notitle|1888|",
        "corpus:nodate|unknown|x",
    }


def test_the_withheld_count_is_reported(served, caplog):
    """Silently dropping them would be the same class of problem as
    counting them."""
    with caplog.at_level("INFO"):
        get_missing_references(min_citations=2, limit=50)
    assert "withheld 2" in caplog.text
    assert "#155" in caplog.text


def test_the_count_respects_the_citation_threshold(served, caplog):
    """The report has to describe the query that was run, not the table."""
    with caplog.at_level("INFO"):
        get_missing_references(min_citations=20, limit=50)
    assert "withheld 1" in caplog.text


def test_nothing_withheld_says_nothing(served, caplog):
    with caplog.at_level("INFO"):
        get_missing_references(min_citations=100, limit=50)
    assert "withheld" not in caplog.text


# --- saying the rest out loud -------------------------------------------


def test_the_tool_declares_itself_best_effort():
    """PLAN.md's other option, and they are complementary: the slice
    removes what is definitely not a lead, and the docstring says the
    remainder still needs verifying. The docstring is the MCP tool
    description, so it is the channel a client actually reads."""
    doc = (get_missing_references.fn if hasattr(get_missing_references, "fn")
           else get_missing_references).__doc__
    assert "Best-effort" in doc
    assert "not proof that a work is absent" in doc
    assert "resolve_reference" in doc
    assert "reference_reconciliation.py" in doc
    assert "#155" in doc


def test_the_signature_is_unchanged():
    """1.0 freezes tool inputs and the return annotation. The slice is
    unconditional rather than a parameter, so no caller changes a call —
    and there is no query the withheld population answers."""
    import inspect

    sig = inspect.signature(
        get_missing_references.fn if hasattr(get_missing_references, "fn")
        else get_missing_references)
    assert list(sig.parameters) == [
        "min_citations", "year_from", "year_to", "limit"]
