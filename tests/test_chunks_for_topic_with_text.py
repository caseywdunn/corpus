"""Tests for the with_text= mode on get_chunks_for_topic (#82).

Pins the shape difference: with_text=True keeps the existing
contract (full chunk text in `text` field); with_text=False drops
`text` and emits `len_chars` instead. The scan-then-drill-down
workflow that motivates the change relies on this exact contract,
so a regression here directly invalidates the token-cost lever.

The semantic-search code path is tested by injecting a fake topic
searcher into the index — we don't need real LanceDB for the shape
contract, just deterministic search results.
"""
from __future__ import annotations

import types
from typing import Dict, List, Tuple

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools.chunks import get_chunks_for_topic


class _FakeTable:
    def __init__(self, results: List[Dict]):
        self._results = results

    def search(self, _qvec):
        return self

    def limit(self, _k: int):
        return self

    def where(self, _clause: str):
        return self

    def to_list(self) -> List[Dict]:
        return self._results


class _FakeEmbedder:
    def embed(self, queries: List[str]) -> List[List[float]]:
        return [[0.0] * 3 for _ in queries]


class _FakeBiblioDB:
    """Minimal biblio_db stub exposing ``.conn.execute(...).fetchall()``
    for the cited_work_ids query in get_chunks_for_topic(with_cites=True)."""

    def __init__(self, cites_by_paper: Dict[str, List[str]]):
        self._cites = cites_by_paper

    class _Conn:
        def __init__(self, outer):
            self._outer = outer

        def execute(self, _sql: str, params: Tuple):
            paper_hash = params[0]
            rows = [(wid,) for wid in self._outer._cites.get(paper_hash, [])]
            return types.SimpleNamespace(fetchall=lambda: rows)

    @property
    def conn(self):
        return self._Conn(self)


def _make_index(results: List[Dict], biblio_db=None):
    """Build a CorpusIndex-shaped stub with a synthetic topic-searcher
    that returns the given LanceDB row records."""
    embedder = _FakeEmbedder()
    table = _FakeTable(results)
    papers = {"aaaaaaaaaaaa": {"title": "Paper A", "year": 2005}}
    return types.SimpleNamespace(
        papers=papers,
        get_topic_searcher=lambda: (embedder, table),
        biblio_db=biblio_db,
    )


@pytest.fixture
def index():
    """One LanceDB-shape row, mentioning paper aaa, chunk c0 with a
    medium-length body text. Score is a cosine *distance* so lower
    is better."""
    results = [{
        "metadata": {
            "pdf_hash": "aaaaaaaaaaaa",
            "title": "Paper A",
            "year": 2005,
            "chunk_id": "c0",
            "section_class": "description",
            "headings": ["Description"],
        },
        "text": "A" * 600,
        "_distance": 0.12,
    }]
    fake = _make_index(results)
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    yield fake
    mcp_app.set_index(original)


def test_default_with_text_true_includes_text(index):
    out = get_chunks_for_topic("query")
    assert "text" in out[0]
    assert "len_chars" not in out[0]
    assert out[0]["text"] == "A" * 600


def test_with_text_false_drops_text_emits_len_chars(index):
    out = get_chunks_for_topic("query", with_text=False)
    assert "text" not in out[0]
    assert out[0]["len_chars"] == 600


def test_with_text_false_preserves_metadata_fields(index):
    """The whole point of with_text=False is metadata-only scan —
    everything except text must still be there for the LLM to pick
    chunks for drill-down."""
    out = get_chunks_for_topic("query", with_text=False)
    row = out[0]
    for key in ("paper_hash", "paper_title", "paper_year",
                "chunk_id", "section_class", "headings", "score"):
        assert key in row, f"missing {key} in metadata-only row"
    assert row["chunk_id"] == "c0"
    assert row["section_class"] == "description"


def test_with_text_false_token_cut_vs_with_text_true(index):
    """Concretely demonstrates the cost lever. Compare serialized
    JSON length — with_text=False should be dramatically shorter."""
    import json
    with_text = json.dumps(get_chunks_for_topic("query"))
    no_text = json.dumps(get_chunks_for_topic("query", with_text=False))
    # 600-char text should dominate the row size; metadata-only
    # should be a small fraction.
    assert len(no_text) < len(with_text) // 3, (
        f"with_text=False not significantly shorter: "
        f"{len(no_text)} vs {len(with_text)}"
    )


# ---------------------------------------------------------------------------
# with_cites= (#88)
# ---------------------------------------------------------------------------


def test_with_cites_false_omits_cited_work_ids(index):
    out = get_chunks_for_topic("query")
    assert "cited_work_ids" not in out[0]


def test_with_cites_true_attaches_parent_paper_citations(monkeypatch):
    results = [{
        "metadata": {"pdf_hash": "aaaaaaaaaaaa", "title": "Paper A",
                     "year": 2005, "chunk_id": "c0",
                     "section_class": "description", "headings": []},
        "text": "body",
        "_distance": 0.1,
    }]
    biblio = _FakeBiblioDB({"aaaaaaaaaaaa": [
        "corpus:dunn|2005|x", "doi:10.1/abc",
    ]})
    fake = _make_index(results, biblio_db=biblio)
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    try:
        out = get_chunks_for_topic("query", with_cites=True)
        assert out[0]["cited_work_ids"] == ["corpus:dunn|2005|x", "doi:10.1/abc"]
    finally:
        mcp_app.set_index(original)


def test_with_cites_true_empty_when_no_biblio_db(index):
    # The `index` fixture has biblio_db=None — with_cites must degrade
    # to an empty list rather than raising.
    out = get_chunks_for_topic("query", with_cites=True)
    assert out[0]["cited_work_ids"] == []
