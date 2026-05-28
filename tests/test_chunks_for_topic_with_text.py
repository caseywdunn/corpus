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


def _make_index(results: List[Dict]):
    """Build a CorpusIndex-shaped stub with a synthetic topic-searcher
    that returns the given LanceDB row records."""
    embedder = _FakeEmbedder()
    table = _FakeTable(results)
    papers = {"aaaaaaaaaaaa": {"title": "Paper A", "year": 2005}}
    return types.SimpleNamespace(
        papers=papers,
        get_topic_searcher=lambda: (embedder, table),
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
