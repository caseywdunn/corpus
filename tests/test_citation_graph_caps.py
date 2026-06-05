"""Breadth + total-edge caps on get_citation_graph (#87).

A hub work (cited by hundreds) must not return a runaway graph.
``max_edges_per_node`` caps per-node fan-out — ranking survivors by
``cited_by_count`` so the most-cited edges are kept — and
``max_total_edges`` caps the whole walk. ``truncated`` reports whether
either fired. Defaults are generous, so the change is additive.
"""
from __future__ import annotations

import types

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools.bibliography import get_citation_graph


class _FakeBiblio:
    def __init__(self, citing_edges, counts):
        self._citing = citing_edges      # wid -> [neighbour work_ids]
        self._counts = counts            # wid -> cited_by_count

    def get_work(self, wid):
        return {"work_id": wid, "title": wid} if wid in self._counts or True else None

    def get_work_by_corpus_hash(self, _h):
        return None

    def get_authors(self, _wid):
        return []

    def citation_count(self, wid):
        return self._counts.get(wid, 0)

    def citing(self, wid):
        return [{"work_id": n} for n in self._citing.get(wid, [])]

    def cited_by(self, _wid):
        return []


def _set_index(biblio):
    fake = types.SimpleNamespace(biblio_db=biblio)
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    return original


@pytest.fixture(autouse=True)
def _restore():
    original = mcp_app._INDEX
    yield
    mcp_app.set_index(original)


def test_no_truncation_under_caps():
    biblio = _FakeBiblio(
        {"root": ["a", "b", "c"]},
        {"a": 10, "b": 5, "c": 1},
    )
    _set_index(biblio)
    out = get_citation_graph(work_id="root", direction="citing", depth=1,
                             max_edges_per_node=10)
    assert out["truncated"] is False
    assert {e["work_id"] for e in out["citing"]} == {"a", "b", "c"}


def test_per_node_cap_keeps_most_cited():
    biblio = _FakeBiblio(
        {"root": ["a", "b", "c"]},
        {"a": 10, "b": 5, "c": 1},
    )
    _set_index(biblio)
    out = get_citation_graph(work_id="root", direction="citing", depth=1,
                             max_edges_per_node=2)
    assert out["truncated"] is True
    # ranked by cited_by_count desc → a (10), b (5); c (1) dropped.
    assert [e["work_id"] for e in out["citing"]] == ["a", "b"]


def test_total_edge_cap_stops_walk():
    biblio = _FakeBiblio(
        {"root": ["a", "b", "c", "d", "e"]},
        {k: 1 for k in "abcde"},
    )
    _set_index(biblio)
    out = get_citation_graph(work_id="root", direction="citing", depth=1,
                             max_edges_per_node=50, max_total_edges=3)
    assert out["truncated"] is True
    assert len(out["citing"]) == 3
