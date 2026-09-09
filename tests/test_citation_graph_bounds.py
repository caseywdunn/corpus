"""`get_citation_graph` reports what it returned, not just whether it cut
(#166).

`truncated` reported only whether *this tool* dropped edges, so it read
`false` — accurately — while the client failed to deliver the payload. A
caller checking it alone concluded it had received everything.

Measured on the 1,775-document reference corpus: **72 works return
150-500 edges with `truncated: false`, up to 145 kB**; 14 exceed the
default 500-edge cap and are truncated honestly; and the largest
bibliography is 2,277 edges / 625 kB uncapped. The reported client
failure was at ~55 kB, so this is not one hub paper — it is 72 of them,
at up to 2.6x that size.

The server cannot know an individual client's transport limit, which is
why the counts matter more than the ceiling's value. Input schema is
untouched: the 1.0 freeze pins tool inputs, and these are response
fields.
"""
from __future__ import annotations

import types

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools import bibliography as biblio_tools
from mcpsrv.tools.bibliography import get_citation_graph


class _FakeBiblio:
    def __init__(self, cited_by=(), citing=()):
        self._cited_by = list(cited_by)
        self._citing = list(citing)

    def get_work(self, work_id):
        return {"work_id": work_id, "title": "Root", "year": 1965}

    def get_authors(self, _work_id):
        return [{"surname": "Totton"}]

    def citation_count(self, _work_id):
        return 1

    def cited_by(self, work_id):
        return [dict(r) for r in self._cited_by] if work_id == "root" else []

    def citing(self, work_id):
        return [dict(r) for r in self._citing] if work_id == "root" else []


def _edges(n, prefix="w", pad=""):
    return [{"work_id": f"{prefix}{i}", "title": f"Cited work {i}{pad}",
             "year": 1900 + (i % 100)} for i in range(n)]


@pytest.fixture
def served(request):
    cited_by, citing = getattr(request, "param", (_edges(10), []))
    original = mcp_app._INDEX
    mcp_app.set_index(types.SimpleNamespace(
        biblio_db=_FakeBiblio(cited_by, citing)))
    yield
    mcp_app.set_index(original)


# --- the honest counts --------------------------------------------------


@pytest.mark.parametrize("served", [(_edges(210), [])], indirect=True)
def test_a_complete_answer_reports_that_it_is_complete(served):
    """The reported case: 210 edges, nothing dropped. `truncated: false`
    was true and uninformative; the counts and the size are what let a
    caller near its own limit see the payload coming."""
    out = get_citation_graph(work_id="root", direction="cited_by")
    assert out["truncated"] is False
    assert "truncated_reason" not in out
    assert out["edges_available"] == {"cited_by": 210}
    assert out["edges_returned"] == {"cited_by": 210}
    assert out["response_bytes"] > 0


@pytest.mark.parametrize("served", [(_edges(2277), [])], indirect=True)
def test_a_capped_answer_says_how_much_it_did_not_return(served):
    """`edges_available` is the number a reader of a depth-1 graph is
    actually asking about: how long is this paper's bibliography."""
    out = get_citation_graph(work_id="root", direction="cited_by")
    assert out["truncated"] is True
    assert out["edges_available"] == {"cited_by": 2277}
    assert out["edges_returned"]["cited_by"] == 500
    assert "max_edges_per_node" in out["truncated_reason"]


@pytest.mark.parametrize("served", [(_edges(20), _edges(5, "c"))], indirect=True)
def test_both_directions_are_counted_separately(served):
    out = get_citation_graph(work_id="root", direction="both")
    assert out["edges_available"] == {"citing": 5, "cited_by": 20}
    assert out["edges_returned"] == {"citing": 5, "cited_by": 20}


@pytest.mark.parametrize("served", [(_edges(20), _edges(5, "c"))], indirect=True)
def test_a_single_direction_reports_only_that_direction(served):
    out = get_citation_graph(work_id="root", direction="citing")
    assert set(out["edges_available"]) == {"citing"}


# --- the ceiling --------------------------------------------------------


@pytest.mark.parametrize("served", [(_edges(2000, pad="x" * 400), [])],
                         indirect=True)
def test_the_byte_ceiling_bounds_the_payload_and_admits_it(served, monkeypatch):
    monkeypatch.setattr(biblio_tools, "CITATION_GRAPH_MAX_BYTES", 20_000)
    out = get_citation_graph(work_id="root", direction="cited_by",
                             max_edges_per_node=2000, max_total_edges=2000)
    assert out["response_bytes"] <= 20_000
    assert out["truncated"] is True
    assert "response_bytes" in out["truncated_reason"]
    # And it still says how much there was.
    assert out["edges_available"] == {"cited_by": 2000}
    assert out["edges_returned"]["cited_by"] < 2000


@pytest.mark.parametrize("served",
                         [(_edges(400, pad="x" * 200), _edges(3, "c"))],
                         indirect=True)
def test_trimming_for_size_leaves_both_directions_present(served, monkeypatch):
    """Trim the larger side first, or a `direction="both"` call can come
    back with one side silently empty — which looks like "this work cites
    nothing"."""
    monkeypatch.setattr(biblio_tools, "CITATION_GRAPH_MAX_BYTES", 15_000)
    out = get_citation_graph(work_id="root", direction="both",
                             max_edges_per_node=400, max_total_edges=400)
    assert out["response_bytes"] <= 15_000
    assert out["edges_returned"]["citing"] == 3
    assert 0 < out["edges_returned"]["cited_by"] < 400


@pytest.mark.parametrize("served", [(_edges(3), [])], indirect=True)
def test_a_small_graph_is_untouched(served):
    """No truncation noise on the common case."""
    out = get_citation_graph(work_id="root", direction="cited_by")
    assert out["truncated"] is False
    assert len(out["cited_by"]) == 3
    assert "truncated_reason" not in out


def test_the_ceiling_counts_its_own_reporting(monkeypatch):
    """A ceiling that does not include the fields it adds is not a
    ceiling: the first cut trimmed to the limit and then stamped
    `edges_available`, `edges_returned`, `truncated_reason` and
    `response_bytes` on top, coming back 29 bytes over."""
    import types

    from mcpsrv import app as mcp_app
    original = mcp_app._INDEX
    mcp_app.set_index(types.SimpleNamespace(
        biblio_db=_FakeBiblio(_edges(300, pad="y" * 300), [])))
    try:
        monkeypatch.setattr(biblio_tools, "CITATION_GRAPH_MAX_BYTES", 12_000)
        out = get_citation_graph(work_id="root", direction="cited_by",
                                 max_edges_per_node=300, max_total_edges=300)
        import json
        real = len(json.dumps(out, default=str).encode("utf-8"))
        assert out["response_bytes"] == real, "reported size must be the real one"
        assert real <= 12_000
    finally:
        mcp_app.set_index(original)


def test_the_ceiling_is_operator_overridable():
    """A client with a stricter transport limit needs a lever, and the
    server cannot infer one."""
    import inspect
    src = inspect.getsource(biblio_tools)
    assert "CORPUS_CITATION_GRAPH_MAX_BYTES" in src


# --- the frozen input surface ------------------------------------------


def test_the_tool_signature_is_unchanged():
    """1.0 freezes tool *inputs*; all of the above is response shape, so
    no client's call has to change."""
    import inspect
    sig = inspect.signature(get_citation_graph.fn
                            if hasattr(get_citation_graph, "fn")
                            else get_citation_graph)
    assert list(sig.parameters) == [
        "work_id", "paper_hash", "direction", "depth",
        "max_edges_per_node", "max_total_edges",
    ]
