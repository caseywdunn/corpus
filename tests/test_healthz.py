"""Health-check + degraded-capability tests (#91).

Three layers:

* :meth:`CorpusIndex.capabilities` classifies each backing index as
  ``ok`` / ``absent`` / ``degraded`` — the single source of truth.
* ``_HealthzASGI`` turns that into a JSON body + a 200/503 status so a
  load balancer pulls a degraded server out of rotation.
* ``get_chunks_for_topic`` raises a hard MCP error (not an empty
  ``[{error: …}]`` row) when its backing semantic index is degraded.
"""
from __future__ import annotations

import asyncio
import json
import types
from pathlib import Path

import pytest

from mcpsrv.indexes import CorpusIndex
from mcpsrv.transport import _HealthzASGI


# ---------------------------------------------------------------------------
# capabilities() classification
# ---------------------------------------------------------------------------


def _bare_index(tmp_path: Path) -> CorpusIndex:
    """A CorpusIndex over an empty output dir — no documents, no
    vector store, no structured DBs. Avoids load() (which needs a
    documents/ tree) by constructing directly."""
    return CorpusIndex(tmp_path)


def test_structured_only_corpus_is_healthy(tmp_path):
    """No vector store + no structured DBs is a valid deploy: every
    optional capability is ``absent`` (not ``degraded``), so the server
    is healthy. The one exception is ``papers`` — see below."""
    idx = _bare_index(tmp_path)
    idx.papers = {"deadbeef0000": {"title": "x"}}  # one paper → core ok
    caps = idx.capabilities()

    assert caps["topic_search"]["state"] == "absent"
    assert caps["taxonomy"]["state"] == "absent"
    assert caps["bibliography"]["state"] == "absent"
    assert caps["taxon_mentions"]["state"] == "absent"
    assert caps["papers"]["state"] == "ok"
    assert idx.is_healthy() is True


def test_empty_paper_index_is_degraded(tmp_path):
    """An empty paper index is the v0.5 zero-yield failure mode — a
    broken bundle, not a valid structured-only corpus. Degraded."""
    idx = _bare_index(tmp_path)
    assert idx.papers == {}
    caps = idx.capabilities()
    assert caps["papers"]["state"] == "degraded"
    assert idx.is_healthy() is False


def test_structured_dbs_report_ok_when_loaded(tmp_path):
    idx = _bare_index(tmp_path)
    idx.papers = {"deadbeef0000": {"title": "x"}}
    idx.taxonomy_db = object()
    idx.biblio_db = object()
    idx.taxon_mention_db = object()
    caps = idx.capabilities()
    assert caps["taxonomy"]["state"] == "ok"
    assert caps["bibliography"]["state"] == "ok"
    assert caps["taxon_mentions"]["state"] == "ok"


def test_recorded_degraded_status_takes_precedence(tmp_path):
    """Once a real query records a degraded reason (e.g. an embedding-
    model dim mismatch only detectable by loading the model),
    capabilities() must report it rather than re-running the cheap
    structural probe that can't see it."""
    idx = _bare_index(tmp_path)
    idx.papers = {"deadbeef0000": {"title": "x"}}
    idx._topic_search_status = (
        "degraded: embedding-model dim mismatch (model emits 768, index expects 1024)"
    )
    caps = idx.capabilities()
    assert caps["topic_search"]["state"] == "degraded"
    assert "dim mismatch" in caps["topic_search"]["detail"]
    # detail is path-free (safe for unauthenticated /healthz)
    assert "/" not in caps["topic_search"]["detail"]
    assert idx.is_healthy() is False


def test_missing_vector_store_probes_absent_without_loading_embedder(tmp_path):
    """The cheap structural probe must not force the ~600 MB embedder
    load: a missing vector_db/ resolves to ``absent`` with the embedder
    untouched."""
    idx = _bare_index(tmp_path)
    idx.papers = {"deadbeef0000": {"title": "x"}}
    assert idx.capabilities()["topic_search"]["state"] == "absent"
    assert idx._embedder is None  # never loaded


# ---------------------------------------------------------------------------
# _HealthzASGI status code + JSON body
# ---------------------------------------------------------------------------


def _probe_healthz(idx) -> tuple[int, dict]:
    """Drive _HealthzASGI through a GET /healthz scope, capturing the
    response status + decoded JSON body."""
    app = _HealthzASGI(app=_unreachable_app, idx=idx)
    sent: list[dict] = []

    async def send(msg):
        sent.append(msg)

    async def receive():  # pragma: no cover - not used by the healthz path
        return {"type": "http.request"}

    scope = {"type": "http", "method": "GET", "path": "/healthz"}
    asyncio.run(app(scope, receive, send))

    start = next(m for m in sent if m["type"] == "http.response.start")
    body_msg = next(m for m in sent if m["type"] == "http.response.body")
    return start["status"], json.loads(body_msg["body"])


async def _unreachable_app(scope, receive, send):  # pragma: no cover
    raise AssertionError("/healthz must short-circuit before delegating")


def test_healthz_returns_200_when_healthy(tmp_path):
    idx = _bare_index(tmp_path)
    idx.papers = {"deadbeef0000": {"title": "x"}}
    status, body = _probe_healthz(idx)
    assert status == 200
    assert body["status"] == "ok"
    assert body["capabilities"]["papers"]["state"] == "ok"
    assert "tools" in body and "total_calls" in body["tools"]


def test_healthz_returns_503_when_degraded(tmp_path):
    idx = _bare_index(tmp_path)
    idx.papers = {"deadbeef0000": {"title": "x"}}
    idx._topic_search_status = "degraded: cannot open LanceDB index"
    status, body = _probe_healthz(idx)
    assert status == 503
    assert body["status"] == "degraded"
    assert body["capabilities"]["topic_search"]["state"] == "degraded"


def test_healthz_delegates_non_probe_paths(tmp_path):
    """A non-/healthz request must fall through to the wrapped app."""
    idx = _bare_index(tmp_path)
    delegated = {"hit": False}

    async def downstream(scope, receive, send):
        delegated["hit"] = True

    app = _HealthzASGI(app=downstream, idx=idx)
    asyncio.run(app({"type": "http", "method": "GET", "path": "/sse"},
                    None, None))
    assert delegated["hit"] is True


# ---------------------------------------------------------------------------
# get_chunks_for_topic: hard error on degraded, friendly row on absent
# ---------------------------------------------------------------------------


def _inject(idx):
    from mcpsrv import app as mcp_app
    original = mcp_app._INDEX
    mcp_app.set_index(idx)
    return original


def test_topic_search_absent_returns_build_index_row(tmp_path):
    """Structured-only corpus: a clear, non-raising guidance row."""
    from mcpsrv import app as mcp_app
    from mcpsrv.tools.chunks import get_chunks_for_topic

    idx = _bare_index(tmp_path)
    idx.papers = {"deadbeef0000": {"title": "x"}}
    original = _inject(idx)
    try:
        out = get_chunks_for_topic(query="anything")
        assert isinstance(out, list)
        assert "no LanceDB index available" in out[0]["error"]
    finally:
        mcp_app.set_index(original)


def test_topic_search_degraded_raises_hard_error(tmp_path):
    """Operational fault: refuse to serve with a raised error rather
    than returning empty rows that read as 'no matches'."""
    from mcpsrv import app as mcp_app
    from mcpsrv.tools.chunks import get_chunks_for_topic

    idx = _bare_index(tmp_path)
    idx.papers = {"deadbeef0000": {"title": "x"}}

    # Mimic a real degraded probe: get_topic_searcher records the
    # reason and returns (None, None). (Setting the status directly
    # wouldn't survive, since get_chunks_for_topic re-probes first.)
    def _degraded_searcher():
        idx._topic_search_status = "degraded: embedding backend unavailable"
        return None, None

    idx.get_topic_searcher = _degraded_searcher
    original = _inject(idx)
    try:
        with pytest.raises(RuntimeError, match="degraded and refusing to serve"):
            get_chunks_for_topic(query="anything")
    finally:
        mcp_app.set_index(original)
