"""Phase-0 per-tool instrumentation shim (mcpsrv/app.py).

Verifies the dispatch-layer interceptor records call counts, error
counts, and latency without touching FastMCP's schema generation, and
that the args summary never leaks argument *values*. Feeds #90 (run
log) and #91 (/healthz).
"""
from __future__ import annotations

import asyncio

import pytest

from mcpsrv.app import TOOL_STATS, _args_summary, _ToolStats, mcp


def test_tool_stats_math():
    s = _ToolStats()
    s.record("t", 10.0, ok=True)
    s.record("t", 20.0, ok=False)
    snap = s.snapshot()
    assert snap["t"] == {"calls": 2, "errors": 1, "avg_latency_ms": 15.0}
    assert s.total_calls == 2
    assert s.total_errors == 1


def test_args_summary_is_keys_only():
    # Sorted key names only — never values, so logs can't leak secrets.
    assert _args_summary({"b": 2, "a": 1}) == "a,b"
    assert _args_summary(None) == ""
    assert _args_summary("garbage") == ""
    assert "hunter2" not in _args_summary({"token": "hunter2"})


def test_call_tool_is_instrumented():
    # The interceptor must be installed on the live FastMCP instance,
    # and installation must be idempotent under re-import.
    assert getattr(mcp._tool_manager, "_corpus_instrumented", False) is True


def test_dispatch_error_is_recorded_and_reraised():
    name = "does_not_exist_xyz"
    before_calls = TOOL_STATS.calls.get(name, 0)
    before_errors = TOOL_STATS.errors.get(name, 0)
    with pytest.raises(Exception):
        asyncio.run(mcp._tool_manager.call_tool(name, {}))
    assert TOOL_STATS.calls[name] == before_calls + 1
    assert TOOL_STATS.errors[name] == before_errors + 1
