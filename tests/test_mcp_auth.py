"""Unit tests for the bearer-auth ASGI middleware in mcp_server.py.

The middleware is the only piece of the HTTP transport that enforces
access control, so test it directly at the ASGI protocol level.  We
don't spin up a real server — just feed the middleware hand-rolled
scope/receive/send callables.
"""

from __future__ import annotations

import asyncio
import importlib.util
import sys
from pathlib import Path
from typing import List

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

# Load mcp_server.py as a module — can't import normally because it
# pulls in the FastMCP dep which isn't always present in a dev env.
# Skip gracefully if that dep is missing; the middleware itself is
# pure asyncio code that doesn't need mcp installed to be tested.
_spec = importlib.util.spec_from_file_location(
    "mcp_server", REPO_ROOT / "mcp_server.py"
)
try:
    _mcp_server = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_mcp_server)
    _HAS_MCP = True
except ModuleNotFoundError:
    _HAS_MCP = False


pytestmark = pytest.mark.skipif(
    not _HAS_MCP,
    reason="mcp_server imports the mcp sdk; skip when not installed",
)


# ── Capture-only ASGI app ───────────────────────────────────────────

class _CaptureApp:
    """An ASGI app that records whether it got invoked, so tests can
    tell whether the middleware passed the request through or shorted."""

    def __init__(self):
        self.called = False
        self.scope = None

    async def __call__(self, scope, receive, send):
        self.called = True
        self.scope = scope
        await send({
            "type": "http.response.start",
            "status": 200,
            "headers": [(b"content-type", b"text/plain")],
        })
        await send({"type": "http.response.body", "body": b"ok"})


async def _send_collector():
    """Helper that records all send() messages for inspection."""
    messages: List[dict] = []

    async def send(msg):
        messages.append(msg)

    return messages, send


def _run(coro):
    return asyncio.new_event_loop().run_until_complete(coro)


# ── Tests ───────────────────────────────────────────────────────────

def test_missing_auth_header_rejects_with_401():
    inner = _CaptureApp()
    mw = _mcp_server._BearerAuthASGI(inner, token="secret")

    async def run():
        scope = {"type": "http", "headers": []}
        messages = []
        async def send(m): messages.append(m)
        async def receive(): return {"type": "http.request", "body": b""}
        await mw(scope, receive, send)
        return messages

    messages = _run(run())
    assert inner.called is False
    assert messages[0]["status"] == 401
    # body describes the reason
    body = messages[1]["body"].decode()
    assert "missing bearer" in body


def test_wrong_token_rejects_with_401():
    inner = _CaptureApp()
    mw = _mcp_server._BearerAuthASGI(inner, token="secret")

    async def run():
        scope = {"type": "http", "headers": [(b"authorization", b"Bearer wrong")]}
        messages = []
        async def send(m): messages.append(m)
        async def receive(): return {"type": "http.request", "body": b""}
        await mw(scope, receive, send)
        return messages

    messages = _run(run())
    assert inner.called is False
    assert messages[0]["status"] == 401
    body = messages[1]["body"].decode()
    assert "invalid bearer" in body


def test_correct_token_passes_through():
    inner = _CaptureApp()
    mw = _mcp_server._BearerAuthASGI(inner, token="secret")

    async def run():
        scope = {
            "type": "http",
            "headers": [(b"authorization", b"Bearer secret")],
        }
        messages = []
        async def send(m): messages.append(m)
        async def receive(): return {"type": "http.request", "body": b""}
        await mw(scope, receive, send)
        return messages

    messages = _run(run())
    assert inner.called is True
    assert messages[0]["status"] == 200


def test_non_http_scope_passes_through_unauthenticated():
    """Lifespan + websocket scopes bypass auth — FastMCP relies on
    lifespan startup events and we shouldn't block those with token
    checks the client has no way to supply.
    """
    inner = _CaptureApp()
    mw = _mcp_server._BearerAuthASGI(inner, token="secret")

    async def run():
        scope = {"type": "lifespan"}
        async def send(m): pass
        async def receive(): return {"type": "lifespan.startup"}
        await mw(scope, receive, send)

    _run(run())
    assert inner.called is True


def test_token_comparison_is_constant_time():
    """Regression: use hmac.compare_digest, not plain ==.  Guard
    against a well-meaning refactor re-introducing the timing leak.
    """
    import inspect
    src = inspect.getsource(_mcp_server._BearerAuthASGI.__call__)
    assert "compare_digest" in src, (
        "token comparison must use hmac.compare_digest for constant-time "
        "behavior — do not switch to == without thinking about timing attacks"
    )
