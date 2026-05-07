"""Transport + auth for the MCP server.

* :class:`_BearerAuthASGI` — middleware that gates SSE traffic on a
  shared bearer token (dev_docs/PLAN.md §10).
* :func:`_load_auth_token` — resolves the active token from
  ``--auth-token-file``, ``CORPUS_MCP_TOKEN``, or ``/etc/corpus/token``.
* :func:`_run_sse` — wires the FastMCP SSE app behind the auth
  middleware and serves it on host:port via uvicorn.
"""
from __future__ import annotations

import hmac
import logging
import os
from pathlib import Path
from typing import Optional

from .app import mcp

logger = logging.getLogger(__name__)


class _BearerAuthASGI:
    """ASGI middleware: 401 unless ``Authorization: Bearer <token>`` matches.

    Wraps any ASGI app (FastMCP's ``sse_app()``) at construction time —
    avoids the post-init middleware caveats on finalized Starlette apps.
    Lifespan and websocket scopes pass through unauthenticated; only
    ``http`` requests are gated.

    ``hmac.compare_digest`` makes the comparison constant-time.  For a
    shared-secret setup with ~20 trusted collaborators this is cheap
    insurance against timing attacks.
    """

    def __init__(self, app, token: str):
        self.app = app
        self._token_bytes = token.encode("utf-8")

    async def __call__(self, scope, receive, send):
        if scope.get("type") != "http":
            await self.app(scope, receive, send)
            return
        raw = dict(scope.get("headers", [])).get(b"authorization", b"")
        auth = raw.decode("latin-1", errors="replace")
        if not auth.startswith("Bearer "):
            await _send_401(send, "missing bearer token")
            return
        offered = auth[7:].encode("utf-8")
        if not hmac.compare_digest(offered, self._token_bytes):
            await _send_401(send, "invalid bearer token")
            return
        await self.app(scope, receive, send)


async def _send_401(send, reason: str) -> None:
    body = f"Unauthorized: {reason}\n".encode("utf-8")
    await send({
        "type": "http.response.start",
        "status": 401,
        "headers": [
            (b"content-type", b"text/plain; charset=utf-8"),
            (b"www-authenticate", b'Bearer realm="corpus-mcp"'),
        ],
    })
    await send({"type": "http.response.body", "body": body})


def _load_auth_token(cli_token_file: Optional[Path]) -> Optional[str]:
    """Resolve the bearer-auth token.

    Search order:
      1. ``--auth-token-file <path>`` — read, strip.
      2. ``CORPUS_MCP_TOKEN`` env var.
      3. None (server runs open; main() logs a big warning).

    ``--auth-token`` is deliberately *not* a CLI flag — secrets on the
    command line leak via ``ps``/proc/<pid>/cmdline.
    """
    if cli_token_file is not None:
        return cli_token_file.read_text().strip()
    env = os.environ.get("CORPUS_MCP_TOKEN", "").strip()
    return env or None


def _run_sse(host: str, port: int, token: Optional[str]) -> None:
    """Serve the FastMCP app over SSE on ``host:port``, optionally
    with bearer-token auth.  Requires uvicorn (pulled in transitively
    by the ``mcp`` package for its HTTP transports)."""
    try:
        import uvicorn
    except ImportError as e:
        raise RuntimeError(
            "uvicorn is required for --transport sse. "
            "It should come in via `pip install mcp`; if you're seeing "
            "this, your mcp install may be incomplete."
        ) from e

    app = mcp.sse_app()
    if token:
        app = _BearerAuthASGI(app, token)
        logger.info("Bearer-token auth enabled")
    else:
        logger.warning(
            "*** Running WITHOUT auth: anyone who can reach %s:%d can call "
            "every MCP tool. Set CORPUS_MCP_TOKEN or --auth-token-file "
            "before exposing this beyond localhost. ***",
            host, port,
        )
    logger.info("Serving SSE on http://%s:%d", host, port)
    uvicorn.run(app, host=host, port=port, log_level="info")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


