"""Transport + auth for the MCP server.

* :class:`_BearerAuthASGI` — middleware that gates SSE traffic on a
  shared bearer token (dev_docs/PLAN.md §10).
* :func:`_load_auth_token` — resolves the active token from
  ``--auth-token-file``, ``CORPUS_MCP_TOKEN``, or ``/etc/corpus/token``.
* :func:`_run_sse` — wires the FastMCP SSE app behind the auth
  middleware and serves it on host:port via uvicorn, with the figure
  HTTP route from :mod:`.figure_http` mounted alongside.
* :func:`start_stdio_figure_server` — spawn a daemon-thread uvicorn
  serving only the figure HTTP route on a loopback port, so a stdio
  MCP server can still expose figure bytes via ``curl -o`` (#69).
"""
from __future__ import annotations

import hmac
import json
import logging
import os
import secrets
import socket
import threading
import time
from pathlib import Path
from typing import Optional, Tuple

from .app import mcp
from .figure_http import make_figure_app

logger = logging.getLogger(__name__)


class _HealthzASGI:
    """ASGI middleware: short-circuit ``GET /healthz`` with a JSON
    capability report, delegate everything else.

    Mounted *outside* :class:`_BearerAuthASGI` so the probe is reachable
    without a bearer token — uptime monitors / ALB readiness probes
    should not need the shared secret. The body reports per-capability
    state (#91); detail strings are kept filesystem-path-free so this
    stays safe to expose unauthenticated.

    Status code is **200** when every capability is ``ok`` or ``absent``
    and **503** when any capability is ``degraded`` — so a load balancer
    pulls a server whose backing index broke (e.g. an embedding-model
    dim mismatch) out of rotation instead of serving silently-empty
    semantic searches.
    """

    def __init__(self, app, idx=None):
        self.app = app
        self.idx = idx

    def _report(self) -> Tuple[int, dict]:
        from .app import TOOL_STATS
        if self.idx is None:
            # No index wired (shouldn't happen in SSE serve) — process
            # is up but can't report capabilities.
            return 200, {"status": "ok", "capabilities": {}, "tools": {}}
        caps = self.idx.capabilities()
        healthy = not any(c["state"] == "degraded" for c in caps.values())
        body = {
            "status": "ok" if healthy else "degraded",
            "capabilities": caps,
            "tools": {
                "total_calls": TOOL_STATS.total_calls,
                "total_errors": TOOL_STATS.total_errors,
            },
        }
        return (200 if healthy else 503), body

    async def __call__(self, scope, receive, send):
        if (
            scope.get("type") == "http"
            and scope.get("method") == "GET"
            and scope.get("path") == "/healthz"
        ):
            status, body = self._report()
            payload = (json.dumps(body) + "\n").encode("utf-8")
            await send({
                "type": "http.response.start",
                "status": status,
                "headers": [(b"content-type", b"application/json; charset=utf-8")],
            })
            await send({"type": "http.response.body", "body": payload})
            return
        await self.app(scope, receive, send)


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
      1. ``--auth-token-file <path>`` — read, strip; raise on empty.
      2. ``CORPUS_MCP_TOKEN`` env var.
      3. None (server runs open; main() logs a big warning).

    Raises ``ValueError`` when the token file exists but resolves to
    an empty string. An empty file is ambiguous — the operator clearly
    intended to set a token but ended up with no auth at all, which
    would otherwise silently open the SSE server.

    ``--auth-token`` is deliberately *not* a CLI flag — secrets on the
    command line leak via ``ps``/proc/<pid>/cmdline.
    """
    if cli_token_file is not None:
        token = cli_token_file.read_text().strip()
        if not token:
            raise ValueError(
                f"--auth-token-file {cli_token_file} is empty (no non-"
                "whitespace bytes). Refusing to start: an empty token "
                "file is almost certainly a configuration mistake. Put "
                "a secret in the file, or unset --auth-token-file to "
                "run open intentionally."
            )
        return token
    env = os.environ.get("CORPUS_MCP_TOKEN", "").strip()
    return env or None


class _RouteMuxASGI:
    """ASGI multiplexer: route ``/figures/...`` to the figure HTTP app,
    everything else to the SSE app. Both apps are pre-wrapped with the
    auth middleware so this layer is purely path-based dispatch.
    """

    def __init__(self, sse_app, figure_app):
        self.sse_app = sse_app
        self.figure_app = figure_app

    async def __call__(self, scope, receive, send):
        if (
            scope.get("type") == "http"
            and scope.get("path", "").startswith("/figures/")
        ):
            await self.figure_app(scope, receive, send)
            return
        await self.sse_app(scope, receive, send)


def _run_sse(
    host: str,
    port: int,
    token: Optional[str],
    idx=None,
    default_profile: Optional[str] = None,
) -> None:
    """Serve the FastMCP app over SSE on ``host:port``, optionally
    with bearer-token auth.  Requires uvicorn (pulled in transitively
    by the ``mcp`` package for its HTTP transports).

    ``idx`` is the live :class:`CorpusIndex`; if supplied, the figure
    HTTP route from :mod:`.figure_http` is mounted alongside ``/sse``
    so ``get_figure_url`` can hand operators a working URL (#69).
    """
    try:
        import uvicorn
    except ImportError as e:
        raise RuntimeError(
            "uvicorn is required for --transport sse. "
            "It should come in via `pip install mcp`; if you're seeing "
            "this, your mcp install may be incomplete."
        ) from e

    sse_app = mcp.sse_app()
    if idx is not None:
        figure_app = make_figure_app(idx, default_profile=default_profile)
        app = _RouteMuxASGI(sse_app, figure_app)
    else:
        app = sse_app
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
    # Healthz wraps the auth layer so probes don't need the bearer token.
    app = _HealthzASGI(app, idx=idx)
    logger.info("Serving SSE on http://%s:%d (healthz: GET /healthz)", host, port)
    if idx is not None:
        logger.info(
            "Figure HTTP route mounted at GET /figures/<paper_hash>/<figure_id> "
            "(used by get_figure_url; same bearer-auth as /sse)"
        )
    try:
        uvicorn.run(app, host=host, port=port, log_level="info")
    except OSError as e:
        # #47: surface common startup failures with actionable messages
        # rather than letting uvicorn's traceback bubble up.
        if e.errno in (48, 98):  # EADDRINUSE on macOS / Linux
            logger.error(
                "port %s:%d already in use. Check `lsof -i :%d` (or "
                "`ss -ltnp 'sport = :%d'`) and either stop the conflicting "
                "process or pass a different --port. Run `corpus serve --check` "
                "to pre-flight before binding.",
                host, port, port, port,
            )
            raise SystemExit(3) from e
        raise


def start_stdio_figure_server(
    idx,
    default_profile: Optional[str] = None,
    host: str = "127.0.0.1",
) -> Tuple[str, int, str]:
    """Spin up a daemon-thread uvicorn serving only the figure HTTP
    route (#69), on a loopback ephemeral port. Returns the
    ``(host, port, token)`` tuple so the caller can stash it on the
    index for ``get_figure_url`` to build URLs.

    Used by stdio MCP mode, where there's no existing HTTP endpoint.
    The token is a one-shot 32-byte hex string generated here; the
    parent stdio MCP session is the only thing that knows it.
    """
    try:
        import uvicorn
    except ImportError as e:
        raise RuntimeError("uvicorn required for figure HTTP route") from e

    # Pre-bind a socket so we know the port before uvicorn starts; pass
    # it to uvicorn via the ``fd=`` mechanism so there's no race window.
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.bind((host, 0))
    port = sock.getsockname()[1]
    sock.set_inheritable(True)

    token = secrets.token_hex(32)
    figure_app = make_figure_app(idx, default_profile=default_profile)
    app = _BearerAuthASGI(figure_app, token)

    config = uvicorn.Config(
        app, fd=sock.fileno(),
        log_level="warning",   # keep uvicorn quiet — main logger handles status
        access_log=False,
    )
    server = uvicorn.Server(config)

    def _run():
        try:
            server.run()
        except Exception as e:
            logger.error("figure HTTP server crashed: %s", e)

    thread = threading.Thread(target=_run, name="corpus-figure-http", daemon=True)
    thread.start()

    # Wait briefly for uvicorn to flip server.started before returning,
    # so the first get_figure_url call doesn't race the bind.
    deadline = time.monotonic() + 5.0
    while not server.started and time.monotonic() < deadline:
        time.sleep(0.05)
    if not server.started:
        logger.warning(
            "figure HTTP server did not signal `started` within 5s; "
            "subsequent figure URL fetches may briefly 502 until it's up"
        )

    logger.info(
        "Figure HTTP route mounted at http://%s:%d/figures/<paper_hash>/<figure_id> "
        "(stdio side-car; bearer-token gated)",
        host, port,
    )
    return host, port, token


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


