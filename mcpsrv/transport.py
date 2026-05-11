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
    """ASGI middleware: short-circuit ``GET /healthz`` with a plain
    200, delegate everything else.

    Mounted *outside* :class:`_BearerAuthASGI` so the probe is reachable
    without a bearer token.  ``/healthz`` reveals only that the process
    is up, which is already inferable from a successful TCP connect,
    and uptime monitors / nginx readiness probes should not need the
    shared secret.
    """

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if (
            scope.get("type") == "http"
            and scope.get("method") == "GET"
            and scope.get("path") == "/healthz"
        ):
            await send({
                "type": "http.response.start",
                "status": 200,
                "headers": [(b"content-type", b"text/plain; charset=utf-8")],
            })
            await send({"type": "http.response.body", "body": b"ok\n"})
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
    allow_unpublishable: bool = False,
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
        figure_app = make_figure_app(idx, allow_unpublishable=allow_unpublishable)
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
    app = _HealthzASGI(app)
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
    allow_unpublishable: bool = False,
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
    figure_app = make_figure_app(idx, allow_unpublishable=allow_unpublishable)
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


