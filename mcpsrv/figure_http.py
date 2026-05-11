"""Out-of-band HTTP route for serving figure PNG bytes (#69).

Background: MCP tool results carry image bytes through the SDK's
``Image()`` content type, which clients render inline for the human
reader but do *not* expose to the model as raw/base64 bytes. So the
model has no way to write a figure to disk via the standard
``Write``/``Bash`` tools — a regression from v0.1, where
``get_figure``'s ``image_path`` was an absolute filesystem path that
local stdio clients could ``Read`` directly. v0.2 scrubbed the
absolute path (sensible for SSE deploys, where server-side paths
aren't reachable from a remote client), but left local stdio without
a replacement.

This module provides a transport-agnostic alternative: an HTTP route
the server mounts so the model can ``curl -o /tmp/fig.png <url>``
via Bash. The bytes flow over HTTP outside the MCP JSON-RPC channel,
so they don't burn context tokens — the tool result is just a short
URL string (~50 tokens regardless of figure size). Same shape on
local stdio (loopback HTTP server bound to 127.0.0.1) and SSE/AWS
(uvicorn already running for the SSE endpoint).

Honored gates:

* ``_BearerAuthASGI`` (same middleware that guards ``/sse``).
* ``#51`` publishable check (refuses figures the parent work isn't
  licensed-cleared for, unless the server was started with
  ``--allow-unpublishable``).
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# Same constraint as the MCP tool: only allow short hex hashes
# (12 chars, matching ``short_hash``) and figure ids matching the
# usual ``<paper_hash>__fig-NN`` etc. shape. Defense in depth — the
# Path() construction below also rejects ``..`` traversal.
_HASH_RE = re.compile(r"^[a-f0-9]{12}$")
_FIGURE_ID_RE = re.compile(r"^[A-Za-z0-9._\-]+$")


def make_figure_app(idx, allow_unpublishable: bool = False):
    """Build an ASGI handler that serves
    ``GET /figures/<paper_hash>/<figure_id>?label=...``.

    ``idx`` is the live :class:`mcpsrv.indexes.CorpusIndex` so we can
    resolve ``paper_hash → hash_dir`` and probe the bibliographic
    authority for license fields.
    """
    # Lazy import — keeps the load-time cost off `import mcpsrv` for
    # non-serve uses (tests, bundle distillation).
    from .tools.figures import _license_metadata_for_paper

    async def app(scope, receive, send):
        if scope.get("type") != "http":
            await _send_text(send, 404, "not found")
            return
        method = scope.get("method")
        if method not in ("GET", "HEAD"):
            await _send_text(send, 405, f"method {method} not allowed")
            return

        path = scope.get("path", "")
        if not path.startswith("/figures/"):
            await _send_text(send, 404, "not found")
            return
        tail = path[len("/figures/"):]
        parts = tail.split("/")
        if len(parts) != 2:
            await _send_text(send, 400, "expected /figures/<paper_hash>/<figure_id>")
            return
        paper_hash, figure_id = parts
        if not _HASH_RE.match(paper_hash):
            await _send_text(send, 400, "malformed paper_hash")
            return
        if not _FIGURE_ID_RE.match(figure_id):
            await _send_text(send, 400, "malformed figure_id")
            return

        # Optional ``?label=<panel>`` to fetch a sub-panel crop.
        label: Optional[str] = None
        qs = scope.get("query_string", b"").decode("latin-1", errors="replace")
        if qs:
            for piece in qs.split("&"):
                if piece.startswith("label="):
                    label = piece[len("label="):]
                    if not _FIGURE_ID_RE.match(label):
                        await _send_text(send, 400, "malformed label")
                        return

        p = idx.papers.get(paper_hash)
        if not p:
            await _send_text(send, 404, f"no such paper_hash: {paper_hash}")
            return

        # Publishable gate (#51) — mirror the get_figure_image policy.
        if not allow_unpublishable:
            lic = _license_metadata_for_paper(paper_hash)
            if not lic.get("publishable"):
                body = {
                    "error": "figure not publishable",
                    "license": lic.get("license") or "unknown",
                    "license_source": lic.get("license_source"),
                    "hint": (
                        "restart the server with --allow-unpublishable for "
                        "local rights-holder cases, or read get_figure() for "
                        "the raw license fields if your jurisdiction differs"
                    ),
                }
                await _send_json(send, 403, body)
                return

        hash_dir = Path(p["hash_dir"])
        figs_path = hash_dir / "figures.json"
        try:
            figs = json.loads(figs_path.read_text(encoding="utf-8")) if figs_path.exists() else {}
        except Exception as e:
            await _send_text(send, 500, f"failed to load figures.json: {e}")
            return
        fig = next(
            (f for f in figs.get("figures", []) or [] if f.get("figure_id") == figure_id),
            None,
        )
        if fig is None:
            await _send_text(send, 404, f"no such figure_id {figure_id!r}")
            return

        whole = hash_dir / "figures" / (fig.get("filename") or "")
        # Defense in depth: confirm the resolved file is inside the
        # corpuscle's figures dir before sending it. Filters out any
        # traversal that snuck past the regex above.
        try:
            whole = whole.resolve(strict=True)
            figures_root = (hash_dir / "figures").resolve(strict=True)
        except FileNotFoundError:
            await _send_text(send, 404, "figure file missing on disk")
            return
        if figures_root not in whole.parents and whole != figures_root:
            await _send_text(send, 400, "path traversal blocked")
            return

        target = whole
        if label is not None:
            # Look for an existing panel crop; if absent, fall through
            # to the whole figure (matches get_figure_image's fallback).
            crops_dir = hash_dir / "figures" / "crops"
            crop_path = crops_dir / f"{figure_id}__{label}.png"
            if crop_path.exists():
                target = crop_path.resolve()

        try:
            data = target.read_bytes()
        except OSError as e:
            await _send_text(send, 500, f"could not read figure: {e}")
            return

        await _send_bytes(send, 200, data, content_type="image/png", head_only=(method == "HEAD"))

    return app


# ---------------------------------------------------------------------------
# Tiny ASGI response helpers — avoid pulling in Starlette just for this.
# ---------------------------------------------------------------------------


async def _send_text(send, status: int, msg: str) -> None:
    body = (msg + "\n").encode("utf-8")
    await send({
        "type": "http.response.start",
        "status": status,
        "headers": [
            (b"content-type", b"text/plain; charset=utf-8"),
            (b"content-length", str(len(body)).encode("ascii")),
        ],
    })
    await send({"type": "http.response.body", "body": body})


async def _send_json(send, status: int, obj) -> None:
    body = (json.dumps(obj) + "\n").encode("utf-8")
    await send({
        "type": "http.response.start",
        "status": status,
        "headers": [
            (b"content-type", b"application/json; charset=utf-8"),
            (b"content-length", str(len(body)).encode("ascii")),
        ],
    })
    await send({"type": "http.response.body", "body": body})


async def _send_bytes(
    send, status: int, body: bytes, content_type: str = "application/octet-stream",
    head_only: bool = False,
) -> None:
    await send({
        "type": "http.response.start",
        "status": status,
        "headers": [
            (b"content-type", content_type.encode("ascii")),
            (b"content-length", str(len(body)).encode("ascii")),
        ],
    })
    await send({"type": "http.response.body", "body": b"" if head_only else body})
