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

* ``_BearerAuthASGI`` accepts expiring figure-only capabilities or an
  operator's bearer credential. Figure capabilities cannot authorize MCP.
* ``#51`` / ``#101`` figure-licensing check keyed to the request
  ``?profile=`` (refuses figures the parent work isn't licensed-cleared
  for under a *strict* profile; the permissive ``report`` default
  allows them). The profile is encoded into the URL by
  ``get_figure_url`` so the fetch enforces the same policy.
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Optional
from urllib.parse import parse_qs

logger = logging.getLogger(__name__)


# Same constraint as the MCP tool: only allow short hex hashes
# (12 chars, matching ``short_hash``) and figure ids matching the
# usual ``<paper_hash>__fig-NN`` etc. shape. Defense in depth — the
# Path() construction below also rejects ``..`` traversal.
_HASH_RE = re.compile(r"^[a-f0-9]{12}$")
_FIGURE_ID_RE = re.compile(r"^[A-Za-z0-9._\-]+$")


def make_figure_app(idx, default_profile: Optional[str] = None):
    """Build an ASGI handler that serves
    ``GET /figures/<paper_hash>/<figure_id>?profile=...&label=...``.

    ``idx`` is the live :class:`mcpsrv.indexes.CorpusIndex` so we can
    resolve ``paper_hash → hash_dir`` and probe the bibliographic
    authority for license fields.

    The figure-licensing gate (#101) is keyed to the request's
    ``?profile=`` (the value :func:`get_figure_url` encodes into the URL
    it hands out), falling back to ``default_profile`` — the server's
    ``--default-profile`` — when the query omits it. A strict client's
    URL therefore stays strict on fetch; only an unprofiled URL falls
    back to the server default.
    """
    # Lazy import — keeps the load-time cost off `import mcpsrv` for
    # non-serve uses (tests, bundle distillation).
    from .tools.figures import _license_metadata_for_paper
    from .profiles import get_profile, resolve_profile, unknown_profile_error

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

        # Optional ``?label=<panel>`` (sub-panel crop) and ``?profile=``
        # (the policy the URL was issued under).
        #
        # parse_qs rather than a hand-rolled split (#154): the previous
        # loop never URL-decoded, so a label containing anything
        # percent-encoded (``?label=A%20B``) failed the charset check or
        # silently missed its crop file and fell back to the whole figure.
        # It also mishandled repeated parameters.
        label: Optional[str] = None
        req_profile: Optional[str] = None
        qs = scope.get("query_string", b"").decode("latin-1", errors="replace")
        if qs:
            params = parse_qs(qs, keep_blank_values=True)
            if "label" in params:
                label = params["label"][-1]
                if not _FIGURE_ID_RE.match(label):
                    await _send_text(send, 400, "malformed label")
                    return
            if "profile" in params:
                req_profile = params["profile"][-1]
                # Validate up front (#154). resolve_profile's contract says
                # callers do this: unknown names fall *through* to the
                # server default, so a typo'd ``?profile=manuscipt`` would
                # silently be served under a different policy than the
                # client asked for, with no error. The MCP tools validate;
                # this route must agree with them.
                if get_profile(req_profile) is None:
                    await _send_json(send, 400, {
                        **unknown_profile_error(req_profile),
                        "hint": (
                            "the profile in a figure URL must be one of the "
                            "built-ins; omit it to use the server default"
                        ),
                    })
                    return

        p = idx.papers.get(paper_hash)
        if not p:
            await _send_text(send, 404, f"no such paper_hash: {paper_hash}")
            return

        # Figure-licensing gate (#101) — keyed to the request profile,
        # falling back to the server default. Mirrors get_figure_url.
        active = resolve_profile(req_profile, default_profile)
        if active.figure_licensing == "strict":
            lic = _license_metadata_for_paper(paper_hash)
            if not lic.get("publishable"):
                body = {
                    "error": "figure withheld under a strict profile",
                    "profile": active.name,
                    # #154 §2 — say which state caused it. `no_record` is
                    # an absence of evidence, not a refusal.
                    "publication_clearance": lic.get("publication_clearance"),
                    "license": lic.get("license") or "none recorded",
                    "license_source": lic.get("license_source"),
                    "hint": (
                        "this URL was issued under a strict profile; request "
                        "with profile=report for in-chat display, or read "
                        "get_figure() for the raw license fields"
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

        data = None
        if label is not None:
            roi = next((r for r in fig.get("rois") or [] if r.get("label") == label and r.get("roi_px")), None)
            if roi is not None:
                from .figure_cache import crop_figure
                try:
                    _, data = crop_figure(idx, whole, roi["roi_px"])
                except (OSError, ValueError, TypeError, OverflowError) as exc:
                    await _send_text(send, 500, f"could not crop figure: {exc}")
                    return

        try:
            if data is None:
                data = whole.read_bytes()
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
            (b"cache-control", b"private, no-store"),
            (b"referrer-policy", b"no-referrer"),
        ],
    })
    await send({"type": "http.response.body", "body": b"" if head_only else body})
