"""Module-level FastMCP instance, index accessor, and shared helpers.

Every tools module imports ``mcp`` from here and registers via the
``@mcp.tool()`` decorator. The package's ``__init__`` triggers those
imports so all tools land on the same FastMCP instance.

``_INDEX`` is populated at startup by :func:`mcpsrv.main.main`. Tools
call :func:`_need_index` rather than touching ``_INDEX`` directly so
the not-initialized error is uniform.
"""
from __future__ import annotations

import json
import logging
import time
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:
    # Forward reference for _INDEX's type annotation. Circular at
    # runtime (indexes.py imports from this module via mcp), so
    # gated behind TYPE_CHECKING — tools/test_no_undefined_names.py
    # requires the name be importable in some context.
    from .indexes import CorpusIndex  # noqa: F401

from dotenv import load_dotenv
load_dotenv()  # Ensure ANTHROPIC_API_KEY etc. are visible to the stdio subprocess
                # Claude Code launches — no shell inheritance to rely on.

from mcp.server.fastmcp import FastMCP

logger = logging.getLogger(__name__)


def _load_json(path: Path, default: Any = None) -> Any:
    """Best-effort JSON load. Missing or malformed files return ``default``
    — the index is lossy by design so a single broken per-paper file
    doesn't prevent the server from starting."""
    try:
        if path.exists() and path.stat().st_size > 0:
            with path.open(encoding="utf-8") as f:
                return json.load(f)
    except Exception as e:
        logger.warning("Could not read %s: %s", path, e)
    return default


# Safety cap on every paginated MCP tool. Existing per-tool defaults
# (50, 100, …) live below this; the cap just prevents pathological
# limit=10_000 calls from spilling massive payloads into chat context.
MAX_LIMIT = 500


def _validated_limit(limit: int, *, max_value: int = MAX_LIMIT) -> int:
    """Validate + clamp a tool's ``limit`` parameter (#86).

    Returns the clamped int. Raises ``ValueError`` on ``limit < 1`` so
    a caller passing ``limit=0`` gets a clear error instead of an
    unbounded response. Several tools used to treat ``limit=0`` as
    "unlimited" via ``rows[:limit] if limit else rows`` — a footgun
    the moment a client mistakes 0 for "no rows please". Callers
    convert the exception into their usual error-response shape::

        try:
            n = _validated_limit(limit)
        except ValueError as e:
            return [{"error": str(e)}]
    """
    n = int(limit)
    if n < 1:
        raise ValueError(
            f"limit must be >= 1 (got {n}); pass a positive integer or "
            "omit the parameter to use the tool's default"
        )
    return min(n, max_value)


# ---------------------------------------------------------------------------
# Per-tool instrumentation (Phase 0 — feeds #90 run-log + #91 /healthz)
# ---------------------------------------------------------------------------
#
# A single dispatch-layer interceptor records, per tool, the call count,
# transport-level error count, and cumulative latency. We hook
# ``ToolManager.call_tool`` rather than wrapping each ``@mcp.tool()``
# function: FastMCP builds every tool's JSON-schema from the *function
# signature* (``Tool.from_function``), so wrapping the functions would
# risk corrupting that introspection. ``FastMCP.call_tool`` looks up
# ``self._tool_manager.call_tool`` dynamically on every call
# (server.py), so replacing that bound method takes effect for both the
# stdio and SSE transports.
#
# "error" here means the dispatch raised (unknown tool, validation
# failure, or an unhandled exception in the tool body). Tools that
# *return* an ``{"error": ...}`` row are a normal response at this layer
# and are not counted — the Phase-3 uniform error-row shape standardizes
# those payloads separately.


class _ToolStats:
    """Process-wide per-tool counters. Read via :meth:`snapshot`."""

    def __init__(self) -> None:
        self.calls: dict[str, int] = defaultdict(int)
        self.errors: dict[str, int] = defaultdict(int)
        self._latency_ms: dict[str, float] = defaultdict(float)

    def record(self, name: str, latency_ms: float, ok: bool) -> None:
        self.calls[name] += 1
        if not ok:
            self.errors[name] += 1
        self._latency_ms[name] += latency_ms

    @property
    def total_calls(self) -> int:
        return sum(self.calls.values())

    @property
    def total_errors(self) -> int:
        return sum(self.errors.values())

    def snapshot(self) -> dict[str, dict[str, Any]]:
        """Per-tool ``{calls, errors, avg_latency_ms}``, sorted by name."""
        return {
            name: {
                "calls": self.calls[name],
                "errors": self.errors[name],
                "avg_latency_ms": (
                    round(self._latency_ms[name] / self.calls[name], 1)
                    if self.calls[name]
                    else 0.0
                ),
            }
            for name in sorted(self.calls)
        }


TOOL_STATS = _ToolStats()


def _args_summary(arguments: Any) -> str:
    """Comma-joined sorted arg names — never values, to keep logs safe."""
    if isinstance(arguments, dict):
        return ",".join(sorted(str(k) for k in arguments))
    return ""


def _install_call_instrumentation(mcp_instance: "FastMCP") -> None:
    """Wrap ``_tool_manager.call_tool`` to feed :data:`TOOL_STATS`."""
    tm = mcp_instance._tool_manager
    if getattr(tm, "_corpus_instrumented", False):
        return  # idempotent — guard against double-import
    orig_call_tool = tm.call_tool

    async def _instrumented_call_tool(name, arguments=None, *args, **kwargs):
        start = time.perf_counter()
        ok = True
        try:
            return await orig_call_tool(name, arguments, *args, **kwargs)
        except Exception:
            ok = False
            raise
        finally:
            latency_ms = (time.perf_counter() - start) * 1000.0
            TOOL_STATS.record(name, latency_ms, ok)
            logger.info(
                "mcp_tool name=%s args=[%s] latency_ms=%.1f status=%s",
                name, _args_summary(arguments), latency_ms,
                "ok" if ok else "error",
            )

    tm.call_tool = _instrumented_call_tool  # type: ignore[method-assign]
    tm._corpus_instrumented = True  # type: ignore[attr-defined]


# Module-level FastMCP instance — every tools module decorates this one.
mcp = FastMCP("corpus")
_install_call_instrumentation(mcp)

# Set by :func:`mcpsrv.main.main` after the corpus index is built.
# Forward-declared here so tools can reference ``_INDEX`` and
# ``_need_index`` without a circular import on the indexes module.
_INDEX: Optional["CorpusIndex"] = None  # type: ignore  # noqa: F821


# ---------------------------------------------------------------------------
# Uniform tool error payload (Phase 3 freeze gate — #88 §3)
# ---------------------------------------------------------------------------
#
# Every tool error *return* (as opposed to a raised hard error) carries a
# human-readable ``error`` message plus a stable machine-readable
# ``code``. 1.0 freezes this shape, so clients can branch on ``code``
# without string-matching the message. The canonical codes:
#
#   not_found        a named entity (paper/chunk/taxon/work/figure) doesn't exist
#   ambiguous        a lookup matched more than one candidate
#   invalid_argument a parameter was missing, malformed, or out of range
#   not_configured   the tool's backing index/DB isn't part of this corpus
#   no_results       a valid query simply produced nothing
#   unavailable      a required capability is degraded / an upstream call failed
#   empty_item       a batch element was empty/blank
#   forbidden        a policy gate refused (e.g. figure licensing / profile)

ERROR_CODES = frozenset({
    "not_found", "ambiguous", "invalid_argument", "not_configured",
    "no_results", "unavailable", "empty_item", "forbidden",
})


def error(message: str, code: str, **extra: Any) -> dict:
    """Build a uniform tool error payload ``{"error", "code", **extra}``.

    ``extra`` carries any tool-specific context (e.g. ``matches``,
    ``supported_styles``, ``queried``) alongside the frozen two keys.
    """
    return {"error": message, "code": code, **extra}


def _need_index():
    """Return the active CorpusIndex; raise if main() hasn't run yet."""
    if _INDEX is None:
        raise RuntimeError(
            "CorpusIndex not initialized; start the server via main()"
        )
    return _INDEX


def set_index(index) -> None:
    """Setter used by :func:`mcpsrv.main.main` after building the index."""
    global _INDEX
    _INDEX = index
