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


# Module-level FastMCP instance — every tools module decorates this one.
mcp = FastMCP("corpus")

# Set by :func:`mcpsrv.main.main` after the corpus index is built.
# Forward-declared here so tools can reference ``_INDEX`` and
# ``_need_index`` without a circular import on the indexes module.
_INDEX: Optional["CorpusIndex"] = None  # type: ignore  # noqa: F821


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
