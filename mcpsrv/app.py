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
from typing import Any, Optional

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
