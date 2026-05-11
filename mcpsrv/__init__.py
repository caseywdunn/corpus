"""Internal package backing the corpus MCP server CLI (#15).

The root-level ``mcp_server.py`` is a thin shim that delegates here so
``python mcp_server.py …`` invocations keep working. The real
implementation is split per-concern across submodules:

* :mod:`mcpsrv.app` — FastMCP instance + index accessor + ``_load_json``
* :mod:`mcpsrv.indexes` — ``CorpusIndex``, ``TaxonMentionDB``, ``BiblioAuthority``
* :mod:`mcpsrv.tools.papers` / ``.taxonomy`` / ``.figures`` / ``.chunks``
  / ``.bibliography`` — the 27 ``@mcp.tool()``-decorated functions, grouped
  by concern
* :mod:`mcpsrv.transport` — bearer-auth ASGI middleware + SSE runner
* :mod:`mcpsrv.main` — argparse + index construction + transport dispatch
"""
# Honor CORPUS_LOG_LEVEL env var before any submodule's basicConfig
# runs (see pipeline/__init__.py for the rationale).
import logging as _logging  # noqa: E402
import os as _os  # noqa: E402
_log_level = _os.environ.get("CORPUS_LOG_LEVEL", "").upper()
if _log_level in {"WARNING", "INFO", "DEBUG"}:
    _logging.basicConfig(
        level=getattr(_logging, _log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

from .main import main  # noqa: E402

__all__ = ["main"]
