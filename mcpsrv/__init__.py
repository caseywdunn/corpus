"""Internal package backing the corpus MCP server CLI (#15).

The root-level ``mcp_server.py`` is a thin shim that delegates here so
``python mcp_server.py …`` invocations keep working. The real
implementation is split per-concern across submodules:

* :mod:`mcpsrv.app` — FastMCP instance + index accessor + ``_load_json``
* :mod:`mcpsrv.indexes` — ``CorpusIndex``, ``TaxonMentionDB``, ``BiblioAuthority``
* :mod:`mcpsrv.tools.papers` / ``.taxonomy`` / ``.figures`` / ``.chunks``
  / ``.bibliography`` / ``.lexicon`` — the ``@mcp.tool()``-decorated
  functions, grouped by concern (see dev_docs/MCP_TOOLS.md for the catalog)
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

# Lazy access for `from mcpsrv import main` — eagerly importing
# `mcpsrv.main` here would put it in sys.modules at package-import time,
# which trips a runpy warning when `python -m mcpsrv.main` is then run
# (the corpus `serve` verb does exactly that).
def __getattr__(name):  # noqa: E402
    if name == "main":
        from .main import main
        globals()["main"] = main
        return main
    raise AttributeError(f"module 'mcpsrv' has no attribute {name!r}")


__all__ = ["main"]
