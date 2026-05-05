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
from .main import main

__all__ = ["main"]
