"""Tool surface for the corpus MCP server (#15).

Each submodule defines a thematic group of ``@mcp.tool()``-decorated
functions. Importing this package triggers the imports of every
submodule, which fires the decorators and registers the tools on the
shared FastMCP instance in :mod:`mcpsrv.app`.

Adding a new tool: create or extend the appropriate submodule, then
import it here so the decorator runs at startup.
"""
from . import (  # noqa: F401  (side-effect imports register MCP tools)
    papers,
    taxonomy,
    figures,
    chunks,
    bibliography,
    lexicon,
)
