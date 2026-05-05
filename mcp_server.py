#!/usr/bin/env python3
"""MCP server CLI shim (#15).

The implementation moved to the ``mcpsrv/`` package. This file stays
at the repo root so existing ``python mcp_server.py …`` invocations
keep working unchanged. See ``mcpsrv/__init__.py`` for the new layout.

Usage:

    python mcp_server.py /path/to/output_dir
    python mcp_server.py /path/to/output_dir --transport sse --port 8080

MCP clients (Claude Desktop, Claude Code, Cursor) typically launch
this script via a config file. See README.md for the snippets.
"""
import sys

from mcpsrv import main

if __name__ == "__main__":
    sys.exit(main())
