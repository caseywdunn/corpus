"""`corpus` CLI entry point.

Stub for v0.3 packaging (#58). The unified subcommand router
(run/check/status/serve/init/bib) lands in #60. This module exists
so `pip install -e .` can wire up the `[project.scripts] corpus`
entry point and `corpus --help` runs.
"""
from __future__ import annotations

import sys

from .version import __version__


_USAGE = """\
usage: corpus [-V|--version] [-h|--help] <command> [<args>]

corpus — MCP server for taxonomic literature corpora.

The unified subcommand surface (run, check, status, serve, init, bib,
completion) lands incrementally across v0.3.0 — see issues #58–#68
and dev_docs/PLAN.md. Until then, use the legacy root-level scripts
(process_corpus.py, mcp_server.py, etc.) as documented in AGENTS.md.
"""


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if any(a in ("-V", "--version") for a in args):
        print(f"corpus {__version__}")
        return 0
    if not args or any(a in ("-h", "--help") for a in args):
        print(_USAGE)
        return 0
    print(_USAGE, file=sys.stderr)
    print(f"corpus: unknown command {args[0]!r}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
