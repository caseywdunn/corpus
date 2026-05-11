"""Entry-point glue: argparse → build CorpusIndex → mcp.run()/SSE.

Imports :mod:`mcpsrv.tools` for its side effects (decorator
registration). Don't remove that import even though it looks unused —
removing it deletes every tool from the running surface.

#47: ``--check`` runs a serve-time pre-flight (bundle present, token
readable, port free, cross-paper DBs loadable) and exits without
binding the port. Distinct from the pipeline-side ``corpus check``
(#62), which validates the *next-run* surface; this one validates
the *next-serve* surface.
"""
from __future__ import annotations

import argparse
import logging
import socket
import sqlite3
import sys
from pathlib import Path
from typing import List, Optional

from .app import _load_json, mcp, set_index
from .indexes import BiblioAuthority, CorpusIndex, TaxonMentionDB
from . import tools as _tools  # noqa: F401  (registers @mcp.tool() decorators)
from .transport import _load_auth_token, _run_sse

from pipeline.taxa import TaxonomyDB

logger = logging.getLogger(__name__)


# Re-exported here so external callers can import without pulling all
# of pipeline.cli; these mirror the v0.3 exit-code convention (#61).
EXIT_OK = 0
EXIT_GENERIC = 1
EXIT_CONFIG_ERROR = 2
EXIT_PRECONDITION = 3


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Root of the processed corpus (contains documents/)",
    )
    parser.add_argument(
        "--taxonomy-db",
        type=Path,
        default=None,
        help="Override path to the Darwin Core taxonomy SQLite "
             "(default: <output_dir>/taxonomy.sqlite). "
             "Build with: python -m pipeline.taxonomy_ingest --source <dwc|dwca|worms> ...",
    )
    parser.add_argument(
        "--biblio-sqlite",
        type=Path,
        default=None,
        help="Override path to the bibliographic authority SQLite "
             "(default: <output_dir>/biblio_authority.sqlite)",
    )
    parser.add_argument(
        "--taxon-mention-sqlite",
        type=Path,
        default=None,
        help="Override path to the taxon mention SQLite "
             "(default: <output_dir>/taxon_mentions.sqlite)",
    )
    parser.add_argument(
        "--instructions",
        type=Path,
        default=None,
        help="Optional markdown file with corpuscle-specific nudges, "
             "prepended to the server's packaged defaults "
             "(mcpsrv/default_instructions.md). The combined text is "
             "returned to MCP clients in InitializeResult.instructions, "
             "which well-behaved clients (Claude Desktop, Claude Code) "
             "inject into the LLM's context at session start. "
             "Default: <output_dir>/instructions.md if it exists.",
    )
    parser.add_argument(
        "--embedding-model",
        default=None,
        help="Override the default HuggingFace embedding model id "
             "(default: BAAI/bge-m3). Must emit vectors of the same "
             "dim as the LanceDB index.",
    )
    parser.add_argument(
        "--transport",
        choices=["stdio", "sse"],
        default="stdio",
        help="MCP transport. stdio (default) is for local MCP clients "
             "that launch this process themselves (Claude Desktop, Claude "
             "Code, Cursor). sse serves over HTTP/Server-Sent-Events for "
             "remote deployments (dev_docs/PLAN.md §10).",
    )
    parser.add_argument(
        "--host", default="127.0.0.1",
        help="Interface to bind when --transport sse (default: 127.0.0.1)",
    )
    parser.add_argument(
        "--port", type=int, default=8080,
        help="Port to bind when --transport sse (default: 8080)",
    )
    parser.add_argument(
        "--auth-token-file", type=Path, default=None,
        help="Path to a file whose contents are the bearer-auth token. "
             "Alternative: set CORPUS_MCP_TOKEN. CLI-literal tokens are "
             "deliberately not supported — they leak via `ps`.",
    )
    parser.add_argument(
        "--name", default=None,
        help="Server name reported to MCP clients (e.g. corpus:my_group). "
             "Default: 'corpus:<basename of output_dir>'. Lets clients tell "
             "multiple corpuscle deployments apart in their server list.",
    )
    parser.add_argument(
        "--check", action="store_true",
        help="Run a serve-time pre-flight (bundle present, token readable, "
             "port free, cross-paper DBs loadable) and exit without binding "
             "the port. Exits 0 on green, 3 on precondition failure. (#47)",
    )
    parser.add_argument(
        "--allow-unpublishable", action="store_true",
        help="Bypass the #51 publishable gate on get_figure_image — return "
             "image bytes regardless of license / age. Use only for local "
             "rights-holder cases (you own the figures, or you're operating "
             "under fair use). Never set on a public-facing deploy.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    if args.check:
        return _serve_check(args)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        stream=sys.stderr,  # stdout is the MCP transport
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    # Corpuscle-aware server name (#17). FastMCP exposes ``name`` as a
    # read-only property backed by the lower-level Server's ``name``
    # attribute, which is plain mutable state — we can rename after
    # the module-level FastMCP is constructed and the change takes
    # effect on subsequent ``initialize`` calls.
    server_name = args.name or f"corpus:{args.output_dir.name}"
    mcp._mcp_server.name = server_name
    logger.info("MCP server name: %s", server_name)

    # MCP instructions surfaced to clients in InitializeResult.instructions.
    # Always serve a packaged default (mcpsrv/default_instructions.md);
    # if a corpuscle-specific instructions.md is present, prepend it so
    # corpus-specific nudges land before the universal defaults.
    # FastMCP exposes ``instructions`` as a read-only property backed by
    # ``_mcp_server.instructions``; that backing attribute is plain
    # mutable state, so setting it here takes effect on the next
    # ``initialize`` call.
    parts: list[str] = []
    corpuscle_instructions = args.instructions or (args.output_dir / "instructions.md")
    if corpuscle_instructions.exists():
        try:
            parts.append(corpuscle_instructions.read_text(encoding="utf-8").strip())
        except Exception as e:
            logger.warning(
                "Could not read corpuscle instructions %s: %s",
                corpuscle_instructions, e,
            )
    default_instructions = Path(__file__).parent / "default_instructions.md"
    try:
        parts.append(default_instructions.read_text(encoding="utf-8").strip())
    except Exception as e:
        logger.warning(
            "Could not read default instructions %s: %s",
            default_instructions, e,
        )
    if parts:
        text = "\n\n".join(p for p in parts if p)
        mcp._mcp_server.instructions = text
        logger.info(
            "Loaded MCP instructions (%d chars from %d source%s)",
            len(text), len(parts), "" if len(parts) == 1 else "s",
        )

    taxonomy_path = args.taxonomy_db or (args.output_dir / "taxonomy.sqlite")
    taxonomy: Optional[TaxonomyDB] = None
    if taxonomy_path.exists():
        try:
            taxonomy = TaxonomyDB(taxonomy_path)
            logger.info("Taxonomy snapshot loaded from %s (%d names)",
                        taxonomy_path, len(taxonomy.name_set()))
        except Exception as e:
            logger.warning("Could not open taxonomy DB at %s: %s", taxonomy_path, e)
    else:
        logger.warning(
            "Taxonomy snapshot not found at %s — search_taxon and taxon tools "
            "will return an error. Build it: "
            "python -m pipeline.taxonomy_ingest --source <dwc|dwca|worms> ...",
            taxonomy_path,
        )

    biblio_path = args.biblio_sqlite or (args.output_dir / "biblio_authority.sqlite")
    biblio: Optional[BiblioAuthority] = None
    if biblio_path.exists():
        try:
            biblio = BiblioAuthority(biblio_path)
            n_works = biblio.conn.execute("SELECT COUNT(*) FROM works").fetchone()[0]
            logger.info("Bibliographic authority loaded from %s (%d works)",
                        biblio_path, n_works)
        except Exception as e:
            logger.warning("Could not open biblio authority at %s: %s", biblio_path, e)
    else:
        logger.warning(
            "Bibliographic authority not found at %s — citation graph tools "
            "will return an error. Build it: python build_biblio_authority.py",
            biblio_path,
        )

    taxon_mention_path = args.taxon_mention_sqlite or (args.output_dir / "taxon_mentions.sqlite")
    taxon_mention_db: Optional[TaxonMentionDB] = None
    if taxon_mention_path.exists():
        try:
            taxon_mention_db = TaxonMentionDB(taxon_mention_path)
            stats = taxon_mention_db.stats()
            logger.info(
                "Taxon mention DB loaded from %s (%s mentions, %s taxa)",
                taxon_mention_path,
                stats.get("total_mentions", "?"),
                stats.get("unique_taxon_ids", "?"),
            )
        except Exception as e:
            logger.warning("Could not open taxon mention DB at %s: %s",
                           taxon_mention_path, e)
    else:
        logger.info(
            "Taxon mention DB not found at %s — get_taxon_mentions will "
            "fall back to per-paper taxa.json scanning. Build it: "
            "python build_taxon_mentions.py",
            taxon_mention_path,
        )

    index = CorpusIndex(
        args.output_dir,
        taxonomy_db=taxonomy,
        biblio_db=biblio,
        taxon_mention_db=taxon_mention_db,
        embedding_model=args.embedding_model,
    )
    # #51 — server-level allow-unpublishable override; carried on the
    # index since that's what every figure tool already has access to.
    index.allow_unpublishable = bool(args.allow_unpublishable)
    if index.allow_unpublishable:
        logger.warning(
            "*** --allow-unpublishable enabled: get_figure_image returns "
            "bytes regardless of license. Use only for local rights-holder "
            "cases; never on a public-facing deploy. ***",
        )
    n = index.load()
    set_index(index)
    logger.info("Serving corpus: %s (%d papers)", args.output_dir, n)

    if args.transport == "stdio":
        mcp.run()
    elif args.transport == "sse":
        token = _load_auth_token(args.auth_token_file)
        _run_sse(args.host, args.port, token)
    else:  # argparse's choices= should prevent this
        raise ValueError(f"unknown transport: {args.transport!r}")
    return 0


def _serve_check(args: argparse.Namespace) -> int:
    """Serve-time pre-flight (#47).

    Distinct from the pipeline-side ``corpus check`` (#62): this one
    validates the surface needed to *serve* the bundle — bundle layout,
    cross-paper DBs loadable, auth token readable, port free. Doesn't
    bind the port; doesn't load any heavy index. Returns ``EXIT_OK`` on
    green, ``EXIT_PRECONDITION`` on any failure.
    """
    # Local import to avoid the heavy console/rich pull-in when the
    # mcpsrv package is imported for non-CLI use.
    from pipeline.console import print_status

    failures: List[str] = []

    # 1. Bundle layout (output_dir + documents/ + at least one summary.json)
    if not args.output_dir.exists():
        print_status(f"output_dir not found: {args.output_dir}", status="fail")
        failures.append(f"output_dir does not exist: {args.output_dir}")
    else:
        documents_dir = args.output_dir / "documents"
        if not documents_dir.is_dir():
            print_status(
                f"documents/ subdir missing under {args.output_dir} — not a "
                "corpus output tree (or bundle layout is wrong)",
                status="fail",
            )
            failures.append("missing documents/")
        else:
            n_docs = sum(1 for p in documents_dir.iterdir() if p.is_dir())
            n_with_summary = sum(
                1 for p in documents_dir.iterdir()
                if p.is_dir() and (p / "summary.json").is_file()
            )
            if n_docs == 0:
                print_status("documents/ is empty (no per-paper subdirs)",
                             status="warn")
            elif n_with_summary == 0:
                print_status(
                    f"documents/ has {n_docs} hash dirs but none carry "
                    "summary.json — bundle distillation may be incomplete",
                    status="fail",
                )
                failures.append("no summary.json found")
            else:
                print_status(
                    f"bundle: {n_docs} hash dirs ({n_with_summary} with summary.json)",
                    status="ok",
                )

        manifest = args.output_dir / "bundle_manifest.json"
        serve_manifest = args.output_dir / "_serve" / "bundle_manifest.json"
        if manifest.is_file():
            # The output_dir we were pointed at is itself a distilled
            # served bundle (e.g. shipped from another host).
            print_status(
                "bundle_manifest.json present — serving a distilled bundle",
                status="ok",
            )
        elif serve_manifest.is_file():
            # We are serving the build bundle, but `corpus run` did
            # also produce a distilled bundle next to it. Both options
            # are valid; this is the expected local-build state.
            print_status(
                "build bundle has no top-level manifest (expected); a "
                "distilled bundle is also available at ./_serve/ "
                "(for shipping the corpus to a different host)",
                status="ok",
            )
        else:
            # No manifest anywhere — `corpus run` either skipped or
            # failed the distillation step.
            print_status(
                "no bundle_manifest.json here and no ./_serve/ distilled "
                "bundle either — `corpus run` should have produced one. "
                "Re-run, or distill manually with `python -m mcpsrv.bundle "
                "<output> <serve_dir> --version vX.Y.Z`",
                status="warn",
            )

    # 2. Cross-paper DBs loadable
    for label, override, default_name in [
        ("taxonomy", args.taxonomy_db, "taxonomy.sqlite"),
        ("biblio_authority", args.biblio_sqlite, "biblio_authority.sqlite"),
        ("taxon_mentions", args.taxon_mention_sqlite, "taxon_mentions.sqlite"),
    ]:
        path = override or (args.output_dir / default_name)
        if not path.exists():
            print_status(
                f"{label}: not present at {path} — affected MCP tools "
                "will return errors at query time",
                status="warn",
            )
            continue
        try:
            conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
            try:
                conn.execute("SELECT 1").fetchone()
            finally:
                conn.close()
            print_status(f"{label}: loadable ({path.name})", status="ok")
        except sqlite3.Error as e:
            print_status(
                f"{label}: file present but unreadable as SQLite ({path}): {e}",
                status="fail",
            )
            failures.append(f"{label} unreadable")

    # 3. Vector DB
    vector_db = args.output_dir / "vector_db" / "lancedb"
    if vector_db.is_dir():
        print_status(f"vector_db/lancedb: present", status="ok")
    else:
        print_status(
            "vector_db/lancedb: absent — semantic-search tools will return "
            "an error. Re-run `corpus run` after `pipeline.embed` populates it.",
            status="warn",
        )

    # 4. Instructions readable (if specified or default exists)
    instructions = args.instructions or (args.output_dir / "instructions.md")
    if instructions.is_file():
        try:
            instructions.read_text(encoding="utf-8")
            print_status(
                f"instructions: readable ({instructions.name})", status="ok",
            )
        except OSError as e:
            print_status(
                f"instructions: file at {instructions} unreadable ({e})",
                status="fail",
            )
            failures.append("instructions unreadable")

    # 5. Auth token (only when explicitly requested)
    import os
    requested_auth = args.auth_token_file is not None or os.environ.get("CORPUS_MCP_TOKEN")
    if args.auth_token_file is not None:
        if not args.auth_token_file.exists():
            print_status(
                f"auth-token-file: {args.auth_token_file} not found",
                status="fail",
            )
            failures.append("auth-token-file missing")
        else:
            try:
                tok = args.auth_token_file.read_text().strip()
                if not tok:
                    print_status(
                        f"auth-token-file: {args.auth_token_file} is empty",
                        status="fail",
                    )
                    failures.append("auth-token-file empty")
                else:
                    print_status(
                        f"auth-token-file: readable ({len(tok)} chars)",
                        status="ok",
                    )
            except OSError as e:
                print_status(
                    f"auth-token-file: unreadable ({e})", status="fail",
                )
                failures.append("auth-token-file unreadable")
    elif args.transport == "sse" and not requested_auth:
        print_status(
            "auth: no --auth-token-file and CORPUS_MCP_TOKEN unset — "
            "server will run open. Fine for localhost; never safe for "
            "public-facing deploys.",
            status="warn",
        )

    # 6. Port free (only when --transport sse)
    if args.transport == "sse":
        if _port_in_use(args.host, args.port):
            print_status(
                f"port: {args.host}:{args.port} already in use. "
                f"Check `lsof -i :{args.port}` (or `ss -ltnp 'sport = :{args.port}'`) "
                f"and either stop the conflicting process or pass "
                f"a different `--port`.",
                status="fail",
            )
            failures.append(f"port {args.port} in use")
        else:
            print_status(f"port: {args.host}:{args.port} free", status="ok")

    print()
    if failures:
        print_status(
            f"{len(failures)} precondition(s) failed:", status="fail",
        )
        for f in failures:
            print(f"  - {f}")
        return EXIT_PRECONDITION
    print_status("ready: `corpus serve` should start cleanly.", status="ok")
    return EXIT_OK


def _port_in_use(host: str, port: int) -> bool:
    """Probe whether ``host:port`` is bindable. Returns True on EADDRINUSE."""
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            s.bind((host, port))
        return False
    except OSError:
        return True


if __name__ == "__main__":
    sys.exit(main())
