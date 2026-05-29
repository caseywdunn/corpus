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
        "--default-profile", choices=["report", "manuscript", "presentation"],
        default="report",
        help="Output profile (#101) applied to figure/citation calls that "
             "omit profile=. Selection is really a per-call client property; "
             "this is only the server fallback. 'report' (default) serves "
             "figures permissively (in-chat display = fair use); 'manuscript'/"
             "'presentation' gate on the publishable flag. Publication-bound "
             "clients should pass profile='manuscript' per call regardless.",
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
    # #101 — server-level output-profile fallback for calls that omit
    # profile=. Carried on the index, which every figure tool can read.
    # The real selection is per-call; this is only the default.
    from .profiles import get_profile
    index.default_profile = args.default_profile
    if get_profile(args.default_profile).figure_licensing == "permissive":
        logger.warning(
            "*** default profile %r serves figures permissively: a call "
            "without profile= returns image bytes/URLs regardless of "
            "license. This is the intended chat default; publication-bound "
            "clients must pass profile='manuscript' per call. ***",
            args.default_profile,
        )
    n = index.load()
    set_index(index)
    logger.info("Serving corpus: %s (%d papers)", args.output_dir, n)

    # #91 — surface capability health at startup. Degraded capabilities
    # (index present but unserveable) get a loud warning; they also make
    # /healthz return 503 and the backed tools raise hard errors rather
    # than returning empty rows. ``absent`` is normal (structured-only
    # corpus) and logged at info level only.
    caps = index.capabilities()
    degraded = {k: c["detail"] for k, c in caps.items() if c["state"] == "degraded"}
    if degraded:
        logger.warning(
            "*** Serving with DEGRADED capabilities: %s. /healthz will "
            "report 503 and the backed tools will refuse to serve. ***",
            ", ".join(f"{k} ({d})" for k, d in degraded.items()),
        )
    else:
        logger.info(
            "Capabilities: %s",
            ", ".join(f"{k}={c['state']}" for k, c in caps.items()),
        )

    if args.transport == "stdio":
        # #69 — start the figure HTTP side-car so get_figure_url can
        # return a working URL on local stdio MCP. Daemon thread; dies
        # when this process exits.
        from .transport import start_stdio_figure_server
        try:
            fig_host, fig_port, fig_token = start_stdio_figure_server(
                index, default_profile=index.default_profile,
            )
            index.figure_url_base = f"http://{fig_host}:{fig_port}"
            index.figure_auth_token = fig_token
        except Exception as e:
            logger.warning(
                "could not start figure HTTP side-car (%s); get_figure_url "
                "will return an error but other tools work normally", e,
            )
            index.figure_url_base = None
            index.figure_auth_token = None
        mcp.run()
    elif args.transport == "sse":
        try:
            token = _load_auth_token(args.auth_token_file)
        except ValueError as e:
            logger.error("auth-token config error: %s", e)
            return 3
        # #69 — build a public URL prefix the tool can hand out. We
        # don't try to be clever about the operator's reverse-proxy
        # hostname; bind host + bind port is the most we know, and
        # operators behind nginx/ALB can rewrite client-side.
        index.figure_url_base = f"http://{args.host}:{args.port}"
        index.figure_auth_token = token
        _run_sse(
            args.host, args.port, token,
            idx=index,
            default_profile=index.default_profile,
        )
    else:  # argparse's choices= should prevent this
        raise ValueError(f"unknown transport: {args.transport!r}")
    return 0


def _check_one_bundle(
    label: str,
    path: Path,
    args: argparse.Namespace,
    failures: List[str],
) -> None:
    """Validate one bundle directory (build or distilled). Each line is
    prefixed with ``[<label>]`` so the dual-bundle output stays
    attributable when both checks fire.

    ``label`` is the human-readable bundle name; ``"distilled bundle"``
    triggers manifest-required semantics, anything else (typically
    ``"build bundle"``) treats the top-level manifest as expected-absent.
    """
    from pipeline.console import print_status

    is_distilled = label == "distilled bundle"

    # documents/ — presence + summary.json count
    documents_dir = path / "documents"
    if not documents_dir.is_dir():
        print_status(
            f"[{label}] documents/ subdir missing under {path}",
            status="fail",
        )
        failures.append(f"[{label}] missing documents/")
    else:
        n_docs = sum(1 for p in documents_dir.iterdir() if p.is_dir())
        n_with_summary = sum(
            1 for p in documents_dir.iterdir()
            if p.is_dir() and (p / "summary.json").is_file()
        )
        if n_docs == 0:
            print_status(
                f"[{label}] documents/ is empty (no per-paper subdirs)",
                status="warn",
            )
        elif n_with_summary == 0:
            print_status(
                f"[{label}] documents/ has {n_docs} hash dirs but none "
                "carry summary.json — bundle layout broken",
                status="fail",
            )
            failures.append(f"[{label}] no summary.json")
        else:
            print_status(
                f"[{label}] {n_docs} hash dirs "
                f"({n_with_summary} with summary.json)",
                status="ok",
            )

    # bundle_manifest.json — required for distilled, not for build
    manifest = path / "bundle_manifest.json"
    if is_distilled:
        if manifest.is_file():
            print_status(
                f"[{label}] bundle_manifest.json present",
                status="ok",
            )
        else:
            print_status(
                f"[{label}] bundle_manifest.json missing — distillation "
                "incomplete. Re-run `corpus run` (it produces _serve/ "
                "automatically) or `python -m mcpsrv.bundle <output> "
                "<serve_dir> --version vX.Y.Z`",
                status="warn",
            )

    # Cross-paper DBs loadable. Per-component overrides apply only to
    # the build (primary) bundle, since ``args`` carries a single
    # override path each — the distilled bundle always uses the default
    # filename relative to its own root.
    for db_label, override, default_name in [
        ("taxonomy", args.taxonomy_db, "taxonomy.sqlite"),
        ("biblio_authority", args.biblio_sqlite, "biblio_authority.sqlite"),
        ("taxon_mentions", args.taxon_mention_sqlite, "taxon_mentions.sqlite"),
    ]:
        if override is not None and not is_distilled:
            db_path = override
        else:
            db_path = path / default_name
        if not db_path.exists():
            print_status(
                f"[{label}] {db_label}: not present at {db_path} — "
                "affected MCP tools will return errors at query time",
                status="warn",
            )
            continue
        try:
            conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
            try:
                conn.execute("SELECT 1").fetchone()
            finally:
                conn.close()
            print_status(
                f"[{label}] {db_label}: loadable ({db_path.name})",
                status="ok",
            )
        except sqlite3.Error as e:
            print_status(
                f"[{label}] {db_label}: file present but unreadable as "
                f"SQLite ({db_path}): {e}",
                status="fail",
            )
            failures.append(f"[{label}] {db_label} unreadable")

    # Vector DB
    vector_db = path / "vector_db" / "lancedb"
    if vector_db.is_dir():
        print_status(f"[{label}] vector_db/lancedb: present", status="ok")
    else:
        print_status(
            f"[{label}] vector_db/lancedb: absent — semantic-search "
            "tools will return an error",
            status="warn",
        )


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

    # 1. Bundle layout — `corpus run` produces two bundles: the build
    # bundle at <output_dir>/ (everything `corpus run` writes) and the
    # distilled bundle at <output_dir>/_serve/ (subset shipped to remote
    # hosts; carries bundle_manifest.json). Validate whichever ones
    # exist, and tell the user when one is missing.
    if not args.output_dir.exists():
        print_status(f"output_dir not found: {args.output_dir}", status="fail")
        failures.append(f"output_dir does not exist: {args.output_dir}")
    else:
        build_dir = args.output_dir
        serve_dir = build_dir / "_serve"

        # Heuristic: if the user pointed straight at a distilled bundle
        # (basename `_serve` or a top-level manifest exists), there is
        # no nested `_serve/` to also check; treat as single-bundle.
        is_already_serve = (
            build_dir.name == "_serve"
            or (build_dir / "bundle_manifest.json").is_file()
        )

        if is_already_serve:
            _check_one_bundle("distilled bundle", build_dir, args, failures)
        else:
            _check_one_bundle("build bundle", build_dir, args, failures)
            if serve_dir.is_dir():
                _check_one_bundle("distilled bundle", serve_dir, args, failures)
            else:
                print_status(
                    "distilled bundle (./_serve/) not found — local serving "
                    "works off the build bundle, but remote deploys ship the "
                    "_serve/ tree. `corpus run` produces both; re-run if the "
                    "missing _serve is unexpected, or build it manually with "
                    "`python -m mcpsrv.bundle <output> <serve_dir> "
                    "--version vX.Y.Z`.",
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
