"""Entry-point glue: argparse → build CorpusIndex → mcp.run()/SSE.

Imports :mod:`mcpsrv.tools` for its side effects (decorator
registration). Don't remove that import even though it looks unused —
removing it deletes every tool from the running surface.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

from .app import _load_json, mcp, set_index
from .indexes import BiblioAuthority, CorpusIndex, TaxonMentionDB
from . import tools as _tools  # noqa: F401  (registers @mcp.tool() decorators)
from .transport import _load_auth_token, _run_sse

from taxa import TaxonomyDB

logger = logging.getLogger(__name__)


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
             "Build with: python ingest_taxonomy.py --source <dwc|dwca|worms> ...",
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
        help="Markdown file whose contents are returned to MCP clients in "
             "the InitializeResult.instructions field, so well-behaved "
             "clients (Claude Desktop, Claude Code) inject it into the LLM's "
             "context at session start.  Use it for per-corpus nudges — e.g. "
             "'Velella is not a siphonophore'.  "
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
             "remote deployments (PLAN.md §10).",
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
        help="Server name reported to MCP clients (e.g. corpus:siphonophores). "
             "Default: 'corpus:<basename of output_dir>'. Lets clients tell "
             "multiple corpuscle deployments apart in their server list.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

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

    # Optional per-corpus instructions surfaced to clients in
    # InitializeResult.instructions.  FastMCP exposes ``instructions`` as
    # a read-only property backed by ``_mcp_server.instructions``; that
    # backing attribute is plain mutable state on the lower-level Server,
    # so setting it here (after the module-level FastMCP is constructed)
    # takes effect on the next ``initialize`` call.
    instructions_path = args.instructions or (args.output_dir / "instructions.md")
    if instructions_path.exists():
        try:
            text = instructions_path.read_text(encoding="utf-8").strip()
            mcp._mcp_server.instructions = text
            logger.info(
                "Loaded MCP instructions from %s (%d chars)",
                instructions_path, len(text),
            )
        except Exception as e:
            logger.warning(
                "Could not read instructions %s: %s", instructions_path, e,
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
            "python ingest_taxonomy.py --source <dwc|dwca|worms> ...",
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


if __name__ == "__main__":
    sys.exit(main())
