#!/usr/bin/env python3
"""Embed per-paper chunks into a LanceDB vector index.

Phase E: defaults to a local sentence-transformers backend (BGE-M3,
1024-dim) running on the best available accelerator (CUDA → MPS → CPU).
The OpenAI backend is still available behind ``--backend openai`` for
comparison, but is no longer the default — see PLAN.md §7.

Reads ``<output_dir>/documents/<HASH>/chunks.json`` for each
already-processed PDF, batches the chunk text through the embedding
backend, and writes records into ``<output_dir>/vector_db/lancedb/``.
A per-hash marker file at ``<output_dir>/vector_db/<HASH>_embedded.done``
records the backend + model + dim used, so ``--resume`` knows what to
skip and ``--rebuild`` knows when to drop and re-embed.

Usage:
    python embed_chunks.py demo_output                         # local default
    python embed_chunks.py demo_output --backend openai        # legacy path
    python embed_chunks.py demo_output --model BAAI/bge-m3
    python embed_chunks.py demo_output --pdf-hash af043530e5dd # one doc
    python embed_chunks.py demo_output --resume                # skip embedded
    python embed_chunks.py demo_output --rebuild               # drop the table first
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Type

import lancedb
from lancedb.pydantic import LanceModel, Vector

from dotenv import load_dotenv
load_dotenv()

from embeddings import EmbeddingBackend, EmbeddingError, get_embedder

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# LanceDB schema (built dynamically per backend dim)
# ---------------------------------------------------------------------------


class ChunkMetadata(LanceModel):
    pdf_hash: str
    filename: str = ""
    title: str = ""
    authors: List[str] = []
    year: Optional[int] = None
    journal: str = ""
    doi: str = ""
    total_pages: Optional[int] = None
    chunk_id: str = ""
    relative_paths: List[str] = []
    section_class: Optional[str] = None
    headings: List[str] = []


def make_chunk_model(dim: int) -> Type[LanceModel]:
    """Build a LanceDB model class with a Vector of the given dim.

    Different embedding backends have different output dimensions
    (text-embedding-3-small=1536, BGE-M3=1024, MiniLM-L6=384). Rather
    than hard-code one schema, we generate the model class at runtime
    so the table created on disk matches what the backend produces.

    The class is given a stable ``__qualname__`` so LanceDB's
    pydantic-derived schema stays deterministic across runs.
    """
    cls = type(
        f"DocumentChunk{dim}",
        (LanceModel,),
        {
            "__annotations__": {
                "text": str,
                "vector": Vector(dim),
                "metadata": ChunkMetadata,
            },
            "text": "",
            "metadata": ChunkMetadata(pdf_hash=""),
        },
    )
    return cls


# ---------------------------------------------------------------------------
# Per-document embedding
# ---------------------------------------------------------------------------


def load_document_data(hash_dir: Path) -> Dict:
    """Load summary + chunks + metadata for one paper."""
    data: Dict = {}
    for fname in ("summary.json", "chunks.json", "metadata.json"):
        f = hash_dir / fname
        if not f.exists():
            continue
        try:
            with f.open(encoding="utf-8") as fh:
                data[fname.removesuffix(".json")] = json.load(fh)
        except Exception as e:
            logger.warning("Could not read %s: %s", f, e)
    return data


def _build_records(
    hash_dir: Path,
    chunks: List[Dict],
    vectors: List[List[float]],
    chunk_model: Type[LanceModel],
    summary: Dict,
    metadata: Dict,
    chunks_meta: Dict,
) -> List[LanceModel]:
    """Stitch chunk text + computed vector + paper metadata into LanceDB rows."""
    records = []
    relative_paths = summary.get("relative_paths", []) or []
    filename = relative_paths[0] if relative_paths else ""
    # Authors are stored in metadata.json as list of dicts; flatten to a
    # list of "Forename Surname" strings for the LanceDB record (the
    # schema keeps it simple — joining now beats joining at query time).
    authors_list = []
    for a in metadata.get("authors", []) or []:
        name = f"{(a.get('forename','') or '').strip()} {(a.get('surname','') or '').strip()}".strip()
        if name:
            authors_list.append(name)

    for chunk, vec in zip(chunks, vectors):
        cm = ChunkMetadata(
            pdf_hash=hash_dir.name,
            filename=filename,
            title=metadata.get("title", "") or "",
            authors=authors_list,
            year=metadata.get("year"),
            journal=metadata.get("journal", "") or "",
            doi=metadata.get("doi", "") or "",
            total_pages=(chunks_meta.get("metadata") or {}).get("total_pages"),
            chunk_id=chunk.get("chunk_id", ""),
            relative_paths=relative_paths,
            section_class=chunk.get("section_class"),
            headings=chunk.get("headings", []) or [],
        )
        records.append(chunk_model(
            text=chunk.get("text", "") or "",
            vector=vec,
            metadata=cm,
        ))
    return records


def embed_document(
    hash_dir: Path,
    table,
    chunk_model: Type[LanceModel],
    embedder: EmbeddingBackend,
) -> Optional[int]:
    """Embed one paper's chunks and add them to the table.

    Returns the number of chunks added, or None if the paper had no
    chunks to embed. Raises :class:`EmbeddingError` on failure — the
    caller decides whether to skip the doc or abort.
    """
    pdf_hash = hash_dir.name
    doc = load_document_data(hash_dir)
    chunks_data = doc.get("chunks") or {}
    chunks = chunks_data.get("chunks") or []
    if not chunks:
        logger.info("%s: no chunks; skipping", pdf_hash)
        return None

    texts = [c.get("text", "") or "" for c in chunks]
    logger.info("%s: embedding %d chunks", pdf_hash, len(texts))
    vectors = embedder.embed(texts)
    if len(vectors) != len(texts):
        raise EmbeddingError(
            f"{pdf_hash}: embedder returned {len(vectors)} vectors for "
            f"{len(texts)} chunks"
        )

    records = _build_records(
        hash_dir, chunks, vectors, chunk_model,
        summary=doc.get("summary") or {},
        metadata=doc.get("metadata") or {},
        chunks_meta=chunks_data,
    )
    table.add(records)
    return len(records)


# ---------------------------------------------------------------------------
# Marker file (records what we used for this paper, for resume + audit)
# ---------------------------------------------------------------------------


def _write_marker(marker_file: Path, *, pdf_hash: str, count: int,
                  relative_paths: List[str], embedder: EmbeddingBackend) -> None:
    marker_file.write_text(json.dumps({
        "pdf_hash": pdf_hash,
        "chunks_count": count,
        "relative_paths": relative_paths,
        "embedding_backend": embedder.__class__.__name__,
        "embedding_model": embedder.model_name,
        "embedding_dim": embedder.dim,
        "status": "completed",
    }, indent=2))


def _marker_matches_backend(marker_file: Path, embedder: EmbeddingBackend) -> bool:
    """Return True iff the existing marker was written with the same model
    + dim. Resume should re-embed if the backend changed underneath us."""
    try:
        m = json.loads(marker_file.read_text())
    except Exception:
        return False
    return (
        m.get("embedding_model") == embedder.model_name
        and m.get("embedding_dim") == embedder.dim
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--backend", choices=["local", "openai"], default="local",
        help="Embedding backend (default: local sentence-transformers)",
    )
    parser.add_argument(
        "--model", default=None,
        help="Override the per-backend default model "
             "(local default BAAI/bge-m3, openai default text-embedding-3-small)",
    )
    parser.add_argument(
        "--device", default=None,
        help="Force a torch device for the local backend (cuda|mps|cpu); "
             "default autodetects CUDA → MPS → CPU",
    )
    parser.add_argument("--pdf-hash", help="Process only this 12-char hash")
    parser.add_argument("--resume", action="store_true",
                        help="Skip docs whose marker matches this backend+model")
    parser.add_argument(
        "--rebuild", action="store_true",
        help="Drop the existing LanceDB table before embedding "
             "(needed if you switch backends and don't want stale rows)",
    )
    parser.add_argument(
        "--table-name", default="document_chunks",
        help="LanceDB table name (default: document_chunks)",
    )
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    output_dir = args.output_dir.resolve()
    documents_dir = output_dir / "documents"
    vector_db_dir = output_dir / "vector_db"
    if not documents_dir.exists():
        logger.error("Documents directory %s does not exist", documents_dir)
        return 1
    vector_db_dir.mkdir(exist_ok=True)

    # Build the embedder up front so we know the dim before opening the
    # LanceDB table.
    embedder_kwargs: Dict = {}
    if args.device:
        embedder_kwargs["device"] = args.device
    embedder = get_embedder(args.backend, args.model, **embedder_kwargs)
    logger.info(
        "Embedding backend=%s model=%s dim=%d",
        args.backend, embedder.model_name, embedder.dim,
    )

    chunk_model = make_chunk_model(embedder.dim)

    db = lancedb.connect(str(vector_db_dir / "lancedb"))
    if args.rebuild and args.table_name in db.table_names():
        logger.warning("Dropping existing table %r per --rebuild", args.table_name)
        db.drop_table(args.table_name)

    if args.table_name in db.table_names():
        table = db.open_table(args.table_name)
        # Sanity check: the table on disk has a fixed dim. If the user
        # switched to a different model without --rebuild, fail early.
        existing_dim = table.schema.field("vector").type.list_size
        if existing_dim != embedder.dim:
            logger.error(
                "Table %r has dim=%d but the chosen backend produces dim=%d. "
                "Use --rebuild to drop and recreate, or pick a matching model.",
                args.table_name, existing_dim, embedder.dim,
            )
            return 2
    else:
        table = db.create_table(args.table_name, schema=chunk_model.to_arrow_schema())
        logger.info("Created LanceDB table %r (dim=%d)", args.table_name, embedder.dim)

    # Determine which docs to process
    if args.pdf_hash:
        hash_dirs = [documents_dir / args.pdf_hash]
        if not hash_dirs[0].exists():
            logger.error("PDF hash directory %s not found", args.pdf_hash)
            return 1
    else:
        hash_dirs = sorted(d for d in documents_dir.iterdir() if d.is_dir())
    logger.info("Found %d document(s) to process", len(hash_dirs))

    n_ok = n_skip = n_fail = 0
    for hash_dir in hash_dirs:
        pdf_hash = hash_dir.name
        marker = vector_db_dir / f"{pdf_hash}_embedded.done"

        if args.resume and marker.exists() and _marker_matches_backend(marker, embedder):
            logger.info("Skipping %s (already embedded with %s/%s)",
                        pdf_hash, embedder.__class__.__name__, embedder.model_name)
            n_skip += 1
            continue

        try:
            count = embed_document(hash_dir, table, chunk_model, embedder)
        except EmbeddingError as e:
            logger.error("%s: embedding failed: %s — skipping", pdf_hash, e)
            n_fail += 1
            continue
        except Exception as e:
            logger.exception("%s: unexpected error: %s — skipping", pdf_hash, e)
            n_fail += 1
            continue

        if count is None:
            n_skip += 1
            continue

        # Marker reflects what we just wrote
        summary = load_document_data(hash_dir).get("summary") or {}
        _write_marker(
            marker,
            pdf_hash=pdf_hash, count=count,
            relative_paths=summary.get("relative_paths", []) or [],
            embedder=embedder,
        )
        n_ok += 1

    logger.info("Done. embedded=%d, skipped=%d, failed=%d", n_ok, n_skip, n_fail)
    return 0 if n_fail == 0 else 3


if __name__ == "__main__":
    sys.exit(main())
