#!/usr/bin/env python3
"""Embed per-paper chunks into a LanceDB vector index.

Uses a local sentence-transformers backend (BGE-M3, 1024-dim) on the
best available accelerator (CUDA → MPS → CPU).

Reads ``<output_dir>/documents/<HASH>/chunks.json`` for each
already-processed PDF, batches the chunk text through the embedding
backend, and writes records into ``<output_dir>/vector_db/lancedb/``.
A per-hash marker file at ``<output_dir>/vector_db/<HASH>_embedded.done``
records the input fingerprint, model, dimension and committed generation.
``--resume`` verifies these against the current artifacts and table rows.

Usage:
    python -m pipeline.embed demo_output                         # default
    python -m pipeline.embed demo_output --model BAAI/bge-m3
    python -m pipeline.embed demo_output --pdf-hash af043530e5dd # one doc
    python -m pipeline.embed demo_output --resume                # skip embedded
    python -m pipeline.embed demo_output --rebuild               # drop the table first
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import tempfile
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Type

import lancedb
import pyarrow as pa
from lancedb.pydantic import LanceModel, Vector

from dotenv import load_dotenv
load_dotenv()

from pipeline.embeddings import (
    EmbeddingBackend,
    EmbeddingError,
    get_embedder,
    lancedb_table_names,
    _DEFAULT_LOCAL_MODEL,
)
from pipeline.embedding_state import (
    MARKER_VERSION, input_fingerprint, load_document_data, marker_problem,
    read_marker, record_payloads, row_census,
)
from pipeline.version import __version__

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
                "embedding_generation": Optional[str],
            },
            "text": "",
            "metadata": ChunkMetadata(pdf_hash=""),
            "embedding_generation": None,
        },
    )
    return cls


# ---------------------------------------------------------------------------
# Per-document embedding
# ---------------------------------------------------------------------------


def embed_document(
    hash_dir: Path,
    table,
    chunk_model: Type[LanceModel],
    embedder: EmbeddingBackend,
    *,
    marker_file: Optional[Path] = None,
) -> int:
    """Atomically replace one paper, then publish its completion evidence.

    A fresh generation matches no existing row: one merge inserts all incoming
    chunks and deletes every old row for this hash, including legacy duplicates.
    Failure before the commit leaves the old rows; failure after it leaves a
    stale/missing marker that cannot pass resume. Build directories are single
    writer: do not edit extraction inputs while embedding or bundling.
    """
    pdf_hash = hash_dir.name
    doc = load_document_data(hash_dir)
    fingerprint = input_fingerprint(pdf_hash, doc)
    payloads = record_payloads(pdf_hash, doc)
    texts = [row["text"] for row in payloads]
    logger.info("%s: embedding %d chunks", pdf_hash, len(texts))
    vectors = embedder.embed(texts) if texts else []
    if len(vectors) != len(texts):
        raise EmbeddingError(
            f"{pdf_hash}: embedder returned {len(vectors)} vectors for "
            f"{len(texts)} chunks"
        )

    generation = uuid.uuid4().hex
    records = [chunk_model(**row, vector=vector, embedding_generation=generation)
               .model_dump() for row, vector in zip(payloads, vectors)]
    if fingerprint != input_fingerprint(pdf_hash, load_document_data(hash_dir)):
        raise EmbeddingError(f"{pdf_hash}: inputs changed during embedding; retry")
    # Nullable migration preserves legacy rows until their replacements commit.
    if "embedding_generation" not in table.schema.names:
        table.add_columns(pa.field("embedding_generation", pa.string(), nullable=True))
    predicate = "metadata.pdf_hash = '" + pdf_hash.replace("'", "''") + "'"
    if records:
        rows = pa.Table.from_pylist(records, schema=table.schema)
        (table.merge_insert("embedding_generation")
         .when_not_matched_insert_all()
         .when_not_matched_by_source_delete(predicate)
         .execute(rows))
    else:
        # An explicitly empty chunk list is valid and removes previous rows.
        table.delete(predicate)
    marker_file = marker_file or hash_dir.parent.parent / "vector_db" / f"{pdf_hash}_embedded.done"
    _write_marker(
        marker_file, pdf_hash=pdf_hash, count=len(records),
        relative_paths=doc["summary"].get("relative_paths") or [],
        embedder=embedder, fingerprint=fingerprint, generation=generation,
        table_name=table.name,
    )
    return len(records)


# ---------------------------------------------------------------------------
# Marker file (records what we used for this paper, for resume + audit)
# ---------------------------------------------------------------------------


def _write_marker(marker_file: Path, *, pdf_hash: str, count: int,
                  relative_paths: List[str], embedder: EmbeddingBackend,
                  fingerprint: str, generation: str, table_name: str) -> None:
    data = {
        "marker_version": MARKER_VERSION,
        "pipeline_version": __version__,
        "pdf_hash": pdf_hash,
        "chunks_count": count,
        "relative_paths": relative_paths,
        "embedding_backend": embedder.__class__.__name__,
        "embedding_model": embedder.model_name,
        "embedding_dim": embedder.dim,
        "embedding_table": table_name,
        "input_fingerprint": fingerprint,
        "embedding_generation": generation,
        "status": "completed",
    }
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8",
                                         dir=marker_file.parent, delete=False) as fh:
            temporary = Path(fh.name)
            json.dump(data, fh, indent=2)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(temporary, marker_file)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--model", default=None,
        help="Override the default HuggingFace model id (default: BAAI/bge-m3)",
    )
    parser.add_argument(
        "--device", default=None,
        help="Force a torch device (cuda|mps|cpu); "
             "default autodetects CUDA → MPS → CPU",
    )
    parser.add_argument("--pdf-hash", help="Process only this 12-char hash")
    parser.add_argument("--resume", action="store_true",
                        help="Skip docs with verified current inputs and committed rows")
    parser.add_argument(
        "--rebuild", action="store_true",
        help="Drop the existing LanceDB table before embedding "
             "(needed if you switch models and don't want stale rows)",
    )
    parser.add_argument(
        "--table-name", default="document_chunks",
        help="LanceDB table name (default: document_chunks)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Report which hashes would be embedded without loading the model "
             "or writing to LanceDB.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()
    if args.rebuild and args.pdf_hash:
        parser.error("--rebuild replaces the whole index; it cannot be used with --pdf-hash")
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    output_dir = args.output_dir.resolve()
    documents_dir = output_dir / "documents"
    vector_db_dir = output_dir / "vector_db"
    if not documents_dir.exists():
        if args.dry_run:
            # Nothing has been extracted yet. A dry-run reports the plan
            # for the corpuscle as it stands; "0 to embed" is that plan,
            # not a failure.
            logger.info(
                "Dry-run: %s does not exist yet; the extract step creates "
                "it on a real run. Would embed 0 hash(es).", documents_dir,
            )
            return 0
        logger.error("Documents directory %s does not exist", documents_dir)
        return 1

    if args.pdf_hash:
        if args.pdf_hash in {".", ".."} or Path(args.pdf_hash).name != args.pdf_hash:
            parser.error("--pdf-hash must be a document directory name")
        hash_dirs = [documents_dir / args.pdf_hash]
        if not hash_dirs[0].is_dir():
            logger.error("PDF hash directory %s not found", args.pdf_hash)
            return 1
    else:
        hash_dirs = sorted(d for d in documents_dir.iterdir() if d.is_dir())

    if args.dry_run:
        census = {}
        dim = None
        index_path = vector_db_dir / "lancedb"
        if args.resume and not args.rebuild and index_path.exists():
            db = lancedb.connect(str(index_path))
            if args.table_name in lancedb_table_names(db):
                table = db.open_table(args.table_name)
                census = row_census(table)
                dim = table.schema.field("vector").type.list_size
        n_would_embed = n_skip = 0
        for hash_dir in hash_dirs:
            marker = read_marker(vector_db_dir / f"{hash_dir.name}_embedded.done")
            try:
                current = dim is not None and marker_problem(
                    hash_dir, marker, census,
                    model=args.model or _DEFAULT_LOCAL_MODEL, dim=dim,
                    table_name=args.table_name,
                ) is None
            except (OSError, ValueError, TypeError, AttributeError):
                current = False
            if args.resume and not args.rebuild and current:
                n_skip += 1
            else:
                n_would_embed += 1
        logger.info(
            "Dry-run: would embed %d hash(es); %d verified current. "
            "No model loaded, no LanceDB writes.",
            n_would_embed, n_skip,
        )
        return 0

    vector_db_dir.mkdir(exist_ok=True)

    # Build the embedder up front so we know the dim before opening the
    # LanceDB table.
    embedder_kwargs: Dict = {}
    if args.device:
        embedder_kwargs["device"] = args.device
    embedder = get_embedder(args.model, **embedder_kwargs)
    logger.info("Embedding model=%s dim=%d", embedder.model_name, embedder.dim)

    chunk_model = make_chunk_model(embedder.dim)

    db = lancedb.connect(str(vector_db_dir / "lancedb"))
    table_names = lancedb_table_names(db)
    if args.rebuild and args.table_name in table_names:
        logger.warning("Dropping existing table %r per --rebuild", args.table_name)
        db.drop_table(args.table_name)
        table_names = lancedb_table_names(db)

    if args.table_name in table_names:
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

    census = row_census(table)
    # Dimension alone cannot detect two incompatible models of the same size.
    # A partial model migration would make nearest-neighbour search meaningless.
    for pdf_hash in census:
        recorded_model = read_marker(vector_db_dir / f"{pdf_hash}_embedded.done").get("embedding_model")
        if recorded_model and recorded_model != embedder.model_name:
            logger.error("Index contains model %s; use --rebuild to switch to %s",
                         recorded_model, embedder.model_name)
            return 2
    logger.info("Found %d document(s) to process", len(hash_dirs))

    n_ok = n_skip = n_fail = 0
    for hash_dir in hash_dirs:
        pdf_hash = hash_dir.name
        marker = vector_db_dir / f"{pdf_hash}_embedded.done"

        try:
            if args.resume and not args.rebuild and marker_problem(
                    hash_dir, read_marker(marker), census,
                    model=embedder.model_name, dim=embedder.dim,
                    table_name=args.table_name) is None:
                logger.info("Skipping %s (verified current with %s)", pdf_hash, embedder.model_name)
                n_skip += 1
                continue
            embed_document(hash_dir, table, chunk_model, embedder, marker_file=marker)
        except EmbeddingError as e:
            logger.error("%s: embedding failed: %s — skipping", pdf_hash, e)
            n_fail += 1
            continue
        except Exception as e:
            logger.exception("%s: unexpected error: %s — skipping", pdf_hash, e)
            n_fail += 1
            continue

        n_ok += 1

    logger.info("Done. embedded=%d, skipped=%d, failed=%d", n_ok, n_skip, n_fail)
    return 0 if n_fail == 0 else 3


if __name__ == "__main__":
    sys.exit(main())
