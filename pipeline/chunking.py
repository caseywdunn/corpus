"""Text chunking + vector-db ingestion marker.

* :func:`chunk_text` — hybrid (docling) chunker on the serialized
  DoclingDocument, with a naive char-window fallback when docling
  isn't available. Output: ``chunks.json`` with one record per chunk
  carrying text, heading trail, section_class, and captions.
* :func:`ingest_to_vector_db` — writes a per-paper completion marker
  used by Stage 2's resume logic.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Optional

from . import stamp_artifact
from .config import CONFIG, classify_section

logger = logging.getLogger(__name__)


def chunk_text(
    text_file: Path,
    metadata_file: Optional[Path] = None,
    chunks_output: Optional[Path] = None,
    docling_doc_file: Optional[Path] = None,
):
    """Chunk the document into structurally-aware pieces.

    When a serialized :class:`DoclingDocument` is available (``docling_doc.json``
    produced by :func:`extract_docling_content`), we drive docling's
    :class:`HybridChunker` — tokenizer-aware, respects headings, table and
    figure captions. Each chunk carries its heading trail and a derived
    ``section_class`` (see :func:`classify_section`).

    When no DoclingDocument is available (e.g., docling import failed
    upstream), we fall back to the prior naive character-window behavior so
    downstream stages still run.

    ``metadata_file`` is accepted for backward compatibility with callers
    that still pass it but is unused — chunking has no dependency on
    Grobid output, decoupling Stage 1 from the metadata stage (#28).
    """
    if chunks_output is None:
        raise TypeError("chunk_text: chunks_output is required")
    del metadata_file  # explicitly unused; preserved in signature for compat

    # Resolve default docling_doc_file relative to text_file's directory.
    if docling_doc_file is None:
        docling_doc_file = text_file.parent / "docling_doc.json"

    chunks: List[dict] = []
    chunker_name = "hybrid_chunker"

    if docling_doc_file.exists() and docling_doc_file.stat().st_size > 0:
        try:
            from docling.chunking import HybridChunker
            from docling_core.types.doc import DoclingDocument

            dl_doc = DoclingDocument.load_from_json(docling_doc_file)
            chunker = HybridChunker()
            for i, c in enumerate(chunker.chunk(dl_doc=dl_doc)):
                headings = list(getattr(c.meta, "headings", []) or [])
                captions = list(getattr(c.meta, "captions", []) or [])
                chunks.append(
                    {
                        "chunk_id": f"chunk_{i}",
                        "text": c.text,
                        "headings": headings,
                        "section_class": classify_section(headings),
                        "captions": captions,
                    }
                )
            logger.info("HybridChunker produced %d chunks", len(chunks))
        except Exception as e:
            logger.warning(
                "HybridChunker failed (%s); falling back to naive char chunker", e
            )
            chunks = []
            chunker_name = "naive_char_window"

    if not chunks:
        # Naive fallback.
        chunker_name = "naive_char_window"
        with open(text_file, "r", encoding="utf-8") as f:
            text_data = json.load(f)
        text = text_data.get("text", "")

        chunk_size = int(CONFIG.get("chunking", {}).get("max_tokens", 1000))
        if chunk_size <= 0:
            chunk_size = 1000

        for i in range(0, len(text), chunk_size):
            piece = text[i : i + chunk_size]
            chunks.append(
                {
                    "chunk_id": f"chunk_{len(chunks)}",
                    "text": piece,
                    "headings": [],
                    "section_class": None,
                    "captions": [],
                    "start_char": i,
                    "end_char": min(i + chunk_size, len(text)),
                }
            )
        logger.info("Naive chunker produced %d chunks", len(chunks))

    chunks_data = {
        "chunker": chunker_name,
        "total_chunks": len(chunks),
        "chunks": chunks,
    }

    with open(chunks_output, "w", encoding="utf-8") as f:
        json.dump(stamp_artifact(chunks_data), f, indent=2, ensure_ascii=False)


def ingest_to_vector_db(chunks_file: Path, vector_db_dir: Path, pdf_hash: str):
    """Ingest chunks into vector database."""
    # Placeholder for vector database ingestion
    # This would typically involve embedding the text and storing in a vector database
    
    # Create a simple marker file for now
    ingestion_marker = vector_db_dir / f"{pdf_hash}_embedded.done"
    
    with open(ingestion_marker, 'w') as f:
        json.dump(stamp_artifact({
            "pdf_hash": pdf_hash,
            "chunks_file": str(chunks_file),
            "ingestion_timestamp": str(Path(chunks_file).stat().st_mtime),
            "status": "completed",
        }), f, indent=2)
