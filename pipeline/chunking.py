"""Text chunking (embedding completion belongs exclusively to Stage 2).

* :func:`chunk_text` — hybrid (docling) chunker on the serialized
  DoclingDocument, with a naive char-window fallback when docling
  isn't available. Output: ``chunks.json`` with one record per chunk
  carrying text, heading trail, section_class, and captions.
"""
from __future__ import annotations

import json
import logging
from itertools import groupby
from pathlib import Path
from typing import List, Optional

from . import stamp_artifact
from .config import CONFIG, classify_section
from .treatment_context import TREATMENT_CONTEXT_POLICY

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
            from .treatment_context import (chunk_source_context, context_merge_key,
                                            materialize_treatment_context)
            from .table_structure import SourceTableSerializerProvider, table_chunk_metadata
            context_by_ref = materialize_treatment_context(dl_doc)
            # First obtain item-sized chunks, then apply the pinned Docling
            # token-aware merge within source treatment/section boundaries.
            chunker = HybridChunker(merge_peers=False, serializer_provider=SourceTableSerializerProvider())
            # HybridChunker deliberately over-feeds its tokenizer while
            # measuring where to split, so transformers prints
            # "Token indices sequence length is longer than the specified
            # maximum sequence length for this model (9926 > 512).
            # Running this sequence through the model will result in
            # indexing errors" — straight to the console, with no
            # timestamp or module prefix, unlike every other line in
            # run.log. It is harmless and expected, and it reads like a
            # crash. Quiet it for the duration of the chunk walk only.
            _tok_log = logging.getLogger("transformers.tokenization_utils_base")
            _prev_level = _tok_log.level
            _tok_log.setLevel(logging.ERROR)
            try:
                item_chunks = list(chunker.chunk(dl_doc=dl_doc))
                chunk_iter = []
                for _, group in groupby(item_chunks, key=lambda c: context_merge_key(c, context_by_ref)):
                    chunk_iter.extend(chunker._merge_chunks_with_matching_metadata(list(group)))
            finally:
                _tok_log.setLevel(_prev_level)
            for i, c in enumerate(chunk_iter):
                headings = list(getattr(c.meta, "headings", []) or [])
                captions = list(getattr(c.meta, "captions", []) or [])
                context = chunk_source_context(c.meta.doc_items, context_by_ref)
                table_context = table_chunk_metadata(c.meta.doc_items, c.text, dl_doc)
                if table_context:
                    context["tables"] = table_context
                section_class = classify_section(headings)
                if c.meta.doc_items and all(getattr(item.label, "value", str(item.label)) in
                                             {"caption", "picture", "table"} for item in c.meta.doc_items):
                    section_class = None
                if context["section_type"] == "diagnosis" or (context["section_type"] or "").startswith("description"):
                    section_class = "description"
                chunks.append(
                    {
                        "chunk_id": f"chunk_{i}",
                        "text": c.text,
                        "headings": headings,
                        "section_class": section_class,
                        "captions": captions,
                        **context,
                    }
                )
            logger.info("HybridChunker produced %d chunks", len(chunks))
        except Exception as e:
            # This branch means docling gave us a real document and the
            # *chunker* failed — usually its tokenizer
            # (sentence-transformers/all-MiniLM-L6-v2) not being in the
            # HuggingFace cache under HF_HUB_OFFLINE=1. The naive
            # fallback still produces chunks, so the run exits 0 and
            # nothing downstream complains, while retrieval quality
            # collapses: a 2-page paper chunks to 1 window instead of 16.
            # Log it at ERROR with the remedy — a whole-corpus silent
            # degradation is the failure mode #139 was about.
            logger.error(
                "HybridChunker failed (%s); falling back to the naive char "
                "chunker for this paper. Retrieval quality will be much "
                "worse — chunks stop respecting headings, tables and "
                "captions. If this fires for every paper, the chunker's "
                "tokenizer is missing from the HuggingFace cache: run "
                "`corpus prefetch` on a host with network access.", e,
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
                    "treatment_context": {"status": "unknown", "name": None},
                    "section_type": None,
                    "source_items": [],
                    "start_char": i,
                    "end_char": min(i + chunk_size, len(text)),
                }
            )
        logger.info("Naive chunker produced %d chunks", len(chunks))

    chunks_data = {
        "chunker": chunker_name,
        "treatment_context_policy": TREATMENT_CONTEXT_POLICY,
        "total_chunks": len(chunks),
        "chunks": chunks,
    }

    with open(chunks_output, "w", encoding="utf-8") as f:
        json.dump(stamp_artifact(chunks_data), f, indent=2, ensure_ascii=False)
