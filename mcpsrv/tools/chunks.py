"""Chunk-level MCP tools.

Surfaces: get_chunks_by_section (Grobid-section-typed retrieval),
get_chunks_for_topic (semantic search via the LanceDB vector index),
translate_chunk (Claude-driven on-demand translation, currently the
only LLM-call tool in the surface).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from ..app import _load_json, _need_index, mcp


@mcp.tool()
def get_chunks_by_section(
    paper_hash: str,
    section_class: Optional[str] = None,
    limit: int = 50,
) -> List[Dict]:
    """Chunks of a paper filtered by section class.

    ``section_class`` is one of the canonical values assigned by the
    pipeline: ``abstract``, ``introduction``, ``methods``, ``results``,
    ``description``, ``discussion``, ``conclusion``, ``acknowledgements``,
    ``references``, ``appendix``. Pass ``None`` to return all chunks
    (up to ``limit``).
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return [{"error": f"no such paper_hash: {paper_hash}"}]
    chunks = _load_json(Path(p["hash_dir"]) / "chunks.json", default={}) or {}
    rows: List[Dict] = []
    for c in chunks.get("chunks", []) or []:
        if section_class is not None and c.get("section_class") != section_class:
            continue
        rows.append({
            "paper_hash": paper_hash,
            "chunk_id": c.get("chunk_id"),
            "section_class": c.get("section_class"),
            "headings": c.get("headings") or [],
            "text": c.get("text"),
        })
    return rows[: int(limit)] if limit else rows



@mcp.tool()
def translate_chunk(
    paper_hash: str,
    chunk_id: str,
    target_language: str = "en",
    model: str = "claude-haiku-4-5-20251001",
) -> Dict:
    """Translate one chunk to the target language (default English) via
    the Anthropic Claude API. Cached server-side per chunk so repeats
    are free. If the chunk is already in the target language, returns
    the original with ``translation_needed: False``. Requires
    ``ANTHROPIC_API_KEY``.
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return {"error": f"no such paper_hash: {paper_hash}"}
    hash_dir = Path(p["hash_dir"])

    # Find the chunk text.
    chunks = _load_json(hash_dir / "chunks.json", default={}) or {}
    chunk = next((c for c in chunks.get("chunks", []) or [] if c.get("chunk_id") == chunk_id), None)
    if chunk is None:
        return {"error": f"no such chunk_id {chunk_id!r} in paper {paper_hash}"}
    original = chunk.get("text") or ""

    # If the paper's detected_language already matches, no translation
    # needed. Short-circuit so the user doesn't pay a Claude call to get
    # the same bytes back.
    scan = _load_json(hash_dir / "scan_detection.json", default={}) or {}
    src_lang = scan.get("detected_language")
    if src_lang and src_lang.lower() == target_language.lower():
        return {
            "paper_hash": paper_hash,
            "chunk_id": chunk_id,
            "source_language": src_lang,
            "target_language": target_language,
            "translation_needed": False,
            "original": original,
            "translation": original,
            "cached": False,
        }

    # Check the on-disk cache first.
    cache_file = hash_dir / f"translated_{target_language}.json"
    cache = _load_json(cache_file, default={}) or {}
    if chunk_id in cache:
        hit = cache[chunk_id]
        return {
            "paper_hash": paper_hash,
            "chunk_id": chunk_id,
            "source_language": src_lang,
            "target_language": target_language,
            "translation_needed": True,
            "original": original,
            "translation": hit.get("translation", ""),
            "cached": True,
            "model": hit.get("model"),
        }

    # Translate via Claude.
    try:
        import anthropic
    except ImportError:
        return {"error": "anthropic package not installed (pip install anthropic)"}
    import os
    if not os.environ.get("ANTHROPIC_API_KEY"):
        return {"error": "ANTHROPIC_API_KEY not set in environment"}

    client = anthropic.Anthropic()
    # Short prompt because Claude handles translation well without
    # elaborate instruction. Preserve scientific names (Latin binomials,
    # taxonomic authorities) and quoted terms verbatim — these carry
    # technical meaning that paraphrase would degrade.
    system_prompt = (
        f"Translate the following scientific text into {target_language}. "
        "Preserve any Latin scientific names, taxonomic authorities "
        "(e.g., 'Eschscholtz, 1829'), and terms in quotation marks "
        "exactly as written. Return only the translation — no preface, "
        "no commentary, no source-text repetition."
    )
    try:
        response = client.messages.create(
            model=model,
            max_tokens=max(512, min(4096, len(original))),
            system=system_prompt,
            messages=[{"role": "user", "content": original}],
        )
    except Exception as e:
        return {"error": f"Claude API call failed: {e}"}

    # Extract the text block from the response.
    translated_text = ""
    for block in response.content:
        if getattr(block, "type", None) == "text":
            translated_text += block.text or ""
    translated_text = translated_text.strip()

    # Persist to cache.
    cache[chunk_id] = {
        "translation": translated_text,
        "model": model,
        "source_language": src_lang,
    }
    with cache_file.open("w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2, ensure_ascii=False)

    return {
        "paper_hash": paper_hash,
        "chunk_id": chunk_id,
        "source_language": src_lang,
        "target_language": target_language,
        "translation_needed": True,
        "original": original,
        "translation": translated_text,
        "cached": False,
        "model": model,
    }



@mcp.tool()
def get_chunks_for_topic(
    query: str,
    k: int = 20,
    paper_hash: Optional[str] = None,
) -> List[Dict]:
    """Semantic search over chunks via the LanceDB vector index.

    Returns the top-``k`` chunks most similar to ``query`` by cosine
    similarity. Use this for "how is X discussed in the corpus" style
    questions (dev_docs/PLAN.md §8 Q6, Q8) where the match criterion is
    semantic, not a literal taxon or lexicon term — for those use
    ``get_chunks_for_taxon`` / ``get_figures_for_lexicon_term`` instead,
    which return exhaustive matches.

    Pass ``paper_hash`` to constrain the search to a single paper.

    Returns ``[{error: ...}]`` if no LanceDB index exists yet — run
    ``python embed_chunks.py <output_dir>`` to build one.
    """
    idx = _need_index()
    embedder, table = idx.get_topic_searcher()
    if embedder is None or table is None:
        return [{
            "error": "no LanceDB index available; run "
                     "`python embed_chunks.py <output_dir>` to build one"
        }]
    try:
        qvec = embedder.embed([query])[0]
    except EmbeddingError as e:
        return [{"error": f"embedding the query failed: {e}"}]
    search = table.search(qvec).limit(int(k))
    if paper_hash:
        search = search.where(f"metadata.pdf_hash = '{paper_hash}'")
    results = search.to_list()
    out: List[Dict] = []
    for r in results:
        m = r.get("metadata") or {}
        out.append({
            "paper_hash": m.get("pdf_hash"),
            "paper_title": m.get("title"),
            "paper_year": m.get("year"),
            "chunk_id": m.get("chunk_id"),
            "section_class": m.get("section_class"),
            "headings": m.get("headings") or [],
            "text": r.get("text"),
            "score": r.get("_distance"),  # LanceDB returns cosine distance
        })
    return out


