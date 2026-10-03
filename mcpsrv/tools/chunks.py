"""Chunk-level MCP tools.

Surfaces: get_chunks (batched chunk fetch for one paper),
get_chunks_by_section (stored section/treatment retrieval), and
get_chunks_for_topic (semantic search via the LanceDB vector index).

The server has no LLM-call tools — it is a deterministic retrieval
layer (#124). The former server-side `translate_chunk` was removed
pre-1.0; an MCP client translates retrieved chunk text itself.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from pipeline.embeddings import EmbeddingError

from ..app import _load_json, _need_index, _validate_collection, _validated_limit, error, mcp
from ..chunk_context import ContextProjection


def _validate_context_filters(treatment_name, section_type):
    for name, value in (("treatment_name", treatment_name), ("section_type", section_type)):
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ValueError(f"{name} must be a non-empty string or null")


def _context_is_materialized(chunk):
    context = chunk.get("treatment_context")
    return (isinstance(context, dict) and context.get("status") in {"resolved", "unknown"}
            and "name" in context and "section_type" in chunk
            and isinstance(chunk.get("source_items"), list))


def _context_filter_error(artifact, selected, treatment_name, section_type):
    if treatment_name is None and section_type is None:
        return None
    if not artifact.get("treatment_context_policy") or not all(
        _context_is_materialized(chunk) for chunk in selected
    ):
        return error(
            "This paper lacks materialized treatment/section context. Rebuild its "
            "extraction and chunks, then rebuild annotations, embeddings and the served "
            "bundle before using treatment_name or section_type filters.",
            "rebuild_required", capability="treatment_context",
        )
    return None


def _matches_context(chunk, treatment_name, section_type):
    context = chunk.get("treatment_context") or {}
    if treatment_name is not None and (
        context.get("status") != "resolved" or context.get("name") != treatment_name
    ):
        return False
    return section_type is None or chunk.get("section_type") == section_type



@mcp.tool()
def get_chunks(
    paper_hash: str,
    chunk_ids: Optional[List[str]] = None,
    with_text: bool = True,
    treatment_name: Optional[str] = None,
    section_type: Optional[str] = None,
) -> List[Dict]:
    """Batched chunk fetch — drill-down pair to every ``*_dossier``
    tool (#76). Pick chunk_ids from a dossier's chunk_index, fetch
    text in one call instead of N× ``get_chunk``.

    ``chunk_ids=None`` returns all chunks in paper order. An explicit
    list returns just those (unknown IDs silently skipped).
    ``with_text=False`` omits chunk prose and returns metadata:
    chunk_id, section_class, headings, len_chars, figure_refs, and stored
    source/treatment context. Source prose in context uses bounded preview
    objects (32 characters without text, 96 with text); ``context_projection``
    reports evidence counts and truncation. New context fields have an 8 KiB
    per-row cap and optional evidence arrays share 64 KiB per response.
    Existing row/text semantics are unchanged by these context-only limits.

    Optional ``treatment_name`` matches only an exact stored resolved treatment
    name; it is separate from literal taxon mentions. ``section_type`` matches
    the stored subtype (for example ``diagnosis``). Both filters intersect with
    ``chunk_ids``. Use ``get_chunks_by_section(..., section_type="diagnosis",
    treatment_name=..., limit=...)`` for a bounded discovery route.

    Context filters require a rebuilt artifact and return ``rebuild_required``
    for legacy data. Ordinary legacy fetches remain readable and report
    ``treatment_context.status="unavailable"``. An examined but unassigned
    treatment has ``status="unknown"``. No context is inferred while serving.

    Returns rows with ``chunk_id``, ``section_class``, ``headings``,
    ``figure_refs``, stored context/evidence and ``text`` or ``len_chars``,
    in paper order. ``[{error: ...}]`` on unknown
    paper_hash.
    """
    try:
        _validate_collection(chunk_ids, "chunk_ids")
        _validate_context_filters(treatment_name, section_type)
    except ValueError as exc:
        return [error(str(exc), "invalid_argument")]
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return [error(f"no such paper_hash: {paper_hash}", "not_found")]
    chunks_data = _load_json(Path(p["hash_dir"]) / "chunks.json", default={}) or {}
    all_chunks = chunks_data.get("chunks", []) or []

    if chunk_ids is None:
        selected = all_chunks
    else:
        wanted = set(chunk_ids)
        by_id = {c.get("chunk_id"): c for c in all_chunks}
        # Preserve paper-order rather than caller-order; the caller
        # is usually grabbing a set, not a sequence.
        selected = [by_id[cid] for cid in by_id if cid in wanted]

    context_error = _context_filter_error(chunks_data, selected, treatment_name, section_type)
    if context_error:
        return [context_error]
    out: List[Dict] = []
    context_projection = ContextProjection(with_text=with_text)
    for c in selected:
        if not _matches_context(c, treatment_name, section_type):
            continue
        text = c.get("text") or ""
        row: Dict[str, object] = {
            "chunk_id": c.get("chunk_id"),
            "section_class": c.get("section_class"),
            "headings": c.get("headings") or [],
            "figure_refs": c.get("figure_refs") or [],
            **context_projection.project(c),
        }
        if with_text:
            row["text"] = text
        else:
            row["len_chars"] = len(text)
        out.append(row)
    return out


@mcp.tool()
def get_chunks_by_section(
    paper_hash: str,
    section_class: Optional[str] = None,
    limit: int = 50,
    with_text: bool = True,
    treatment_name: Optional[str] = None,
    section_type: Optional[str] = None,
) -> List[Dict]:
    """Chunks of a paper filtered by section class and stored treatment context.

    ``section_class`` is one of the canonical values assigned by the
    pipeline: ``abstract``, ``introduction``, ``methods``, ``results``,
    ``description``, ``discussion``, ``conclusion``, ``acknowledgements``,
    ``references``, ``appendix``. Pass ``None`` to return all chunks
    (up to ``limit``).

    ``section_type="diagnosis"`` selects materialized diagnostic passages;
    ``treatment_name`` optionally restricts them to an exact stored resolved
    treatment name. Both intersect with the existing ``section_class`` filter.
    All filtering happens before the row limit, in paper order. This retrieves
    name-free treatment text without redefining literal taxon mentions.
    Legacy context filters return ``rebuild_required``; ordinary section-class
    retrieval remains readable with unavailable treatment context.

    Responses include stored ``treatment_context``, ``section_type``, and
    ``source_items``, plus ``text_integrity``, ``tables`` and ``key_branches``
    where the producer materialized them. Unknown context stays explicit.
    Context has an 8 KiB per-row cap; optional evidence arrays share 64 KiB
    per response. ``context_projection`` reports omissions/counts. Source prose
    uses previews (32 characters without text, 96 with text), full character
    counts and original source offsets; full evidence stays in build artifacts.

    ``with_text=False`` (#84) drops the chunk text and adds
    ``len_chars`` — same scan-then-drill-down pattern as
    ``get_chunks_for_topic`` (#82). Pair with
    ``get_chunks(paper_hash, chunk_ids=[...])`` to fetch text for the
    chunks the caller cares about.
    """
    try:
        n = _validated_limit(limit)
        _validate_context_filters(treatment_name, section_type)
    except ValueError as e:
        return [error(str(e), "invalid_argument")]
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return [error(f"no such paper_hash: {paper_hash}", "not_found")]
    chunks = _load_json(Path(p["hash_dir"]) / "chunks.json", default={}) or {}
    candidates = chunks.get("chunks", []) or []
    context_error = _context_filter_error(chunks, candidates, treatment_name, section_type)
    if context_error:
        return [context_error]
    rows: List[Dict] = []
    context_projection = ContextProjection(with_text=with_text)
    for c in candidates:
        if section_class is not None and c.get("section_class") != section_class:
            continue
        if not _matches_context(c, treatment_name, section_type):
            continue
        text = c.get("text") or ""
        row: Dict = {
            "paper_hash": paper_hash,
            "chunk_id": c.get("chunk_id"),
            "section_class": c.get("section_class"),
            "headings": c.get("headings") or [],
            **context_projection.project(c),
        }
        if with_text:
            row["text"] = text
        else:
            row["len_chars"] = len(text)
        rows.append(row)
        if len(rows) >= n:
            break
    return rows



@mcp.tool()
def get_chunks_for_topic(
    query: str,
    k: int = 20,
    paper_hash: Optional[str] = None,
    with_text: bool = True,
    with_cites: bool = False,
) -> List[Dict]:
    """Semantic top-``k`` chunks via the LanceDB cosine index. Use
    for "how is X discussed in the corpus" — for literal taxa /
    lexicon terms use the exhaustive ``get_chunks_for_taxon`` /
    ``get_figures_for_lexicon_term`` instead.

    ``score`` is a cosine **distance**, not a similarity: **lower is
    closer**, rows come back already sorted ascending, and values run
    0–2 (so a ``score`` above 1.0 is normal, not a bug). Do not sort
    descending and do not threshold as if it were a similarity — that
    inverts the ranking.

    ``with_text=False`` (#82) emits a metadata-only shape
    (~80 chars/row vs ~600 with full text) suitable for a
    scan-then-drill-down workflow: scan the result list, pick the
    relevant chunk_ids, fetch full text via
    ``get_chunks(paper_hash, chunk_ids=[...])``. Cuts a typical
    k=10 search's body-text cost by ~45%.

    ``with_cites=True`` (#88) attaches ``cited_work_ids`` to each
    row — the distinct works the chunk's *parent paper* cites in-text,
    resolved via the bibliographic authority (deduplicated per paper).
    Feed them straight to ``format_citations`` to build a reference
    list. Empty list when the authority DB isn't loaded. Citation spans
    are tracked per paper, not per chunk, so every chunk from the same
    paper carries the same list.

    Pass ``paper_hash`` to constrain to one paper. Returns
    ``[{error: ...}]`` if no LanceDB index exists yet — build with
    ``corpus run --only embed``.
    """
    idx = _need_index()
    embedder, table = idx.get_topic_searcher()
    if embedder is None or table is None:
        # #91 — distinguish a legitimately structured-only corpus
        # (return the build-an-index guidance) from an operational
        # fault (raise a hard MCP error so a degraded server refuses to
        # serve silently-empty results that read as "no matches").
        cap = idx.capabilities()["topic_search"]
        if cap["state"] == "degraded":
            raise RuntimeError(
                f"topic search is degraded and refusing to serve: "
                f"{cap['detail']}. Check the server logs and GET /healthz."
            )
        return [error(
            "no LanceDB index available; run "
            "`corpus run --only embed` to build one",
            "not_configured",
        )]
    try:
        qvec = embedder.embed([query])[0]
    except EmbeddingError as e:
        return [error(f"embedding the query failed: {e}", "unavailable")]
    search = table.search(qvec).limit(int(k))
    if paper_hash:
        search = search.where(f"metadata.pdf_hash = '{paper_hash}'")
    results = search.to_list()

    # #88 — per-paper in-text citation targets, loaded once per paper.
    _cites_cache: Dict[str, List[str]] = {}
    biblio = getattr(idx, "biblio_db", None)

    def _cited_work_ids(h: str) -> List[str]:
        if h not in _cites_cache:
            if biblio is None:
                _cites_cache[h] = []
            else:
                rows = biblio.conn.execute(
                    "SELECT DISTINCT cited_work_id FROM citations "
                    "WHERE citing_corpus_hash = ? AND cited_work_id IS NOT NULL "
                    "ORDER BY cited_work_id",
                    (h,),
                ).fetchall()
                _cites_cache[h] = [row[0] for row in rows]
        return _cites_cache[h]

    out: List[Dict] = []
    for r in results:
        m = r.get("metadata") or {}
        h = m.get("pdf_hash")
        # Defense in depth (#54): served bundles drop LanceDB rows for
        # skipped papers at distill time, but a build bundle queried
        # directly might still have them. Filter against idx.papers,
        # which only contains papers whose per-paper artifacts shipped.
        if h and h not in idx.papers:
            continue
        text = r.get("text") or ""
        row: Dict = {
            "paper_hash": h,
            "paper_title": m.get("title"),
            "paper_year": m.get("year"),
            "chunk_id": m.get("chunk_id"),
            "section_class": m.get("section_class"),
            "headings": m.get("headings") or [],
            "score": r.get("_distance"),  # LanceDB returns cosine distance
        }
        if with_text:
            row["text"] = text
        else:
            row["len_chars"] = len(text)
        if with_cites:
            row["cited_work_ids"] = _cited_work_ids(h) if h else []
        out.append(row)
    return out

