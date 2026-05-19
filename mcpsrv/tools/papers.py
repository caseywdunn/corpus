"""Paper-level MCP tools: bundle_info, corpus_summary, list_papers, get_paper, get_chunk."""
from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional

from pipeline.version import __version__
from ..app import _load_json, _need_index, mcp


@mcp.tool()
def bundle_info() -> Dict:
    """Server identity + bundle metadata: version, creation
    timestamp, pipeline SHA, embedding model, paper/chunk/figure
    counts. Returns ``bundle_version: null`` when the server is
    backed by a local build output rather than a packaged bundle.
    """
    idx = _need_index()
    if idx.bundle_manifest is None:
        return {
            "server_name": mcp._mcp_server.name,
            "server_version": __version__,
            "bundle_version": None,
            "note": "no bundle_manifest.json — this server is backed "
                    "by a local build output, not a versioned served "
                    "bundle.  run `python -m mcpsrv.bundle <output_dir> "
                    "<serve_dir> --version vX.Y.Z` to produce one.",
        }
    return {
        "server_name": mcp._mcp_server.name,
        "server_version": __version__,
        **dict(idx.bundle_manifest),
    }


# #76 — defaults sized for the cache-friendly dossier-tools work.
# These bound the typical payload size; the user can override
# top_taxa and per-category top_terms via parameters.
_CORPUS_SUMMARY_TOP_TAXA_DEFAULT = 25
_CORPUS_SUMMARY_TOP_TERMS_PER_CATEGORY_DEFAULT = 15


@mcp.tool()
def corpus_summary(
    top_taxa: int = _CORPUS_SUMMARY_TOP_TAXA_DEFAULT,
    top_terms_per_category: int = _CORPUS_SUMMARY_TOP_TERMS_PER_CATEGORY_DEFAULT,
) -> Dict[str, Any]:
    """One-call corpus orientation: paper counts by decade, lexicon
    coverage per category, top taxa, figure totals, bundle identity
    (#76). Supersedes the N× ``list_papers`` pagination pattern. Pair
    with ``list_papers`` for per-paper detail or ``get_taxon_dossier``
    for one-taxon drill-down. Fixed-shape ~2–5 k tokens.

    Counts sorted desc, ties by name asc. Caps: ``top_taxa`` (default
    25), ``top_terms_per_category`` (default 15); 0 omits the list.

    Returns ``{n_papers, year_range, by_decade, lexicon_categories,
    lexicon_coverage: {category: {n_terms_hit, n_papers_with_hits,
    n_mentions_total, top_terms: [{term, n_papers, n_mentions}]}},
    n_figures_total, n_unique_taxa, top_taxa: [{taxon_id, name,
    rank?, n_papers, n_mentions}], bundle_version, server_version}``.
    """
    idx = _need_index()

    # Year coverage. None years are skipped from both min/max and
    # the decade histogram.
    years = [p.get("year") for p in idx.papers.values() if p.get("year")]
    if years:
        year_range: Optional[Dict[str, int]] = {
            "min": min(years), "max": max(years),
        }
        decades = Counter((y // 10) * 10 for y in years)
        by_decade = [
            {"decade": d, "n_papers": n}
            for d, n in sorted(decades.items())
        ]
    else:
        year_range = None
        by_decade = []

    # Lexicon coverage per category. lexicon_to_papers is
    # {category: {term: [paper_hash, ...]}}; lexicon_mention_counts is
    # {category: {term: {paper_hash: n_mentions}}}.
    lexicon_categories = sorted(idx.lexicon_to_papers.keys())
    lexicon_coverage: Dict[str, Dict[str, Any]] = {}
    for category in lexicon_categories:
        terms = idx.lexicon_to_papers[category]
        mention_counts = idx.lexicon_mention_counts.get(category, {})
        papers_with_hits: set = set()
        n_mentions_total = 0
        term_stats: List[Dict[str, Any]] = []
        for term, paper_hashes in terms.items():
            unique_papers = set(paper_hashes)
            papers_with_hits |= unique_papers
            n_mentions = sum(mention_counts.get(term, {}).values())
            n_mentions_total += n_mentions
            term_stats.append({
                "term": term,
                "n_papers": len(unique_papers),
                "n_mentions": n_mentions,
            })
        term_stats.sort(key=lambda t: (-t["n_mentions"], t["term"]))
        lexicon_coverage[category] = {
            "n_terms_hit": len(terms),
            "n_papers_with_hits": len(papers_with_hits),
            "n_mentions_total": n_mentions_total,
            "top_terms": (
                term_stats[: top_terms_per_category]
                if top_terms_per_category > 0 else []
            ),
        }

    # Taxa: pick the N most-mentioned across the corpus.
    taxon_stats: List[Dict[str, Any]] = []
    for aid, paper_hashes in idx.taxon_to_papers.items():
        unique_papers = set(paper_hashes)
        n_mentions = sum(
            idx.taxon_mention_counts.get(aid, {}).values()
        )
        display = idx.taxon_display.get(aid, {})
        taxon_stats.append({
            "taxon_id": aid,
            "name": display.get("accepted_name") or aid,
            "rank": display.get("rank"),
            "n_papers": len(unique_papers),
            "n_mentions": n_mentions,
        })
    taxon_stats.sort(key=lambda t: (-t["n_mentions"], t["name"]))
    top_taxa_list = taxon_stats[: top_taxa] if top_taxa > 0 else []
    # Drop sparse-encoded null rank from each entry.
    for t in top_taxa_list:
        if t.get("rank") is None:
            t.pop("rank", None)

    n_figures_total = sum(
        p.get("n_figures") or 0 for p in idx.papers.values()
    )

    summary: Dict[str, Any] = {
        "n_papers": len(idx.papers),
        "year_range": year_range,
        "by_decade": by_decade,
        "lexicon_categories": lexicon_categories,
        "lexicon_coverage": lexicon_coverage,
        "n_figures_total": n_figures_total,
        "n_unique_taxa": len(idx.taxon_to_papers),
        "top_taxa": top_taxa_list,
        "bundle_version": (
            idx.bundle_manifest.get("bundle_version")
            if idx.bundle_manifest else None
        ),
        "server_version": __version__,
    }
    return summary


@mcp.tool()
def list_papers(
    year_from: Optional[int] = None,
    year_to: Optional[int] = None,
    limit: int = 100,
    offset: int = 0,
) -> List[Dict]:
    """List papers in the corpus, ordered by year. One row per paper:
    hash, title, year, first author, counts of chunks/figures/taxa
    plus per-category lexicon-term counts.

    Filter by year range with ``year_from`` / ``year_to`` (inclusive).
    Defaults to 100 rows; paginate via ``offset``. Call
    ``get_paper(hash)`` for full metadata (authors, abstract, DOI,
    filename, top taxa, top lexicon terms).
    """
    idx = _need_index()
    rows: List[Dict] = []
    for h, p in sorted(idx.papers.items(), key=lambda kv: (kv[1].get("year") or 0, kv[0])):
        year = p.get("year")
        if year_from is not None and (year is None or year < year_from):
            continue
        if year_to is not None and (year is None or year > year_to):
            continue
        first_author = ""
        if p.get("authors"):
            a0 = p["authors"][0]
            first_author = f"{a0.get('forename','').strip()} {a0.get('surname','').strip()}".strip()
        rows.append({
            "hash": h,
            "title": p.get("title"),
            "year": year,
            "first_author": first_author,
            "n_chunks": p.get("n_chunks", 0),
            "n_figures": p.get("n_figures", 0),
            "n_taxa": p.get("n_taxa", 0),
            "n_lexicon_terms": p.get("n_lexicon_terms") or {},
        })
    return rows[offset: offset + int(limit)]


_NON_LEXICON_FILES = {
    "summary.json", "metadata.json", "references.json",
    "taxa.json", "figures.json", "chunks.json",
    "scan_detection.json", "docling_doc.json",
    "intext_citations.json", "pipeline_state.json",
}


@mcp.tool()
def get_paper(paper_hash: str) -> Dict:
    """Full metadata for one paper: title, authors, year, abstract, DOI,
    plus top-10 taxa and the top-10 terms in each configured lexicon
    category (anatomy, biogeography, …) found on this paper."""
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return {"error": f"no such paper_hash: {paper_hash}"}
    hash_dir = Path(p["hash_dir"])
    taxa = _load_json(hash_dir / "taxa.json", default={}) or {}
    top_lexicon_terms: Dict[str, list] = {}
    for child in hash_dir.glob("*.json"):
        if child.name in _NON_LEXICON_FILES:
            continue
        payload = _load_json(child, default=None)
        if isinstance(payload, dict) and "category" in payload:
            top_lexicon_terms[payload["category"]] = (payload.get("terms") or [])[:10]
    return {
        **{k: v for k, v in p.items() if k != "hash_dir"},
        "top_taxa": (taxa.get("taxa") or [])[:10],
        "top_lexicon_terms": top_lexicon_terms,
    }

@mcp.tool()
def get_papers(
    hashes: List[str],
    fields: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
    """Batched ``get_paper`` for many hashes at once (#76).

    Replaces the N×``get_paper`` pattern when a dossier returns a
    paper list and the LLM wants metadata for each. Server-side this
    is a hash-map lookup — already free.

    ``fields`` (optional) whitelists which keys to return per paper.
    Unknown fields are silently dropped; the standard set is
    ``{hash, title, year, authors, doi, journal, abstract, filename,
    n_chunks, n_figures, n_taxa, n_lexicon_terms, top_taxa,
    top_lexicon_terms}``. Pass ``["hash", "title", "year",
    "first_author"]`` for a thin payload.

    Unknown hashes surface as ``{"hash": <h>, "error": "not_found"}``
    so the caller can locate misses without re-issuing per-hash
    requests. Output is in input order so the caller can zip back to
    its prompt context.
    """
    idx = _need_index()
    out: List[Dict[str, Any]] = []
    for h in hashes:
        p = idx.papers.get(h)
        if p is None:
            out.append({"hash": h, "error": "not_found"})
            continue
        # Materialize the full record (mirrors get_paper). The
        # field-whitelist below trims after, so callers can ask for a
        # cheap projection without us having to short-circuit reads.
        hash_dir = Path(p["hash_dir"])
        taxa = _load_json(hash_dir / "taxa.json", default={}) or {}
        top_lexicon_terms: Dict[str, list] = {}
        for child in hash_dir.glob("*.json"):
            if child.name in _NON_LEXICON_FILES:
                continue
            payload = _load_json(child, default=None)
            if isinstance(payload, dict) and "category" in payload:
                top_lexicon_terms[payload["category"]] = (
                    payload.get("terms") or []
                )[:10]
        # Synthesized "first_author" convenience — saves the caller
        # walking authors[0].surname themselves.
        authors = p.get("authors") or []
        first_author = (
            (authors[0].get("surname") or "") if authors else ""
        ) or None

        record: Dict[str, Any] = {
            **{k: v for k, v in p.items() if k != "hash_dir"},
            "first_author": first_author,
            "top_taxa": (taxa.get("taxa") or [])[:10],
            "top_lexicon_terms": top_lexicon_terms,
        }
        if fields is not None:
            record = {k: v for k, v in record.items() if k in set(fields)}
            # Always carry the hash so callers can locate rows.
            record["hash"] = h
        out.append(record)
    return out


@mcp.tool()
def get_chunk(paper_hash: str, chunk_id: str) -> Dict:
    """One chunk's full record: text, headings, section_class,
    figure_refs."""
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return {"error": f"no such paper_hash: {paper_hash}"}
    chunks = _load_json(Path(p["hash_dir"]) / "chunks.json", default={}) or {}
    for c in chunks.get("chunks", []) or []:
        if c.get("chunk_id") == chunk_id:
            return {**c, "paper_hash": paper_hash, "paper_title": p.get("title")}
    return {"error": f"no such chunk_id {chunk_id!r} in paper {paper_hash}"}


# ---------------------------------------------------------------------------
# HTTP transport + bearer-token auth (dev_docs/PLAN.md §10)
# ---------------------------------------------------------------------------

