"""Paper-level MCP tools: bundle_info, list_papers, get_paper, get_chunk."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

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

