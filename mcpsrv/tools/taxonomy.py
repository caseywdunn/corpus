"""Taxonomy-keyed MCP tools.

Surfaces: search_taxon, get_papers_for_taxon, get_chunks_for_taxon,
get_taxon_mentions, list_valid_species_under, get_papers_by_author.
The first four resolve queries through synonymy via the configured
DwC taxonomy snapshot — callers never need to know whether a given
form is the accepted name, an unaccepted one, or a registered synonym.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional

from ..app import _load_json, _need_index, mcp


# Helper used by the get_taxon_mentions fallback path.
_CHUNK_INDEX_RE_MCP = re.compile(r"chunk_(\d+)")


def parse_chunk_index(chunk_id: str) -> int:
    m = _CHUNK_INDEX_RE_MCP.search(chunk_id or "")
    return int(m.group(1)) if m else -1


@mcp.tool()
def search_taxon(name: str) -> Dict:
    """Resolve a taxon name against the configured DwC taxonomy snapshot.

    Accepts any form the snapshot knows: accepted names, historical
    synonyms, common misspellings registered as ``names`` rows. Returns
    the accepted DwC ``taxonID``, scientific name, authorship, rank, and
    the ``name_type`` of the match (``accepted`` | ``unaccepted`` |
    ``synonym``). Missing names return a structured ``not_found`` result
    rather than an error.
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return {"error": "no taxonomy snapshot configured"}
    hit = idx.taxonomy_db.lookup(name)
    if not hit:
        return {"not_found": True, "queried": name}
    in_corpus = hit["accepted_taxon_id"] in idx.taxon_to_papers
    return {
        **hit,
        "queried": name,
        "in_corpus": in_corpus,
        "mentioning_paper_count": len(idx.taxon_to_papers.get(hit["accepted_taxon_id"], [])),
    }


@mcp.tool()
def get_papers_for_taxon(
    taxon_name: str,
    *,
    min_mentions: int = 1,
) -> List[Dict]:
    """List papers mentioning a taxon, resolved through synonymy.

    Papers are returned ordered by mention count (desc). Synonymy is
    followed to the accepted ``taxonID`` — so a query for
    'Stephanomia amphitridis' returns papers citing its accepted synonym
    *Apolemia uvaria*, without the caller needing to know the synonymy.
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return [{"error": "no taxonomy snapshot configured"}]
    hit = idx.taxonomy_db.lookup(taxon_name)
    if not hit:
        return []
    aid = hit["accepted_taxon_id"]
    hashes = idx.taxon_to_papers.get(aid, [])
    counts = idx.taxon_mention_counts.get(aid, {})
    ordered = sorted(set(hashes), key=lambda h: -counts.get(h, 0))
    out = []
    for h in ordered:
        c = counts.get(h, 0)
        if c < min_mentions:
            continue
        p = idx.papers.get(h)
        if not p:
            continue
        first_author = ""
        if p.get("authors"):
            a0 = p["authors"][0]
            first_author = f"{a0.get('forename','').strip()} {a0.get('surname','').strip()}".strip()
        out.append({
            "hash": h,
            "title": p.get("title"),
            "year": p.get("year"),
            "first_author": first_author,
            "mention_count": c,
        })
    return out


@mcp.tool()
def get_chunks_for_taxon(
    taxon_name: str,
    paper_hash: Optional[str] = None,
    limit: int = 200,
    offset: int = 0,
) -> List[Dict]:
    """Every chunk that mentions the taxon (resolved through synonymy).

    By default scans the whole corpus. Restrict to a single paper with
    ``paper_hash``. Results include the full chunk text so the LLM can
    synthesize directly — don't lower ``limit`` unless you know the
    query's answer fits in fewer chunks; enumerative queries (§8 Q2,
    monographic review of a genus) need the full set.
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return [{"error": "no taxonomy snapshot configured"}]
    hit = idx.taxonomy_db.lookup(taxon_name)
    if not hit:
        return []
    aid = hit["accepted_taxon_id"]
    target_hashes = (
        [paper_hash] if paper_hash else idx.taxon_to_papers.get(aid, [])
    )

    out: List[Dict] = []
    for h in target_hashes:
        p = idx.papers.get(h)
        if not p:
            continue
        taxa = _load_json(Path(p["hash_dir"]) / "taxa.json", default={}) or {}
        chunks_data = _load_json(Path(p["hash_dir"]) / "chunks.json", default={}) or {}
        chunks_by_id = {c.get("chunk_id"): c for c in chunks_data.get("chunks", []) or []}
        matching_chunk_ids: List[str] = []
        # A single chunk may have multiple mentions of the taxon; dedup
        # so each chunk appears once per paper.
        seen = set()
        for m in taxa.get("mentions", []) or []:
            if m.get("accepted_taxon_id") != aid:
                continue
            cid = m.get("chunk_id")
            if cid and cid not in seen:
                matching_chunk_ids.append(cid)
                seen.add(cid)

        for cid in matching_chunk_ids:
            ch = chunks_by_id.get(cid) or {}
            out.append({
                "paper_hash": h,
                "paper_title": p.get("title"),
                "paper_year": p.get("year"),
                "chunk_id": cid,
                "section_class": ch.get("section_class"),
                "headings": ch.get("headings") or [],
                "text": ch.get("text"),
            })

    return out[offset: offset + int(limit)] if limit else out[offset:]


@mcp.tool()
def get_taxon_mentions(
    taxon_name: str,
    paper_hash: Optional[str] = None,
    limit: int = 500,
    offset: int = 0,
) -> List[Dict]:
    """All text-span mentions of a taxon across the corpus, resolved
    through synonymy.

    Returns individual mention records with character offsets into the
    chunk text — finer-grained than ``get_chunks_for_taxon`` which
    returns whole chunks. Each record includes ``char_start``,
    ``char_end``, ``mention_text``, and the chunk/paper context.

    Requires the taxon mention database (built by
    ``build_taxon_mentions.py``). Falls back to per-paper
    ``taxa.json`` scanning if the database is not available.
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return [{"error": "no taxonomy snapshot configured"}]
    hit = idx.taxonomy_db.lookup(taxon_name)
    if not hit:
        return []
    aid = hit["accepted_taxon_id"]

    if idx.taxon_mention_db is not None:
        # Fast path: query the corpus-wide SQLite
        rows = idx.taxon_mention_db.mentions_for_taxon(
            aid, corpus_hash=paper_hash, limit=limit, offset=offset,
        )
        # Enrich with paper-level metadata
        out = []
        for r in rows:
            p = idx.papers.get(r["corpus_hash"], {})
            out.append({
                "paper_hash": r["corpus_hash"],
                "paper_title": p.get("title"),
                "paper_year": p.get("year"),
                "chunk_id": r["chunk_id"],
                "chunk_index": r["chunk_index"],
                "char_start": r["char_start"],
                "char_end": r["char_end"],
                "mention_text": r["mention_text"],
                "matched_name": r["matched_name"],
                "accepted_name": r["accepted_name"],
                # SQLite column is ``taxon_rank`` (build_taxon_mentions.py);
                # output API key stays ``rank`` for client compatibility.
                "rank": r["taxon_rank"],
                "name_type": r["name_type"],
                "method": r["method"],
            })
        return out

    # Fallback: scan per-paper taxa.json files (slower, same result shape)
    target_hashes = (
        [paper_hash] if paper_hash else idx.taxon_to_papers.get(aid, [])
    )
    out: List[Dict] = []
    for h in target_hashes:
        p = idx.papers.get(h)
        if not p:
            continue
        taxa = _load_json(Path(p["hash_dir"]) / "taxa.json", default={}) or {}
        for m in taxa.get("mentions", []) or []:
            if m.get("accepted_taxon_id") != aid:
                continue
            span = m.get("text_span", [0, 0])
            chunk_id = m.get("chunk_id", "")
            out.append({
                "paper_hash": h,
                "paper_title": p.get("title"),
                "paper_year": p.get("year"),
                "chunk_id": chunk_id,
                "chunk_index": parse_chunk_index(chunk_id),
                "char_start": span[0] if len(span) > 0 else 0,
                "char_end": span[1] if len(span) > 1 else 0,
                "mention_text": m.get("matched_text", ""),
                "matched_name": m.get("matched_text", ""),
                "accepted_name": m.get("accepted_name", ""),
                "rank": m.get("rank", ""),
                "name_type": m.get("name_type", ""),
                "method": "regex_taxonomy",
            })
    return out[offset: offset + int(limit)] if limit else out[offset:]

@mcp.tool()
def get_papers_by_author(surname: str) -> List[Dict]:
    """Papers authored by the given surname (case-insensitive).

    Matches the parsed Grobid authors in each paper's metadata.json.
    Useful for "summarize all of author X's comments about Y"
    (§8 Q5) — filter papers here then hand the list to
    get_chunks_for_taxon or get_chunks_by_section.
    """
    idx = _need_index()
    key = (surname or "").strip().lower()
    if not key:
        return []
    hashes = idx.author_to_papers.get(key, [])
    out: List[Dict] = []
    for h in hashes:
        p = idx.papers.get(h)
        if not p:
            continue
        all_authors = ", ".join(
            f"{a.get('forename','').strip()} {a.get('surname','').strip()}".strip()
            for a in p.get("authors", []) or []
        )
        out.append({
            "hash": h,
            "title": p.get("title"),
            "year": p.get("year"),
            "authors": all_authors,
            "filename": p.get("filename"),
        })
    return out


@mcp.tool()
def list_valid_species_under(parent_taxon_name: str) -> List[Dict]:
    """All currently-valid species descending from the given taxon in
    the configured Darwin Core taxonomy snapshot.

    Accepts any rank above species (genus, family, order, …). The result
    is a filtered view of the taxonomy snapshot; it does not consult the
    corpus — pair with :func:`get_papers_for_taxon` for per-species
    corpus coverage.
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return [{"error": "no taxonomy snapshot configured"}]
    hit = idx.taxonomy_db.lookup(parent_taxon_name)
    if not hit:
        return []
    parent_id = hit["accepted_taxon_id"]
    # Walk the parent_name_usage_id tree in the snapshot. BFS — the tree
    # may not be strictly shallow. Only accepted species/subspecies are
    # returned, matching the DwC taxonomicStatus convention.
    conn = idx.taxonomy_db.conn
    frontier = [parent_id]
    descendants: List[str] = []
    seen: set = set()
    while frontier:
        parent = frontier.pop(0)
        if parent in seen:
            continue
        seen.add(parent)
        cur = conn.execute(
            "SELECT taxon_id, taxon_rank, taxonomic_status FROM taxa "
            "WHERE parent_name_usage_id = ?",
            (parent,),
        )
        for row in cur:
            descendants.append(row[0])
            frontier.append(row[0])
    if not descendants:
        return []
    placeholders = ",".join("?" * len(descendants))
    cur = conn.execute(
        f"""
        SELECT taxon_id, scientific_name, scientific_name_authorship, taxon_rank
        FROM taxa
        WHERE taxon_id IN ({placeholders})
          AND taxonomic_status = 'accepted'
          AND lower(taxon_rank) IN ('species', 'subspecies')
        ORDER BY scientific_name
        """,
        descendants,
    )
    out: List[Dict] = []
    for row in cur:
        tid = row[0]
        out.append({
            "accepted_taxon_id": tid,
            "accepted_name": row[1],
            "authorship": row[2],
            "rank": row[3],
            "mentioning_paper_count": len(idx.taxon_to_papers.get(tid, [])),
        })
    return out


