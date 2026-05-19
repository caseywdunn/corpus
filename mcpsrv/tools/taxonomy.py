"""Taxonomy-keyed MCP tools.

Surfaces: search_taxon, get_papers_for_taxon, get_chunks_for_taxon,
get_taxon_mentions, list_valid_species_under, get_papers_by_author.
The first four resolve queries through synonymy via the configured
DwC taxonomy snapshot — callers never need to know whether a given
form is the accepted name, an unaccepted one, or a registered synonym.
"""
from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..app import _load_json, _need_index, mcp

# Filenames at the per-paper root that are NOT lexicon outputs.
# Mirrors the same list in CorpusIndex.load() so dossier readers
# walking *.json files can distinguish lexicon outputs from the
# fixed pipeline artifacts.
_NON_LEXICON_PER_PAPER_FILES = frozenset({
    "summary.json", "metadata.json", "references.json",
    "taxa.json", "figures.json", "chunks.json",
    "scan_detection.json", "docling_doc.json",
    "intext_citations.json", "pipeline_state.json",
})


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
    limit: int = 50,
    offset: int = 0,
) -> List[Dict]:
    """Every chunk that mentions the taxon (resolved through synonymy).

    By default scans the whole corpus and returns the first 50 chunks
    with full text. Restrict to a single paper with ``paper_hash``, or
    paginate via ``offset`` / ``limit`` for enumerative queries
    (monographic reviews of a genus). Each chunk is roughly 100–200
    tokens, so a 50-chunk response is ~5–10k tokens.
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


# #76 — dossier defaults sized for typical 3–10 k-token responses.
_DOSSIER_MAX_PAPERS_DEFAULT = 50
_DOSSIER_MAX_CHUNKS_DEFAULT = 100
_DOSSIER_MAX_FIGURES_DEFAULT = 25
_DOSSIER_LEXICON_TOP_N_DEFAULT = 10
_DOSSIER_COOCCURRING_TOP_N_DEFAULT = 20
_DOSSIER_SECTIONS = frozenset({
    "papers", "chunks", "figures", "lexicon", "cooccurring_taxa",
})


@mcp.tool()
def get_taxon_dossier(
    taxon_name: str,
    include: Optional[List[str]] = None,
    max_papers: int = _DOSSIER_MAX_PAPERS_DEFAULT,
    max_chunks: int = _DOSSIER_MAX_CHUNKS_DEFAULT,
    max_figures: int = _DOSSIER_MAX_FIGURES_DEFAULT,
    top_n_lexicon: int = _DOSSIER_LEXICON_TOP_N_DEFAULT,
    top_n_cooccurring: int = _DOSSIER_COOCCURRING_TOP_N_DEFAULT,
) -> Dict[str, Any]:
    """One-call comprehensive view of a taxon across the corpus (#76).

    Supersedes the chain ``search_taxon`` + ``get_papers_for_taxon``
    + N× ``get_paper`` + N× ``get_chunks_for_taxon`` +
    ``get_figures_for_taxon`` for enumerative taxon-anatomy and
    monographic prompts (~45 round-trips → 1). Typical payload
    3–10 k tokens.

    Pair with ``get_chunks(paper_hash, chunk_ids=[...])`` for
    selective full-text drill-down — ``chunk_index`` returns IDs +
    section + headings, **not text**. Pair with
    ``get_figure_image(paper_hash, figure_id)`` to pull a figure's
    bytes once the figure_index identifies the relevant one.

    ``include`` (default = all five sections) trims the response:
    one or more of ``papers``, ``chunks``, ``figures``, ``lexicon``,
    ``cooccurring_taxa``. Use it to drop sections an enumerative
    prompt doesn't need (e.g. ``include=["papers", "chunks"]`` skips
    figures + lexicon + cooccurring for a text-only review).

    Per-section caps: ``max_papers`` (default 50, sorted by mention
    count desc), ``max_chunks`` (100, paper-then-chunk order),
    ``max_figures`` (25, caption-relevance-first), ``top_n_lexicon``
    (10 per category), ``top_n_cooccurring`` (20 taxa).

    Returned shape (sections present iff requested via ``include``):

        {
          "taxon": {taxon_id, accepted_name, rank, authority,
                    matched_name?},
          "n_papers_mentioning": int,
          "papers": [{hash, title, year, first_author, n_mentions}, ...],
          "chunk_index": [{paper_hash, chunk_id, section_class,
                           headings}, ...],
          "figure_index": [{paper_hash, figure_id, figure_type,
                            page, figure_number, caption_has_taxon}, ...],
          "lexicon_aggregated": {
            category: [{term, n_papers, n_mentions}, ...],
            ...
          },
          "cooccurring_taxa": [{taxon_id, name, rank?,
                                n_shared_papers}, ...],
        }
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return {"error": "no taxonomy snapshot configured"}
    hit = idx.taxonomy_db.lookup(taxon_name)
    if not hit:
        return {"not_found": True, "queried": taxon_name}

    aid = hit["accepted_taxon_id"]
    requested = (
        _DOSSIER_SECTIONS
        if include is None
        else _DOSSIER_SECTIONS & set(include)
    )

    # Always emit taxon metadata + paper-mention count.
    out: Dict[str, Any] = {
        "taxon": {
            "taxon_id": aid,
            "accepted_name": hit.get("accepted_name"),
            "rank": hit.get("rank"),
            "authority": hit.get("authority"),
        },
        "n_papers_mentioning": len(idx.taxon_to_papers.get(aid, [])),
    }
    # matched_name only when it diverges from accepted_name (synonym hit)
    matched_name = hit.get("matched_name")
    if matched_name and matched_name != hit.get("accepted_name"):
        out["taxon"]["matched_name"] = matched_name
    # Drop null rank / authority for sparse encoding.
    for k in ("rank", "authority"):
        if out["taxon"].get(k) is None:
            out["taxon"].pop(k, None)

    paper_hashes = list(idx.taxon_to_papers.get(aid, []))
    if not paper_hashes:
        return out  # taxon known but no corpus papers mention it

    # Sort papers by mention count desc, ties by year desc.
    paper_hashes.sort(
        key=lambda h: (
            -idx.taxon_mention_counts.get(aid, {}).get(h, 0),
            -(idx.papers.get(h, {}).get("year") or 0),
        ),
    )
    in_scope = paper_hashes[: max_papers]

    if "papers" in requested:
        papers_out: List[Dict[str, Any]] = []
        for h in in_scope:
            p = idx.papers.get(h) or {}
            authors = p.get("authors") or []
            first_author = (
                (authors[0].get("surname") or "") if authors else ""
            )
            papers_out.append({
                "hash": h,
                "title": p.get("title"),
                "year": p.get("year"),
                "first_author": first_author or None,
                "n_mentions": idx.taxon_mention_counts.get(aid, {}).get(h, 0),
            })
        # Sparse-encode null first_author.
        for row in papers_out:
            if row.get("first_author") is None:
                row.pop("first_author", None)
        out["papers"] = papers_out

    if "chunks" in requested:
        chunk_index: List[Dict[str, Any]] = []
        for h in in_scope:
            p = idx.papers.get(h)
            if not p:
                continue
            taxa = _load_json(Path(p["hash_dir"]) / "taxa.json", default={}) or {}
            chunks_data = _load_json(
                Path(p["hash_dir"]) / "chunks.json", default={},
            ) or {}
            chunks_by_id = {
                c.get("chunk_id"): c
                for c in chunks_data.get("chunks", []) or []
            }
            seen: set = set()
            for m in taxa.get("mentions", []) or []:
                if m.get("accepted_taxon_id") != aid:
                    continue
                cid = m.get("chunk_id")
                if not cid or cid in seen:
                    continue
                seen.add(cid)
                ch = chunks_by_id.get(cid) or {}
                chunk_index.append({
                    "paper_hash": h,
                    "chunk_id": cid,
                    "section_class": ch.get("section_class"),
                    "headings": ch.get("headings") or [],
                })
                if len(chunk_index) >= max_chunks:
                    break
            if len(chunk_index) >= max_chunks:
                break
        out["chunk_index"] = chunk_index

    if "figures" in requested:
        from .figures import _REAL_FIGURE_TYPES  # local import: avoid cycle

        accepted_low = (hit.get("accepted_name") or "").lower()
        matched_low = (matched_name or "").lower()
        fig_rows: List[Dict[str, Any]] = []
        for h in in_scope:
            p = idx.papers.get(h)
            if not p:
                continue
            figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
            for f in figs.get("figures", []) or []:
                if f.get("figure_type") not in _REAL_FIGURE_TYPES:
                    continue
                caption = (f.get("caption_text") or f.get("caption") or "").lower()
                hit_in_caption = bool(
                    accepted_low and accepted_low in caption
                ) or bool(matched_low and matched_low in caption)
                fig_rows.append({
                    "paper_hash": h,
                    "figure_id": f.get("figure_id"),
                    "figure_type": f.get("figure_type"),
                    "page": f.get("page"),
                    "figure_number": f.get("figure_number"),
                    "caption_has_taxon": hit_in_caption,
                    "_score": (100 if hit_in_caption else 0)
                    + idx.taxon_mention_counts.get(aid, {}).get(h, 0),
                })
        fig_rows.sort(key=lambda r: -r["_score"])
        for r in fig_rows:
            r.pop("_score", None)
        out["figure_index"] = fig_rows[: max_figures]

    if "lexicon" in requested:
        # Aggregate lexicon term counts across papers in scope by
        # walking per-paper <category>.json files. The corpus-wide
        # lexicon_to_papers index is paper-level but doesn't carry
        # mention counts on a per-paper basis sliced by taxon scope —
        # so we accumulate from disk.
        per_category: Dict[str, Dict[str, Dict[str, int]]] = {}
        for h in in_scope:
            p = idx.papers.get(h)
            if not p:
                continue
            hash_dir = Path(p["hash_dir"])
            for child in hash_dir.glob("*.json"):
                if child.name in _NON_LEXICON_PER_PAPER_FILES:
                    continue
                payload = _load_json(child, default=None)
                if not isinstance(payload, dict):
                    continue
                category = payload.get("category")
                if not category:
                    continue
                cat_acc = per_category.setdefault(category, {})
                for term in payload.get("terms", []) or []:
                    canon = term.get("canonical")
                    if not canon:
                        continue
                    rec = cat_acc.setdefault(canon, {"n_papers": 0, "n_mentions": 0})
                    rec["n_papers"] += 1
                    rec["n_mentions"] += term.get("mention_count", 0) or 0
        lexicon_aggregated: Dict[str, List[Dict[str, Any]]] = {}
        for category, terms in per_category.items():
            ranked = sorted(
                ({"term": t, **rec} for t, rec in terms.items()),
                key=lambda r: (-r["n_mentions"], r["term"]),
            )
            lexicon_aggregated[category] = ranked[: top_n_lexicon]
        out["lexicon_aggregated"] = lexicon_aggregated

    if "cooccurring_taxa" in requested:
        # Counter of taxon_id (other than aid) keyed by how many of
        # the in-scope papers mention each one.
        co_papers: Counter = Counter()
        for h in in_scope:
            p = idx.papers.get(h)
            if not p:
                continue
            taxa = _load_json(Path(p["hash_dir"]) / "taxa.json", default={}) or {}
            seen_here: set = set()
            for t in taxa.get("taxa", []) or []:
                other_aid = t.get("accepted_taxon_id")
                if other_aid and other_aid != aid and other_aid not in seen_here:
                    co_papers[other_aid] += 1
                    seen_here.add(other_aid)
        cooccur_rows: List[Dict[str, Any]] = []
        for other_aid, n_shared in co_papers.most_common(top_n_cooccurring):
            display = idx.taxon_display.get(other_aid, {})
            row: Dict[str, Any] = {
                "taxon_id": other_aid,
                "name": display.get("accepted_name") or other_aid,
                "n_shared_papers": n_shared,
            }
            if display.get("rank"):
                row["rank"] = display.get("rank")
            cooccur_rows.append(row)
        out["cooccurring_taxa"] = cooccur_rows

    return out


@mcp.tool()
def get_taxon_lexicon_slice(
    taxon_name: str,
    category: str,
    top_n: int = 25,
    max_paper_examples: int = 5,
) -> Dict[str, Any]:
    """Lexicon coverage for one taxon under one category, joined at
    the chunk level (#76, priority 3/6).

    Answers "for taxon T, what does the corpus say about it under
    lens C?" — supersedes the implicit
    N× ``get_chunks_for_taxon`` + per-chunk lexicon scanning loop in
    p06 (country names per taxon), p10 (chemosensory terms), p12
    (genetic markers), p11 (life-history traits). Single call,
    ~1–3 k tokens typical.

    The join is **same-chunk co-occurrence**: a term counts only when
    it appears in a chunk where the taxon is also mentioned. This is
    tighter than the paper-level rollups in ``corpus_summary`` /
    ``get_taxon_dossier``'s ``lexicon_aggregated`` — chunk co-occurrence
    is what tells the LLM "this term is *applied to* this taxon" vs
    "this term appears somewhere in a paper that also mentions this
    taxon."

    Category-agnostic: returns ``{"error": "unknown_category",
    "available": [...]}`` if the bundle doesn't declare it.

    Returned shape (terms sorted by n_chunks desc, then term asc):

        {
          "taxon": {taxon_id, accepted_name, rank?, authority?},
          "category": str,
          "n_chunks_total": int,        # chunks where taxon mentioned
          "n_papers_total": int,        # papers where taxon mentioned
          "terms": [{
            "term": str,
            "n_chunks": int,            # chunks with taxon + term
            "n_papers": int,
            "paper_examples": [hash, ...],  # up to max_paper_examples
          }, ...]
        }
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return {"error": "no taxonomy snapshot configured"}
    hit = idx.taxonomy_db.lookup(taxon_name)
    if not hit:
        return {"not_found": True, "queried": taxon_name}

    available = sorted(idx.lexicon_to_papers.keys())
    if category not in available:
        return {
            "error": "unknown_category",
            "queried_category": category,
            "available": available,
        }

    aid = hit["accepted_taxon_id"]
    paper_hashes = list(idx.taxon_to_papers.get(aid, []))

    taxon_block: Dict[str, Any] = {
        "taxon_id": aid,
        "accepted_name": hit.get("accepted_name"),
    }
    if hit.get("rank"):
        taxon_block["rank"] = hit["rank"]
    if hit.get("authority"):
        taxon_block["authority"] = hit["authority"]

    if not paper_hashes:
        return {
            "taxon": taxon_block,
            "category": category,
            "n_chunks_total": 0,
            "n_papers_total": 0,
            "terms": [],
        }

    # Walk papers: pair the taxon's chunk_ids per paper with the
    # category's chunk_ids per paper, intersect.
    n_chunks_total = 0
    n_papers_with_taxon = 0
    # canonical → {"chunks": set of (paper_hash, chunk_id),
    #              "papers": set of paper_hash}
    term_aggregate: Dict[str, Dict[str, set]] = {}

    for h in paper_hashes:
        p = idx.papers.get(h)
        if not p:
            continue
        hash_dir = Path(p["hash_dir"])
        taxa = _load_json(hash_dir / "taxa.json", default={}) or {}
        # Set of chunk_ids in this paper where the taxon is mentioned.
        taxon_chunks: set = set()
        for m in taxa.get("mentions", []) or []:
            if m.get("accepted_taxon_id") == aid:
                cid = m.get("chunk_id")
                if cid:
                    taxon_chunks.add(cid)
        if not taxon_chunks:
            continue
        n_chunks_total += len(taxon_chunks)
        n_papers_with_taxon += 1

        lex_payload = _load_json(
            hash_dir / f"{category}.json", default=None,
        )
        if not isinstance(lex_payload, dict):
            continue
        for m in lex_payload.get("mentions", []) or []:
            cid = m.get("chunk_id")
            canonical = m.get("canonical")
            if not cid or not canonical or cid not in taxon_chunks:
                continue
            rec = term_aggregate.setdefault(
                canonical, {"chunks": set(), "papers": set()},
            )
            rec["chunks"].add((h, cid))
            rec["papers"].add(h)

    terms_out = sorted(
        (
            {
                "term": canonical,
                "n_chunks": len(rec["chunks"]),
                "n_papers": len(rec["papers"]),
                "paper_examples": sorted(rec["papers"])[: max_paper_examples],
            }
            for canonical, rec in term_aggregate.items()
        ),
        key=lambda r: (-r["n_chunks"], r["term"]),
    )

    return {
        "taxon": taxon_block,
        "category": category,
        "n_chunks_total": n_chunks_total,
        "n_papers_total": n_papers_with_taxon,
        "terms": terms_out[: top_n],
    }


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


