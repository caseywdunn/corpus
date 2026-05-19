"""Lexicon-centric dossier tools (#76, priority 4/6).

Two tools, both category-agnostic per #76 principle 7 — the bundle
declares its categories, the server reflects them, no name is
hardwired:

- ``lexicon_matrix(category, terms, top_n, paper_hashes, year_from,
  year_to)``: tabular paper × term mention-count grid. Lets an LLM
  render a coverage section (p01) or a keyword-sanity table (p13)
  without reading any chunk text. Typical ~2 k tokens for 100 × 20.

- ``get_lexicon_term_dossier(category, term, max_papers,
  max_chunk_examples)``: per-term cross-corpus view. Replaces the
  "for every paper, fetch top_lexicon_terms and scan for term X"
  pattern in p10 / p13.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from ..app import _load_json, _need_index, mcp


# #76 — bounded budgets for typical ~1–5 k-token responses.
_MATRIX_TOP_N_DEFAULT = 20
_TERM_DOSSIER_MAX_PAPERS_DEFAULT = 50
_TERM_DOSSIER_MAX_CHUNK_EXAMPLES_DEFAULT = 10


def _category_or_error(idx, category: str) -> Optional[Dict[str, Any]]:
    """Validate ``category`` against the bundle's declared list.
    Returns an error dict if invalid, otherwise ``None``."""
    available = sorted(idx.lexicon_to_papers.keys())
    if category not in available:
        return {
            "error": "unknown_category",
            "queried_category": category,
            "available": available,
        }
    return None


def _top_terms_by_mention_count(
    idx, category: str, top_n: int,
) -> List[str]:
    """Pick the most-mentioned terms in the category, summed across
    all papers. Used by lexicon_matrix when the caller doesn't pass
    an explicit ``terms`` list."""
    mention_counts = idx.lexicon_mention_counts.get(category, {})
    totals = [
        (term, sum(per_paper.values()))
        for term, per_paper in mention_counts.items()
    ]
    totals.sort(key=lambda t: (-t[1], t[0]))
    return [t for t, _ in totals[: top_n]]


def _paper_filter(
    idx, paper_hashes: Optional[List[str]],
    year_from: Optional[int], year_to: Optional[int],
) -> List[str]:
    """Resolve the row set for lexicon_matrix. ``paper_hashes`` wins
    if provided; otherwise everything in the year range. Skips
    unknown hashes silently — the caller may be passing a stale
    cached list."""
    if paper_hashes:
        candidates = [h for h in paper_hashes if h in idx.papers]
    else:
        candidates = list(idx.papers.keys())
    if year_from is not None or year_to is not None:
        def _in_range(h: str) -> bool:
            y = idx.papers[h].get("year")
            if y is None:
                return False
            if year_from is not None and y < year_from:
                return False
            if year_to is not None and y > year_to:
                return False
            return True
        candidates = [h for h in candidates if _in_range(h)]
    candidates.sort(key=lambda h: (
        -(idx.papers[h].get("year") or 0), h,
    ))
    return candidates


@mcp.tool()
def lexicon_matrix(
    category: str,
    terms: Optional[List[str]] = None,
    top_n: int = _MATRIX_TOP_N_DEFAULT,
    paper_hashes: Optional[List[str]] = None,
    year_from: Optional[int] = None,
    year_to: Optional[int] = None,
) -> Dict[str, Any]:
    """Paper × term mention-count grid for one lexicon category (#76).

    Replaces ad-hoc rendering of "lexicon coverage" sections (p01) and
    keyword-sanity tables (p13's TaxoDros check) — instead of fetching
    each paper's top_lexicon_terms and joining client-side, the server
    emits the table directly. Typical payload ~2 k tokens for a
    100-paper × 20-term grid.

    Columns: pass ``terms=[...]`` for a specific column set, or
    accept the default — the ``top_n`` most-mentioned terms across
    the whole corpus, sorted by total mention count.

    Rows: pass ``paper_hashes=[...]`` for a specific row set, or
    accept the default — every paper, optionally filtered by
    ``year_from`` / ``year_to`` (inclusive). Rows are sorted by year
    desc, then hash for stable ordering.

    Cells: integer mention counts read from each paper's
    ``<category>.json`` ``terms[].mention_count`` field. Papers with
    no ``<category>.json`` contribute a zero row.

    Category-agnostic. Returns ``{"error": "unknown_category",
    "available": [...]}`` if the bundle doesn't declare it.

    Returned shape::

        {
          "category": str,
          "terms": [canonical, ...],          # column order
          "rows": [{
            "paper_hash": str,
            "title": str,
            "year": int | null,
            "counts": [int, ...],             # parallel to terms
          }, ...]
        }
    """
    idx = _need_index()
    err = _category_or_error(idx, category)
    if err is not None:
        return err

    if terms is None:
        columns = _top_terms_by_mention_count(idx, category, top_n)
    else:
        # Caller-provided columns: preserve their order, dedup,
        # silently drop unknowns. Drift between cached prompts and
        # the bundle's current term set shouldn't crash.
        known = set(idx.lexicon_to_papers.get(category, {}).keys())
        seen: set = set()
        columns = []
        for t in terms:
            if t in known and t not in seen:
                columns.append(t)
                seen.add(t)

    row_hashes = _paper_filter(idx, paper_hashes, year_from, year_to)

    # Per-paper, read the <category>.json once and build a vector
    # parallel to columns.
    rows: List[Dict[str, Any]] = []
    column_index = {t: i for i, t in enumerate(columns)}
    for h in row_hashes:
        p = idx.papers[h]
        payload = _load_json(
            Path(p["hash_dir"]) / f"{category}.json", default=None,
        )
        counts = [0] * len(columns)
        if isinstance(payload, dict):
            for term_entry in payload.get("terms", []) or []:
                canonical = term_entry.get("canonical")
                if canonical in column_index:
                    counts[column_index[canonical]] = (
                        term_entry.get("mention_count", 0) or 0
                    )
        rows.append({
            "paper_hash": h,
            "title": p.get("title"),
            "year": p.get("year"),
            "counts": counts,
        })

    return {
        "category": category,
        "terms": columns,
        "rows": rows,
    }


@mcp.tool()
def get_lexicon_term_dossier(
    category: str,
    term: str,
    max_papers: int = _TERM_DOSSIER_MAX_PAPERS_DEFAULT,
    max_chunk_examples: int = _TERM_DOSSIER_MAX_CHUNK_EXAMPLES_DEFAULT,
) -> Dict[str, Any]:
    """Per-term cross-corpus dossier (#76).

    Replaces the "for every paper, fetch top_lexicon_terms and scan
    for term X" pattern in p10 / p13. Single call, ~1–3 k tokens
    typical.

    Pair with ``get_chunks(paper_hash, chunk_ids=[...])`` for full
    text from the ``chunk_examples`` list.

    Category-agnostic. Returns ``{"error": "unknown_category",
    "available": [...]}`` if the bundle doesn't declare the category,
    or ``{"error": "unknown_term", "category", "queried_term"}`` if
    the term isn't hit anywhere in that category.

    Returned shape::

        {
          "category": str,
          "term": str,
          "description": str,        # from lexicon definition; "" if absent
          "n_papers": int,
          "n_mentions_total": int,
          "papers": [{"hash", "title", "year", "n_mentions"}, ...],
          "chunk_examples": [{
            "paper_hash", "chunk_id", "section_class", "headings"
          }, ...]                    # pair with get_chunks for full text
        }
    """
    idx = _need_index()
    err = _category_or_error(idx, category)
    if err is not None:
        return err

    paper_hashes_for_term = list(
        idx.lexicon_to_papers.get(category, {}).get(term, [])
    )
    if not paper_hashes_for_term:
        return {
            "error": "unknown_term",
            "category": category,
            "queried_term": term,
        }

    mention_counts = idx.lexicon_mention_counts.get(category, {}).get(term, {})

    # Per-paper rollups (sorted by mention count desc, hash asc for tie-break).
    paper_stats: List[Dict[str, Any]] = []
    for h in set(paper_hashes_for_term):
        p = idx.papers.get(h)
        if not p:
            continue
        paper_stats.append({
            "hash": h,
            "title": p.get("title"),
            "year": p.get("year"),
            "n_mentions": mention_counts.get(h, 0),
        })
    paper_stats.sort(key=lambda r: (-r["n_mentions"], r["hash"]))
    n_mentions_total = sum(p["n_mentions"] for p in paper_stats)

    # Chunk examples: walk the top-mention papers' <category>.json
    # files, collect chunk_ids where this term appears, enrich with
    # section_class + headings from chunks.json.
    chunk_examples: List[Dict[str, Any]] = []
    description = ""
    for ps in paper_stats[: max_papers]:
        if len(chunk_examples) >= max_chunk_examples:
            break
        p = idx.papers.get(ps["hash"])
        if not p:
            continue
        hash_dir = Path(p["hash_dir"])
        payload = _load_json(hash_dir / f"{category}.json", default=None)
        if not isinstance(payload, dict):
            continue

        # Pluck description from the first paper that carries it.
        if not description:
            for t in payload.get("terms", []) or []:
                if t.get("canonical") == term:
                    desc = t.get("description") or ""
                    if desc:
                        description = desc
                    break

        chunks_data = _load_json(hash_dir / "chunks.json", default=None) or {}
        chunks_by_id = {
            c.get("chunk_id"): c
            for c in chunks_data.get("chunks", []) or []
        }
        seen_chunks: set = set()
        for m in payload.get("mentions", []) or []:
            if m.get("canonical") != term:
                continue
            cid = m.get("chunk_id")
            if not cid or cid in seen_chunks:
                continue
            seen_chunks.add(cid)
            ch = chunks_by_id.get(cid) or {}
            chunk_examples.append({
                "paper_hash": ps["hash"],
                "chunk_id": cid,
                "section_class": ch.get("section_class"),
                "headings": ch.get("headings") or [],
            })
            if len(chunk_examples) >= max_chunk_examples:
                break

    return {
        "category": category,
        "term": term,
        "description": description,
        "n_papers": len(paper_stats),
        "n_mentions_total": n_mentions_total,
        "papers": paper_stats[: max_papers],
        "chunk_examples": chunk_examples,
    }
