#!/usr/bin/env python3
"""MCP server exposing the corpus output as a queryable agent surface.

Phase G (after Phases A–F + filename-year). See PLAN.md §3 "Endpoint: MCP
server over open-format indices" and §8 for the target queries this
serves.

An MCP server is intentionally thin: every tool is a small function over
the JSON/SQLite artifacts written by ``process_corpus.py``. There is no
hidden state, no hidden DB — the server is a view, not a store. Anyone
else can bypass the server entirely and read the same files directly
(DuckDB over parquet, Jupyter, Snakemake, another lab's pipeline).

Startup ingests per-paper JSON files once into in-memory indexes
(~seconds even at corpus scale), so each tool call is fast. Restart the
server to pick up newly-processed papers.

Usage
-----

Run standalone (for local Claude Desktop / Claude Code / Cursor):

    python mcp_server.py /path/to/demo_output

Or with absolute overrides:

    python mcp_server.py /path/to/demo_output \\
        --worms-sqlite /path/to/worms.sqlite

MCP clients typically discover the server through a config file that
tells them to launch this script. See README.md for the Claude Desktop
/ Claude Code snippets.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional

from mcp.server.fastmcp import FastMCP

from taxa import WormsDB

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# In-memory index built once at startup
# ---------------------------------------------------------------------------


def _load_json(path: Path, default: Any = None) -> Any:
    """Best-effort JSON load. Missing or malformed files return ``default``
    — the index is lossy by design so a single broken per-paper file
    doesn't prevent the server from starting."""
    try:
        if path.exists() and path.stat().st_size > 0:
            with path.open(encoding="utf-8") as f:
                return json.load(f)
    except Exception as e:
        logger.warning("Could not read %s: %s", path, e)
    return default


class CorpusIndex:
    """In-memory view over the per-paper artifacts in ``documents/``.

    Holds just enough to answer tool calls quickly without re-scanning
    the whole documents tree every time. For full chunk text / figure
    captions the tool reopens the specific JSON file — cheap at per-
    paper granularity.

    Designed for a single-process stdio MCP server; don't share across
    threads without external locking.
    """

    def __init__(self, output_dir: Path, worms_db: Optional[WormsDB] = None):
        self.output_dir = Path(output_dir).resolve()
        self.documents_dir = self.output_dir / "documents"
        self.worms_db = worms_db

        # Per-paper header, keyed on the 12-char short hash (directory name).
        self.papers: Dict[str, Dict] = {}
        # Reverse indexes — all point to a list of paper hashes.
        self.taxon_to_papers: Dict[int, List[str]] = defaultdict(list)
        self.taxon_mention_counts: Dict[int, Dict[str, int]] = defaultdict(dict)
        self.anatomy_to_papers: Dict[str, List[str]] = defaultdict(list)
        self.anatomy_mention_counts: Dict[str, Dict[str, int]] = defaultdict(dict)
        self.author_to_papers: Dict[str, List[str]] = defaultdict(list)
        # WoRMS accepted_aphia_id → accepted name (for display when we only
        # have IDs from the reverse index).
        self.taxon_display: Dict[int, Dict] = {}

    def load(self) -> int:
        """Populate the indexes from ``documents/*/``. Returns the number
        of papers successfully indexed."""
        if not self.documents_dir.exists():
            raise FileNotFoundError(
                f"documents/ not found under {self.output_dir}. "
                f"Point the server at a processed corpus output directory."
            )

        count = 0
        for hash_dir in sorted(self.documents_dir.iterdir()):
            if not hash_dir.is_dir():
                continue
            summary = _load_json(hash_dir / "summary.json", default={}) or {}
            metadata = _load_json(hash_dir / "metadata.json", default={}) or {}
            taxa = _load_json(hash_dir / "taxa.json", default={"taxa": []}) or {}
            anatomy = _load_json(hash_dir / "anatomy.json", default={"terms": []}) or {}
            figures = _load_json(hash_dir / "figures.json", default={"figures": []}) or {}
            chunks = _load_json(hash_dir / "chunks.json", default={"chunks": []}) or {}

            paper_hash = hash_dir.name
            self.papers[paper_hash] = {
                "hash": paper_hash,
                "hash_full": summary.get("pdf_hash_full"),
                "relative_paths": summary.get("relative_paths", []) or [],
                "title": metadata.get("title"),
                "year": metadata.get("year"),
                "year_source": metadata.get("year_source"),
                "authors": metadata.get("authors", []) or [],
                "doi": metadata.get("doi"),
                "journal": metadata.get("journal"),
                "abstract": metadata.get("abstract"),
                "filename": metadata.get("filename"),
                "hash_dir": str(hash_dir),
                "n_chunks": chunks.get("total_chunks", len(chunks.get("chunks", []) or [])),
                "n_figures": figures.get("total_figures", len(figures.get("figures", []) or [])),
                "n_taxa": taxa.get("unique_taxa", 0),
                "n_anatomy_terms": anatomy.get("unique_terms", 0),
                "scan_file_type": (
                    _load_json(hash_dir / "scan_detection.json", default={}) or {}
                ).get("file_type"),
            }

            for t in taxa.get("taxa", []) or []:
                aid = t.get("accepted_aphia_id")
                if not aid:
                    continue
                self.taxon_to_papers[aid].append(paper_hash)
                self.taxon_mention_counts[aid][paper_hash] = t.get("mention_count", 0)
                self.taxon_display.setdefault(aid, {
                    "accepted_aphia_id": aid,
                    "accepted_name": t.get("accepted_name"),
                    "authority": t.get("authority"),
                    "rank": t.get("rank"),
                })

            for term in anatomy.get("terms", []) or []:
                canonical = term.get("canonical")
                if not canonical:
                    continue
                self.anatomy_to_papers[canonical].append(paper_hash)
                self.anatomy_mention_counts[canonical][paper_hash] = term.get("mention_count", 0)

            for author in metadata.get("authors", []) or []:
                surname = (author.get("surname") or "").strip().lower()
                if surname:
                    self.author_to_papers[surname].append(paper_hash)

            count += 1

        logger.info(
            "Index built: %d papers, %d taxa, %d anatomy terms, %d authors",
            count, len(self.taxon_to_papers),
            len(self.anatomy_to_papers), len(self.author_to_papers),
        )
        return count


# ---------------------------------------------------------------------------
# Module-level singletons populated in main()
# ---------------------------------------------------------------------------

mcp = FastMCP("corpus")
_INDEX: Optional[CorpusIndex] = None


def _need_index() -> CorpusIndex:
    if _INDEX is None:
        raise RuntimeError(
            "CorpusIndex not initialized; start the server via main()"
        )
    return _INDEX


# ---------------------------------------------------------------------------
# Tools
# ---------------------------------------------------------------------------


@mcp.tool()
def list_papers(
    year_from: Optional[int] = None,
    year_to: Optional[int] = None,
    limit: Optional[int] = None,
) -> List[Dict]:
    """List every paper in the corpus with bibliographic + annotation counts.

    Filter by year range with ``year_from`` / ``year_to`` (inclusive).
    Returns at most ``limit`` rows (default: all papers). One row per
    paper with: hash, title, year, first author surname, counts of
    chunks, figures, taxa, anatomy terms.
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
            "year_source": p.get("year_source"),
            "first_author": first_author,
            "n_authors": len(p.get("authors", []) or []),
            "filename": p.get("filename"),
            "n_chunks": p.get("n_chunks", 0),
            "n_figures": p.get("n_figures", 0),
            "n_taxa": p.get("n_taxa", 0),
            "n_anatomy_terms": p.get("n_anatomy_terms", 0),
            "scan_file_type": p.get("scan_file_type"),
        })
    return rows if limit is None else rows[: int(limit)]


@mcp.tool()
def get_paper(paper_hash: str) -> Dict:
    """Full metadata for one paper: title, authors, year, abstract, DOI,
    plus top-10 taxa and top anatomy terms."""
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return {"error": f"no such paper_hash: {paper_hash}"}
    taxa = _load_json(Path(p["hash_dir"]) / "taxa.json", default={}) or {}
    anat = _load_json(Path(p["hash_dir"]) / "anatomy.json", default={}) or {}
    return {
        **{k: v for k, v in p.items() if k != "hash_dir"},
        "top_taxa": (taxa.get("taxa") or [])[:10],
        "top_anatomy": (anat.get("terms") or [])[:10],
    }


@mcp.tool()
def search_taxon(name: str) -> Dict:
    """Resolve a taxon name against the WoRMS snapshot.

    Accepts any form WoRMS knows: accepted names, historical synonyms,
    common misspellings present in the snapshot. Returns the accepted
    AphiaID, name, authority, rank, and the name_type of the match
    ('accepted' | 'unaccepted' | 'synonym'). Missing names return a
    structured ``not_found`` result rather than an error.
    """
    idx = _need_index()
    if idx.worms_db is None:
        return {"error": "no WoRMS snapshot configured"}
    hit = idx.worms_db.lookup(name)
    if not hit:
        return {"not_found": True, "queried": name}
    in_corpus = hit["accepted_aphia_id"] in idx.taxon_to_papers
    return {
        **hit,
        "queried": name,
        "in_corpus": in_corpus,
        "mentioning_paper_count": len(idx.taxon_to_papers.get(hit["accepted_aphia_id"], [])),
    }


@mcp.tool()
def get_papers_for_taxon(
    taxon_name: str,
    *,
    min_mentions: int = 1,
) -> List[Dict]:
    """List papers mentioning a taxon, resolved through synonymy.

    Papers are returned ordered by mention count (desc). Synonymy is
    followed to the accepted AphiaID — so a query for
    'Stephanomia amphitridis' returns papers citing its accepted synonym
    *Apolemia uvaria*, without the caller needing to know the synonymy.
    """
    idx = _need_index()
    if idx.worms_db is None:
        return [{"error": "no WoRMS snapshot configured"}]
    hit = idx.worms_db.lookup(taxon_name)
    if not hit:
        return []
    aid = hit["accepted_aphia_id"]
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
    if idx.worms_db is None:
        return [{"error": "no WoRMS snapshot configured"}]
    hit = idx.worms_db.lookup(taxon_name)
    if not hit:
        return []
    aid = hit["accepted_aphia_id"]
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
            if m.get("accepted_aphia_id") != aid:
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
def get_figures_for_taxon(
    taxon_name: str,
    paper_hash: Optional[str] = None,
    limit: int = 50,
) -> List[Dict]:
    """Figures from papers that mention the taxon, ranked by caption
    relevance.

    A figure whose caption names the taxon directly scores higher than
    a figure from a paper that merely mentions it elsewhere. Each
    result carries the on-disk figure image path so a downstream
    consumer (or the LLM client) can display it.
    """
    idx = _need_index()
    if idx.worms_db is None:
        return [{"error": "no WoRMS snapshot configured"}]
    hit = idx.worms_db.lookup(taxon_name)
    if not hit:
        return []
    aid = hit["accepted_aphia_id"]
    accepted_name_low = (hit["accepted_name"] or "").lower()
    matched_name_low = (hit["matched_name"] or "").lower()
    target_hashes = (
        [paper_hash] if paper_hash else idx.taxon_to_papers.get(aid, [])
    )

    rows: List[Dict] = []
    for h in target_hashes:
        p = idx.papers.get(h)
        if not p:
            continue
        figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
        for f in figs.get("figures", []) or []:
            caption = (f.get("caption_text") or f.get("caption") or "").lower()
            caption_hit = accepted_name_low in caption or (
                matched_name_low and matched_name_low in caption
            )
            rows.append({
                "paper_hash": h,
                "paper_title": p.get("title"),
                "figure_id": f.get("figure_id"),
                "page": f.get("page"),
                "caption_text": f.get("caption_text") or f.get("caption"),
                "figure_number": f.get("figure_number"),
                "image_path": str(Path(p["hash_dir"]) / "figures" / (f.get("filename") or "")),
                "caption_has_taxon": caption_hit,
                "score": (100 if caption_hit else 0) + idx.taxon_mention_counts.get(aid, {}).get(h, 0),
            })
    rows.sort(key=lambda r: -r["score"])
    return rows[: int(limit)] if limit else rows


@mcp.tool()
def get_figures_for_anatomy(
    anatomy_term: str,
    paper_hash: Optional[str] = None,
    limit: int = 50,
) -> List[Dict]:
    """Figures whose captions mention an anatomy term.

    Takes the canonical anatomy term or any of its configured synonyms
    / translations (see ``resources/anatomy_lexicon.yaml``). Matches the
    term case-insensitively in figure captions and returns the figure
    records with image paths and caption text, ranked by number of
    occurrences in the caption.
    """
    idx = _need_index()
    term_low = (anatomy_term or "").strip().lower()
    if not term_low:
        return []
    target_hashes = (
        [paper_hash] if paper_hash else list(idx.papers.keys())
    )
    rows: List[Dict] = []
    for h in target_hashes:
        p = idx.papers.get(h)
        if not p:
            continue
        figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
        for f in figs.get("figures", []) or []:
            caption = (f.get("caption_text") or f.get("caption") or "")
            occ = caption.lower().count(term_low)
            if occ == 0:
                continue
            rows.append({
                "paper_hash": h,
                "paper_title": p.get("title"),
                "figure_id": f.get("figure_id"),
                "page": f.get("page"),
                "figure_number": f.get("figure_number"),
                "caption_text": caption,
                "image_path": str(Path(p["hash_dir"]) / "figures" / (f.get("filename") or "")),
                "match_count": occ,
            })
    rows.sort(key=lambda r: -r["match_count"])
    return rows[: int(limit)] if limit else rows


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
def get_bibliography(paper_hash: str) -> List[Dict]:
    """Parsed references for one paper (from Grobid TEI)."""
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return [{"error": f"no such paper_hash: {paper_hash}"}]
    refs = _load_json(Path(p["hash_dir"]) / "references.json", default={}) or {}
    return refs.get("references", []) or []


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
    """All currently-valid species that descend from the given taxon in
    the WoRMS snapshot.

    Accepts any rank above species (genus, family, order, …). The result
    is a filtered view of the WoRMS snapshot; it does not consult the
    corpus — pair with :func:`get_papers_for_taxon` for per-species
    corpus coverage.
    """
    idx = _need_index()
    if idx.worms_db is None:
        return [{"error": "no WoRMS snapshot configured"}]
    hit = idx.worms_db.lookup(parent_taxon_name)
    if not hit:
        return []
    parent_aid = hit["accepted_aphia_id"]
    # Walk the parent_id tree in the snapshot. The tree may not be
    # strictly shallow so we BFS. Only include status='accepted'
    # Species (and Subspecies for completeness).
    conn = idx.worms_db.conn
    frontier = [parent_aid]
    descendants: List[int] = []
    seen = set()
    while frontier:
        parent = frontier.pop(0)
        if parent in seen:
            continue
        seen.add(parent)
        cur = conn.execute(
            "SELECT aphia_id, rank, status FROM taxa WHERE parent_id = ?", (parent,)
        )
        for row in cur:
            descendants.append(row[0])
            frontier.append(row[0])
    if not descendants:
        return []
    # Now fetch the species-level records
    placeholders = ",".join("?" * len(descendants))
    cur = conn.execute(
        f"""
        SELECT aphia_id, scientific_name, authority, rank
        FROM taxa
        WHERE aphia_id IN ({placeholders})
          AND status = 'accepted'
          AND rank IN ('Species', 'Subspecies')
        ORDER BY scientific_name
        """,
        descendants,
    )
    out: List[Dict] = []
    for row in cur:
        aid = row[0]
        name = row[1]
        out.append({
            "accepted_aphia_id": aid,
            "accepted_name": name,
            "authority": row[2],
            "rank": row[3],
            "mentioning_paper_count": len(idx.taxon_to_papers.get(aid, [])),
        })
    return out


@mcp.tool()
def get_figure(paper_hash: str, figure_id: str) -> Dict:
    """One figure's full record: caption, page, bbox, image path,
    cross-references."""
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return {"error": f"no such paper_hash: {paper_hash}"}
    figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
    for f in figs.get("figures", []) or []:
        if f.get("figure_id") == figure_id:
            return {
                **f,
                "paper_hash": paper_hash,
                "paper_title": p.get("title"),
                "image_path": str(
                    Path(p["hash_dir"]) / "figures" / (f.get("filename") or "")
                ),
            }
    return {"error": f"no such figure_id {figure_id!r} in paper {paper_hash}"}


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
# Entry point
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Root of the processed corpus (contains documents/)",
    )
    parser.add_argument(
        "--worms-sqlite",
        type=Path,
        default=None,
        help="Override path to the WoRMS SQLite (default: resources/worms_siphonophorae.sqlite under repo root)",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        stream=sys.stderr,  # stdout is the MCP transport
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    worms_path = args.worms_sqlite or (
        Path(__file__).parent / "resources" / "worms_siphonophorae.sqlite"
    )
    worms: Optional[WormsDB] = None
    if worms_path.exists():
        try:
            worms = WormsDB(worms_path)
            logger.info("WoRMS snapshot loaded from %s (%d names)",
                        worms_path, len(worms.name_set()))
        except Exception as e:
            logger.warning("Could not open WoRMS at %s: %s", worms_path, e)
    else:
        logger.warning(
            "WoRMS snapshot not found at %s — search_taxon and taxon tools "
            "will return an error. Build it: python scripts/ingest_worms.py",
            worms_path,
        )

    global _INDEX
    _INDEX = CorpusIndex(args.output_dir, worms_db=worms)
    n = _INDEX.load()
    logger.info("Serving corpus: %s (%d papers)", args.output_dir, n)

    mcp.run()  # default transport is stdio
    return 0


if __name__ == "__main__":
    sys.exit(main())
