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

from dotenv import load_dotenv
load_dotenv()  # Ensure ANTHROPIC_API_KEY etc. are visible to the stdio subprocess
                # Claude Code launches — no shell inheritance to rely on.

from mcp.server.fastmcp import FastMCP

from taxa import WormsDB
from embeddings import EmbeddingBackend, EmbeddingError, get_embedder

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

    def __init__(self, output_dir: Path, worms_db: Optional[WormsDB] = None,
                 embedding_backend: str = "local",
                 embedding_model: Optional[str] = None):
        self.output_dir = Path(output_dir).resolve()
        self.documents_dir = self.output_dir / "documents"
        self.vector_db_dir = self.output_dir / "vector_db" / "lancedb"
        self.worms_db = worms_db
        # Embedder + LanceDB table are loaded lazily on the first
        # get_chunks_for_topic call — they cost ~600 MB resident and
        # several seconds to load, both wasteful for users who only
        # touch the structured tools.
        self._embedding_backend_name = embedding_backend
        self._embedding_model_override = embedding_model
        self._embedder: Optional[EmbeddingBackend] = None
        self._lance_table = None

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

    def get_topic_searcher(self):
        """Return (embedder, table) for semantic chunk search, loading on
        first use. Returns ``(None, None)`` if there is no LanceDB table
        on disk — get_chunks_for_topic surfaces this as a clear error
        rather than crashing.
        """
        if self._lance_table is not None and self._embedder is not None:
            return self._embedder, self._lance_table
        try:
            import lancedb
        except ImportError:
            return None, None
        if not self.vector_db_dir.exists():
            return None, None
        try:
            db = lancedb.connect(str(self.vector_db_dir))
            if "document_chunks" not in db.table_names():
                return None, None
            table = db.open_table("document_chunks")
        except Exception as e:
            logger.warning("Could not open LanceDB at %s: %s", self.vector_db_dir, e)
            return None, None

        # Pin the model to whatever was used to build the table — schema
        # carries the dim and we choose a backend with matching dim.
        existing_dim = table.schema.field("vector").type.list_size
        try:
            embedder = get_embedder(
                self._embedding_backend_name,
                self._embedding_model_override,
            )
        except EmbeddingError as e:
            logger.warning("Could not load embedding backend: %s", e)
            return None, None
        if embedder.dim != existing_dim:
            logger.warning(
                "Embedding-model dim mismatch: backend=%s/%s emits %d, "
                "table expects %d. Topic search disabled until the model "
                "matches the index.",
                self._embedding_backend_name, embedder.model_name,
                embedder.dim, existing_dim,
            )
            return None, None
        self._embedder = embedder
        self._lance_table = table
        return embedder, table

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


_REAL_FIGURE_TYPES = {"figure", "plate", "subpanel"}


@mcp.tool()
def get_figures_for_taxon(
    taxon_name: str,
    paper_hash: Optional[str] = None,
    limit: int = 50,
    include_all: bool = False,
) -> List[Dict]:
    """Figures from papers that mention the taxon, ranked by caption
    relevance.

    A figure whose caption names the taxon directly scores higher than
    a figure from a paper that merely mentions it elsewhere. Each
    result carries the on-disk figure image path.

    By default only returns items classified as ``figure`` or ``plate``
    (skipping journal furniture, subpanels of already-returned figures,
    and unclassifiable graphical elements). Pass ``include_all=True`` to
    see every extracted item including the review bucket.
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
            ftype = f.get("figure_type")
            if not include_all and ftype not in _REAL_FIGURE_TYPES:
                continue
            caption = (f.get("caption_text") or f.get("caption") or "").lower()
            caption_hit = accepted_name_low in caption or (
                matched_name_low and matched_name_low in caption
            )
            rows.append({
                "paper_hash": h,
                "paper_title": p.get("title"),
                "figure_id": f.get("figure_id"),
                "figure_type": ftype,
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
    include_all: bool = False,
) -> List[Dict]:
    """Figures whose captions mention an anatomy term.

    Takes the canonical anatomy term or any of its configured synonyms
    / translations (see ``resources/anatomy_lexicon.yaml``). Matches the
    term case-insensitively in figure captions and returns the figure
    records with image paths and caption text, ranked by number of
    occurrences in the caption.

    By default only real figures and plates are returned; pass
    ``include_all=True`` to include the review bucket.
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
            ftype = f.get("figure_type")
            if not include_all and ftype not in _REAL_FIGURE_TYPES:
                continue
            caption = (f.get("caption_text") or f.get("caption") or "")
            occ = caption.lower().count(term_low)
            if occ == 0:
                continue
            rows.append({
                "paper_hash": h,
                "paper_title": p.get("title"),
                "figure_id": f.get("figure_id"),
                "figure_type": ftype,
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
def translate_chunk(
    paper_hash: str,
    chunk_id: str,
    target_language: str = "en",
    model: str = "claude-haiku-4-5-20251001",
) -> Dict:
    """Translate one chunk to the target language (default English), via
    the Anthropic Claude API.

    On-demand — nothing is translated at ingest time (PLAN.md §3 decision).
    Results are cached at ``<hash_dir>/translated_{target_language}.json``
    keyed by ``chunk_id`` so repeated calls are free.

    Returns the original text, the translation, the source language
    reported by ``scan_detection.json``, and whether the result was
    cached or freshly generated.

    Returns ``{error: …}`` if ``ANTHROPIC_API_KEY`` isn't set or the
    Anthropic SDK isn't installed. If the chunk is already in the target
    language (matching ``scan_detection.detected_language``) the original
    text is returned with ``translation_needed: False``.
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
    questions (PLAN.md §8 Q6, Q8) where the match criterion is
    semantic, not a literal taxon or anatomy term — for those use
    ``get_chunks_for_taxon`` / ``get_figures_for_anatomy`` instead,
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


@mcp.tool()
def list_figure_rois(paper_hash: str, figure_id: str) -> List[Dict]:
    """Return the per-panel / per-subfigure ROIs annotated on a figure.

    ROIs are populated by Pass 2.5 (caption-derived ``panels_from_caption``)
    and Pass 3a (OCR-derived ``rois`` with pixel bboxes). The caption-
    derived list gives labels + descriptions even when Pass 3a wasn't
    run or OCR didn't find the panel; ``rois`` gives pixel coordinates
    when available for :func:`get_figure_roi_image`.
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return [{"error": f"no such paper_hash: {paper_hash}"}]
    figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
    for f in figs.get("figures", []) or []:
        if f.get("figure_id") == figure_id:
            return [{
                "paper_hash": paper_hash,
                "figure_id": figure_id,
                "figure_number": f.get("figure_number"),
                "filename": f.get("filename"),
                "panels_from_caption": f.get("panels_from_caption") or [],
                "panel_count_from_caption": f.get("panel_count_from_caption", 0),
                "rois": f.get("rois") or [],
                "pass3_status": f.get("pass3_status"),
                "image_size_px": f.get("image_size_px"),
            }]
    return [{"error": f"no such figure_id {figure_id!r} in paper {paper_hash}"}]


@mcp.tool()
def get_figure_roi_image(
    paper_hash: str,
    figure_id: str,
    label: str,
) -> Dict:
    """Crop a panel ROI out of a figure image and return the crop's path.

    ``label`` is the panel letter (e.g. ``"A"``, ``"B"``) as stored in
    ``figures.json`` ``rois[*].label``. If Pass 3a found a pixel ROI
    for that label, we crop the figure PNG to that region and cache the
    result under ``<hash_dir>/figures/crops/{figure_id}__{label}.png``
    so repeated asks are free.

    Returns the crop path + caption description. If the label exists in
    ``panels_from_caption`` but Pass 3a couldn't find a pixel ROI for it
    (OCR missed the label), returns the whole figure's image path with
    ``crop: false`` — the LLM can still display the whole figure and
    reason about which region is which panel using the caption.
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return {"error": f"no such paper_hash: {paper_hash}"}
    hash_dir = Path(p["hash_dir"])
    figs = _load_json(hash_dir / "figures.json", default={}) or {}
    fig = next(
        (f for f in figs.get("figures", []) or [] if f.get("figure_id") == figure_id),
        None,
    )
    if fig is None:
        return {"error": f"no such figure_id {figure_id!r} in paper {paper_hash}"}

    whole_image = hash_dir / "figures" / (fig.get("filename") or "")
    caption_text = fig.get("caption_text") or fig.get("caption") or ""
    description_from_caption = next(
        (p["description"] for p in fig.get("panels_from_caption") or []
         if p.get("label") == label),
        None,
    )
    roi_entry = next(
        (r for r in fig.get("rois") or []
         if r.get("label") == label and r.get("roi_px")),
        None,
    )

    if roi_entry is None:
        # No pixel ROI — return whole figure with the caption info so the
        # LLM can display + caption-filter mentally.
        return {
            "paper_hash": paper_hash,
            "figure_id": figure_id,
            "label": label,
            "crop": False,
            "image_path": str(whole_image),
            "caption_text": caption_text,
            "description_from_caption": description_from_caption,
            "reason": "no_pixel_roi — Pass 3a didn't locate this label in the image",
        }

    # Cache crops so a second retrieval is free.
    crops_dir = hash_dir / "figures" / "crops"
    crops_dir.mkdir(parents=True, exist_ok=True)
    crop_path = crops_dir / f"{figure_id}__{label}.png"
    if not crop_path.exists():
        try:
            from PIL import Image
            with Image.open(whole_image) as im:
                x0, y0, x1, y1 = [int(v) for v in roi_entry["roi_px"]]
                # Clamp the ROI to the image bounds defensively — ROI
                # computation can exceed image dims at the edges.
                x0 = max(0, min(x0, im.width - 1))
                y0 = max(0, min(y0, im.height - 1))
                x1 = max(x0 + 1, min(x1, im.width))
                y1 = max(y0 + 1, min(y1, im.height))
                im.crop((x0, y0, x1, y1)).save(str(crop_path))
        except Exception as e:
            return {"error": f"could not crop figure: {e}"}

    return {
        "paper_hash": paper_hash,
        "figure_id": figure_id,
        "label": label,
        "crop": True,
        "image_path": str(crop_path),
        "roi_px": roi_entry.get("roi_px"),
        "ocr_confidence": roi_entry.get("ocr_confidence"),
        "caption_text": caption_text,
        "description_from_caption": description_from_caption,
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
    parser.add_argument(
        "--embedding-backend",
        choices=["local", "openai"],
        default="local",
        help="Embedding backend used by get_chunks_for_topic (default: local). "
             "Must produce vectors of the same dim as the LanceDB index.",
    )
    parser.add_argument(
        "--embedding-model",
        default=None,
        help="Override the per-backend default embedding model "
             "(local default BAAI/bge-m3, openai default text-embedding-3-small).",
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
    _INDEX = CorpusIndex(
        args.output_dir,
        worms_db=worms,
        embedding_backend=args.embedding_backend,
        embedding_model=args.embedding_model,
    )
    n = _INDEX.load()
    logger.info("Serving corpus: %s (%d papers)", args.output_dir, n)

    mcp.run()  # default transport is stdio
    return 0


if __name__ == "__main__":
    sys.exit(main())
