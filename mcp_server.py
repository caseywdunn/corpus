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
        --taxonomy-db /path/to/taxonomy.sqlite

MCP clients typically discover the server through a config file that
tells them to launch this script. See README.md for the Claude Desktop
/ Claude Code snippets.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional

from dotenv import load_dotenv
load_dotenv()  # Ensure ANTHROPIC_API_KEY etc. are visible to the stdio subprocess
                # Claude Code launches — no shell inheritance to rely on.

from mcp.server.fastmcp import FastMCP

import hmac
import os
import sqlite3

from taxa import TaxonomyDB
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

    def __init__(self, output_dir: Path, taxonomy_db: Optional[TaxonomyDB] = None,
                 biblio_db: Optional[BiblioAuthority] = None,
                 taxon_mention_db: Optional[TaxonMentionDB] = None,
                 embedding_model: Optional[str] = None):
        self.output_dir = Path(output_dir).resolve()
        self.documents_dir = self.output_dir / "documents"
        self.vector_db_dir = self.output_dir / "vector_db" / "lancedb"
        self.taxonomy_db = taxonomy_db
        self.biblio_db = biblio_db
        self.taxon_mention_db = taxon_mention_db
        # Embedder + LanceDB table are loaded lazily on the first
        # get_chunks_for_topic call — they cost ~600 MB resident and
        # several seconds to load, both wasteful for users who only
        # touch the structured tools.
        self._embedding_model_override = embedding_model
        self._embedder: Optional[EmbeddingBackend] = None
        self._lance_table = None

        # Per-paper header, keyed on the 12-char short hash (directory name).
        self.papers: Dict[str, Dict] = {}
        # Bundle manifest (PLAN.md §10) if this is a served bundle.  None
        # for build outputs — bundle_info surfaces this distinction so
        # clients can tell a local dev run from a versioned deploy.
        self.bundle_manifest: Optional[Dict] = None
        # Reverse indexes — all point to a list of paper hashes.
        self.taxon_to_papers: Dict[str, List[str]] = defaultdict(list)
        self.taxon_mention_counts: Dict[str, Dict[str, int]] = defaultdict(dict)
        self.anatomy_to_papers: Dict[str, List[str]] = defaultdict(list)
        self.anatomy_mention_counts: Dict[str, Dict[str, int]] = defaultdict(dict)
        self.author_to_papers: Dict[str, List[str]] = defaultdict(list)
        # accepted_taxon_id → accepted name (for display when we only
        # have IDs from the reverse index).
        self.taxon_display: Dict[str, Dict] = {}

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
        # carries the dim and we choose a model with matching dim.
        existing_dim = table.schema.field("vector").type.list_size
        try:
            embedder = get_embedder(self._embedding_model_override)
        except EmbeddingError as e:
            logger.warning("Could not load embedding backend: %s", e)
            return None, None
        if embedder.dim != existing_dim:
            logger.warning(
                "Embedding-model dim mismatch: model=%s emits %d, "
                "table expects %d. Topic search disabled until the model "
                "matches the index.",
                embedder.model_name, embedder.dim, existing_dim,
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

        # PLAN.md §10: served bundles ship a bundle_manifest.json at the
        # root of the output dir.  Build outputs don't.  Read on startup
        # so bundle_info can report it without re-opening files.
        self.bundle_manifest = _load_json(
            self.output_dir / "bundle_manifest.json", default=None
        )
        if self.bundle_manifest:
            logger.info(
                "Bundle manifest loaded: %s (created %s)",
                self.bundle_manifest.get("bundle_version", "?"),
                self.bundle_manifest.get("created_at", "?"),
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
                aid = t.get("accepted_taxon_id")
                if not aid:
                    continue
                self.taxon_to_papers[aid].append(paper_hash)
                self.taxon_mention_counts[aid][paper_hash] = t.get("mention_count", 0)
                self.taxon_display.setdefault(aid, {
                    "accepted_taxon_id": aid,
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
# Bibliographic authority (read-only wrapper)
# ---------------------------------------------------------------------------


class TaxonMentionDB:
    """Read-only wrapper around ``taxon_mentions.sqlite``.

    Provides corpus-wide span-level taxon queries — §12 Layer 2.
    """

    def __init__(self, db_path: Path):
        self.conn = sqlite3.connect(
            f"file:{db_path}?mode=ro", uri=True, check_same_thread=False,
        )
        self.conn.row_factory = sqlite3.Row

    def mentions_for_taxon(self, taxon_id: int, *,
                           corpus_hash: Optional[str] = None,
                           limit: int = 500,
                           offset: int = 0) -> List[Dict]:
        query = """SELECT * FROM taxon_mentions WHERE taxon_id = ?"""
        params: list = [taxon_id]
        if corpus_hash:
            query += " AND corpus_hash = ?"
            params.append(corpus_hash)
        query += " ORDER BY corpus_hash, chunk_index, char_start"
        query += " LIMIT ? OFFSET ?"
        params.extend([limit, offset])
        return [dict(r) for r in self.conn.execute(query, params)]

    def mention_count_for_taxon(self, taxon_id: int) -> int:
        cur = self.conn.execute(
            "SELECT COUNT(*) FROM taxon_mentions WHERE taxon_id = ?",
            (taxon_id,),
        )
        return cur.fetchone()[0]

    def papers_for_taxon_id(self, taxon_id: int) -> List[Dict]:
        """Distinct papers mentioning a taxon, with mention count."""
        cur = self.conn.execute(
            """SELECT corpus_hash, COUNT(*) as mention_count
               FROM taxon_mentions WHERE taxon_id = ?
               GROUP BY corpus_hash
               ORDER BY mention_count DESC""",
            (taxon_id,),
        )
        return [dict(r) for r in cur]

    def stats(self) -> Dict:
        out = {}
        for key in ("total_mentions", "total_papers", "unique_taxon_ids",
                     "unique_accepted_names"):
            cur = self.conn.execute(
                "SELECT value FROM build_meta WHERE key = ?", (key,),
            )
            row = cur.fetchone()
            out[key] = row[0] if row else None
        return out


class BiblioAuthority:
    """Read-only wrapper around ``biblio_authority.sqlite``."""

    def __init__(self, db_path: Path):
        self.conn = sqlite3.connect(
            f"file:{db_path}?mode=ro", uri=True, check_same_thread=False,
        )
        self.conn.row_factory = sqlite3.Row

    def get_work(self, work_id: str) -> Optional[Dict]:
        cur = self.conn.execute("SELECT * FROM works WHERE work_id = ?", (work_id,))
        row = cur.fetchone()
        return dict(row) if row else None

    def get_work_by_corpus_hash(self, corpus_hash: str) -> Optional[Dict]:
        cur = self.conn.execute(
            "SELECT * FROM works WHERE corpus_hash = ?", (corpus_hash,),
        )
        row = cur.fetchone()
        return dict(row) if row else None

    def get_authors(self, work_id: str) -> List[Dict]:
        cur = self.conn.execute(
            "SELECT surname, forename, position FROM work_authors "
            "WHERE work_id = ? ORDER BY position",
            (work_id,),
        )
        return [dict(r) for r in cur]

    def citing(self, work_id: str) -> List[Dict]:
        """Works that cite the given work."""
        cur = self.conn.execute(
            """SELECT w.work_id, w.title, w.year, w.in_corpus, w.corpus_hash,
                      c.match_method, c.match_score
               FROM citations c JOIN works w ON c.citing_work_id = w.work_id
               WHERE c.cited_work_id = ?
               ORDER BY w.year""",
            (work_id,),
        )
        return [dict(r) for r in cur]

    def cited_by(self, work_id: str) -> List[Dict]:
        """Works cited by the given work."""
        cur = self.conn.execute(
            """SELECT w.work_id, w.title, w.year, w.in_corpus, w.corpus_hash,
                      c.match_method, c.match_score
               FROM citations c JOIN works w ON c.cited_work_id = w.work_id
               WHERE c.citing_work_id = ?
               ORDER BY w.year""",
            (work_id,),
        )
        return [dict(r) for r in cur]

    def search_works(self, surname: str, year: Optional[int] = None,
                     title_fragment: Optional[str] = None) -> List[Dict]:
        """Search works by author surname, optional year, optional title."""
        query = """
            SELECT DISTINCT w.work_id, w.title, w.year, w.journal,
                   w.in_corpus, w.corpus_hash, w.guid_type, w.doi
            FROM works w JOIN work_authors wa ON w.work_id = wa.work_id
            WHERE wa.surname_normalized = ?
        """
        params: list = [surname.strip().lower()]
        if year:
            query += " AND w.year = ?"
            params.append(year)
        cur = self.conn.execute(query + " ORDER BY w.year", params)
        results = [dict(r) for r in cur]
        if title_fragment and results:
            frag = title_fragment.lower()
            results = [r for r in results if r.get("title") and frag in r["title"].lower()]
        return results

    def citation_count(self, work_id: str) -> int:
        cur = self.conn.execute(
            "SELECT COUNT(*) FROM citations WHERE cited_work_id = ?", (work_id,),
        )
        return cur.fetchone()[0]

    def taxon_links_for_work(self, work_id: str) -> List[Dict]:
        cur = self.conn.execute(
            "SELECT taxon_id, link_type, confidence FROM taxon_work_links WHERE work_id = ?",
            (work_id,),
        )
        return [dict(r) for r in cur]

    def work_for_taxon(self, taxon_id: int) -> List[Dict]:
        cur = self.conn.execute(
            """SELECT w.work_id, w.title, w.year, w.in_corpus, w.corpus_hash,
                      wl.link_type, wl.confidence
               FROM taxon_work_links wl JOIN works w ON wl.work_id = w.work_id
               WHERE wl.taxon_id = ?""",
            (taxon_id,),
        )
        return [dict(r) for r in cur]


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
def bundle_info() -> Dict:
    """Return metadata about the served bundle the server is backed by.

    For bundles produced by ``package_for_serve.py`` this reports the
    version string, creation timestamp, pipeline git SHA, embedding
    model + dim, paper / chunk / figure counts, and whether PDFs were
    included.  Clients can call this on startup to detect stale
    endpoints ("I queried against v1.0.0 but the server now reports
    v1.1.0") and to cite a corpus version in downstream work.

    For local build outputs (no ``bundle_manifest.json`` at the root),
    returns ``{"bundle_version": null}`` — the server still works but
    there is no stable identifier for what it's serving.
    """
    idx = _need_index()
    if idx.bundle_manifest is None:
        return {
            "bundle_version": None,
            "note": "no bundle_manifest.json — this server is backed "
                    "by a local build output, not a versioned served "
                    "bundle.  run package_for_serve.py to produce one.",
        }
    return dict(idx.bundle_manifest)


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
                "rank": r["rank"],
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


# Need the chunk-index parser for the fallback path
_CHUNK_INDEX_RE_MCP = re.compile(r"chunk_(\d+)")
def parse_chunk_index(chunk_id: str) -> int:
    m = _CHUNK_INDEX_RE_MCP.search(chunk_id or "")
    return int(m.group(1)) if m else -1


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
    if idx.taxonomy_db is None:
        return [{"error": "no taxonomy snapshot configured"}]
    hit = idx.taxonomy_db.lookup(taxon_name)
    if not hit:
        return []
    aid = hit["accepted_taxon_id"]
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
    / translations (defined by the user-supplied ``--anatomy-lexicon``
    YAML at process time; see ``demo/anatomy_lexicon.yaml`` for an
    example). Matches the
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
def get_bibliography(
    paper_hash: str,
    resolved: bool = False,
) -> List[Dict]:
    """Parsed references for one paper (from Grobid TEI).

    With ``resolved=True``, each reference is enriched with its
    bibliographic authority record: ``work_id``, ``in_corpus``,
    ``corpus_hash`` (if in corpus), and ``cited_by_count`` (how many
    corpus papers cite this work). This turns a flat list into a
    navigable bibliography — click through to a cited work that's also
    in the corpus via ``get_paper(corpus_hash)``.

    Requires the bibliographic authority database; falls back to
    unresolved output if unavailable.
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return [{"error": f"no such paper_hash: {paper_hash}"}]
    refs = _load_json(Path(p["hash_dir"]) / "references.json", default={}) or {}
    ref_list = refs.get("references", []) or []
    if not resolved or idx.biblio_db is None:
        return ref_list

    # Enrich each reference with authority DB info
    for ref in ref_list:
        cur = idx.biblio_db.conn.execute(
            """SELECT c.cited_work_id, w.in_corpus, w.corpus_hash,
                      w.title AS resolved_title, w.doi AS resolved_doi,
                      c.match_method, c.match_score
               FROM citations c JOIN works w ON c.cited_work_id = w.work_id
               WHERE c.citing_corpus_hash = ? AND c.grobid_xml_id = ?""",
            (paper_hash, ref.get("xml_id", "")),
        )
        row = cur.fetchone()
        if row:
            ref["work_id"] = row["cited_work_id"]
            ref["in_corpus"] = bool(row["in_corpus"])
            ref["corpus_hash"] = row["corpus_hash"]
            ref["match_method"] = row["match_method"]
            ref["cited_by_count"] = idx.biblio_db.citation_count(row["cited_work_id"])
        else:
            ref["work_id"] = None
            ref["in_corpus"] = False
            ref["corpus_hash"] = None
    return ref_list


@mcp.tool()
def get_citation_graph(
    work_id: Optional[str] = None,
    paper_hash: Optional[str] = None,
    direction: str = "both",
    depth: int = 1,
) -> Dict:
    """Citation graph around a work.

    Accepts either a ``work_id`` (DOI, ``bhl:…``, or ``corpus:…``) or
    a ``paper_hash`` (for corpus papers). ``direction`` controls which
    edges to follow: ``"citing"`` (papers that cite this work),
    ``"cited_by"`` (papers cited by this work), or ``"both"``.
    ``depth > 1`` follows transitive citations (max 3).

    Returns the root work plus the citation edges.
    """
    idx = _need_index()
    if idx.biblio_db is None:
        return {"error": "bibliographic authority database not configured"}
    # Resolve paper_hash to work_id if needed
    if not work_id and paper_hash:
        w = idx.biblio_db.get_work_by_corpus_hash(paper_hash)
        if not w:
            return {"error": f"no work found for paper_hash: {paper_hash}"}
        work_id = w["work_id"]
    if not work_id:
        return {"error": "provide either work_id or paper_hash"}

    root = idx.biblio_db.get_work(work_id)
    if not root:
        return {"error": f"no such work_id: {work_id}"}

    depth = min(int(depth), 3)
    result: Dict[str, Any] = {
        "root": {
            **root,
            "authors": idx.biblio_db.get_authors(work_id),
            "cited_by_count": idx.biblio_db.citation_count(work_id),
        },
    }

    if direction in ("citing", "both"):
        citing = _walk_citations(idx.biblio_db, work_id, "citing", depth)
        result["citing"] = citing
    if direction in ("cited_by", "both"):
        cited_by = _walk_citations(idx.biblio_db, work_id, "cited_by", depth)
        result["cited_by"] = cited_by

    return result


def _walk_citations(
    biblio: BiblioAuthority, work_id: str, direction: str, depth: int,
) -> List[Dict]:
    """BFS citation walk."""
    visited: set = set()
    frontier = [work_id]
    results: List[Dict] = []
    for d in range(depth):
        next_frontier: List[str] = []
        for wid in frontier:
            if wid in visited:
                continue
            visited.add(wid)
            rows = biblio.citing(wid) if direction == "citing" else biblio.cited_by(wid)
            for r in rows:
                r["depth"] = d + 1
                results.append(r)
                if r["work_id"] not in visited:
                    next_frontier.append(r["work_id"])
        frontier = next_frontier
    return results


@mcp.tool()
def resolve_reference(
    query: str,
    author: Optional[str] = None,
    year: Optional[int] = None,
) -> Dict:
    """Resolve a free-text bibliographic reference to a work in the
    authority database.

    Accepts forms like ``"Haeckel 1888"``, ``"Totton 1965 A Synopsis"``,
    or a raw citation string. Optionally pass ``author`` and ``year``
    separately for more precise matching. Returns the matched work with
    all known identifiers, in-corpus status, and citation counts.
    """
    idx = _need_index()
    if idx.biblio_db is None:
        return {"error": "bibliographic authority database not configured"}

    # Parse author and year from query if not provided separately
    if not author:
        import re
        # Try to extract "Author Year" or "Author, Year"
        m = re.match(r"^([A-Za-zÀ-ÿ\-\s]+?)[\s,]+(\d{4})\b", query.strip())
        if m:
            author = m.group(1).strip()
            if not year:
                year = int(m.group(2))

    if not author:
        return {"error": "could not parse author from query; pass author= explicitly"}

    # Title fragment is whatever remains after author+year
    title_frag = query
    if author:
        title_frag = title_frag.replace(author, "", 1).strip()
    if year:
        title_frag = title_frag.replace(str(year), "", 1).strip()
    title_frag = title_frag.strip(" ,;:")

    results = idx.biblio_db.search_works(
        author, year, title_frag if title_frag else None,
    )
    if not results:
        # Broaden: try without title
        results = idx.biblio_db.search_works(author, year)
    if not results:
        return {"not_found": True, "queried": query, "parsed_author": author, "parsed_year": year}
    if len(results) == 1:
        w = results[0]
        w["authors"] = idx.biblio_db.get_authors(w["work_id"])
        w["cited_by_count"] = idx.biblio_db.citation_count(w["work_id"])
        w["taxon_links"] = idx.biblio_db.taxon_links_for_work(w["work_id"])
        return w
    # Multiple matches — return all with citation counts
    for w in results:
        w["cited_by_count"] = idx.biblio_db.citation_count(w["work_id"])
    return {"matches": results, "count": len(results), "queried": query}


@mcp.tool()
def get_missing_references(
    min_citations: int = 2,
    year_from: Optional[int] = None,
    year_to: Optional[int] = None,
    limit: int = 50,
) -> List[Dict]:
    """Works cited by corpus papers that are NOT in the corpus.

    Sorted by citation count (most-cited missing works first). Useful
    for identifying high-impact papers to add to the corpus. Filter by
    year range to focus on a particular era.
    """
    idx = _need_index()
    if idx.biblio_db is None:
        return [{"error": "bibliographic authority database not configured"}]

    query = """
        SELECT w.work_id, w.title, w.year, w.journal, w.doi,
               w.guid_type, COUNT(*) AS cited_by_count
        FROM citations c JOIN works w ON c.cited_work_id = w.work_id
        WHERE w.in_corpus = 0
    """
    params: list = []
    if year_from:
        query += " AND w.year >= ?"
        params.append(year_from)
    if year_to:
        query += " AND w.year <= ?"
        params.append(year_to)
    query += " GROUP BY c.cited_work_id HAVING COUNT(*) >= ?"
    params.append(int(min_citations))
    query += " ORDER BY cited_by_count DESC LIMIT ?"
    params.append(int(limit))

    cur = idx.biblio_db.conn.execute(query, params)
    results = []
    for row in cur:
        r = dict(row)
        r["authors"] = idx.biblio_db.get_authors(r["work_id"])
        results.append(r)
    return results


@mcp.tool()
def get_original_description(taxon_name: str) -> Dict:
    """Find the original description paper for a taxon.

    Resolves the taxon through DwC synonymy (``acceptedNameUsageID``),
    parses ``scientificNameAuthorship``, and looks up the corresponding
    work in the bibliographic authority database. Returns the work
    record, whether it's in the corpus, and if so the ``corpus_hash``
    for direct access via ``get_paper()``.
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return {"error": "no taxonomy snapshot configured"}
    if idx.biblio_db is None:
        return {"error": "bibliographic authority database not configured"}

    hit = idx.taxonomy_db.lookup(taxon_name)
    if not hit:
        return {"not_found": True, "queried": taxon_name}

    aid = hit["accepted_taxon_id"]
    works = idx.biblio_db.work_for_taxon(aid)
    if not works:
        return {
            "taxon": hit,
            "original_description": None,
            "note": "no matching work found in the authority database",
        }

    # Enrich with authors and citation count
    for w in works:
        w["authors"] = idx.biblio_db.get_authors(w["work_id"])
        w["cited_by_count"] = idx.biblio_db.citation_count(w["work_id"])

    return {
        "taxon": hit,
        "original_description": works[0] if len(works) == 1 else None,
        "candidate_works": works if len(works) > 1 else None,
        "work": works[0] if len(works) == 1 else works,
    }


@mcp.tool()
def get_works_by_author(
    surname: str,
    in_corpus_only: bool = False,
    limit: int = 100,
) -> List[Dict]:
    """All works by an author across the full bibliographic authority
    database — not just corpus papers, but also cited references and
    stub works synthesized from taxonomic-authority strings.

    Each result includes ``cited_by_count`` (how many corpus papers
    cite it) and ``in_corpus`` flag. Use ``in_corpus_only=True`` to
    restrict to papers physically in the corpus.
    """
    idx = _need_index()
    if idx.biblio_db is None:
        return [{"error": "bibliographic authority database not configured"}]

    norm = surname.strip().lower()
    query = """
        SELECT DISTINCT w.work_id, w.title, w.year, w.journal, w.doi,
               w.in_corpus, w.corpus_hash, w.guid_type
        FROM works w JOIN work_authors wa ON w.work_id = wa.work_id
        WHERE wa.surname_normalized = ?
    """
    params: list = [norm]
    if in_corpus_only:
        query += " AND w.in_corpus = 1"
    query += " ORDER BY w.year LIMIT ?"
    params.append(int(limit))

    cur = idx.biblio_db.conn.execute(query, params)
    results = []
    for row in cur:
        r = dict(row)
        r["authors"] = idx.biblio_db.get_authors(r["work_id"])
        r["cited_by_count"] = idx.biblio_db.citation_count(r["work_id"])
        results.append(r)
    return results


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
# HTTP transport + bearer-token auth (PLAN.md §10)
# ---------------------------------------------------------------------------


class _BearerAuthASGI:
    """ASGI middleware: 401 unless ``Authorization: Bearer <token>`` matches.

    Wraps any ASGI app (FastMCP's ``sse_app()``) at construction time —
    avoids the post-init middleware caveats on finalized Starlette apps.
    Lifespan and websocket scopes pass through unauthenticated; only
    ``http`` requests are gated.

    ``hmac.compare_digest`` makes the comparison constant-time.  For a
    shared-secret setup with ~20 trusted collaborators this is cheap
    insurance against timing attacks.
    """

    def __init__(self, app, token: str):
        self.app = app
        self._token_bytes = token.encode("utf-8")

    async def __call__(self, scope, receive, send):
        if scope.get("type") != "http":
            await self.app(scope, receive, send)
            return
        raw = dict(scope.get("headers", [])).get(b"authorization", b"")
        auth = raw.decode("latin-1", errors="replace")
        if not auth.startswith("Bearer "):
            await _send_401(send, "missing bearer token")
            return
        offered = auth[7:].encode("utf-8")
        if not hmac.compare_digest(offered, self._token_bytes):
            await _send_401(send, "invalid bearer token")
            return
        await self.app(scope, receive, send)


async def _send_401(send, reason: str) -> None:
    body = f"Unauthorized: {reason}\n".encode("utf-8")
    await send({
        "type": "http.response.start",
        "status": 401,
        "headers": [
            (b"content-type", b"text/plain; charset=utf-8"),
            (b"www-authenticate", b'Bearer realm="corpus-mcp"'),
        ],
    })
    await send({"type": "http.response.body", "body": body})


def _load_auth_token(cli_token_file: Optional[Path]) -> Optional[str]:
    """Resolve the bearer-auth token.

    Search order:
      1. ``--auth-token-file <path>`` — read, strip.
      2. ``CORPUS_MCP_TOKEN`` env var.
      3. None (server runs open; main() logs a big warning).

    ``--auth-token`` is deliberately *not* a CLI flag — secrets on the
    command line leak via ``ps``/proc/<pid>/cmdline.
    """
    if cli_token_file is not None:
        return cli_token_file.read_text().strip()
    env = os.environ.get("CORPUS_MCP_TOKEN", "").strip()
    return env or None


def _run_sse(host: str, port: int, token: Optional[str]) -> None:
    """Serve the FastMCP app over SSE on ``host:port``, optionally
    with bearer-token auth.  Requires uvicorn (pulled in transitively
    by the ``mcp`` package for its HTTP transports)."""
    try:
        import uvicorn
    except ImportError as e:
        raise RuntimeError(
            "uvicorn is required for --transport sse. "
            "It should come in via `pip install mcp`; if you're seeing "
            "this, your mcp install may be incomplete."
        ) from e

    app = mcp.sse_app()
    if token:
        app = _BearerAuthASGI(app, token)
        logger.info("Bearer-token auth enabled")
    else:
        logger.warning(
            "*** Running WITHOUT auth: anyone who can reach %s:%d can call "
            "every MCP tool. Set CORPUS_MCP_TOKEN or --auth-token-file "
            "before exposing this beyond localhost. ***",
            host, port,
        )
    logger.info("Serving SSE on http://%s:%d", host, port)
    uvicorn.run(app, host=host, port=port, log_level="info")


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
        "--taxonomy-db",
        type=Path,
        default=None,
        help="Override path to the Darwin Core taxonomy SQLite "
             "(default: <output_dir>/taxonomy.sqlite). "
             "Build with: python ingest_taxonomy.py --source <dwc|dwca|worms> ...",
    )
    parser.add_argument(
        "--biblio-sqlite",
        type=Path,
        default=None,
        help="Override path to the bibliographic authority SQLite "
             "(default: <output_dir>/biblio_authority.sqlite)",
    )
    parser.add_argument(
        "--taxon-mention-sqlite",
        type=Path,
        default=None,
        help="Override path to the taxon mention SQLite "
             "(default: <output_dir>/taxon_mentions.sqlite)",
    )
    parser.add_argument(
        "--instructions",
        type=Path,
        default=None,
        help="Markdown file whose contents are returned to MCP clients in "
             "the InitializeResult.instructions field, so well-behaved "
             "clients (Claude Desktop, Claude Code) inject it into the LLM's "
             "context at session start.  Use it for per-corpus nudges — e.g. "
             "'Velella is not a siphonophore'.  "
             "Default: <output_dir>/instructions.md if it exists.",
    )
    parser.add_argument(
        "--embedding-model",
        default=None,
        help="Override the default HuggingFace embedding model id "
             "(default: BAAI/bge-m3). Must emit vectors of the same "
             "dim as the LanceDB index.",
    )
    parser.add_argument(
        "--transport",
        choices=["stdio", "sse"],
        default="stdio",
        help="MCP transport. stdio (default) is for local MCP clients "
             "that launch this process themselves (Claude Desktop, Claude "
             "Code, Cursor). sse serves over HTTP/Server-Sent-Events for "
             "remote deployments (PLAN.md §10).",
    )
    parser.add_argument(
        "--host", default="127.0.0.1",
        help="Interface to bind when --transport sse (default: 127.0.0.1)",
    )
    parser.add_argument(
        "--port", type=int, default=8080,
        help="Port to bind when --transport sse (default: 8080)",
    )
    parser.add_argument(
        "--auth-token-file", type=Path, default=None,
        help="Path to a file whose contents are the bearer-auth token. "
             "Alternative: set CORPUS_MCP_TOKEN. CLI-literal tokens are "
             "deliberately not supported — they leak via `ps`.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        stream=sys.stderr,  # stdout is the MCP transport
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    # Optional per-corpus instructions surfaced to clients in
    # InitializeResult.instructions.  FastMCP exposes ``instructions`` as
    # a read-only property backed by ``_mcp_server.instructions``; that
    # backing attribute is plain mutable state on the lower-level Server,
    # so setting it here (after the module-level FastMCP is constructed)
    # takes effect on the next ``initialize`` call.
    instructions_path = args.instructions or (args.output_dir / "instructions.md")
    if instructions_path.exists():
        try:
            text = instructions_path.read_text(encoding="utf-8").strip()
            mcp._mcp_server.instructions = text
            logger.info(
                "Loaded MCP instructions from %s (%d chars)",
                instructions_path, len(text),
            )
        except Exception as e:
            logger.warning(
                "Could not read instructions %s: %s", instructions_path, e,
            )

    taxonomy_path = args.taxonomy_db or (args.output_dir / "taxonomy.sqlite")
    taxonomy: Optional[TaxonomyDB] = None
    if taxonomy_path.exists():
        try:
            taxonomy = TaxonomyDB(taxonomy_path)
            logger.info("Taxonomy snapshot loaded from %s (%d names)",
                        taxonomy_path, len(taxonomy.name_set()))
        except Exception as e:
            logger.warning("Could not open taxonomy DB at %s: %s", taxonomy_path, e)
    else:
        logger.warning(
            "Taxonomy snapshot not found at %s — search_taxon and taxon tools "
            "will return an error. Build it: "
            "python ingest_taxonomy.py --source <dwc|dwca|worms> ...",
            taxonomy_path,
        )

    biblio_path = args.biblio_sqlite or (args.output_dir / "biblio_authority.sqlite")
    biblio: Optional[BiblioAuthority] = None
    if biblio_path.exists():
        try:
            biblio = BiblioAuthority(biblio_path)
            n_works = biblio.conn.execute("SELECT COUNT(*) FROM works").fetchone()[0]
            logger.info("Bibliographic authority loaded from %s (%d works)",
                        biblio_path, n_works)
        except Exception as e:
            logger.warning("Could not open biblio authority at %s: %s", biblio_path, e)
    else:
        logger.warning(
            "Bibliographic authority not found at %s — citation graph tools "
            "will return an error. Build it: python build_biblio_authority.py",
            biblio_path,
        )

    taxon_mention_path = args.taxon_mention_sqlite or (args.output_dir / "taxon_mentions.sqlite")
    taxon_mention_db: Optional[TaxonMentionDB] = None
    if taxon_mention_path.exists():
        try:
            taxon_mention_db = TaxonMentionDB(taxon_mention_path)
            stats = taxon_mention_db.stats()
            logger.info(
                "Taxon mention DB loaded from %s (%s mentions, %s taxa)",
                taxon_mention_path,
                stats.get("total_mentions", "?"),
                stats.get("unique_taxon_ids", "?"),
            )
        except Exception as e:
            logger.warning("Could not open taxon mention DB at %s: %s",
                           taxon_mention_path, e)
    else:
        logger.info(
            "Taxon mention DB not found at %s — get_taxon_mentions will "
            "fall back to per-paper taxa.json scanning. Build it: "
            "python build_taxon_mentions.py",
            taxon_mention_path,
        )

    global _INDEX
    _INDEX = CorpusIndex(
        args.output_dir,
        taxonomy_db=taxonomy,
        biblio_db=biblio,
        taxon_mention_db=taxon_mention_db,
        embedding_model=args.embedding_model,
    )
    n = _INDEX.load()
    logger.info("Serving corpus: %s (%d papers)", args.output_dir, n)

    if args.transport == "stdio":
        mcp.run()
    elif args.transport == "sse":
        token = _load_auth_token(args.auth_token_file)
        _run_sse(args.host, args.port, token)
    else:  # argparse's choices= should prevent this
        raise ValueError(f"unknown transport: {args.transport!r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
