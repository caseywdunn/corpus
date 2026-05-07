"""SQLite + LanceDB wrappers and the corpus-wide in-memory index.

* :class:`CorpusIndex` — built once at startup by :func:`mcpsrv.main.main`.
  Holds per-paper headers + reverse indexes for taxon/lexicon/author
  lookup so tool calls don't re-scan the documents tree.
* :class:`TaxonMentionDB` — read-only wrapper over
  ``taxon_mentions.sqlite`` (built by ``build_taxon_mentions.py``).
* :class:`BiblioAuthority` — read-only wrapper over
  ``biblio_authority.sqlite`` (built by ``build_biblio_authority.py``).

The CorpusIndex defers loading the embedder + LanceDB table until the
first :func:`get_chunks_for_topic` call — both cost ~600 MB resident
and several seconds, wasteful for runs that only touch the structured
tools.
"""
from __future__ import annotations

import logging
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

from pipeline.taxa import TaxonomyDB
from pipeline.embeddings import EmbeddingBackend, EmbeddingError, get_embedder

from .app import _load_json

logger = logging.getLogger(__name__)


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
        # Bundle manifest (dev_docs/PLAN.md §10) if this is a served bundle.  None
        # for build outputs — bundle_info surfaces this distinction so
        # clients can tell a local dev run from a versioned deploy.
        self.bundle_manifest: Optional[Dict] = None
        # Reverse indexes — all point to a list of paper hashes.
        self.taxon_to_papers: Dict[str, List[str]] = defaultdict(list)
        self.taxon_mention_counts: Dict[str, Dict[str, int]] = defaultdict(dict)
        # Lexicon reverse indexes are nested by category so a future
        # `--lexicon` with anatomy + biogeography + methods produces three
        # independent term spaces. ``lexicon_to_papers["anatomy"][term]``
        # is the list of paper hashes whose anatomy.json mentions term.
        self.lexicon_to_papers: Dict[str, Dict[str, List[str]]] = defaultdict(
            lambda: defaultdict(list),
        )
        self.lexicon_mention_counts: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(
            lambda: defaultdict(dict),
        )
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

        # dev_docs/PLAN.md §10: served bundles ship a bundle_manifest.json at the
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
        # Files at the per-paper root that aren't lexicon outputs. Anything
        # else with a top-level ``category`` field is treated as one.
        non_lexicon_files = {
            "summary.json", "metadata.json", "references.json",
            "taxa.json", "figures.json", "chunks.json",
            "scan_detection.json", "docling_doc.json",
            "intext_citations.json", "pipeline_state.json",
        }
        for hash_dir in sorted(self.documents_dir.iterdir()):
            if not hash_dir.is_dir():
                continue
            summary = _load_json(hash_dir / "summary.json", default={}) or {}
            metadata = _load_json(hash_dir / "metadata.json", default={}) or {}
            taxa = _load_json(hash_dir / "taxa.json", default={"taxa": []}) or {}
            figures = _load_json(hash_dir / "figures.json", default={"figures": []}) or {}
            chunks = _load_json(hash_dir / "chunks.json", default={"chunks": []}) or {}

            # Discover lexicon artifacts dynamically. Any *.json in the
            # paper dir with a top-level ``category`` field is one (this is
            # exactly what taxa.extract_lexicon_mentions stamps).
            lexicons_for_paper: Dict[str, Dict] = {}
            for child in hash_dir.glob("*.json"):
                if child.name in non_lexicon_files:
                    continue
                payload = _load_json(child, default=None)
                if isinstance(payload, dict) and "category" in payload:
                    lexicons_for_paper[payload["category"]] = payload

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
                "n_lexicon_terms": {
                    cat: payload.get("unique_terms", 0)
                    for cat, payload in lexicons_for_paper.items()
                },
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

            for category, payload in lexicons_for_paper.items():
                for term in payload.get("terms", []) or []:
                    canonical = term.get("canonical")
                    if not canonical:
                        continue
                    self.lexicon_to_papers[category][canonical].append(paper_hash)
                    self.lexicon_mention_counts[category][canonical][paper_hash] = (
                        term.get("mention_count", 0)
                    )

            for author in metadata.get("authors", []) or []:
                surname = (author.get("surname") or "").strip().lower()
                if surname:
                    self.author_to_papers[surname].append(paper_hash)

            count += 1

        n_lexicon_terms = sum(len(by_term) for by_term in self.lexicon_to_papers.values())
        logger.info(
            "Index built: %d papers, %d taxa, %d lexicon term(s) across %d "
            "categor(ies), %d authors",
            count, len(self.taxon_to_papers),
            n_lexicon_terms, len(self.lexicon_to_papers),
            len(self.author_to_papers),
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
