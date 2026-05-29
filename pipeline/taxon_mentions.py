#!/usr/bin/env python3
"""Build the corpus-wide taxon mention database from per-paper taxa.json.

Reads existing per-paper artifacts only — no re-extraction needed. Walks
``<corpuscle>/documents/*/taxa.json`` and rolls up every mention into a
single SQLite database at ``<corpuscle>/taxon_mentions.sqlite``.

This is §12 Layer 2 from dev_docs/PLAN.md: a cross-paper mention table that enables
span-level taxon queries and (later) taxon × locality joins with Layer 3.

The per-paper taxa.json files remain the source of truth — this database
is a derived index that can be rebuilt at any time.

Idempotency (#30): each paper has a row in the ``papers_processed``
table; subsequent runs skip already-processed papers at the top of the
ingest loop. Re-running is a near-zero-work no-op when the corpus
hasn't changed; adding a new paper picks it up incrementally without
disturbing existing rows. ``--rebuild`` drops everything and
re-ingests from scratch.

Usage:
    python build_taxon_mentions.py /path/to/output
    python build_taxon_mentions.py /path/to/output --rebuild
    python build_taxon_mentions.py /path/to/output -o custom_path.sqlite
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import signal
import sqlite3
import sys
import time
from pathlib import Path
from typing import Optional

logger = logging.getLogger("build_taxon_mentions")

# Default is derived per-corpus from the output_dir positional arg in
# main(); see "corpuscle" layout in README.md.

# ── Schema ───────────────────────────────────────────────────────────

def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript("""
        -- Per-span mention records. ``taxon_id`` matches taxa.taxon_id in
        -- the Darwin Core taxonomy snapshot (TEXT — DwC taxonID is a
        -- string, e.g. "1371" for WoRMS or "urn:lsid:..." for GBIF).
        CREATE TABLE IF NOT EXISTS taxon_mentions (
            mention_id     INTEGER PRIMARY KEY AUTOINCREMENT,
            taxon_id       TEXT,
            matched_name   TEXT NOT NULL,
            accepted_name  TEXT,
            taxon_rank     TEXT,
            corpus_hash    TEXT NOT NULL,
            chunk_id       TEXT NOT NULL,
            chunk_index    INTEGER NOT NULL,
            char_start     INTEGER NOT NULL,
            char_end       INTEGER NOT NULL,
            mention_text   TEXT NOT NULL,
            name_type      TEXT,
            confidence     REAL DEFAULT 1.0,
            method         TEXT NOT NULL DEFAULT 'regex_taxonomy'
        );

        -- Per-paper processing status for resumability.
        -- ``source_mtime`` is the mtime of taxa.json at the time we
        -- ingested it; on subsequent passes we re-ingest if the file
        -- on disk is newer (auto-reconcile after per-paper re-runs).
        -- ``taxonomy_sha`` is the sha256 of the taxonomy.sqlite that
        -- the ingested taxa.json was resolved against, read from the
        -- paper's taxa-stage fingerprint (#95). mtime is unreliable
        -- across HPC array-task nodes with skewed clocks, so a change
        -- of recorded backbone forces a re-ingest even when mtime looks
        -- fresh.
        CREATE TABLE IF NOT EXISTS papers_processed (
            corpus_hash    TEXT PRIMARY KEY,
            n_mentions     INTEGER NOT NULL,
            n_unique_taxa  INTEGER NOT NULL,
            processed_at   REAL NOT NULL,
            source_mtime   REAL,
            taxonomy_sha   TEXT
        );

        CREATE TABLE IF NOT EXISTS build_meta (
            key   TEXT PRIMARY KEY,
            value TEXT
        );

        -- Indexes for the queries we care about
        CREATE INDEX IF NOT EXISTS idx_tm_taxon_id
            ON taxon_mentions(taxon_id);
        CREATE INDEX IF NOT EXISTS idx_tm_corpus_hash
            ON taxon_mentions(corpus_hash);
        CREATE INDEX IF NOT EXISTS idx_tm_chunk
            ON taxon_mentions(corpus_hash, chunk_id);
        CREATE INDEX IF NOT EXISTS idx_tm_accepted_name
            ON taxon_mentions(accepted_name);
    """)
    # Backward-compat ALTER for DBs created before the source_mtime
    # column existed. CREATE TABLE IF NOT EXISTS leaves an older schema
    # untouched, so without this migration the staleness check would
    # silently read NULL on every paper and re-ingest every time.
    have_cols = {row[1] for row in conn.execute("PRAGMA table_info(papers_processed)")}
    if "source_mtime" not in have_cols:
        conn.execute("ALTER TABLE papers_processed ADD COLUMN source_mtime REAL")
    # #95 — same back-compat ALTER for the taxonomy fingerprint column.
    # Older DBs read NULL taxonomy_sha, which the staleness check treats
    # as "ingested against an unknown backbone" so the first fingerprint-
    # aware pass re-ingests them once.
    if "taxonomy_sha" not in have_cols:
        conn.execute("ALTER TABLE papers_processed ADD COLUMN taxonomy_sha TEXT")
    conn.commit()


# ── Chunk ID parsing ────────────────────────────────────────────────

_CHUNK_INDEX_RE = re.compile(r"chunk_(\d+)")


def parse_chunk_index(chunk_id: str) -> int:
    """Extract integer index from chunk_id like 'chunk_5' → 5."""
    m = _CHUNK_INDEX_RE.search(chunk_id or "")
    return int(m.group(1)) if m else -1


def _paper_taxonomy_sha(hash_dir: Path) -> Optional[str]:
    """Return the sha256 of the taxonomy.sqlite that this paper's
    taxa.json was resolved against (#95), read from the taxa-stage
    ``input_fingerprint`` in ``pipeline_state.json``.

    Returns ``None`` when the file is missing/malformed, the taxa stage
    has no recorded fingerprint, or the corpus has no taxonomy — all of
    which mean "no recorded backbone to compare", and the staleness
    check falls back to mtime alone.
    """
    state_path = hash_dir / "pipeline_state.json"
    try:
        state = json.loads(state_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(state, dict):
        return None
    rec = (state.get("stages") or {}).get("taxa_and_lexicon_extraction") or {}
    taxonomy_fp = (rec.get("input_fingerprint") or {}).get("taxonomy") or {}
    sha = taxonomy_fp.get("sha256")
    return sha if isinstance(sha, str) else None


# ── Main build logic ────────────────────────────────────────────────

def ingest_paper(conn: sqlite3.Connection, corpus_hash: str,
                 taxa_data: dict) -> int:
    """Insert all mentions from one paper's taxa.json. Returns mention count."""
    mentions = taxa_data.get("mentions", [])
    if not mentions:
        return 0

    rows = []
    for m in mentions:
        chunk_id = m.get("chunk_id", "")
        span = m.get("text_span", [0, 0])
        rows.append((
            m.get("accepted_taxon_id"),
            m.get("matched_text", ""),
            m.get("accepted_name", ""),
            m.get("rank", ""),
            corpus_hash,
            chunk_id,
            parse_chunk_index(chunk_id),
            span[0] if len(span) > 0 else 0,
            span[1] if len(span) > 1 else 0,
            m.get("matched_text", ""),
            m.get("name_type", ""),
            1.0,
            "regex_taxonomy",
        ))

    conn.executemany(
        """INSERT INTO taxon_mentions
           (taxon_id, matched_name, accepted_name, taxon_rank, corpus_hash,
            chunk_id, chunk_index, char_start, char_end, mention_text,
            name_type, confidence, method)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        rows,
    )
    return len(rows)


def build(conn: sqlite3.Connection, output_dir: Path) -> dict:
    """Walk all papers and ingest their taxa.json mentions."""
    docs_dir = output_dir / "documents"
    if not docs_dir.is_dir():
        logger.error("Documents directory not found: %s", docs_dir)
        return {"papers": 0, "mentions": 0, "errors": 0}

    total_papers = 0
    total_mentions = 0
    skipped = 0
    refreshed = 0
    errors = 0
    lagging = 0
    batch = 0

    # #95 — hash the current taxonomy.sqlite once so we can flag papers
    # whose taxa.json was built against an older backbone (an operator
    # hint to re-run `corpus run`, which regenerates taxa.json before a
    # mention rebuild). None when the corpus has no taxonomy.
    from .stages import _file_sha256
    taxonomy_path = output_dir / "taxonomy.sqlite"
    current_taxonomy_sha: Optional[str] = None
    if taxonomy_path.exists():
        try:
            current_taxonomy_sha = _file_sha256(taxonomy_path)
        except OSError as e:
            logger.warning("Could not hash %s: %s", taxonomy_path, e)

    for hash_dir in sorted(docs_dir.iterdir()):
        if not hash_dir.is_dir():
            continue
        taxa_path = hash_dir / "taxa.json"
        if not taxa_path.exists():
            continue

        corpus_hash = hash_dir.name

        # Staleness check: keep the prior row when the source artifact
        # hasn't changed since we last ingested it; re-ingest otherwise.
        # Older DBs may have NULL ``source_mtime`` — treat that as
        # "ingested at unknown time" so the next pass refreshes them
        # rather than silently skipping forever.
        try:
            current_mtime = taxa_path.stat().st_mtime
        except OSError as e:
            logger.warning("Skipping %s: %s", taxa_path, e)
            errors += 1
            continue
        # #95 — the backbone a paper's taxa.json was resolved against,
        # recorded in its taxa-stage fingerprint. Drives re-ingest
        # independently of mtime.
        paper_taxo_sha = _paper_taxonomy_sha(hash_dir)
        if (
            current_taxonomy_sha is not None
            and paper_taxo_sha is not None
            and paper_taxo_sha != current_taxonomy_sha
        ):
            lagging += 1

        cur = conn.execute(
            "SELECT source_mtime, taxonomy_sha FROM papers_processed "
            "WHERE corpus_hash = ?",
            (corpus_hash,),
        )
        row = cur.fetchone()
        is_refresh = False
        if row is not None:
            stored_mtime, stored_taxo_sha = row
            mtime_fresh = stored_mtime is not None and stored_mtime >= current_mtime
            # Skip only when BOTH signals say "unchanged": taxa.json is
            # no newer than our last ingest AND the recorded backbone is
            # the same one we ingested against. A changed backbone forces
            # a re-ingest even when mtime looks fresh (#95 — the
            # cross-component case mtime alone misses on HPC nodes with
            # skewed clocks).
            taxo_unchanged = paper_taxo_sha == stored_taxo_sha
            if mtime_fresh and taxo_unchanged:
                skipped += 1
                continue
            is_refresh = True
            # Drop the prior rows for this corpus_hash before re-ingest
            # so we don't accumulate duplicate mentions when the source
            # file was regenerated.
            conn.execute(
                "DELETE FROM taxon_mentions WHERE corpus_hash = ?",
                (corpus_hash,),
            )

        try:
            taxa_data = json.loads(taxa_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as e:
            logger.warning("Skipping %s: %s", taxa_path, e)
            errors += 1
            continue

        n_mentions = ingest_paper(conn, corpus_hash, taxa_data)
        n_unique = taxa_data.get("unique_taxa", 0)

        conn.execute(
            """INSERT OR REPLACE INTO papers_processed
               (corpus_hash, n_mentions, n_unique_taxa, processed_at,
                source_mtime, taxonomy_sha)
               VALUES (?, ?, ?, ?, ?, ?)""",
            (corpus_hash, n_mentions, n_unique, time.time(), current_mtime,
             paper_taxo_sha),
        )

        total_papers += 1
        total_mentions += n_mentions
        if is_refresh:
            refreshed += 1
        batch += 1
        if batch >= 50:
            conn.commit()
            batch = 0
            logger.info(
                "Progress: %d papers, %d mentions (%d skipped, %d refreshed)",
                total_papers, total_mentions, skipped, refreshed,
            )

    conn.commit()
    logger.info(
        "Build complete: %d papers, %d mentions, %d skipped, %d refreshed, %d errors",
        total_papers, total_mentions, skipped, refreshed, errors,
    )
    if lagging:
        logger.warning(
            "%d paper(s) have taxa.json resolved against an OLDER taxonomy "
            "than the current %s. Their mentions reflect the stale backbone; "
            "re-run `corpus run` to regenerate taxa.json before rebuilding "
            "mentions.",
            lagging, taxonomy_path.name,
        )
    return {
        "papers": total_papers,
        "mentions": total_mentions,
        "skipped": skipped,
        "refreshed": refreshed,
        "errors": errors,
        "lagging": lagging,
    }


def write_stats(conn: sqlite3.Connection) -> None:
    """Write summary stats to build_meta."""
    for key, query in [
        ("total_mentions", "SELECT COUNT(*) FROM taxon_mentions"),
        ("total_papers", "SELECT COUNT(*) FROM papers_processed"),
        ("unique_taxon_ids", "SELECT COUNT(DISTINCT taxon_id) FROM taxon_mentions WHERE taxon_id IS NOT NULL"),
        ("unique_accepted_names", "SELECT COUNT(DISTINCT accepted_name) FROM taxon_mentions WHERE accepted_name IS NOT NULL AND accepted_name != ''"),
    ]:
        cur = conn.execute(query)
        val = cur.fetchone()[0]
        conn.execute(
            "INSERT OR REPLACE INTO build_meta (key, value) VALUES (?, ?)",
            (key, str(val)),
        )
        logger.info("  %s: %s", key, val)

    conn.execute(
        "INSERT OR REPLACE INTO build_meta (key, value) VALUES (?, ?)",
        ("build_complete_ts", str(time.time())),
    )
    conn.commit()


# ── CLI ──────────────────────────────────────────────────────────────

def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "output_dir", type=Path,
        help="Corpus output directory (contains documents/<hash>/ subdirs)",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=None,
        help="SQLite output path (default: <output_dir>/taxon_mentions.sqlite)",
    )
    parser.add_argument(
        "--rebuild", action="store_true",
        help="Drop and recreate all tables before building",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Report what would be processed without writing to the SQLite.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.output is None:
        args.output = args.output_dir / "taxon_mentions.sqlite"

    if args.dry_run:
        documents_dir = args.output_dir / "documents"
        if not documents_dir.is_dir():
            logger.error("Not a corpus output dir: %s (no documents/)", args.output_dir)
            return 1
        n_papers = sum(1 for d in documents_dir.iterdir() if d.is_dir())
        n_chunks_files = sum(
            1 for d in documents_dir.iterdir()
            if d.is_dir() and (d / "chunks.json").exists()
        )
        action = "rebuild from scratch" if args.rebuild else "incrementally update"
        logger.info(
            "Dry-run: would %s %s. Source: %d hash dirs, %d with chunks.json. "
            "No SQLite writes.",
            action, args.output, n_papers, n_chunks_files,
        )
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(args.output)

    # Graceful Ctrl-C
    def _sigint(signum, frame):  # noqa: ARG001
        logger.info("SIGINT received; committing and exiting")
        conn.commit()
        conn.close()
        sys.exit(130)
    signal.signal(signal.SIGINT, _sigint)

    try:
        if args.rebuild:
            logger.info("Rebuilding: dropping all tables")
            conn.executescript("""
                DROP TABLE IF EXISTS taxon_mentions;
                DROP TABLE IF EXISTS papers_processed;
                DROP TABLE IF EXISTS build_meta;
            """)

        create_schema(conn)

        conn.execute(
            "INSERT OR REPLACE INTO build_meta (key, value) VALUES (?, ?)",
            ("last_build_ts", str(time.time())),
        )
        conn.execute(
            "INSERT OR REPLACE INTO build_meta (key, value) VALUES (?, ?)",
            ("output_dir", str(args.output_dir)),
        )
        conn.commit()

        logger.info("═══ Building taxon mention index ═══")
        stats = build(conn, args.output_dir)

        logger.info("═══ Summary ═══")
        write_stats(conn)

        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
