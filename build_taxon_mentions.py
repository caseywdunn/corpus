#!/usr/bin/env python3
"""Build the corpus-wide taxon mention database from per-paper taxa.json.

Reads existing per-paper artifacts only — no re-extraction needed. Walks
``output/documents/*/taxa.json`` and rolls up every mention into a single
SQLite database at ``resources/taxon_mentions.sqlite``.

This is §12 Layer 2 from PLAN.md: a cross-paper mention table that enables
span-level taxon queries and (later) taxon × locality joins with Layer 3.

The per-paper taxa.json files remain the source of truth — this database
is a derived index that can be rebuilt at any time.

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

DEFAULT_OUTPUT = Path(__file__).resolve().parent / "resources" / "taxon_mentions.sqlite"

# ── Schema ───────────────────────────────────────────────────────────

def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript("""
        -- Per-span mention records
        CREATE TABLE IF NOT EXISTS taxon_mentions (
            mention_id     INTEGER PRIMARY KEY AUTOINCREMENT,
            aphia_id       INTEGER,
            matched_name   TEXT NOT NULL,
            accepted_name  TEXT,
            rank           TEXT,
            corpus_hash    TEXT NOT NULL,
            chunk_id       TEXT NOT NULL,
            chunk_index    INTEGER NOT NULL,
            char_start     INTEGER NOT NULL,
            char_end       INTEGER NOT NULL,
            mention_text   TEXT NOT NULL,
            name_type      TEXT,
            confidence     REAL DEFAULT 1.0,
            method         TEXT NOT NULL DEFAULT 'regex_worms'
        );

        -- Per-paper processing status for resumability
        CREATE TABLE IF NOT EXISTS papers_processed (
            corpus_hash    TEXT PRIMARY KEY,
            n_mentions     INTEGER NOT NULL,
            n_unique_taxa  INTEGER NOT NULL,
            processed_at   REAL NOT NULL
        );

        CREATE TABLE IF NOT EXISTS build_meta (
            key   TEXT PRIMARY KEY,
            value TEXT
        );

        -- Indexes for the queries we care about
        CREATE INDEX IF NOT EXISTS idx_tm_aphia
            ON taxon_mentions(aphia_id);
        CREATE INDEX IF NOT EXISTS idx_tm_corpus_hash
            ON taxon_mentions(corpus_hash);
        CREATE INDEX IF NOT EXISTS idx_tm_chunk
            ON taxon_mentions(corpus_hash, chunk_id);
        CREATE INDEX IF NOT EXISTS idx_tm_accepted_name
            ON taxon_mentions(accepted_name);
    """)
    conn.commit()


# ── Chunk ID parsing ────────────────────────────────────────────────

_CHUNK_INDEX_RE = re.compile(r"chunk_(\d+)")


def parse_chunk_index(chunk_id: str) -> int:
    """Extract integer index from chunk_id like 'chunk_5' → 5."""
    m = _CHUNK_INDEX_RE.search(chunk_id or "")
    return int(m.group(1)) if m else -1


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
            m.get("accepted_aphia_id"),
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
            "regex_worms",
        ))

    conn.executemany(
        """INSERT INTO taxon_mentions
           (aphia_id, matched_name, accepted_name, rank, corpus_hash,
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
    errors = 0
    batch = 0

    for hash_dir in sorted(docs_dir.iterdir()):
        if not hash_dir.is_dir():
            continue
        taxa_path = hash_dir / "taxa.json"
        if not taxa_path.exists():
            continue

        corpus_hash = hash_dir.name

        # Skip if already processed
        cur = conn.execute(
            "SELECT 1 FROM papers_processed WHERE corpus_hash = ?",
            (corpus_hash,),
        )
        if cur.fetchone():
            skipped += 1
            continue

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
               (corpus_hash, n_mentions, n_unique_taxa, processed_at)
               VALUES (?, ?, ?, ?)""",
            (corpus_hash, n_mentions, n_unique, time.time()),
        )

        total_papers += 1
        total_mentions += n_mentions
        batch += 1
        if batch >= 50:
            conn.commit()
            batch = 0
            logger.info(
                "Progress: %d papers, %d mentions (%d skipped)",
                total_papers, total_mentions, skipped,
            )

    conn.commit()
    logger.info(
        "Build complete: %d papers, %d mentions, %d skipped, %d errors",
        total_papers, total_mentions, skipped, errors,
    )
    return {
        "papers": total_papers,
        "mentions": total_mentions,
        "skipped": skipped,
        "errors": errors,
    }


def write_stats(conn: sqlite3.Connection) -> None:
    """Write summary stats to build_meta."""
    for key, query in [
        ("total_mentions", "SELECT COUNT(*) FROM taxon_mentions"),
        ("total_papers", "SELECT COUNT(*) FROM papers_processed"),
        ("unique_aphia_ids", "SELECT COUNT(DISTINCT aphia_id) FROM taxon_mentions WHERE aphia_id IS NOT NULL"),
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
        "-o", "--output", type=Path, default=DEFAULT_OUTPUT,
        help=f"SQLite output path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--rebuild", action="store_true",
        help="Drop and recreate all tables before building",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

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
