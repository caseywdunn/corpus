#!/usr/bin/env python3
"""Build a SQLite snapshot of the WoRMS taxonomy under order Siphonophorae.

Starts at AphiaID 1371 (Siphonophorae) and recursively fetches every
descendant and every synonym, storing the result as a two-table SQLite
that the corpus pipeline consumes at run time.

Why SQLite (not JSON, not the live API):
  * The live API has rate limits and would couple pipeline runs to
    network availability. A snapshot is reproducible, version-pinnable,
    and lets every per-PDF extraction step be deterministic.
  * SQLite supports a case-insensitive name lookup in O(log n) via an
    index on a lowercased column; a JSON blob would force a linear scan.
  * It is trivially shippable to Bouchet with the rest of the repo
    artifacts and loadable read-only in parallel worker processes.

Schema:
  taxa     one row per AphiaID (both accepted and non-accepted)
  names    flat (name, aphia_id, name_type) — many-to-one to taxa

The script is idempotent and resumable: re-running picks up where it
left off by skipping AphiaIDs already in the ``taxa`` table. Interrupt
with SIGINT; the SQLite is committed in small batches.

Usage:
    python ingest_worms.py                       # default output
    python ingest_worms.py -o /path/to/out.sqlite
    python ingest_worms.py --root-aphia 1371     # default is Siphonophorae
"""

from __future__ import annotations

import argparse
import logging
import signal
import sqlite3
import sys
import time
from pathlib import Path
from typing import Iterable, List, Optional

import requests


WORMS_API = "https://www.marinespecies.org/rest"
DEFAULT_ROOT = 1371  # Siphonophorae
DEFAULT_OUTPUT = Path(__file__).resolve().parent.parent / "resources" / "worms_siphonophorae.sqlite"


logger = logging.getLogger("ingest_worms")


class RateLimitedSession:
    """Wrap requests so we never hammer the WoRMS API faster than
    ``min_interval`` seconds per call. The API docs ask for courtesy; a
    half-second pacing keeps us well-behaved without crawling."""

    def __init__(self, min_interval: float = 0.3):
        self.min_interval = min_interval
        self._last_call = 0.0
        self.session = requests.Session()

    def get(self, url: str, **kw) -> requests.Response:
        delay = self.min_interval - (time.monotonic() - self._last_call)
        if delay > 0:
            time.sleep(delay)
        self._last_call = time.monotonic()
        return self.session.get(url, **kw)


def create_schema(conn: sqlite3.Connection) -> None:
    """Create the two-table schema if it doesn't exist."""
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS taxa (
            aphia_id INTEGER PRIMARY KEY,
            scientific_name TEXT NOT NULL,
            authority TEXT,
            rank TEXT,
            status TEXT,
            parent_id INTEGER,
            valid_aphia_id INTEGER,   -- for synonyms, the accepted name's AphiaID
            valid_name TEXT,
            kingdom TEXT, phylum TEXT, class TEXT, "order" TEXT,
            family TEXT, genus TEXT,
            citation TEXT,
            lsid TEXT,
            url TEXT,
            fetched_at REAL NOT NULL
        );

        CREATE TABLE IF NOT EXISTS names (
            name TEXT NOT NULL,
            name_lowercase TEXT NOT NULL,
            aphia_id INTEGER NOT NULL,
            name_type TEXT NOT NULL,  -- 'accepted' | 'unaccepted' | 'synonym'
            FOREIGN KEY (aphia_id) REFERENCES taxa(aphia_id)
        );
        CREATE INDEX IF NOT EXISTS idx_names_lowercase ON names(name_lowercase);
        CREATE INDEX IF NOT EXISTS idx_names_aphia_id ON names(aphia_id);
        CREATE INDEX IF NOT EXISTS idx_taxa_valid ON taxa(valid_aphia_id);
        CREATE INDEX IF NOT EXISTS idx_taxa_parent ON taxa(parent_id);

        CREATE TABLE IF NOT EXISTS meta (
            key TEXT PRIMARY KEY,
            value TEXT
        );
        """
    )
    conn.commit()


def record_exists(conn: sqlite3.Connection, aphia_id: int) -> bool:
    cur = conn.execute("SELECT 1 FROM taxa WHERE aphia_id = ?", (aphia_id,))
    return cur.fetchone() is not None


def insert_taxon(conn: sqlite3.Connection, rec: dict, name_type: str) -> None:
    """Insert a taxon record + its primary name. Synonyms are added later."""
    conn.execute(
        """
        INSERT OR REPLACE INTO taxa (
            aphia_id, scientific_name, authority, rank, status,
            parent_id, valid_aphia_id, valid_name,
            kingdom, phylum, class, "order", family, genus,
            citation, lsid, url, fetched_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            rec.get("AphiaID"),
            rec.get("scientificname") or "",
            rec.get("authority"),
            rec.get("rank"),
            rec.get("status"),
            rec.get("parentNameUsageID"),
            rec.get("valid_AphiaID"),
            rec.get("valid_name"),
            rec.get("kingdom"), rec.get("phylum"), rec.get("class"), rec.get("order"),
            rec.get("family"), rec.get("genus"),
            rec.get("citation"), rec.get("lsid"), rec.get("url"),
            time.time(),
        ),
    )
    # Insert the scientific name into the flat names index. A single row
    # per name/AphiaID/type combination is fine — dedup on insert below.
    name = rec.get("scientificname") or ""
    if name:
        conn.execute(
            "INSERT INTO names (name, name_lowercase, aphia_id, name_type) VALUES (?, ?, ?, ?)",
            (name, name.lower(), rec.get("AphiaID"), name_type),
        )


def insert_synonym_name(conn: sqlite3.Connection, accepted_aphia_id: int, syn_name: str) -> None:
    """Index a synonym name pointing at the accepted taxon."""
    conn.execute(
        "INSERT INTO names (name, name_lowercase, aphia_id, name_type) VALUES (?, ?, ?, ?)",
        (syn_name, syn_name.lower(), accepted_aphia_id, "synonym"),
    )


def fetch_record(sess: RateLimitedSession, aphia_id: int) -> Optional[dict]:
    r = sess.get(f"{WORMS_API}/AphiaRecordByAphiaID/{aphia_id}")
    if r.status_code == 204:
        return None
    r.raise_for_status()
    return r.json()


def fetch_children(sess: RateLimitedSession, aphia_id: int) -> List[dict]:
    """Page through /AphiaChildrenByAphiaID. The endpoint returns up to 50
    records per call and asks for an ``offset`` 1-indexed; we stop when a
    page comes back short."""
    out: List[dict] = []
    offset = 1
    while True:
        r = sess.get(
            f"{WORMS_API}/AphiaChildrenByAphiaID/{aphia_id}",
            params={"offset": offset, "marine_only": "false"},
        )
        if r.status_code == 204:
            break
        r.raise_for_status()
        page = r.json() or []
        if not page:
            break
        out.extend(page)
        if len(page) < 50:
            break
        offset += 50
    return out


def fetch_synonyms(sess: RateLimitedSession, aphia_id: int) -> List[dict]:
    """Page through /AphiaSynonymsByAphiaID. Same paging contract as
    children. Returns the TaxonRecord for each synonym (so we also learn
    their AphiaIDs and authorities, not just names)."""
    out: List[dict] = []
    offset = 1
    while True:
        r = sess.get(
            f"{WORMS_API}/AphiaSynonymsByAphiaID/{aphia_id}",
            params={"offset": offset},
        )
        if r.status_code == 204:
            break
        r.raise_for_status()
        page = r.json() or []
        if not page:
            break
        out.extend(page)
        if len(page) < 50:
            break
        offset += 50
    return out


def walk(sess: RateLimitedSession, conn: sqlite3.Connection, root_id: int) -> None:
    """Breadth-first walk from ``root_id`` storing every descendant.

    Skips AphiaIDs already in the local DB so a killed run can resume.
    Commits every 50 taxa so a crash leaves progress intact.
    """
    to_visit: List[int] = [root_id]
    batch_count = 0
    visited_total = 0

    while to_visit:
        aphia_id = to_visit.pop(0)
        if record_exists(conn, aphia_id):
            continue
        try:
            rec = fetch_record(sess, aphia_id)
        except Exception as e:
            logger.warning("record fetch failed for %s: %s", aphia_id, e)
            continue
        if not rec:
            continue

        # Insert taxon itself (accepted or unaccepted based on its status).
        name_type = "accepted" if rec.get("status") == "accepted" else "unaccepted"
        insert_taxon(conn, rec, name_type=name_type)

        # For accepted taxa, also fetch + store synonyms (name records
        # only; the synonyms themselves get their own taxa row if/when we
        # encounter them as siblings).
        if rec.get("status") == "accepted":
            try:
                syns = fetch_synonyms(sess, aphia_id)
            except Exception as e:
                logger.warning("synonyms fetch failed for %s: %s", aphia_id, e)
                syns = []
            for syn in syns:
                syn_name = syn.get("scientificname") or ""
                syn_id = syn.get("AphiaID")
                if not syn_id:
                    continue
                # Record the synonym taxon if we don't have it; it's common
                # for the same unaccepted name to be a synonym of multiple
                # accepted names, so INSERT OR REPLACE is safe.
                if not record_exists(conn, syn_id):
                    insert_taxon(conn, syn, name_type="synonym")
                # Also register the synonym name against the *accepted*
                # taxon so name-based lookup can resolve directly.
                insert_synonym_name(conn, aphia_id, syn_name)

        # Enqueue children for recursion.
        try:
            children = fetch_children(sess, aphia_id)
        except Exception as e:
            logger.warning("children fetch failed for %s: %s", aphia_id, e)
            children = []
        for child in children:
            cid = child.get("AphiaID")
            if cid and not record_exists(conn, cid):
                to_visit.append(cid)

        visited_total += 1
        batch_count += 1
        if batch_count >= 50:
            conn.commit()
            batch_count = 0
            logger.info(
                "Checkpoint: %d taxa, queue=%d (last=%s %s)",
                visited_total, len(to_visit),
                rec.get("rank"), rec.get("scientificname"),
            )

    conn.commit()
    logger.info("Walk complete: %d taxa fetched from root %d", visited_total, root_id)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o", "--output", type=Path, default=DEFAULT_OUTPUT,
        help=f"SQLite output path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--root-aphia", type=int, default=DEFAULT_ROOT,
        help=f"WoRMS AphiaID to start walking from (default: {DEFAULT_ROOT} = Siphonophorae)",
    )
    parser.add_argument(
        "--min-interval", type=float, default=0.3,
        help="Minimum seconds between API calls (default: 0.3; be courteous)",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(args.output)
    try:
        create_schema(conn)
        conn.execute(
            "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
            ("root_aphia_id", str(args.root_aphia)),
        )
        conn.execute(
            "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
            ("last_ingest_ts", str(time.time())),
        )
        conn.commit()

        # Graceful Ctrl-C: commit and exit.
        def _sigint(signum, frame):  # noqa: ARG001
            logger.info("SIGINT received; committing and exiting")
            conn.commit()
            conn.close()
            sys.exit(130)

        signal.signal(signal.SIGINT, _sigint)

        sess = RateLimitedSession(min_interval=args.min_interval)
        walk(sess, conn, args.root_aphia)

        # Summary stats
        cur = conn.execute("SELECT COUNT(*) FROM taxa")
        n_taxa = cur.fetchone()[0]
        cur = conn.execute("SELECT COUNT(*) FROM names")
        n_names = cur.fetchone()[0]
        cur = conn.execute(
            "SELECT rank, COUNT(*) FROM taxa GROUP BY rank ORDER BY COUNT(*) DESC LIMIT 10"
        )
        ranks = list(cur)
        logger.info("Done. %d taxa, %d names in %s", n_taxa, n_names, args.output)
        logger.info("Rank distribution (top 10): %s", ranks)
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
