#!/usr/bin/env python3
"""Merge DOI-keyed ghost rows with their corpus_key-keyed siblings.

A corpus paper with a DOI often ends up split across two ``works``
rows in the biblio authority DB:

  * a ``doi:<doi>`` row that accumulates citations from references
    that carried the DOI through Grobid, and
  * a ``corpus:<surname>|<year>|<title>`` row — either the in-corpus
    paper itself (Phase 1) or a ghost built by Phase 2 from
    references that didn't carry a DOI.

Both rows describe the same work, so ``get_missing_references`` and
the citation graph look split.  Example from the live DB:

    corpus:dunn|2005|molecular phylogenetics ...    in_corpus=1   10 cites
    10.1080/10635150500354837                       in_corpus=0   59 cites

After this pass the corpus row absorbs the DOI row's 59 cites and
carries the DOI in its ``doi`` column — one canonical row per work.

Algorithm (runs per DOI ghost):

  1. Candidates: all ``corpus_key`` rows with the same first-author
     surname and year as the DOI row.
  2. Score each candidate's title against the DOI row's title via
     the dual-metric gate from #2's cascade: ``token_set_ratio ≥ 85
     AND ratio ≥ 60``.  If exactly one candidate clears the bar,
     it's the unambiguous match.  If zero match, leave the pair
     alone.  If ≥ 2 match, skip (ambiguous — needs a human call).
  3. Pick the survivor:
     - If exactly one side is ``in_corpus=1``, that row survives and
       keeps its work_id.  The DOI moves into its ``doi`` field.
     - If both are ghosts, the DOI row survives (DOI is the more
       canonical identifier; its work_id is globally unique).
     - If both are in_corpus=1 — rare and surprising — skip.
  4. Redirect incoming + outgoing citations onto the survivor, copy
     authors and aliases the survivor doesn't already carry, drop
     the non-survivor row.

Idempotent: once a work lives under a single row with a populated
``doi`` field, no DOI ghost remains to pair with.

Usage:
    python unify_doi_corpus_key.py
    python unify_doi_corpus_key.py --dry-run
    python unify_doi_corpus_key.py --min-set-score 90
"""

from __future__ import annotations

import argparse
import logging
import re
import sqlite3
import sys
import time
import unicodedata
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from rapidfuzz import fuzz
    _HAS_RAPIDFUZZ = True
except ImportError:
    fuzz = None  # type: ignore
    _HAS_RAPIDFUZZ = False

logger = logging.getLogger("unify_doi")

# Default is derived per-corpus from the output_dir positional arg in
# main(); see "corpuscle" layout in README.md.


# ── Normalization (kept in sync with build_biblio_authority.py) ─────

def normalize_for_key(s: str) -> str:
    s = unicodedata.normalize("NFD", s)
    s = "".join(c for c in s if unicodedata.category(c) != "Mn")
    s = s.lower().strip()
    s = re.sub(r"[^\w\s]", "", s)
    s = re.sub(r"\s+", " ", s)
    return s


# ── Candidate lookup ────────────────────────────────────────────────

def find_corpus_key_siblings(conn: sqlite3.Connection, doi_row: sqlite3.Row
                             ) -> List[sqlite3.Row]:
    """Return corpus_key rows in the same (first-author surname, year)
    family as ``doi_row``."""
    surname = doi_row["surname"]
    year = doi_row["year"]
    if not surname or year is None:
        return []
    cur = conn.execute(
        """SELECT w.work_id, w.title, w.in_corpus, w.corpus_hash, w.doi,
                  wa.surname_normalized AS surname, w.year,
                  (SELECT COUNT(*) FROM citations c WHERE c.cited_work_id = w.work_id) AS cites
           FROM works w
           JOIN work_authors wa ON w.work_id = wa.work_id
           WHERE wa.position = 0
             AND wa.surname_normalized = ?
             AND w.year = ?
             AND w.guid_type = 'corpus_key'""",
        (surname, year),
    )
    return list(cur)


def score_title(a: str, b: str) -> Tuple[int, int]:
    """Return (token_set_ratio, partial_ratio) for two titles.

    ``partial_ratio`` (instead of plain ``ratio``) handles the common
    case where the DOI-keyed row carries the full published title
    with subtitle and the corpus_key row has the shorter Grobid-
    extracted variant.  The canonical "Molecular Phylogenetics of
    the Siphonophora" scores partial_ratio=100 against the longer
    "... (Cnidaria), with implications for the evolution of
    functional specialization" — a plain ratio gives ~53 because it
    penalizes length mismatch, letting obvious duplicates slip past.
    """
    if not a or not b:
        return 0, 0
    na, nb = normalize_for_key(a), normalize_for_key(b)
    return int(fuzz.token_set_ratio(na, nb)), int(fuzz.partial_ratio(na, nb))


def pick_match(doi_row: sqlite3.Row,
               candidates: List[sqlite3.Row],
               min_set_score: int,
               min_partial_ratio: int) -> Tuple[Optional[sqlite3.Row], str, int]:
    """Return (match, status, n_qualifying) where status is
    'matched' | 'no_match' | 'ambiguous'.

    Tie-break: when multiple candidates clear the gate, accept the
    pick if exactly one of them is ``in_corpus=1`` (the real paper
    beats its own OCR ghosts that dedup left behind below the 95-set
    threshold).  Otherwise — multiple in-corpus candidates, or zero
    in-corpus among qualifiers — ambiguous and skip.
    """
    qualifying = []
    for c in candidates:
        set_s, partial_s = score_title(doi_row["title"], c["title"])
        if set_s >= min_set_score and partial_s >= min_partial_ratio:
            qualifying.append((c, set_s, partial_s))
    if not qualifying:
        return None, "no_match", 0
    if len(qualifying) == 1:
        return qualifying[0][0], "matched", 1
    in_corpus_hits = [q for q in qualifying if q[0]["in_corpus"] == 1]
    if len(in_corpus_hits) == 1:
        return in_corpus_hits[0][0], "matched", len(qualifying)
    return qualifying[0][0], "ambiguous", len(qualifying)


# ── Merge target selection ──────────────────────────────────────────

def pick_survivor(doi_row: sqlite3.Row, corpus_row: sqlite3.Row) -> Tuple[str, str, str]:
    """Return (survivor_id, dup_id, reason).  Raises if both are
    in_corpus=1 (rare; caller skips)."""
    if doi_row["in_corpus"] == 1 and corpus_row["in_corpus"] == 1:
        raise ValueError("both in_corpus — not auto-mergeable")
    if corpus_row["in_corpus"] == 1:
        return corpus_row["work_id"], doi_row["work_id"], "corpus_paper_wins"
    if doi_row["in_corpus"] == 1:
        # Unusual: DOI work_id on an in-corpus row.  Keep it.
        return doi_row["work_id"], corpus_row["work_id"], "doi_in_corpus_wins"
    # Both ghosts — DOI is the more canonical identifier.
    return doi_row["work_id"], corpus_row["work_id"], "doi_wins_ghost_ghost"


# ── Merge ───────────────────────────────────────────────────────────

def merge_duplicate(conn: sqlite3.Connection,
                    survivor_id: str,
                    dup_id: str,
                    doi_value: Optional[str]) -> None:
    """Redirect edges, copy metadata, set DOI on survivor, delete dup.

    Near-identical to ``dedup_ghost_works.merge_ghost_into_canonical``
    but additionally writes the DOI onto the survivor's ``doi``
    column so downstream ``lookup_by_doi`` still resolves."""
    if survivor_id == dup_id:
        return

    # 1. Redirect incoming citations (dup was cited -> survivor is cited).
    conn.execute(
        "UPDATE OR IGNORE citations SET cited_work_id = ? WHERE cited_work_id = ?",
        (survivor_id, dup_id),
    )
    conn.execute("DELETE FROM citations WHERE cited_work_id = ?", (dup_id,))

    # 2. Redirect outgoing citations (survivor now does the citing).
    conn.execute(
        "UPDATE OR IGNORE citations SET citing_work_id = ? WHERE citing_work_id = ?",
        (survivor_id, dup_id),
    )
    conn.execute("DELETE FROM citations WHERE citing_work_id = ?", (dup_id,))

    # 3. Copy authors the survivor doesn't already have.
    cur = conn.execute(
        "SELECT position, surname, surname_normalized, forename "
        "FROM work_authors WHERE work_id = ? ORDER BY position",
        (dup_id,),
    )
    dup_authors = list(cur)
    cur = conn.execute(
        "SELECT MAX(position) FROM work_authors WHERE work_id = ?", (survivor_id,),
    )
    row = cur.fetchone()
    next_pos = (row[0] + 1) if row and row[0] is not None else 0
    for _, surname, surname_norm, forename in dup_authors:
        dup_check = conn.execute(
            "SELECT 1 FROM work_authors WHERE work_id = ? AND surname_normalized = ?",
            (survivor_id, surname_norm),
        ).fetchone()
        if dup_check:
            continue
        conn.execute(
            """INSERT OR IGNORE INTO work_authors
               (work_id, position, surname, surname_normalized, forename)
               VALUES (?, ?, ?, ?, ?)""",
            (survivor_id, next_pos, surname, surname_norm, forename),
        )
        next_pos += 1

    # 4. Copy aliases the survivor doesn't already have.
    conn.execute(
        """INSERT OR IGNORE INTO work_aliases (alias_key, work_id)
           SELECT alias_key, ? FROM work_aliases WHERE work_id = ?""",
        (survivor_id, dup_id),
    )

    # 5. Propagate the DOI onto the survivor if it doesn't carry one.
    if doi_value:
        cur = conn.execute(
            "SELECT doi FROM works WHERE work_id = ?", (survivor_id,)
        )
        row = cur.fetchone()
        current = row[0] if row else None
        if not current:
            conn.execute(
                "UPDATE works SET doi = ?, updated_at = ? WHERE work_id = ?",
                (doi_value, time.time(), survivor_id),
            )

    # 6. Drop dependents and the dup row itself.
    conn.execute("DELETE FROM work_aliases WHERE work_id = ?", (dup_id,))
    conn.execute("DELETE FROM work_authors WHERE work_id = ?", (dup_id,))
    conn.execute("DELETE FROM works WHERE work_id = ?", (dup_id,))


# ── Main loop ───────────────────────────────────────────────────────

def load_doi_ghosts(conn: sqlite3.Connection) -> List[sqlite3.Row]:
    """Every DOI-keyed ghost that has a first-author + year we can
    use to look up corpus_key siblings."""
    cur = conn.execute(
        """SELECT w.work_id, w.title, w.year, w.in_corpus, w.doi,
                  wa.surname_normalized AS surname,
                  (SELECT COUNT(*) FROM citations c WHERE c.cited_work_id = w.work_id) AS cites
           FROM works w
           JOIN work_authors wa ON w.work_id = wa.work_id
           WHERE wa.position = 0
             AND w.guid_type = 'doi'
             AND w.in_corpus = 0
             AND w.year IS NOT NULL
             AND LENGTH(wa.surname_normalized) >= 3"""
    )
    return list(cur)


def unify(conn: sqlite3.Connection,
          min_set_score: int = 85,
          min_partial_ratio: int = 75,
          dry_run: bool = False) -> Dict[str, int]:
    conn.row_factory = sqlite3.Row
    doi_ghosts = load_doi_ghosts(conn)
    logger.info("Scanning %d DOI ghost rows", len(doi_ghosts))

    counts = {"merged": 0, "no_match": 0, "ambiguous": 0,
              "both_in_corpus": 0, "no_title": 0}

    for doi_row in doi_ghosts:
        if not doi_row["title"]:
            counts["no_title"] += 1
            continue
        candidates = find_corpus_key_siblings(conn, doi_row)
        if not candidates:
            counts["no_match"] += 1
            continue
        match, status, n_qual = pick_match(
            doi_row, candidates, min_set_score, min_partial_ratio,
        )
        if status == "ambiguous":
            counts["ambiguous"] += 1
            logger.info(
                "ambiguous %s (cites=%d, %d qualifying of %d candidates)",
                doi_row["work_id"][:60], doi_row["cites"],
                n_qual, len(candidates),
            )
            continue
        if status == "no_match" or match is None:
            counts["no_match"] += 1
            continue

        try:
            survivor_id, dup_id, reason = pick_survivor(doi_row, match)
        except ValueError:
            counts["both_in_corpus"] += 1
            logger.warning(
                "both in_corpus — skipping: %s vs %s",
                doi_row["work_id"][:50], match["work_id"][:50],
            )
            continue

        counts["merged"] += 1
        logger.info(
            "merge [%s] %s (%d cites) + %s (%d cites) -> %s",
            reason,
            doi_row["work_id"][:50], doi_row["cites"],
            match["work_id"][:50], match["cites"],
            survivor_id[:50],
        )
        if not dry_run:
            merge_duplicate(conn, survivor_id, dup_id, doi_row["doi"])

    if not dry_run:
        conn.commit()
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("output_dir", type=Path,
                        help="Corpus output directory (used to locate biblio_authority.sqlite)")
    parser.add_argument("-d", "--db", type=Path, default=None,
                        help="Biblio authority SQLite "
                             "(default: <output_dir>/biblio_authority.sqlite)")
    parser.add_argument("--min-set-score", type=int, default=85,
                        help="token_set_ratio threshold for merge (default: 85)")
    parser.add_argument("--min-partial-ratio", type=int, default=75,
                        dest="min_partial_ratio",
                        help="partial_ratio threshold for merge (default: 75)")
    parser.add_argument("--dry-run", action="store_true",
                        help="scan and log would-be merges without writing")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if not _HAS_RAPIDFUZZ:
        sys.stderr.write(
            "rapidfuzz is required. Install with `pip install rapidfuzz`.\n"
        )
        return 1
    if args.db is None:
        args.db = args.output_dir / "biblio_authority.sqlite"

    if not args.db.exists():
        logger.error("DB not found: %s", args.db)
        return 1

    conn = sqlite3.connect(args.db)
    try:
        counts = unify(
            conn,
            min_set_score=args.min_set_score,
            min_partial_ratio=args.min_partial_ratio,
            dry_run=args.dry_run,
        )
    finally:
        conn.close()

    logger.info("═══ DOI/corpus_key unify %s ═══",
                "dry-run" if args.dry_run else "complete")
    for k in ("merged", "ambiguous", "no_match", "no_title", "both_in_corpus"):
        logger.info("  %-16s %d", k, counts.get(k, 0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
