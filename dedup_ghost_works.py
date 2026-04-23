#!/usr/bin/env python3
"""Dedupe OCR-variant ghost rows within (surname, year) families.

The Phase 2 cascade in ``build_biblio_authority.py`` (after #2) correctly
creates a new ghost row whenever it sees a materially-different title in
some paper's ``references.json``. But Grobid's reference parser emits
many near-duplicate titles for the same underlying work — the Haeckel
1888 Challenger report alone has ~27 ghosts in the live DB (canonical
"Report on the Siphonophorae" plus "Siphonophorae of the Challenger",
"Rep. sci. Res. H.M.S. Challenger Zool.", empty-title rows, etc.).

The match-time picker in ``reconcile_corpus_to_biblio.py`` already
handles this at the single-PDF level. This script does the same at the
ghost-family level: for each (normalized_surname, year) with ≥ 2
``in_corpus=0`` rows, identify the canonical and merge OCR/fragment
variants into it.

Algorithm per family:

  1. Skip families where the surname is < 4 chars (Grobid-misparse
     single-letter tokens) or the year is missing — nothing good comes
     of collapsing the "p, 1900" family.
  2. Never touch DOI-keyed rows (``guid_type='doi'``) — the DOI-vs-
     corpus_key dedup is its own concern (see issue #4).
  3. Pick the canonical: the row with the most incoming citations that
     also has a substantive title (≥ 20 alpha chars). If no row has a
     substantive title, the family can't be disambiguated — skip.
  4. For every other row in the family:
     - Empty or near-empty title (< 20 alpha chars, typically journal
       fragments like "Rep. sci. Res." or "Discovery Rep") → merge.
     - ``token_set_ratio ≥ 95`` AND ``ratio ≥ 40`` against the canonical
       title → merge. (The threshold is tighter than #2's cascade
       because we're asking a different question: given same surname
       and year, is this candidate a superset-with-extra-punctuation
       of the canonical, or a genuinely different work? Real OCR
       variants score set=100 because the canonical's tokens are all
       present; distinct works like Haeckel 1888's *Prodromus
       Medusarum* vs. his Challenger report score set≤50.)
     - Otherwise → keep separate. Genuinely different works by the same
       (author, year) (Totton 1965's Synopsis vs. Lensia paper; Haeckel
       1888's Challenger report vs. Prodromus Medusarum) score below
       the threshold and survive.

Merge = redirect citation edges from the dup onto the canonical, copy
authors/aliases the canonical doesn't already carry, delete the dup.
Idempotent — re-runs find nothing to do.

Usage:
    python dedup_ghost_works.py
    python dedup_ghost_works.py --dry-run
    python dedup_ghost_works.py --min-set-score 90 --min-ratio 65
"""

from __future__ import annotations

import argparse
import logging
import re
import sqlite3
import sys
import unicodedata
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from rapidfuzz import fuzz
    _HAS_RAPIDFUZZ = True
except ImportError:
    fuzz = None  # type: ignore
    _HAS_RAPIDFUZZ = False

logger = logging.getLogger("dedup_ghosts")

DEFAULT_DB = Path(__file__).resolve().parent / "resources" / "biblio_authority.sqlite"


# ── Normalization (kept in sync with build_biblio_authority.py) ─────

def normalize_for_key(s: str) -> str:
    s = unicodedata.normalize("NFD", s)
    s = "".join(c for c in s if unicodedata.category(c) != "Mn")
    s = s.lower().strip()
    s = re.sub(r"[^\w\s]", "", s)
    s = re.sub(r"\s+", " ", s)
    return s


def _alpha_len(s: str) -> int:
    """Count of alphabetic characters in ``s`` — a better length signal
    than raw string length when titles are mostly publisher abbreviations
    and dots."""
    return sum(1 for c in s if c.isalpha())


# ── Family enumeration ──────────────────────────────────────────────

def load_families(conn: sqlite3.Connection,
                  min_surname_len: int = 4) -> Dict[Tuple[str, int], List[sqlite3.Row]]:
    """Return {(surname_normalized, year): [rows]} for every (surname,
    year) family where there's at least one ``in_corpus=0`` ghost to
    consider merging.

    Includes ``in_corpus=1`` rows in the family so they can be picked
    as the canonical target — a ghost OCR variant of Totton 1965's
    *Synopsis* should flow into the corpus-paper row for that work,
    not into a sibling ghost. Corpus-paper rows are never merged out
    of their position (the caller enforces this).

    Only position=0 authors. DOI-keyed rows are excluded — issue #4
    handles DOI-vs-corpus_key merges separately.
    """
    cur = conn.execute(
        """SELECT wa.surname_normalized AS surname, w.year AS year,
                  w.work_id, w.title, w.guid_type, w.in_corpus,
                  (SELECT COUNT(*) FROM citations c WHERE c.cited_work_id = w.work_id) AS cites
           FROM works w
           JOIN work_authors wa ON w.work_id = wa.work_id
           WHERE wa.position = 0
             AND w.guid_type != 'doi'
             AND w.year IS NOT NULL
             AND LENGTH(wa.surname_normalized) >= ?""",
        (min_surname_len,),
    )
    families: Dict[Tuple[str, int], List[sqlite3.Row]] = {}
    for row in cur:
        key = (row["surname"], row["year"])
        families.setdefault(key, []).append(row)
    # Keep only families that have at least one ghost to consider merging
    # AND more than one row total (else there's nothing to dedup).
    return {
        k: v for k, v in families.items()
        if len(v) >= 2 and any(r["in_corpus"] == 0 for r in v)
    }


# ── Canonical picker ────────────────────────────────────────────────

def pick_canonical(family: List[sqlite3.Row],
                   min_alpha: int = 20) -> Optional[sqlite3.Row]:
    """Return the highest-cites row in the family whose title has at
    least ``min_alpha`` alphabetic chars — i.e., the row most likely
    to represent the real underlying work. Returns None if no row in
    the family clears the bar (family can't be disambiguated)."""
    substantive = [r for r in family if r["title"] and _alpha_len(r["title"]) >= min_alpha]
    if not substantive:
        return None
    return max(substantive, key=lambda r: r["cites"])


# ── Merge decision ──────────────────────────────────────────────────

def should_merge(candidate: sqlite3.Row, canonical: sqlite3.Row,
                 min_set_score: int = 95, min_ratio_score: int = 40,
                 fragment_alpha_threshold: int = 20) -> Tuple[bool, str]:
    """Return (merge, reason) for a candidate row against the canonical.

    A candidate is merged when it looks like a duplicate of the
    canonical rather than a distinct work that happens to share
    surname+year.
    """
    cand_title = candidate["title"] or ""
    cand_alpha = _alpha_len(cand_title)
    if cand_alpha == 0:
        return True, "empty_title"
    if cand_alpha < fragment_alpha_threshold:
        return True, "fragment"

    norm_cand = normalize_for_key(cand_title)
    norm_canon = normalize_for_key(canonical["title"])
    set_score = int(fuzz.token_set_ratio(norm_cand, norm_canon))
    ratio_score = int(fuzz.ratio(norm_cand, norm_canon))
    if set_score >= min_set_score and ratio_score >= min_ratio_score:
        return True, f"fuzzy({set_score}/{ratio_score})"
    return False, f"distinct({set_score}/{ratio_score})"


# ── Merge ───────────────────────────────────────────────────────────

def merge_ghost_into_canonical(conn: sqlite3.Connection,
                               dup_work_id: str,
                               canonical_work_id: str) -> None:
    """Redirect citations, copy authors+aliases, delete the dup row.

    Unlike ``reconcile_corpus_to_biblio.merge_phase1_into_ghost`` this
    merges two ghost rows (both ``in_corpus=0``) so we don't touch
    ``in_corpus`` or ``corpus_hash`` on the survivor.
    """
    if dup_work_id == canonical_work_id:
        return

    # 1. Redirect incoming citations (dup was cited -> canonical is cited).
    conn.execute(
        "UPDATE OR IGNORE citations SET cited_work_id = ? WHERE cited_work_id = ?",
        (canonical_work_id, dup_work_id),
    )
    conn.execute("DELETE FROM citations WHERE cited_work_id = ?", (dup_work_id,))

    # 2. Redirect outgoing citations (rare for ghosts but defensive).
    conn.execute(
        "UPDATE OR IGNORE citations SET citing_work_id = ? WHERE citing_work_id = ?",
        (canonical_work_id, dup_work_id),
    )
    conn.execute("DELETE FROM citations WHERE citing_work_id = ?", (dup_work_id,))

    # 3. Copy authors the canonical doesn't already have (appending).
    cur = conn.execute(
        "SELECT position, surname, surname_normalized, forename "
        "FROM work_authors WHERE work_id = ? ORDER BY position",
        (dup_work_id,),
    )
    dup_authors = list(cur)
    cur = conn.execute(
        "SELECT MAX(position) FROM work_authors WHERE work_id = ?", (canonical_work_id,),
    )
    row = cur.fetchone()
    next_pos = (row[0] + 1) if row and row[0] is not None else 0
    for _, surname, surname_norm, forename in dup_authors:
        dup_check = conn.execute(
            "SELECT 1 FROM work_authors WHERE work_id = ? AND surname_normalized = ?",
            (canonical_work_id, surname_norm),
        ).fetchone()
        if dup_check:
            continue
        conn.execute(
            """INSERT OR IGNORE INTO work_authors
               (work_id, position, surname, surname_normalized, forename)
               VALUES (?, ?, ?, ?, ?)""",
            (canonical_work_id, next_pos, surname, surname_norm, forename),
        )
        next_pos += 1

    # 4. Copy aliases that aren't already on the canonical.
    conn.execute(
        """INSERT OR IGNORE INTO work_aliases (alias_key, work_id)
           SELECT alias_key, ? FROM work_aliases WHERE work_id = ?""",
        (canonical_work_id, dup_work_id),
    )

    # 5. Drop dependent rows then the dup itself.
    conn.execute("DELETE FROM work_aliases WHERE work_id = ?", (dup_work_id,))
    conn.execute("DELETE FROM work_authors WHERE work_id = ?", (dup_work_id,))
    conn.execute("DELETE FROM works WHERE work_id = ?", (dup_work_id,))


# ── Main loop ───────────────────────────────────────────────────────

def dedup(conn: sqlite3.Connection,
          min_set_score: int = 95,
          min_ratio_score: int = 40,
          fragment_alpha_threshold: int = 20,
          canonical_alpha_threshold: int = 20,
          min_surname_len: int = 4,
          dry_run: bool = False) -> Dict[str, int]:
    conn.row_factory = sqlite3.Row
    families = load_families(conn, min_surname_len=min_surname_len)
    logger.info("Scanning %d (surname, year) families with ≥ 2 ghosts",
                len(families))

    counts = {"merged": 0, "kept_separate": 0,
              "skipped_no_canonical": 0, "families_touched": 0}

    for (surname, year), rows in sorted(families.items(),
                                        key=lambda kv: (-len(kv[1]), kv[0])):
        canonical = pick_canonical(rows, min_alpha=canonical_alpha_threshold)
        if canonical is None:
            counts["skipped_no_canonical"] += 1
            logger.debug("skip %s %d: no substantive-title canonical",
                         surname, year)
            continue

        family_merged = 0
        for r in rows:
            if r["work_id"] == canonical["work_id"]:
                continue
            if r["in_corpus"] == 1:
                # Never merge a corpus-paper row out of its position —
                # that would break reconcile_corpus_to_biblio's
                # invariant that an in_corpus=1 row stays where it is.
                continue
            merge, reason = should_merge(
                r, canonical,
                min_set_score=min_set_score,
                min_ratio_score=min_ratio_score,
                fragment_alpha_threshold=fragment_alpha_threshold,
            )
            if merge:
                counts["merged"] += 1
                family_merged += 1
                logger.info(
                    "merge %s %d: %s -> %s [%s]",
                    surname, year, r["work_id"][:60],
                    canonical["work_id"][:60], reason,
                )
                if not dry_run:
                    merge_ghost_into_canonical(
                        conn, r["work_id"], canonical["work_id"],
                    )
            else:
                counts["kept_separate"] += 1
                logger.debug(
                    "keep %s %d: %s != %s [%s]",
                    surname, year, r["work_id"][:60],
                    canonical["work_id"][:60], reason,
                )
        if family_merged:
            counts["families_touched"] += 1

    if not dry_run:
        conn.commit()
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-d", "--db", type=Path, default=DEFAULT_DB,
                        help=f"Biblio authority SQLite (default: {DEFAULT_DB})")
    parser.add_argument("--min-set-score", type=int, default=95,
                        help="token_set_ratio threshold for fuzzy merge (default: 95)")
    parser.add_argument("--min-ratio", type=int, default=40, dest="min_ratio_score",
                        help="ratio threshold for fuzzy merge (default: 40)")
    parser.add_argument("--fragment-alpha", type=int, default=20,
                        help="candidates with fewer than this many alpha chars are "
                             "merged as fragments (default: 20)")
    parser.add_argument("--min-surname-len", type=int, default=4,
                        help="skip families whose surname has fewer than this many "
                             "chars (Grobid-misparse guard; default: 4)")
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
    if not args.db.exists():
        logger.error("DB not found: %s", args.db)
        return 1

    conn = sqlite3.connect(args.db)
    try:
        counts = dedup(
            conn,
            min_set_score=args.min_set_score,
            min_ratio_score=args.min_ratio_score,
            fragment_alpha_threshold=args.fragment_alpha,
            min_surname_len=args.min_surname_len,
            dry_run=args.dry_run,
        )
    finally:
        conn.close()

    logger.info("═══ Ghost dedup %s ═══",
                "dry-run" if args.dry_run else "complete")
    for k in ("merged", "kept_separate", "families_touched",
              "skipped_no_canonical"):
        logger.info("  %-22s %d", k, counts.get(k, 0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
