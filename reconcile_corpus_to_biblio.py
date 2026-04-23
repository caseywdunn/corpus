#!/usr/bin/env python3
"""Reconcile corpus papers with ghost cited-reference rows.

``build_biblio_authority.py`` Phase 1 mints a work_id for each corpus
PDF from its Grobid-extracted metadata. On scanned old monographs
Grobid routinely misparses the header (running banners → authors,
taxonomic body text → authors, empty titles, wrong years), so the
Phase-1 work_id for those PDFs is nonsense.

Phase 2 parses references out of other papers' bibliographies, which
works well; the cited-reference rows for those same monographs carry
correct surname/year/title. But the Phase-1 key and the Phase-2 key
never collide, so the most-cited papers in the corpus (Totton 1965,
Totton 1954, Haeckel 1888, Moser 1925, …) show up as ``in_corpus=0``
ghost rows with high ``cited_by_count``, while the actual PDFs sit
under garbage keys with zero incoming citations.

This script merges each affected in-corpus row onto the correct ghost.

Algorithm for each corpus paper with zero incoming citations:

  1. Derive ``(surname, year)`` from the filename in ``metadata.json``
     using a simple regex (``Totton1965a.pdf`` → Totton, 1965). Fall
     back to Grobid's parsed year if the filename has none.
  2. Find ghost candidates: ``in_corpus=0`` rows whose position-0 author
     surname matches and whose ``year`` matches.
  3. Disambiguate via ``rapidfuzz.fuzz.partial_ratio`` between each
     candidate's ``title`` and the first 3 kB of extracted text
     from ``text.json`` — the title-page region, before the
     references section. This handles cases like Totton 1965
     where he has two papers that year, and avoids false positives
     where a paper self-cites in its bibliography.
  4. Require a confident match (score ≥ ``--min-score``, default 80)
     and, when there are multiple candidates, a clear margin over the
     runner-up (≥ ``--margin``, default 15).
  5. Merge: redirect outgoing citations from the Phase-1 row onto the
     ghost, copy authors/aliases, flip the ghost's ``in_corpus=1`` +
     ``corpus_hash``, then delete the Phase-1 row.

Idempotent: once a row is merged, the ghost is ``in_corpus=1`` and
therefore excluded from the candidate set, and the merged-into corpus
paper now has incoming citations and is excluded from the scan.

Usage:
    python reconcile_corpus_to_biblio.py /path/to/output
    python reconcile_corpus_to_biblio.py /path/to/output --dry-run
    python reconcile_corpus_to_biblio.py /path/to/output --min-score 85
"""

from __future__ import annotations

import argparse
import json
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


def _require_rapidfuzz() -> None:
    if not _HAS_RAPIDFUZZ:
        sys.stderr.write(
            "rapidfuzz is required for reconciliation. Install with "
            "`pip install rapidfuzz` or `conda install -c conda-forge rapidfuzz`.\n"
        )
        sys.exit(1)

logger = logging.getLogger("reconcile")

DEFAULT_DB = Path(__file__).resolve().parent.parent / "resources" / "biblio_authority.sqlite"


# ── Normalization (kept in sync with build_biblio_authority.py) ─────

def normalize_for_key(s: str) -> str:
    s = unicodedata.normalize("NFD", s)
    s = "".join(c for c in s if unicodedata.category(c) != "Mn")
    s = s.lower().strip()
    s = re.sub(r"[^\w\s]", "", s)
    s = re.sub(r"\s+", " ", s)
    return s


# ── Filename parsing ────────────────────────────────────────────────

# Captures a leading alphabetic surname and the first 4-digit year.
# Examples that match:
#   Totton1965a.pdf           -> Totton, 1965
#   Totton_1954.pdf           -> Totton, 1954
#   Dunn_Pugh_Haddock2005.pdf -> Dunn, 2005
#   Dunn-etal2005_Marrus.pdf  -> Dunn, 2005
#   Haeckel1888.pdf           -> Haeckel, 1888
_FILENAME_RE = re.compile(r"^([A-Za-zÀ-ÿ]+)[^0-9]*?(\d{4})")


def derive_from_filename(filename: str) -> Tuple[str, Optional[int]]:
    """Derive (surname, year) from a filename. Returns ("", None) on failure.

    Only plausible years (1600–2100) are accepted; otherwise year is None.
    """
    if not filename:
        return "", None
    stem = Path(filename).stem
    m = _FILENAME_RE.match(stem)
    if not m:
        return "", None
    surname = m.group(1)
    year = int(m.group(2))
    if year < 1600 or year > 2100:
        return surname, None
    return surname, year


# ── Text loading ────────────────────────────────────────────────────

def load_first_pages_text(doc_dir: Path, max_chars: int = 3000) -> str:
    """Return the first ``max_chars`` of the paper's extracted plain text."""
    text_path = doc_dir / "text.json"
    if not text_path.exists():
        return ""
    try:
        data = json.loads(text_path.read_text())
    except (json.JSONDecodeError, OSError):
        return ""
    text = data.get("text", "") if isinstance(data, dict) else ""
    return text[:max_chars] if text else ""


def load_filename_and_year(doc_dir: Path) -> Tuple[str, Optional[int]]:
    meta_path = doc_dir / "metadata.json"
    if not meta_path.exists():
        return "", None
    try:
        meta = json.loads(meta_path.read_text())
    except (json.JSONDecodeError, OSError):
        return "", None
    return meta.get("filename", "") or "", meta.get("year")


# ── Candidate lookup ────────────────────────────────────────────────

def find_candidates(conn: sqlite3.Connection, surname: str,
                    year: int) -> List[sqlite3.Row]:
    """Return ghost works (in_corpus=0) whose position-0 author surname
    and year match. The position=0 constraint matches how cited-reference
    keys are minted and is how corpus papers are typically cited.

    ``incoming_cites`` is included so the picker can break fuzzy-title
    ties by preferring the most-cited variant — the canonical ghost
    when Phase 2 has produced many OCR-variant rows for the same work.
    """
    norm = normalize_for_key(surname)
    cur = conn.execute(
        """SELECT w.work_id, w.title, w.year, w.journal,
               (SELECT COUNT(*) FROM citations c
                 WHERE c.cited_work_id = w.work_id) AS incoming_cites
           FROM works w
           JOIN work_authors wa ON w.work_id = wa.work_id
           WHERE wa.surname_normalized = ? AND w.year = ?
             AND wa.position = 0
             AND w.in_corpus = 0""",
        (norm, year),
    )
    return list(cur)


def score_candidates(candidates: List[sqlite3.Row],
                     first_pages_text: str) -> List[Tuple[str, int, str, int]]:
    """Return [(work_id, score, title, incoming_cites), ...] sorted best-first.

    Uses ``partial_ratio`` so a short title scores against its best-matching
    substring in the opening pages rather than being penalized for the
    page text being much larger. Secondary sort by ``incoming_cites``
    descending so the canonical variant (most-cited) outranks obvious
    near-duplicates when the title score ties.
    """
    norm_pages = normalize_for_key(first_pages_text)
    scored: List[Tuple[str, int, str, int]] = []
    for cand in candidates:
        cand_title = cand["title"] or ""
        incoming = cand["incoming_cites"] if "incoming_cites" in cand.keys() else 0
        if not cand_title or not norm_pages:
            scored.append((cand["work_id"], 0, cand_title, incoming))
            continue
        norm_title = normalize_for_key(cand_title)
        if not norm_title:
            scored.append((cand["work_id"], 0, cand_title, incoming))
            continue
        score = int(fuzz.partial_ratio(norm_title, norm_pages))
        scored.append((cand["work_id"], score, cand_title, incoming))
    # Rank by (title_score desc, incoming_cites desc) so ties break on
    # the canonical (most-cited) variant.
    scored.sort(key=lambda t: (t[1], t[3]), reverse=True)
    return scored


def pick_winner(scored: List[Tuple[str, int, str, int]],
                min_score: int, margin: int,
                cites_dominance: int = 2) -> Tuple[str, Optional[str], int]:
    """Return (status, winner_work_id_or_None, score).

    status ∈ {"matched", "low_score", "ambiguous", "no_candidates"}

    Phase 2 produces multiple OCR-variant ghosts for the same work
    (e.g., 18 Totton 1965 rows, 14 Totton 1954 rows). Within each
    variant family the canonical one has dramatically more incoming
    citations, while a shorter-title sibling may score slightly
    higher on ``partial_ratio`` because its full length fits as a
    substring of the PDF text.

    Therefore: take all candidates whose title score is within
    ``margin`` of the top, then pick the most-cited among them. That
    picks the canonical ghost over same-score shorter variants. If
    the top citation count doesn't dominate the runner-up by at least
    ``cites_dominance`` ×, return ``ambiguous``.
    """
    if not scored:
        return ("no_candidates", None, 0)
    top_score = scored[0][1]
    if top_score < min_score:
        return ("low_score", scored[0][0], top_score)

    # Top tier = candidates within ``margin`` of the best score.
    tier = [c for c in scored if top_score - c[1] < margin and c[1] >= min_score]
    # Re-rank by incoming cites (canonical is usually the most-cited).
    tier.sort(key=lambda c: c[3], reverse=True)

    winner = tier[0]
    if len(tier) == 1:
        return ("matched", winner[0], winner[1])
    runner = tier[1]
    # Require citation dominance; otherwise truly ambiguous.
    if winner[3] > 0 and winner[3] >= max(cites_dominance, 1) * max(runner[3], 1):
        return ("matched", winner[0], winner[1])
    return ("ambiguous", winner[0], winner[1])


# ── Merge ───────────────────────────────────────────────────────────

def merge_phase1_into_ghost(conn: sqlite3.Connection,
                            phase1_work_id: str,
                            ghost_work_id: str,
                            corpus_hash: str) -> None:
    """Merge a Phase-1 corpus row onto an existing ghost cited-reference row.

    Assumes (and is protected by our selection predicate) that the
    Phase-1 row has zero incoming citations. Redirects outgoing
    citations, copies authors/aliases that aren't already on the
    ghost, flips the ghost to ``in_corpus=1``, then deletes the
    Phase-1 row.
    """
    if phase1_work_id == ghost_work_id:
        # Already reconciled on a prior run; nothing to do.
        return

    # 1. Redirect outgoing citations: the corpus paper's references.
    # Use OR IGNORE to survive any (rare) PK collisions, then delete
    # leftover rows still pointing at the phase1 id.
    conn.execute(
        "UPDATE OR IGNORE citations SET citing_work_id = ? WHERE citing_work_id = ?",
        (ghost_work_id, phase1_work_id),
    )
    conn.execute(
        "DELETE FROM citations WHERE citing_work_id = ?", (phase1_work_id,),
    )

    # 2. Copy authors the ghost doesn't already have, appending at the end.
    cur = conn.execute(
        "SELECT position, surname, surname_normalized, forename "
        "FROM work_authors WHERE work_id = ? ORDER BY position",
        (phase1_work_id,),
    )
    phase1_authors = list(cur)
    cur = conn.execute(
        "SELECT MAX(position) FROM work_authors WHERE work_id = ?", (ghost_work_id,),
    )
    row = cur.fetchone()
    next_pos = (row[0] + 1) if row and row[0] is not None else 0
    for _, surname, surname_norm, forename in phase1_authors:
        dup = conn.execute(
            "SELECT 1 FROM work_authors WHERE work_id = ? AND surname_normalized = ?",
            (ghost_work_id, surname_norm),
        ).fetchone()
        if dup:
            continue
        conn.execute(
            """INSERT OR IGNORE INTO work_authors
               (work_id, position, surname, surname_normalized, forename)
               VALUES (?, ?, ?, ?, ?)""",
            (ghost_work_id, next_pos, surname, surname_norm, forename),
        )
        next_pos += 1

    # 3. Copy aliases that aren't already on the ghost.
    conn.execute(
        """INSERT OR IGNORE INTO work_aliases (alias_key, work_id)
           SELECT alias_key, ? FROM work_aliases WHERE work_id = ?""",
        (ghost_work_id, phase1_work_id),
    )

    # 4. Upgrade the ghost row to be the corpus paper.
    conn.execute(
        """UPDATE works
           SET in_corpus = 1, corpus_hash = ?, source = 'corpus_paper',
               updated_at = ?
           WHERE work_id = ?""",
        (corpus_hash, time.time(), ghost_work_id),
    )

    # 5. Clean up the now-orphan Phase-1 row.
    conn.execute("DELETE FROM work_aliases WHERE work_id = ?", (phase1_work_id,))
    conn.execute("DELETE FROM work_authors WHERE work_id = ?", (phase1_work_id,))
    conn.execute("DELETE FROM works WHERE work_id = ?", (phase1_work_id,))


# ── Main reconcile loop ─────────────────────────────────────────────

def _unreconciled_corpus_papers(conn: sqlite3.Connection) -> List[sqlite3.Row]:
    cur = conn.execute(
        """SELECT w.work_id, w.corpus_hash, w.title, w.year
           FROM works w
           WHERE w.in_corpus = 1
             AND w.corpus_hash IS NOT NULL
             AND NOT EXISTS (
               SELECT 1 FROM citations c WHERE c.cited_work_id = w.work_id
             )"""
    )
    return list(cur)


def reconcile(conn: sqlite3.Connection, output_dir: Path,
              min_score: int = 80, margin: int = 15,
              max_chars: int = 3000,
              dry_run: bool = False) -> Dict[str, int]:
    docs_dir = output_dir / "documents"
    if not docs_dir.is_dir():
        logger.error("Documents dir not found: %s", docs_dir)
        return {}

    conn.row_factory = sqlite3.Row
    targets = _unreconciled_corpus_papers(conn)
    logger.info("%d corpus papers have zero incoming citations; scanning.",
                len(targets))

    counts = {"matched": 0, "low_score": 0, "ambiguous": 0,
              "no_candidates": 0, "no_filename": 0, "missing_text": 0}

    for row in targets:
        phase1_id = row["work_id"]
        corpus_hash = row["corpus_hash"]
        doc_dir = docs_dir / corpus_hash
        if not doc_dir.is_dir():
            counts["missing_text"] += 1
            continue

        filename, grobid_year = load_filename_and_year(doc_dir)
        surname, year = derive_from_filename(filename)
        if not surname:
            counts["no_filename"] += 1
            logger.debug("no filename-derived surname for %s (%s)",
                         corpus_hash, filename or "<empty>")
            continue
        if not year:
            year = grobid_year
        if not year:
            counts["no_filename"] += 1
            logger.debug("no year for %s (filename=%s grobid_year=%s)",
                         corpus_hash, filename, grobid_year)
            continue

        candidates = find_candidates(conn, surname, year)
        if not candidates:
            counts["no_candidates"] += 1
            logger.debug("no ghost candidates for %s (%s %d)",
                         corpus_hash, surname, year)
            continue

        first_pages = load_first_pages_text(doc_dir, max_chars=max_chars)
        scored = score_candidates(candidates, first_pages)
        status, winner_id, score = pick_winner(scored, min_score, margin)

        if status != "matched":
            counts[status] += 1
            logger.info(
                "%s %s %d: %s (best=%s score=%d, %d candidates)",
                corpus_hash, surname, year, status,
                (winner_id or "")[:60], score, len(scored),
            )
            continue

        counts["matched"] += 1
        logger.info(
            "%s %s %d: match %s (score=%d)",
            corpus_hash, surname, year, winner_id[:80], score,
        )
        if not dry_run:
            merge_phase1_into_ghost(conn, phase1_id, winner_id, corpus_hash)

    if not dry_run:
        conn.commit()
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("output_dir", type=Path,
                        help="Corpus output directory (contains documents/<hash>/ subdirs)")
    parser.add_argument("-d", "--db", type=Path, default=DEFAULT_DB,
                        help=f"Biblio authority SQLite (default: {DEFAULT_DB})")
    parser.add_argument("--min-score", type=int, default=80,
                        help="Minimum rapidfuzz.partial_ratio score to accept a match (default: 80)")
    parser.add_argument("--margin", type=int, default=15,
                        help="Required score margin between winner and runner-up "
                             "when there are multiple candidates (default: 15)")
    parser.add_argument("--max-chars", type=int, default=3000,
                        help="How many chars of text.json to search for titles. Keep small "
                             "(default: 3000) — large windows can match a paper's own "
                             "self-citations in its reference list and produce false positives.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Scan and log would-be merges without writing")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    _require_rapidfuzz()

    if not args.db.exists():
        logger.error("DB not found: %s. Run build_biblio_authority.py first.", args.db)
        return 1
    if not args.output_dir.is_dir():
        logger.error("Output dir not found: %s", args.output_dir)
        return 1

    conn = sqlite3.connect(args.db)
    try:
        counts = reconcile(
            conn, args.output_dir,
            min_score=args.min_score, margin=args.margin,
            max_chars=args.max_chars, dry_run=args.dry_run,
        )
    finally:
        conn.close()

    logger.info("═══ Reconciliation %s ═══", "dry-run" if args.dry_run else "complete")
    for k in ("matched", "ambiguous", "low_score", "no_candidates", "no_filename", "missing_text"):
        logger.info("  %-14s %d", k, counts.get(k, 0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
