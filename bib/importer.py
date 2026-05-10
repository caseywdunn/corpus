#!/usr/bin/env python3
"""Apply hand-edited BibTeX back into biblio_authority.sqlite (#26 import half).

Reads a BibTeX file (typically the output of ``bib_export.py``,
hand-edited by the user) and applies each entry's fields onto the
matching ``works`` row. Authors are replaced wholesale per entry.

The .bib file is the source of truth for user edits — there is no
separate override table. The expected workflow is:

    # 1. Build the authority DB from per-paper Grobid output
    python build_biblio_authority.py output

    # 2. Export current state to BibTeX
    python bib_export.py output -o my_corpus.bib

    # 3. Hand-edit my_corpus.bib (fix titles, dates, authors)

    # 4. Apply edits back
    python bib_import.py output my_corpus.bib

    # 5. (later) After a rebuild, re-apply by re-running step 4 on the
    #    same .bib. Edits are not lost as long as the .bib is preserved.

Identity-matching priority for each BibTeX entry:

    1. ``corpus_hash`` field (the round-trip key emitted by bib_export)
    2. ``doi`` field
    3. otherwise → warn, skip the entry

Behavior:

* Atomic — all updates apply in a single transaction. A parse error
  or any per-entry failure rolls back the whole import.
* ``--dry-run`` reports the diff (changed records + fields) without
  writing.
* Reconcile note — if the user changes author surname or year for a
  corpus paper, ``reconcile_corpus_to_biblio.py`` may want to re-run
  to update the ghost-merge. The import logs a hint when that's
  likely needed.

Out of scope (tracked in #26 follow-ups):

* Lossless ``other_fields`` blob for non-standard BibTeX keys.
* Cite-key generation for new ghost works.
* MCP tool wrapper around export.
"""
from __future__ import annotations

import argparse
import logging
import sqlite3
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .parser import _split_authors, _strip_outer_braces, parse_bibtex

logger = logging.getLogger("bib_import")


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------


def _normalize_doi(s: str) -> str:
    s = s.strip().lower()
    for prefix in (
        "https://doi.org/", "http://doi.org/",
        "https://dx.doi.org/", "http://dx.doi.org/",
        "doi:",
    ):
        if s.startswith(prefix):
            s = s[len(prefix):]
    return s


def find_matching_work_id(
    conn: sqlite3.Connection,
    entry: Dict,
) -> Tuple[Optional[str], str]:
    """Return ``(work_id, match_method)`` for a parsed BibTeX entry.

    Match priority:
      1. corpus_hash field — exact match on ``works.corpus_hash``
      2. DOI field — exact match on ``works.doi`` (after normalization)

    On no match: ``(None, "no_match")``.
    """
    corpus_hash = _strip_outer_braces(entry.get("corpus_hash", "") or "")
    if corpus_hash:
        row = conn.execute(
            "SELECT work_id FROM works WHERE corpus_hash = ?",
            (corpus_hash,),
        ).fetchone()
        if row:
            return row[0], "corpus_hash"

    doi_raw = _strip_outer_braces(entry.get("doi", "") or "")
    if doi_raw:
        doi = _normalize_doi(doi_raw)
        row = conn.execute(
            "SELECT work_id FROM works WHERE doi = ?",
            (doi,),
        ).fetchone()
        if row:
            return row[0], "doi"

    return None, "no_match"


# ---------------------------------------------------------------------------
# Per-entry diff + apply
# ---------------------------------------------------------------------------


# Fields the import touches on the ``works`` row. Author updates go via
# work_authors. Non-listed BibTeX fields (e.g. abstract, volume) are
# silently dropped.
#
# v0.3 (#51 + #54) adds round-trip support for license / licenseurl /
# serve / servereason. The BibTeX field name → works column map covers
# the rename: e.g. ``licenseurl`` (BibTeX flat) → ``license_url`` (SQL).
_WORK_FIELDS = (
    "title", "year", "journal", "doi",
    "license", "license_url", "serve", "serve_reason",
)
_BIBTEX_TO_WORKS = {
    "licenseurl":   "license_url",   # #51 — flat BibTeX → snake_case SQL
    "servereason":  "serve_reason",  # #54
}
_WORKS_TO_BIBTEX = {v: k for k, v in _BIBTEX_TO_WORKS.items()}


def _entry_value(entry: Dict, field: str) -> Optional[str]:
    """Pull a field from a parsed entry, with brace stripping.

    Resolves the BibTeX → works column rename (#51 + #54): callers ask
    for the works-column name (``license_url``, ``serve_reason``); we
    look first under the snake_case key, fall back to the flat BibTeX
    name (``licenseurl``, ``servereason``).
    """
    bib_key = _WORKS_TO_BIBTEX.get(field, field)
    raw = entry.get(field)
    if raw is None or raw == "":
        raw = entry.get(bib_key)
    if raw is None or raw == "":
        return None
    return _strip_outer_braces(str(raw))


def diff_entry_against_work(
    conn: sqlite3.Connection,
    work_id: str,
    entry: Dict,
) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
    """Return ``{field: (db_value, new_value)}`` for fields that differ.

    Only fields present in the entry are considered. Fields whose new
    value equals the DB value are excluded — the diff is what would
    *actually* change.
    """
    cols = ", ".join(_WORK_FIELDS)
    cur = conn.execute(
        f"SELECT {cols} FROM works WHERE work_id = ?",
        (work_id,),
    ).fetchone()
    if cur is None:
        return {}
    db_values: Dict[str, Optional[str]] = {}
    for i, field in enumerate(_WORK_FIELDS):
        v = cur[i]
        if v is None:
            db_values[field] = None
        elif field in {"year", "serve"}:
            db_values[field] = str(v)
        else:
            db_values[field] = v

    changes: Dict[str, Tuple[Optional[str], Optional[str]]] = {}
    for field in _WORK_FIELDS:
        new_val = _entry_value(entry, field)
        if new_val is None:
            continue
        if field == "year":
            try:
                # Compare as ints when both parse — guards against "2010" vs "2010 "
                if str(int(new_val)) == (db_values[field] or ""):
                    continue
            except ValueError:
                pass
        if field == "serve":
            # Normalize {true, 1, yes} vs {false, 0, no}
            new_val = "1" if new_val.strip().lower() in {"1", "true", "yes"} else "0"
        if new_val != (db_values[field] or ""):
            changes[field] = (db_values[field], new_val)

    # Author diff. We compare normalized "Surname, Forename" pairs.
    new_authors = _split_authors(entry.get("author", "") or "")
    db_authors = conn.execute(
        "SELECT surname, COALESCE(forename, '') FROM work_authors "
        "WHERE work_id = ? ORDER BY position",
        (work_id,),
    ).fetchall()
    new_pairs = [(a["surname"], a["forename"]) for a in new_authors]
    db_pairs = [(s, f) for (s, f) in db_authors]
    if new_authors and new_pairs != db_pairs:
        changes["authors"] = (
            "; ".join(f"{s}, {f}" for s, f in db_pairs) or None,
            "; ".join(f"{s}, {f}" for s, f in new_pairs),
        )

    return changes


def apply_entry(
    conn: sqlite3.Connection,
    work_id: str,
    entry: Dict,
) -> int:
    """Apply the entry's fields to ``work_id``. Returns number of fields changed.

    Authors are replaced wholesale (delete + reinsert) when the entry's
    author field differs from the DB's. Caller commits the transaction.
    """
    from .authority import normalize_for_key

    changes = diff_entry_against_work(conn, work_id, entry)
    if not changes:
        return 0

    # Update works.* for the simple fields
    sets: List[str] = []
    params: List = []
    for field in _WORK_FIELDS:
        if field not in changes:
            continue
        new_val: Optional[str] = changes[field][1]
        if field == "year":
            try:
                params.append(int(new_val) if new_val else None)
            except ValueError:
                params.append(None)
        elif field == "doi":
            params.append(_normalize_doi(new_val) if new_val else None)
        elif field == "serve":
            params.append(int(new_val))
        else:
            params.append(new_val)
        sets.append(f"{field} = ?")
    if sets:
        sets.append("updated_at = ?")
        params.append(time.time())
        params.append(work_id)
        conn.execute(
            f"UPDATE works SET {', '.join(sets)} WHERE work_id = ?",
            tuple(params),
        )

    # Replace author list when changed
    if "authors" in changes:
        new_authors = _split_authors(entry.get("author", "") or "")
        conn.execute("DELETE FROM work_authors WHERE work_id = ?", (work_id,))
        for i, a in enumerate(new_authors):
            conn.execute(
                """INSERT INTO work_authors
                   (work_id, position, surname, surname_normalized, forename)
                   VALUES (?, ?, ?, ?, ?)""",
                (
                    work_id, i,
                    a["surname"], normalize_for_key(a["surname"]),
                    a["forename"] or None,
                ),
            )

    return len(changes)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


# #67: sidecar path for `corpus bib import` invoked before the first
# `corpus run`. Picked up at the end of bib.authority.main() so the
# operator can curate ahead of the first build.
PENDING_SIDECAR = ".pending_bib_import.bib"


def stage_pending_bib(output_dir: Path, bib_path: Path) -> Path:
    """Copy ``bib_path`` to ``<output_dir>/.pending_bib_import.bib`` so
    the next ``corpus run`` (specifically bib.authority's tail) applies
    it automatically. Used by `corpus bib import` when the authority DB
    doesn't exist yet (#67)."""
    import shutil
    output_dir.mkdir(parents=True, exist_ok=True)
    dst = output_dir / PENDING_SIDECAR
    shutil.copy(bib_path, dst)
    return dst


def apply_pending_bib(output_dir: Path, db_path: Path) -> Optional[Dict[str, int]]:
    """If a sidecar from #67 is present, apply it against ``db_path`` and
    remove it. Returns the import counters if a sidecar was applied,
    None otherwise.
    """
    sidecar = output_dir / PENDING_SIDECAR
    if not sidecar.exists():
        return None
    if not db_path.exists():
        return None  # nothing to apply against; defer to a later run
    logger.info("Applying staged bib overrides from %s (#67)", sidecar)
    counters = import_bibtex(db_path, sidecar, dry_run=False)
    sidecar.unlink()
    logger.info("Removed sidecar %s after applying %d update(s)",
                sidecar, counters.get("changed", 0))
    return counters


def import_bibtex(
    db_path: Path,
    bib_path: Path,
    *,
    dry_run: bool = False,
) -> Dict[str, int]:
    """Apply ``bib_path`` to ``db_path``. Returns counters for the run."""
    text = bib_path.read_text(encoding="utf-8")
    entries = parse_bibtex(text)
    logger.info("Parsed %d BibTeX entries from %s", len(entries), bib_path)

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    counters = {
        "entries": len(entries),
        "matched_corpus_hash": 0,
        "matched_doi": 0,
        "no_match": 0,
        "no_changes": 0,
        "changed": 0,
        "fields_updated": 0,
        "authors_changed": 0,
    }

    # Run inside one transaction so a partial failure rolls back.
    try:
        for entry in entries:
            cite_key = entry.get("_key", "?")
            work_id, method = find_matching_work_id(conn, entry)
            if work_id is None:
                counters["no_match"] += 1
                logger.warning(
                    "%s: no match (no corpus_hash or DOI) — skipping",
                    cite_key,
                )
                continue
            counters[f"matched_{method}"] += 1

            changes = diff_entry_against_work(conn, work_id, entry)
            if not changes:
                counters["no_changes"] += 1
                logger.debug("%s: no changes (%s)", cite_key, work_id)
                continue

            counters["changed"] += 1
            counters["fields_updated"] += sum(
                1 for f in changes if f != "authors"
            )
            if "authors" in changes:
                counters["authors_changed"] += 1

            for field, (old, new) in changes.items():
                logger.info(
                    "%s [%s]: %s: %r → %r",
                    cite_key, work_id, field, old, new,
                )

            if not dry_run:
                apply_entry(conn, work_id, entry)

        if dry_run:
            conn.rollback()
            logger.info("Dry-run: no changes written")
        else:
            conn.commit()
            logger.info(
                "Committed: %d works updated (%d fields, %d author lists)",
                counters["changed"],
                counters["fields_updated"],
                counters["authors_changed"],
            )
    except Exception:
        conn.rollback()
        logger.error("Import failed; rolled back transaction")
        raise
    finally:
        conn.close()

    # Reconcile-impact hint
    if not dry_run and counters["changed"] > 0:
        logger.info(
            "Note: if author surnames or years changed for in-corpus papers, "
            "consider re-running reconcile_corpus_to_biblio.py — the "
            "(surname, year) match key may surface different ghosts now."
        )

    return counters


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "output_dir", type=Path,
        help="Corpus output directory containing biblio_authority.sqlite.",
    )
    parser.add_argument(
        "bib_file", type=Path,
        help="Hand-edited BibTeX file produced by bib_export.py.",
    )
    parser.add_argument(
        "--db", type=Path, default=None,
        help="Override the SQLite path "
             "(default: <output_dir>/biblio_authority.sqlite).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Report which records and fields would change without writing.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if not args.bib_file.exists():
        logger.error("BibTeX file %s not found", args.bib_file)
        return 1

    db_path = args.db or (args.output_dir / "biblio_authority.sqlite")
    if not db_path.exists():
        # #67: stage for the next `corpus run` rather than erroring,
        # so operators can curate the bibliography ahead of the first
        # build. The sidecar gets picked up + removed at the tail of
        # bib.authority.main().
        if args.dry_run:
            logger.info(
                "biblio_authority.sqlite does not exist yet; would stage "
                "%s for application on the next `corpus run` (dry-run, "
                "no copy made).", args.bib_file,
            )
            return 0
        sidecar = stage_pending_bib(args.output_dir, args.bib_file)
        logger.info(
            "biblio_authority.sqlite does not exist yet; staged %s → %s. "
            "It will be applied automatically on the next `corpus run`.",
            args.bib_file, sidecar,
        )
        return 0

    counters = import_bibtex(db_path, args.bib_file, dry_run=args.dry_run)

    logger.info(
        "Summary: %d entries (%d matched corpus_hash, %d matched DOI, "
        "%d no-match, %d no-change, %d updated)",
        counters["entries"],
        counters["matched_corpus_hash"],
        counters["matched_doi"],
        counters["no_match"],
        counters["no_changes"],
        counters["changed"],
    )
    return 0 if counters["no_match"] == 0 else 0  # warn-but-don't-fail on no-match


if __name__ == "__main__":
    sys.exit(main())
