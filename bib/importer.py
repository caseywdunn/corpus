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

logger = logging.getLogger("corpus.bib.import")


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
      1. corpus_hash field — document membership (legacy scalar fallback)
      2. DOI field — exact match on ``works.doi`` (after normalization)
      3. work_id field — exact match on ``works.work_id`` (#100)

    On no match: ``(None, "no_match")``.
    """
    corpus_hash = _strip_outer_braces(entry.get("corpus_hash", "") or "")
    if corpus_hash:
        from .documents import find_work
        work_id = find_work(conn, corpus_hash)
        if work_id:
            return work_id, "corpus_hash"

    doi_raw = _strip_outer_braces(entry.get("doi", "") or "")
    if doi_raw:
        doi = _normalize_doi(doi_raw)
        row = conn.execute(
            "SELECT work_id FROM works WHERE doi = ?",
            (doi,),
        ).fetchone()
        if row:
            return row[0], "doi"

    # 3. work_id field — exact match on ``works.work_id`` (#100). The
    # exporter writes ``work_id`` into every entry, so a hand-curated
    # cited reference with neither corpus_hash nor DOI (common for
    # pre-DOI literature) can still be matched back to its row and
    # stamped as bib provenance.
    work_id_field = _strip_outer_braces(entry.get("work_id", "") or "")
    if work_id_field:
        row = conn.execute(
            "SELECT work_id FROM works WHERE work_id = ?",
            (work_id_field,),
        ).fetchone()
        if row:
            return row[0], "work_id"

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
    "ocrlang",   # #176 — flat BibTeX name, no rename needed
    "ocrmode",   # #186 — force | redo | skip-text
    "doclang", "pagemap",   # #214 — likewise flat, and read by nothing
    "keeppages",            # #188 — flat too, and this one acts
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


def _document_target(conn, work_id, entry):
    """Resolve local policy edits without guessing among several PDFs."""
    from .documents import DOCUMENT_FIELDS, has_memberships
    sha = _strip_outer_braces(entry.get("corpus_hash", "") or "")
    if not has_memberships(conn):
        return sha
    members = [r[0] for r in conn.execute(
        "SELECT corpus_hash FROM work_documents WHERE work_id=?", (work_id,))]
    local_fields = any(_entry_value(entry, key) is not None for key in DOCUMENT_FIELDS)
    if sha:
        if members and sha not in members and local_fields:
            raise ValueError("Document-local curation corpus_hash is not a member of the matched work")
        return sha
    if len(members) > 1 and local_fields:
        raise ValueError("Document-local curation fields require corpus_hash for a multi-document work")
    return members[0] if len(members) == 1 else ""


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
    have = {row[1] for row in conn.execute("PRAGMA table_info(works)")}
    # Dry-run against an older authority DB must remain read-only. Missing
    # nullable fields compare as NULL; a real import migrates before calling
    # this helper.
    cols = ", ".join(
        field if field in have else f"NULL AS {field}"
        for field in _WORK_FIELDS
    )
    cur = conn.execute(
        f"SELECT {cols} FROM works WHERE work_id = ?",
        (work_id,),
    ).fetchone()
    if cur is None:
        return {}
    from .documents import document_fields, document_metadata
    corpus_hash = _document_target(conn, work_id, entry)
    meta = document_metadata(conn, corpus_hash) if corpus_hash else None
    local_values = document_fields(meta) if meta is not None else {}
    db_values: Dict[str, Optional[str]] = {}
    for i, field in enumerate(_WORK_FIELDS):
        v = local_values.get(field) if field in local_values else cur[i]
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

    now = time.time()

    changes = diff_entry_against_work(conn, work_id, entry)
    if not changes:
        # No field-level diff, but the entry IS present in the user's
        # authoritative .bib — stamp bib provenance anyway (#100).
        # Previously we returned early with bib_imported_at left NULL, so
        # re-importing an unchanged .bib left every reference classified
        # as grobid_reconciled/unresolved and warned by format_citation.
        conn.execute(
            "UPDATE works SET updated_at = ?, bib_imported_at = ? "
            "WHERE work_id = ?",
            (now, now, work_id),
        )
        return 0

    # Update works.* for the simple fields, plus always touch
    # updated_at + bib_imported_at (#79) when any change applies — so an
    # author-only change still records that the importer touched this
    # row, and the format_citation MCP tool can distinguish "from user
    # .bib" from "grobid-reconciled" by checking bib_imported_at IS NOT
    # NULL.
    sets: List[str] = []
    params: List = []
    from .documents import DOCUMENT_FIELDS, document_metadata, update_document_fields
    corpus_hash = _document_target(conn, work_id, entry)
    local = document_metadata(conn, corpus_hash) if corpus_hash else None
    local_changes = {}
    for field in _WORK_FIELDS:
        if field not in changes:
            continue
        new_val: Optional[str] = changes[field][1]
        if local is not None and field in DOCUMENT_FIELDS:
            local_changes[field] = int(new_val) if field == "serve" else new_val
            continue
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
    sets.extend(["updated_at = ?", "bib_imported_at = ?"])
    params.extend([now, now])
    params.append(work_id)
    conn.execute(
        f"UPDATE works SET {', '.join(sets)} WHERE work_id = ?",
        tuple(params),
    )
    if local_changes:
        if "license" in local_changes:
            from .authority import derive_publishable
            publishable, source = derive_publishable(local_changes["license"], local.get("year"))
            local_changes.update(publishable=publishable, license_source=source)
        update_document_fields(conn, corpus_hash, local_changes)
    if "year" in changes:
        from .documents import has_memberships
        from .authority import derive_publishable
        if has_memberships(conn):
            year = conn.execute("SELECT year FROM works WHERE work_id=?", (work_id,)).fetchone()[0]
            for (sha,) in conn.execute("SELECT corpus_hash FROM work_documents WHERE work_id=?", (work_id,)).fetchall():
                meta = document_metadata(conn, sha)
                publishable, source = derive_publishable(meta.get("license"), year)
                update_document_fields(conn, sha, {"year": year, "publishable": publishable, "license_source": source})

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
        "matched_work_id": 0,
        "no_match": 0,
        "no_changes": 0,
        "changed": 0,
        "fields_updated": 0,
        "authors_changed": 0,
    }

    # Run inside one transaction so a partial failure rolls back.
    try:
        # Bring older authority DBs up to the current nullable-column shape
        # before a real import, so adding ocrmode does not require a preceding
        # corpus rebuild (#186). Dry-run stays strictly read-only;
        # diff_entry_against_work represents an absent nullable column as NULL.
        if not dry_run:
            from .authority import _migrate_works_columns
            _migrate_works_columns(conn)
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
                # #142 — still call apply_entry. It has a no-diff branch
                # that stamps bib_imported_at, which is #100's documented
                # intent ("the entry IS present in the user's
                # authoritative .bib"), but this loop used to `continue`
                # here and so never reached it — making that branch
                # unreachable dead code.
                #
                # The consequence was severe and invisible: a full
                # round-trip re-import stamped only the handful of entries
                # that happened to have field edits in the reported build,
                # so CorpusIndex.provenance() kept
                # returning grobid_reconciled for the rest and
                # format_citations kept emitting "generated via
                # reconciliation, check if correct" on works the user had
                # curated by hand (#152). And because an authority-DB
                # rebuild clears bib_imported_at and forces a re-import,
                # curation was defeated again on every rebuild.
                if not dry_run:
                    apply_entry(conn, work_id, entry)
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
    # No-match entries are warned about above but are not a failure: a
    # curated .bib legitimately contains works that aren't in this
    # corpuscle. (This was written as `0 if no_match == 0 else 0` — both
    # branches returning 0 — which read as though a non-zero exit had been
    # intended and lost. It hadn't; simplified to say so plainly.)
    return 0


if __name__ == "__main__":
    sys.exit(main())
