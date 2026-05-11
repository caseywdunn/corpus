#!/usr/bin/env python3
"""Export biblio_authority.sqlite as BibTeX (#26 export half).

Produces a BibTeX file from the works + work_authors tables in
``<output_dir>/biblio_authority.sqlite``. Default scope: papers
in the corpus (``in_corpus = 1``); pass ``--all-works`` to also
emit cited-but-not-in-corpus rows.

Cite key: derived from first-author surname + year, with a
disambiguator letter (``Smith2010``, ``Smith2010a``, ...) when
multiple works collide. The same cite key is produced on every
run so subsequent imports / re-exports are stable. Future work
(import half) will persist user-edited cite keys verbatim.

Usage::

    python bib_export.py /path/to/output                  # to stdout
    python bib_export.py /path/to/output -o corpus.bib    # to a file
    python bib_export.py /path/to/output --all-works      # incl. cited

The export half is the read-only side of round-trip curation. The
import half — match BibTeX entries back, persist user edits across
``build_biblio_authority.py --rebuild`` — has design subtleties
(matching strategy, conflict resolution, override layering) that
warrant their own session; tracked separately.
"""
from __future__ import annotations

import argparse
import logging
import re
import sqlite3
import sys
import unicodedata
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

logger = logging.getLogger("bib_export")


# ---------------------------------------------------------------------------
# Cite key generation
# ---------------------------------------------------------------------------

_CITE_KEY_BAD_CHARS = re.compile(r"[^A-Za-z0-9]+")


def _ascii_surname(surname: str) -> str:
    """Strip diacritics + non-alphanumerics. ``Pugé`` → ``Puge``."""
    if not surname:
        return ""
    s = unicodedata.normalize("NFD", surname)
    s = "".join(c for c in s if unicodedata.category(c) != "Mn")
    s = _CITE_KEY_BAD_CHARS.sub("", s)
    return s


def _alpha_suffix(n: int) -> str:
    """0 → '', 1 → 'a', 2 → 'b', …, 26 → 'z', 27 → 'aa', …."""
    if n <= 0:
        return ""
    out = ""
    n -= 1
    while True:
        out = chr(ord("a") + n % 26) + out
        n //= 26
        if n == 0:
            return out
        n -= 1


def assign_cite_keys(
    works: Iterable[Tuple[str, Optional[int], str]],
) -> Dict[str, str]:
    """Map each work_id to a deterministic cite key.

    ``works`` is an iterable of ``(work_id, year, first_author_surname)``,
    yielded in stable order (e.g., ORDER BY year, work_id). The first
    work to claim ``Smith2010`` gets it; subsequent collisions become
    ``Smith2010a``, ``Smith2010b``, …. Stable as long as the input
    order is.
    """
    keys: Dict[str, str] = {}
    seen: Dict[str, int] = defaultdict(int)
    for work_id, year, surname in works:
        base_surname = _ascii_surname(surname or "")
        base_year = str(year) if year is not None else "nodate"
        if not base_surname:
            base_surname = "Anon"
        base = f"{base_surname}{base_year}"
        # First occurrence: bare base. Subsequent: append a suffix.
        n = seen[base]
        seen[base] = n + 1
        keys[work_id] = base + _alpha_suffix(n)
    return keys


# ---------------------------------------------------------------------------
# BibTeX rendering
# ---------------------------------------------------------------------------

# Characters that need escaping in BibTeX field values. We keep
# diacritics as-is (modern BibTeX with biber/utf8 handles them); only
# brace-balance hazards and explicit BibTeX delimiters get escaped.
_BIBTEX_FIELD_ESCAPES = (
    ("{", r"\{"),
    ("}", r"\}"),
)


def _escape(value: str) -> str:
    if value is None:
        return ""
    out = str(value)
    for src, dst in _BIBTEX_FIELD_ESCAPES:
        out = out.replace(src, dst)
    return out


def render_entry(
    cite_key: str,
    entry_type: str,
    fields: Dict[str, Optional[str]],
) -> str:
    """One BibTeX entry. Empty / None fields are omitted.

    Field order is the order of the input dict; stable across runs
    on Python ≥3.7.
    """
    lines = [f"@{entry_type}{{{cite_key},"]
    for k, v in fields.items():
        if v is None or v == "":
            continue
        lines.append(f"  {k} = {{{_escape(v)}}},")
    # Drop trailing comma on the last field for stylistic taste —
    # tools tolerate either, but consistent with what bibtex emits.
    if lines[-1].endswith(","):
        lines[-1] = lines[-1][:-1]
    lines.append("}")
    return "\n".join(lines)


def _join_authors(rows: List[Tuple[str, str]]) -> str:
    """[(surname, forename), ...] → 'Surname1, Forename1 and Surname2, Forename2'."""
    parts: List[str] = []
    for surname, forename in rows:
        if forename:
            parts.append(f"{surname}, {forename}")
        else:
            parts.append(surname)
    return " and ".join(parts)


# ---------------------------------------------------------------------------
# Main export
# ---------------------------------------------------------------------------


def export_bibtex(
    db_path: Path,
    *,
    corpus_only: bool = True,
    documents_dir: Optional[Path] = None,
) -> str:
    """Walk the works table and produce BibTeX.

    ``documents_dir`` (optional) lets the export attach ``file = {...}``
    entries pointing at the per-paper artifact directory, matching the
    convention used by the existing ``--bib`` flag on
    ``process_corpus.py``. When None, the file field is omitted.
    """
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    where = "WHERE in_corpus = 1" if corpus_only else ""
    works_sql = f"""
        SELECT work_id, title, year, journal, doi, corpus_hash, in_corpus,
               license, license_url, serve, serve_reason
        FROM works
        {where}
        ORDER BY year, work_id
    """
    rows = list(conn.execute(works_sql))

    # First-author surname per work — used for cite-key generation.
    authors_by_work: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    for r in conn.execute(
        "SELECT work_id, surname, COALESCE(forename, '') AS forename "
        "FROM work_authors ORDER BY work_id, position"
    ):
        authors_by_work[r["work_id"]].append((r["surname"], r["forename"]))

    cite_keys = assign_cite_keys(
        (
            r["work_id"],
            r["year"],
            (authors_by_work[r["work_id"]][0][0] if authors_by_work[r["work_id"]] else ""),
        )
        for r in rows
    )

    out_chunks: List[str] = []
    for r in rows:
        work_id = r["work_id"]
        authors_str = _join_authors(authors_by_work.get(work_id, []))
        # Standard BibTeX entry type — corpus-paper filenames typically
        # live under documents/<HASH>/, so a "misc"/"article" choice
        # is mostly aesthetic. Use "article" for in-corpus rows that
        # have a journal, "misc" otherwise.
        entry_type = "article" if (r["journal"] and r["in_corpus"]) else "misc"

        file_field: Optional[str] = None
        if documents_dir is not None and r["corpus_hash"]:
            # Per-paper artifact path. Matches the `file = {...}`
            # convention used by `--bib` import on process_corpus.py.
            hash_dir = documents_dir / r["corpus_hash"]
            file_field = str(hash_dir / "processed.pdf")

        fields: Dict[str, Optional[str]] = {
            "author": authors_str or None,
            "title": r["title"],
            "year": str(r["year"]) if r["year"] is not None else None,
            "journal": r["journal"],
            "doi": r["doi"],
            "file": file_field,
            # #51 — round-trip license metadata so operators can curate
            # in BibTeX. licenseurl is the flat BibTeX spelling; the
            # importer maps it back to the snake_case license_url column.
            "license": r["license"],
            "licenseurl": r["license_url"],
            # #54 — round-trip the deploy-time skip flag.
            "serve": ("false" if (r["serve"] is not None and not r["serve"]) else None),
            "servereason": r["serve_reason"],
            # corpus_hash is opaque but lets future bib_import match
            # entries back to per-paper artifacts even after a rename.
            "corpus_hash": r["corpus_hash"],
            "work_id": work_id,
        }
        out_chunks.append(render_entry(cite_keys[work_id], entry_type, fields))

    conn.close()
    return "\n\n".join(out_chunks) + ("\n" if out_chunks else "")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


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
        "-o", "--output", type=Path, default=None,
        help="Write to this file instead of stdout.",
    )
    parser.add_argument(
        "--db", type=Path, default=None,
        help="Override the SQLite path (default: <output_dir>/biblio_authority.sqlite).",
    )
    parser.add_argument(
        "--all-works", action="store_true",
        help="Include cited-but-not-in-corpus works as well as corpus papers.",
    )
    parser.add_argument(
        "--no-file-field", action="store_true",
        help="Omit `file = {...}` fields. Default: emit them pointing at "
             "<output_dir>/documents/<HASH>/processed.pdf so the round-trip "
             "matches the existing --bib convention.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    db_path = args.db or (args.output_dir / "biblio_authority.sqlite")
    if not db_path.exists():
        logger.error("biblio_authority.sqlite not found at %s. "
                     "Run `corpus run` (or `python -m bib.authority`) first.", db_path)
        return 1

    documents_dir = None if args.no_file_field else (args.output_dir / "documents")

    text = export_bibtex(
        db_path,
        corpus_only=not args.all_works,
        documents_dir=documents_dir,
    )

    if args.output is not None:
        args.output.write_text(text)
        logger.info("Wrote %d entry(ies) to %s",
                    text.count("@") if text else 0, args.output)
    else:
        sys.stdout.write(text)

    return 0


if __name__ == "__main__":
    sys.exit(main())
