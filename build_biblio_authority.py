#!/usr/bin/env python3
"""Build the bibliographic authority database from per-paper artifacts.

Reads metadata.json and references.json from the corpus output directory,
deduplicates cited references across papers, builds a citation graph,
and links Darwin Core taxa to their original-description works via
``scientificNameAuthorship`` parsing.

The result is a single SQLite database that serves as the corpus-level
source of truth for all bibliographic entities — whether or not the
corresponding PDF is physically in the corpus.

GUID priority: DOI > BHL Part/Item ID > normalized citation key.

Three phases:
  1. Seed from corpus papers (metadata.json)
  2. Ingest cited references + build citation graph (references.json)
  3. Link taxonomic-authority strings to works

After a build completes, ``reconcile_corpus_to_biblio.py`` can
be run to merge corpus papers whose Grobid-seeded work_id received no
incoming citations onto matching ghost cited-reference rows.

Idempotent and resumable: re-running merges new papers into the
existing database without duplicating existing records.

Usage:
    python build_biblio_authority.py /path/to/output
    python build_biblio_authority.py /path/to/output --enrich-bhl --bhl-api-key YOUR_KEY
    python build_biblio_authority.py /path/to/output --rebuild
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import signal
import sqlite3
import sys
import time
import unicodedata
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import requests
    from rapidfuzz import fuzz
    _HAS_BHL_DEPS = True
except ImportError:
    _HAS_BHL_DEPS = False

logger = logging.getLogger("build_biblio")

# Defaults are derived per-corpus from the output_dir positional arg in
# main(); see "corpuscle" layout in README.md.

# ── Normalization ────────────────────────────────────────────────────

def normalize_for_key(s: str) -> str:
    """Lowercase, strip diacritics, remove punctuation, collapse whitespace."""
    s = unicodedata.normalize("NFD", s)
    s = "".join(c for c in s if unicodedata.category(c) != "Mn")
    s = s.lower().strip()
    s = re.sub(r"[^\w\s]", "", s)
    s = re.sub(r"\s+", " ", s)
    return s


def normalize_doi(doi: str) -> str:
    """Normalize a DOI to its bare lowercase form."""
    doi = doi.strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "http://dx.doi.org/",
                   "https://dx.doi.org/", "doi:"):
        if doi.startswith(prefix):
            doi = doi[len(prefix):]
    return doi


def make_alias_key(surname: str, year: Optional[int], title: str) -> str:
    """Deterministic key for dedup matching."""
    norm_surname = normalize_for_key(surname)
    year_str = str(year) if year else "unknown"
    norm_title = normalize_for_key(title)[:40]
    return f"{norm_surname}|{year_str}|{norm_title}"


def make_corpus_guid(surname: str, year: Optional[int], title: str) -> str:
    """Generate a corpus: GUID from normalized fields."""
    return "corpus:" + make_alias_key(surname, year, title)


# ── Author parsing ───────────────────────────────────────────────────

# Name particles that are part of the surname, not the forename.
_PARTICLES = {"van", "von", "de", "di", "le", "la", "du", "del", "der", "den", "dos"}


def extract_surname_from_ref_author(author_str: str) -> str:
    """Extract surname from a reference author string like 'F Johnson'
    or 'O Norden Andersen' or 'van Riemsdijk'.

    Heuristic: the surname is the last token, unless preceded by a name
    particle (van, von, de, etc.), in which case particle + last token.
    """
    parts = author_str.strip().split()
    if not parts:
        return ""
    if len(parts) == 1:
        return parts[0]
    # Check if the second-to-last token is a particle
    surname_parts = [parts[-1]]
    i = len(parts) - 2
    while i >= 0 and parts[i].lower() in _PARTICLES:
        surname_parts.insert(0, parts[i])
        i -= 1
    return " ".join(surname_parts)


def extract_forename_from_ref_author(author_str: str) -> str:
    """Extract the forename/initial portion from a ref author string."""
    parts = author_str.strip().split()
    if len(parts) <= 1:
        return ""
    # Find where the surname starts (last token, possibly with particles)
    surname_start = len(parts) - 1
    while surname_start >= 1 and parts[surname_start - 1].lower() in _PARTICLES:
        surname_start -= 1
    return " ".join(parts[:surname_start])


# ── Authority string parsing ─────────────────────────────────────────

_AUTHORITY_RE = re.compile(
    r"^\(?"               # optional opening paren (genus transfer)
    r"(.+?)"              # author(s)
    r",\s*(\d{4})"        # comma + year
    r"\)?"                # optional closing paren
    r"$"
)


def parse_authority(authority: str) -> Optional[Tuple[List[str], int]]:
    """Parse a DwC scientificNameAuthorship string into (author_surnames, year).

    Handles: 'Eschscholtz, 1829', '(Huxley, 1859)', 'Quoy & Gaimard, 1833',
    'L. Agassiz, 1862', 'Lens & van Riemsdijk, 1908'.

    Returns None if unparseable.
    """
    if not authority:
        return None
    m = _AUTHORITY_RE.match(authority.strip())
    if not m:
        return None
    authors_str = m.group(1).strip()
    year = int(m.group(2))
    # Split on ' & ' or ' and '
    raw_authors = re.split(r"\s*&\s*|\s+and\s+", authors_str)
    surnames = []
    for a in raw_authors:
        a = a.strip()
        if not a:
            continue
        # Strip leading initials: "L. Agassiz" -> "Agassiz", "M. Sars" -> "Sars"
        # But keep "van Riemsdijk" as-is
        parts = a.split()
        # Drop parts that look like initials (single letter, optionally with period)
        name_parts = []
        for p in parts:
            if re.match(r"^[A-Z]\.?$", p):
                continue  # initial, skip
            name_parts.append(p)
        if name_parts:
            surnames.append(" ".join(name_parts))
        elif parts:
            # All parts were initials? Use the last one stripped of period
            surnames.append(parts[-1].rstrip("."))
    return surnames, year


# ── Schema ───────────────────────────────────────────────────────────

def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS works (
            work_id        TEXT PRIMARY KEY,
            guid_type      TEXT NOT NULL,
            title          TEXT,
            year           INTEGER,
            journal        TEXT,
            doi            TEXT,
            bhl_item_id    TEXT,
            bhl_part_id    TEXT,
            openalex_id    TEXT,
            corpus_hash    TEXT,
            in_corpus      INTEGER NOT NULL DEFAULT 0,
            source         TEXT NOT NULL,
            confidence     REAL DEFAULT 1.0,
            created_at     REAL NOT NULL,
            updated_at     REAL NOT NULL
        );

        CREATE TABLE IF NOT EXISTS work_authors (
            work_id            TEXT NOT NULL REFERENCES works(work_id),
            position           INTEGER NOT NULL,
            surname            TEXT NOT NULL,
            surname_normalized TEXT NOT NULL,
            forename           TEXT,
            PRIMARY KEY (work_id, position)
        );

        CREATE TABLE IF NOT EXISTS citations (
            citing_work_id     TEXT NOT NULL REFERENCES works(work_id),
            cited_work_id      TEXT NOT NULL REFERENCES works(work_id),
            citing_corpus_hash TEXT NOT NULL,
            grobid_xml_id      TEXT,
            raw_citation       TEXT,
            match_method       TEXT NOT NULL,
            match_score        REAL DEFAULT 1.0,
            PRIMARY KEY (citing_work_id, cited_work_id, citing_corpus_hash)
        );

        CREATE TABLE IF NOT EXISTS work_aliases (
            alias_key  TEXT NOT NULL,
            work_id    TEXT NOT NULL REFERENCES works(work_id),
            PRIMARY KEY (alias_key, work_id)
        );

        -- Links between original-description works and Darwin Core taxa.
        -- ``taxon_id`` is the DwC taxonID from the configured taxonomy
        -- snapshot (TEXT — DwC identifiers are strings).
        CREATE TABLE IF NOT EXISTS taxon_work_links (
            taxon_id   TEXT NOT NULL,
            work_id    TEXT NOT NULL REFERENCES works(work_id),
            link_type  TEXT NOT NULL,
            confidence REAL DEFAULT 1.0,
            PRIMARY KEY (taxon_id, work_id, link_type)
        );

        CREATE TABLE IF NOT EXISTS bhl_lookups (
            work_id      TEXT NOT NULL,
            query        TEXT NOT NULL,
            status       TEXT NOT NULL,
            error_msg    TEXT,
            attempted_at REAL NOT NULL,
            PRIMARY KEY (work_id)
        );

        CREATE TABLE IF NOT EXISTS build_meta (
            key   TEXT PRIMARY KEY,
            value TEXT
        );

        CREATE INDEX IF NOT EXISTS idx_works_doi ON works(doi) WHERE doi IS NOT NULL;
        CREATE INDEX IF NOT EXISTS idx_works_corpus_hash ON works(corpus_hash) WHERE corpus_hash IS NOT NULL;
        CREATE INDEX IF NOT EXISTS idx_works_year ON works(year);
        CREATE INDEX IF NOT EXISTS idx_works_in_corpus ON works(in_corpus);
        CREATE INDEX IF NOT EXISTS idx_work_authors_surname ON work_authors(surname_normalized);
        CREATE INDEX IF NOT EXISTS idx_citations_cited ON citations(cited_work_id);
        CREATE INDEX IF NOT EXISTS idx_citations_citing ON citations(citing_work_id);
        CREATE INDEX IF NOT EXISTS idx_taxon_work_links_work ON taxon_work_links(work_id);
    """)
    conn.commit()


# ── Work insertion helpers ───────────────────────────────────────────

def insert_work(conn: sqlite3.Connection, work_id: str, guid_type: str,
                title: str, year: Optional[int], journal: str, doi: str,
                corpus_hash: Optional[str], in_corpus: bool, source: str,
                confidence: float = 1.0) -> bool:
    """Insert a work. Returns True if inserted, False if already exists."""
    now = time.time()
    try:
        conn.execute(
            """INSERT INTO works (work_id, guid_type, title, year, journal, doi,
               corpus_hash, in_corpus, source, confidence, created_at, updated_at)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (work_id, guid_type, title, year, journal, doi or None,
             corpus_hash, int(in_corpus), source, confidence, now, now),
        )
        return True
    except sqlite3.IntegrityError:
        return False


def insert_authors(conn: sqlite3.Connection, work_id: str,
                   authors: List[Tuple[str, str]]) -> None:
    """Insert author list. authors is [(surname, forename), ...]."""
    for i, (surname, forename) in enumerate(authors):
        try:
            conn.execute(
                """INSERT OR IGNORE INTO work_authors
                   (work_id, position, surname, surname_normalized, forename)
                   VALUES (?, ?, ?, ?, ?)""",
                (work_id, i, surname, normalize_for_key(surname), forename),
            )
        except sqlite3.IntegrityError:
            pass


def insert_alias(conn: sqlite3.Connection, alias_key: str, work_id: str) -> None:
    try:
        conn.execute(
            "INSERT OR IGNORE INTO work_aliases (alias_key, work_id) VALUES (?, ?)",
            (alias_key, work_id),
        )
    except sqlite3.IntegrityError:
        pass


def insert_citation(conn: sqlite3.Connection, citing_work_id: str,
                    cited_work_id: str, citing_corpus_hash: str,
                    grobid_xml_id: str, raw_citation: str,
                    match_method: str, match_score: float) -> None:
    try:
        conn.execute(
            """INSERT OR IGNORE INTO citations
               (citing_work_id, cited_work_id, citing_corpus_hash,
                grobid_xml_id, raw_citation, match_method, match_score)
               VALUES (?, ?, ?, ?, ?, ?, ?)""",
            (citing_work_id, cited_work_id, citing_corpus_hash,
             grobid_xml_id, raw_citation, match_method, match_score),
        )
    except sqlite3.IntegrityError:
        pass


# ── Lookup helpers ───────────────────────────────────────────────────

def lookup_by_doi(conn: sqlite3.Connection, doi: str) -> Optional[str]:
    """Return work_id for a given normalized DOI, or None."""
    cur = conn.execute("SELECT work_id FROM works WHERE doi = ?", (doi,))
    row = cur.fetchone()
    return row[0] if row else None


def lookup_by_alias(conn: sqlite3.Connection, alias_key: str) -> Optional[str]:
    """Return work_id for a given alias key, or None."""
    cur = conn.execute("SELECT work_id FROM work_aliases WHERE alias_key = ?",
                       (alias_key,))
    row = cur.fetchone()
    return row[0] if row else None


def _first_author_candidates(conn: sqlite3.Connection, surname: str,
                             year: Optional[int]) -> List[Tuple[str, str]]:
    """Return [(work_id, title), ...] whose first author matches surname
    (and year, if given). First-author-only so multi-author ghosts don't
    get resurrected via a later-position surname collision.
    """
    norm_surname = normalize_for_key(surname)
    if year:
        cur = conn.execute(
            """SELECT DISTINCT wa.work_id, w.title
               FROM work_authors wa JOIN works w ON wa.work_id = w.work_id
               WHERE wa.surname_normalized = ? AND w.year = ? AND wa.position = 0""",
            (norm_surname, year),
        )
    else:
        cur = conn.execute(
            """SELECT DISTINCT wa.work_id, w.title
               FROM work_authors wa JOIN works w ON wa.work_id = w.work_id
               WHERE wa.surname_normalized = ? AND wa.position = 0""",
            (norm_surname,),
        )
    return [(r[0], r[1] or "") for r in cur.fetchall()]


def fuzzy_match(conn: sqlite3.Connection, surname: str, year: Optional[int],
                title: str) -> Optional[Tuple[str, float]]:
    """Find a matching work by author surname + year, ranked by title similarity.

    Returns (work_id, score) or None. Thin wrapper around
    ``fuzzy_match_with_score`` that applies the same confidence
    thresholds as the resolution cascade.
    """
    result = fuzzy_match_with_score(conn, surname, year, title)
    if result is None:
        return None
    work_id, set_score, ratio_score = result
    if set_score >= 85 and ratio_score >= 60:
        return work_id, set_score / 100.0
    return None


def fuzzy_match_with_score(conn: sqlite3.Connection, surname: str,
                           year: Optional[int],
                           title: str) -> Optional[Tuple[str, int, int]]:
    """Return the best (work_id, token_set_ratio, ratio) candidate, or None.

    Two scores are needed because ``token_set_ratio`` alone misrouted
    references in the live corpus (issue #2). It scores the Totton 1965
    *Lensia* paper title vs *A synopsis of the Siphonophora* at 82 —
    above the old 80 cutoff — because "a of the siphonophora" is a
    shared substring. Levenshtein ``ratio`` scores the same pair at 43,
    which is a truer signal of dissimilarity. The cascade uses both.

    Ranked by ``token_set_ratio`` (primary), breaking ties by ``ratio``
    to prefer candidates that also match by straight Levenshtein.
    """
    if not _HAS_BHL_DEPS:
        return None
    if not surname or not title:
        return None
    candidates = _first_author_candidates(conn, surname, year)
    if not candidates:
        return None
    norm_title = normalize_for_key(title)
    best = None
    for cand_id, cand_title in candidates:
        if not cand_title:
            continue
        norm_cand = normalize_for_key(cand_title)
        set_s = int(fuzz.token_set_ratio(norm_title, norm_cand))
        ratio_s = int(fuzz.ratio(norm_title, norm_cand))
        if best is None or (set_s, ratio_s) > (best[1], best[2]):
            best = (cand_id, set_s, ratio_s)
    return best


def author_year_match(conn: sqlite3.Connection, surname: str,
                      year: Optional[int]) -> Optional[str]:
    """Last-resort match: author surname + year only, if exactly one candidate."""
    if not surname or not year:
        return None
    candidates = _first_author_candidates(conn, surname, year)
    if len(candidates) == 1:
        return candidates[0][0]
    return None


# ── Phase 1: Seed from corpus papers ────────────────────────────────

def phase1_corpus_papers(conn: sqlite3.Connection, output_dir: Path) -> int:
    """Walk output/documents/*/metadata.json and seed works table."""
    docs_dir = output_dir / "documents"
    if not docs_dir.is_dir():
        logger.error("Documents directory not found: %s", docs_dir)
        return 0

    count = 0
    batch = 0
    for hash_dir in sorted(docs_dir.iterdir()):
        if not hash_dir.is_dir():
            continue
        meta_path = hash_dir / "metadata.json"
        if not meta_path.exists():
            continue

        corpus_hash = hash_dir.name
        # Skip if already seeded
        cur = conn.execute("SELECT 1 FROM works WHERE corpus_hash = ?", (corpus_hash,))
        if cur.fetchone():
            continue

        try:
            meta = json.loads(meta_path.read_text())
        except (json.JSONDecodeError, OSError) as e:
            logger.warning("Skipping %s: %s", meta_path, e)
            continue

        title = meta.get("title", "") or ""
        year = meta.get("year")
        journal = meta.get("journal", "") or ""
        doi_raw = meta.get("doi", "") or ""
        authors_raw = meta.get("authors", [])

        # Parse authors
        authors: List[Tuple[str, str]] = []
        first_surname = ""
        for a in authors_raw:
            surname = a.get("surname", "")
            forename = a.get("forename", "")
            if surname:
                authors.append((surname, forename))
                if not first_surname:
                    first_surname = surname

        # Determine work_id
        doi = normalize_doi(doi_raw) if doi_raw else ""
        if doi:
            work_id = doi
            guid_type = "doi"
        elif first_surname and title:
            work_id = make_corpus_guid(first_surname, year, title)
            guid_type = "corpus_key"
        elif first_surname:
            # No title — use filename as title stand-in
            filename = meta.get("filename", corpus_hash)
            work_id = make_corpus_guid(first_surname, year, filename)
            guid_type = "corpus_key"
        else:
            # No author, no DOI — use corpus hash directly
            work_id = f"corpus:{corpus_hash}"
            guid_type = "corpus_key"

        inserted = insert_work(conn, work_id, guid_type, title, year, journal,
                               doi, corpus_hash, in_corpus=True,
                               source="corpus_paper")
        if inserted:
            insert_authors(conn, work_id, authors)
            # Register alias for dedup matching
            if first_surname:
                alias = make_alias_key(first_surname, year, title or meta.get("filename", ""))
                insert_alias(conn, alias, work_id)
            count += 1
            batch += 1
            if batch >= 100:
                conn.commit()
                batch = 0
                logger.info("Phase 1 progress: %d corpus papers seeded", count)

    conn.commit()
    logger.info("Phase 1 complete: %d corpus papers seeded", count)
    return count


# ── Phase 2: Ingest references + citation graph ─────────────────────

def _resolve_reference(conn: sqlite3.Connection, ref: dict,
                       enrich_bhl: bool = False,
                       bhl_api_key: str = "",
                       bhl_max_year: Optional[int] = None) -> Tuple[str, str, float]:
    """Resolve a single reference dict to a work_id.

    Returns (work_id, match_method, match_score).
    Creates the work if no match found.
    """
    title = ref.get("title", "") or ""
    year = ref.get("year")
    journal = ref.get("journal", "") or ""
    doi_raw = ref.get("doi", "") or ""
    authors_raw = ref.get("authors", [])
    raw = ref.get("raw", "") or ""

    # Parse first author surname from ref format ("F Johnson")
    first_surname = ""
    if authors_raw:
        first_surname = extract_surname_from_ref_author(authors_raw[0])

    # ── Cascade step 1: DOI exact ──────────────────────────────────
    doi = normalize_doi(doi_raw) if doi_raw else ""
    if doi:
        existing = lookup_by_doi(conn, doi)
        if existing:
            return existing, "doi_exact", 1.0
        # DOI not yet in DB — create the work with DOI as ID
        work_id = doi
        insert_work(conn, work_id, "doi", title, year, journal, doi,
                     corpus_hash=None, in_corpus=False, source="cited_reference")
        authors = []
        for a in authors_raw:
            sn = extract_surname_from_ref_author(a)
            fn = extract_forename_from_ref_author(a)
            if sn:
                authors.append((sn, fn))
        insert_authors(conn, work_id, authors)
        if first_surname:
            alias = make_alias_key(first_surname, year, title)
            insert_alias(conn, alias, work_id)
        return work_id, "doi_exact", 1.0

    # ── Cascade step 2: Alias key exact ────────────────────────────
    if first_surname and (title or raw):
        alias = make_alias_key(first_surname, year, title or raw)
        existing = lookup_by_alias(conn, alias)
        if existing:
            return existing, "alias_exact", 1.0

    # ── Cascade step 3: BHL lookup (optional) ──────────────────────
    if enrich_bhl and first_surname and title:
        # Compute candidate work_id for resume checking
        candidate_id = make_corpus_guid(first_surname, year, title or raw)
        prev = conn.execute(
            "SELECT status FROM bhl_lookups WHERE work_id = ?", (candidate_id,)
        ).fetchone()
        if prev:
            prev_status = prev[0]
            if prev_status == "found":
                pass  # fall through — alias lookup below will catch it
            elif prev_status == "not_found":
                pass  # skip BHL, fall through to fuzzy match
            elif prev_status == "error":
                prev = None  # retry errors

        if not prev:
            bhl_status, bhl_match, bhl_error = _bhl_lookup(
                first_surname, year, title, api_key=bhl_api_key,
                max_year=bhl_max_year,
            )
            if bhl_status != "skipped":
                now = time.time()
                query_str = f"{first_surname} {year} {title[:60]}"
                conn.execute(
                    """INSERT OR REPLACE INTO bhl_lookups
                       (work_id, query, status, error_msg, attempted_at)
                       VALUES (?, ?, ?, ?, ?)""",
                    (candidate_id, query_str, bhl_status, bhl_error, now),
                )
            if bhl_status == "found":
                bhl_work_id, bhl_item_id, bhl_part_id = bhl_match
                existing = lookup_by_alias(conn, bhl_work_id)
                if existing:
                    return existing, "bhl_lookup", 0.9
                work_id = bhl_work_id
                inserted = insert_work(conn, work_id, "bhl", title, year, journal, doi,
                                       corpus_hash=None, in_corpus=False,
                                       source="cited_reference", confidence=0.9)
                if inserted:
                    authors = []
                    for a in authors_raw:
                        sn = extract_surname_from_ref_author(a)
                        fn = extract_forename_from_ref_author(a)
                        if sn:
                            authors.append((sn, fn))
                    insert_authors(conn, work_id, authors)
                    conn.execute(
                        "UPDATE works SET bhl_item_id = ?, bhl_part_id = ? WHERE work_id = ?",
                        (bhl_item_id, bhl_part_id, work_id),
                    )
                if first_surname:
                    alias = make_alias_key(first_surname, year, title)
                    insert_alias(conn, alias, work_id)
                return work_id, "bhl_lookup", 0.9

    # ── Cascade step 4 & 5: Fuzzy title / author+year fallback ─────
    # When the ref has a title, fuzzy-match it against the candidate
    # set using two metrics (see ``fuzzy_match_with_score`` for why).
    # Low scores mean the candidates are different works — do NOT
    # fall through to author_year_match and do NOT cache via alias.
    # (Single-metric fallthrough + unconditional alias caching is
    # what misrouted 15 corpus citations of Totton 1965's Lensia
    # paper onto the Synopsis row — see issue #2.)
    if first_surname and title:
        fuzzy_result = fuzzy_match_with_score(conn, first_surname, year, title)
        if fuzzy_result is not None:
            matched_id, set_score, ratio_score = fuzzy_result
            if set_score >= 85 and ratio_score >= 60:
                # Confident match. Cache so follow-up refs with the
                # same title resolve via alias_exact. Lower-scoring
                # candidates are treated as different works — better
                # a new ghost we can merge later via reconciliation
                # than a misroute that amplifies via the alias cache.
                alias = make_alias_key(first_surname, year, title)
                insert_alias(conn, alias, matched_id)
                return matched_id, "title_fuzzy", set_score / 100.0
            # Below threshold: fall through to new-work creation.

    # author_year fallback is safe only when we have NO title to
    # evaluate against the candidate. With a title we've already
    # decided above; without one we're truly blind and accept a
    # unique candidate. Never cache an alias on this path.
    if first_surname and year and not title:
        matched_id = author_year_match(conn, first_surname, year)
        if matched_id:
            return matched_id, "author_year_only", 0.6

    # ── No match — create new work ─────────────────────────────────
    if first_surname:
        work_id = make_corpus_guid(first_surname, year, title or raw)
    else:
        # No author at all — use a hash of the raw citation or title
        fallback = title or raw or f"unknown_{time.time_ns()}"
        work_id = f"corpus:unknown|{normalize_for_key(fallback)[:60]}"

    insert_work(conn, work_id, "corpus_key", title, year, journal, doi,
                corpus_hash=None, in_corpus=False, source="cited_reference")
    authors = []
    for a in authors_raw:
        sn = extract_surname_from_ref_author(a)
        fn = extract_forename_from_ref_author(a)
        if sn:
            authors.append((sn, fn))
    insert_authors(conn, work_id, authors)
    if first_surname:
        alias = make_alias_key(first_surname, year, title or raw)
        insert_alias(conn, alias, work_id)
    return work_id, "new", 1.0


# In-memory memoization of BHL search results for the current build run.
# Key is the query string; value is (results, error). Siphonophore lit
# is author-heavy (Haeckel 1888, Totton 1954, Lesson 1843) and the same
# (author, year) appears across many cited references with slightly
# different title formatting — without this, stage 2 does the same
# broad query dozens of times.
_BHL_QUERY_CACHE: Dict[str, Tuple[Optional[list], Optional[str]]] = {}


def _bhl_search(query: str, api_key: str) -> Tuple[Optional[list], Optional[str]]:
    """Call BHL PublicationSearch API. Returns (results_list, error_string)."""
    cached = _BHL_QUERY_CACHE.get(query)
    if cached is not None:
        return cached
    try:
        for attempt in range(3):
            if attempt > 0:
                time.sleep(2)  # back off on retry
            else:
                time.sleep(1)  # BHL API rate limit
            r = requests.get(
                "https://www.biodiversitylibrary.org/api3",
                params={
                    "op": "PublicationSearch",
                    "searchterm": query,
                    "searchtype": "C",
                    "page": 1,
                    "apikey": api_key,
                    "format": "json",
                },
                timeout=15,
            )
            if r.status_code == 500 and attempt < 2:
                time.sleep(2)
                continue
            r.raise_for_status()
            break
        result = (r.json().get("Result", []), None)
    except Exception as e:
        result = (None, str(e))
    _BHL_QUERY_CACHE[query] = result
    return result


def _bhl_match_results(results: list, title: str,
                       n_candidates: int = 10,
                       threshold: int = 75) -> Optional[Tuple[str, str, str]]:
    """Score BHL results against a reference title.

    Returns (work_id, bhl_item_id, bhl_part_id) or None.
    """
    if not results:
        return None
    norm_title = normalize_for_key(title)
    for item in results[:n_candidates]:
        bhl_title = item.get("Title", "")
        score = fuzz.token_set_ratio(norm_title, normalize_for_key(bhl_title))
        if score >= threshold:
            part_id = str(item.get("PartID", "")) if item.get("PartID") else None
            bhl_item_id_val = str(item.get("ItemID", "")) if item.get("ItemID") else None
            if part_id:
                work_id = f"bhl:part/{part_id}"
            elif bhl_item_id_val:
                work_id = f"bhl:item/{bhl_item_id_val}"
            else:
                continue
            return work_id, bhl_item_id_val or "", part_id or ""
    return None


def _bhl_lookup(surname: str, year: Optional[int],
                title: str,
                api_key: str = "",
                max_year: Optional[int] = None) -> Tuple[str, Optional[Tuple[str, str, str]], Optional[str]]:
    """Query BHL API for a matching publication.

    Two-stage search: first a narrow query (author + year + title prefix),
    then a broad fallback (author + year only) with more candidates checked.

    ``max_year`` skips refs newer than BHL's useful coverage window —
    modern papers without a parsed DOI almost never match BHL (its
    strength is pre-1960 natural history) and just burn rate-limited
    calls.

    Returns (status, match, error_msg):
      - ("found",     (work_id, bhl_item_id, bhl_part_id), None)
      - ("not_found", None, None)
      - ("error",     None, error_string)
      - ("skipped",   None, reason)  — missing prereqs (no year/key/deps)
    """
    if not _HAS_BHL_DEPS:
        return ("skipped", None, "rapidfuzz or requests not installed")

    if not year:
        return ("skipped", None, "no year")

    if max_year is not None and year > max_year:
        return ("skipped", None, f"year {year} > max_year {max_year}")

    if not api_key:
        return ("skipped", None, "no API key")

    # Strip chars that cause BHL API 500 errors
    clean_title = re.sub(r'[(){}\[\]"\']', '', title)

    # ── Stage 1: narrow query (author + year + title prefix) ──────
    narrow_query = f"{surname} {year} {clean_title[:60]}"
    results, error = _bhl_search(narrow_query, api_key)
    if error:
        logger.warning("BHL lookup failed for '%s': %s", narrow_query, error)
        return ("error", None, error)

    match = _bhl_match_results(results, title, n_candidates=10)
    if match:
        return ("found", match, None)

    # ── Stage 2: broad query (author + year only) ─────────────────
    broad_query = f"{surname} {year}"
    results, error = _bhl_search(broad_query, api_key)
    if error:
        logger.warning("BHL broad lookup failed for '%s': %s", broad_query, error)
        return ("error", None, error)

    match = _bhl_match_results(results, title, n_candidates=20)
    if match:
        return ("found", match, None)

    return ("not_found", None, None)


def phase2_references(conn: sqlite3.Connection, output_dir: Path,
                      enrich_bhl: bool = False,
                      bhl_api_key: str = "",
                      bhl_max_year: Optional[int] = None) -> Tuple[int, int]:
    """Walk references.json files and build citation graph.

    Returns (n_citations, n_new_works).
    """
    docs_dir = output_dir / "documents"
    n_citations = 0
    n_new_works = 0
    n_bhl_found = 0
    n_bhl_not_found = 0
    n_bhl_error = 0
    n_papers = 0
    batch = 0

    for hash_dir in sorted(docs_dir.iterdir()):
        if not hash_dir.is_dir():
            continue
        refs_path = hash_dir / "references.json"
        if not refs_path.exists():
            continue

        corpus_hash = hash_dir.name

        # Find the citing work_id for this corpus paper
        cur = conn.execute("SELECT work_id FROM works WHERE corpus_hash = ?",
                           (corpus_hash,))
        row = cur.fetchone()
        if not row:
            continue
        citing_work_id = row[0]

        try:
            refs_data = json.loads(refs_path.read_text())
        except (json.JSONDecodeError, OSError) as e:
            logger.warning("Skipping %s: %s", refs_path, e)
            continue

        refs = refs_data.get("references", [])
        for ref in refs:
            cited_work_id, match_method, match_score = _resolve_reference(
                conn, ref, enrich_bhl=enrich_bhl,
                bhl_api_key=bhl_api_key,
                bhl_max_year=bhl_max_year,
            )
            if match_method == "new":
                n_new_works += 1
                match_method = "new_work"

            insert_citation(
                conn, citing_work_id, cited_work_id, corpus_hash,
                grobid_xml_id=ref.get("xml_id", ""),
                raw_citation=ref.get("raw", ""),
                match_method=match_method,
                match_score=match_score,
            )
            n_citations += 1

        n_papers += 1
        batch += 1
        if batch >= 10:
            conn.commit()
            batch = 0
            logger.info("Phase 2 progress: %d papers, %d citations, %d new works",
                        n_papers, n_citations, n_new_works)

    conn.commit()
    logger.info("Phase 2 complete: %d papers, %d citations, %d new works created",
                n_papers, n_citations, n_new_works)
    return n_citations, n_new_works


# ── Phase 3: Link taxonomic-authority strings to works ─────────────

def phase3_authority_links(conn: sqlite3.Connection, taxonomy_path: Path) -> int:
    """Link taxa to original-description works via DwC scientificNameAuthorship.

    Reads ``taxon_id`` + ``scientific_name_authorship`` from the configured
    Darwin Core taxonomy snapshot, parses each authority string into
    (surnames, year), and matches against ``works`` by author+year. When
    no work matches, a stub is inserted so the link still resolves.
    """
    if not taxonomy_path.exists():
        logger.warning("Taxonomy database not found: %s", taxonomy_path)
        return 0

    tx_conn = sqlite3.connect(f"file:{taxonomy_path}?mode=ro", uri=True)
    tx_conn.row_factory = sqlite3.Row

    cur = tx_conn.execute(
        "SELECT taxon_id, scientific_name, scientific_name_authorship "
        "FROM taxa "
        "WHERE scientific_name_authorship IS NOT NULL "
        "  AND scientific_name_authorship != ''"
    )

    n_linked = 0
    n_stubs = 0
    batch = 0

    for row in cur:
        taxon_id = row["taxon_id"]
        authority = row["scientific_name_authorship"]

        # Skip if already linked
        check = conn.execute(
            "SELECT 1 FROM taxon_work_links WHERE taxon_id = ?", (taxon_id,)
        )
        if check.fetchone():
            continue

        parsed = parse_authority(authority)
        if not parsed:
            logger.debug("Could not parse authority '%s' for %s",
                         authority, row["scientific_name"])
            continue

        surnames, year = parsed
        if not surnames:
            continue

        first_surname = surnames[0]

        # Try to find a matching work
        # We don't have a title from the authority string, so use author+year
        matched_id = author_year_match(conn, first_surname, year)
        confidence = 0.7

        # If multiple authors in the authority, try matching with all of them
        if not matched_id and len(surnames) > 1:
            norm_surname = normalize_for_key(first_surname)
            candidate_cur = conn.execute(
                """SELECT DISTINCT wa.work_id
                   FROM work_authors wa JOIN works w ON wa.work_id = w.work_id
                   WHERE wa.surname_normalized = ? AND w.year = ? AND wa.position = 0""",
                (norm_surname, year),
            )
            candidates = candidate_cur.fetchall()
            if len(candidates) > 1:
                # Disambiguate by checking second author
                for cand_row in candidates:
                    cand_id = cand_row[0]
                    auth2_cur = conn.execute(
                        "SELECT surname_normalized FROM work_authors "
                        "WHERE work_id = ? AND position = 1",
                        (cand_id,),
                    )
                    auth2 = auth2_cur.fetchone()
                    if auth2 and len(surnames) > 1:
                        if auth2[0] == normalize_for_key(surnames[1]):
                            matched_id = cand_id
                            confidence = 0.85
                            break

        if matched_id:
            try:
                conn.execute(
                    """INSERT OR IGNORE INTO taxon_work_links
                       (taxon_id, work_id, link_type, confidence)
                       VALUES (?, ?, ?, ?)""",
                    (taxon_id, matched_id, "authority_match", confidence),
                )
                n_linked += 1
            except sqlite3.IntegrityError:
                pass
        else:
            # Create a stub work from the authority string
            work_id = make_corpus_guid(first_surname, year, "")
            inserted = insert_work(
                conn, work_id, "corpus_key", title="", year=year, journal="",
                doi="", corpus_hash=None, in_corpus=False,
                source="taxon_authority", confidence=0.5,
            )
            if inserted:
                authors = [(s, "") for s in surnames]
                insert_authors(conn, work_id, authors)
                n_stubs += 1
            try:
                conn.execute(
                    """INSERT OR IGNORE INTO taxon_work_links
                       (taxon_id, work_id, link_type, confidence)
                       VALUES (?, ?, ?, ?)""",
                    (taxon_id, work_id, "authority_match", 0.5),
                )
                n_linked += 1
            except sqlite3.IntegrityError:
                pass

        batch += 1
        if batch >= 100:
            conn.commit()
            batch = 0

    conn.commit()
    tx_conn.close()
    logger.info("Phase 3 complete: %d taxa linked, %d stub works created",
                n_linked, n_stubs)
    return n_linked


# ── Main ─────────────────────────────────────────────────────────────

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
        help="SQLite output path (default: <output_dir>/biblio_authority.sqlite)",
    )
    parser.add_argument(
        "--taxonomy-db", type=Path, default=None,
        help="Darwin Core taxonomy SQLite path "
             "(default: <output_dir>/taxonomy.sqlite)",
    )
    parser.add_argument(
        "--enrich-bhl", action="store_true",
        help="Query BHL API for references without DOI (requires network + API key)",
    )
    parser.add_argument(
        "--bhl-api-key", type=str,
        default=os.environ.get("BHL_API_KEY", ""),
        help="BHL API key (or set BHL_API_KEY env var). Free at biodiversitylibrary.org/account",
    )
    parser.add_argument(
        "--bhl-max-year", type=int, default=1960,
        help="Skip BHL lookups for refs newer than this (BHL coverage is mostly "
             "pre-1960; modern refs without a DOI rarely match). Default: 1960.",
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

    if args.output is None:
        args.output = args.output_dir / "biblio_authority.sqlite"
    if args.taxonomy_db is None:
        args.taxonomy_db = args.output_dir / "taxonomy.sqlite"

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
            # Preserve bhl_lookups across rebuilds — it's the cache of
            # rate-limited BHL API calls. Dropping it forces every
            # lookup to be re-paid the next time ``--enrich-bhl`` is
            # used, which takes hours and is unrelated to the schema
            # churn --rebuild is meant to address.
            logger.info("Rebuilding: dropping all tables except bhl_lookups")
            conn.executescript("""
                DROP TABLE IF EXISTS taxon_work_links;
                DROP TABLE IF EXISTS citations;
                DROP TABLE IF EXISTS work_aliases;
                DROP TABLE IF EXISTS work_authors;
                DROP TABLE IF EXISTS works;
                DROP TABLE IF EXISTS build_meta;
            """)

        create_schema(conn)

        # Record build metadata
        conn.execute(
            "INSERT OR REPLACE INTO build_meta (key, value) VALUES (?, ?)",
            ("last_build_ts", str(time.time())),
        )
        conn.execute(
            "INSERT OR REPLACE INTO build_meta (key, value) VALUES (?, ?)",
            ("output_dir", str(args.output_dir)),
        )
        conn.commit()

        # Phase 1
        logger.info("═══ Phase 1: Seeding from corpus papers ═══")
        n_papers = phase1_corpus_papers(conn, args.output_dir)

        # Phase 2
        logger.info("═══ Phase 2: Ingesting references + citation graph ═══")
        if args.enrich_bhl:
            if not _HAS_BHL_DEPS:
                logger.warning(
                    "--enrich-bhl set but rapidfuzz/requests are not installed. "
                    "Run `pip install rapidfuzz requests`. BHL lookups will be "
                    "silently skipped."
                )
            elif not args.bhl_api_key:
                logger.warning(
                    "--enrich-bhl set but no API key provided. "
                    "Set BHL_API_KEY env var or pass --bhl-api-key. "
                    "Get a free key at https://www.biodiversitylibrary.org/account/. "
                    "BHL lookups will be skipped."
                )
            else:
                logger.info(
                    "BHL enrichment enabled (max_year=%d, key set)",
                    args.bhl_max_year,
                )
        n_citations, n_new = phase2_references(
            conn, args.output_dir, enrich_bhl=args.enrich_bhl,
            bhl_api_key=args.bhl_api_key,
            bhl_max_year=args.bhl_max_year,
        )

        # Phase 3
        logger.info("═══ Phase 3: Linking taxonomic authorities ═══")
        n_linked = phase3_authority_links(conn, args.taxonomy_db)

        # Summary
        stats = {}
        for table in ("works", "work_authors", "citations", "work_aliases", "taxon_work_links"):
            cur = conn.execute(f"SELECT COUNT(*) FROM {table}")  # noqa: S608
            stats[table] = cur.fetchone()[0]

        cur = conn.execute("SELECT COUNT(*) FROM works WHERE in_corpus = 1")
        n_in_corpus = cur.fetchone()[0]
        cur = conn.execute("SELECT COUNT(*) FROM works WHERE in_corpus = 0")
        n_not_in_corpus = cur.fetchone()[0]
        cur = conn.execute(
            "SELECT guid_type, COUNT(*) FROM works GROUP BY guid_type ORDER BY COUNT(*) DESC"
        )
        guid_dist = list(cur)
        cur = conn.execute(
            "SELECT match_method, COUNT(*) FROM citations GROUP BY match_method ORDER BY COUNT(*) DESC"
        )
        match_dist = list(cur)

        logger.info("═══ Build complete ═══")
        logger.info("  Works:        %d (%d in corpus, %d cited-only)",
                     stats["works"], n_in_corpus, n_not_in_corpus)
        logger.info("  Authors:      %d", stats["work_authors"])
        logger.info("  Citations:    %d", stats["citations"])
        logger.info("  Aliases:      %d", stats["work_aliases"])
        logger.info("  Taxon links:  %d", stats["taxon_work_links"])
        logger.info("  GUID types:   %s", guid_dist)
        logger.info("  Match methods: %s", match_dist)
        cur = conn.execute(
            "SELECT status, COUNT(*) FROM bhl_lookups GROUP BY status ORDER BY COUNT(*) DESC"
        )
        bhl_dist = list(cur)
        if bhl_dist:
            logger.info("  BHL lookups:  %s", bhl_dist)

        conn.execute(
            "INSERT OR REPLACE INTO build_meta (key, value) VALUES (?, ?)",
            ("build_complete_ts", str(time.time())),
        )
        conn.execute(
            "INSERT OR REPLACE INTO build_meta (key, value) VALUES (?, ?)",
            ("total_works", str(stats["works"])),
        )
        conn.execute(
            "INSERT OR REPLACE INTO build_meta (key, value) VALUES (?, ?)",
            ("total_citations", str(stats["citations"])),
        )
        conn.commit()
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
