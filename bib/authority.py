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

After a build completes, ``bib.reconcile`` can merge corpus papers whose
Grobid-seeded work_id received no incoming citations onto matching cited-work
rows. The raw reference observations remain independent of those canonical
work choices.

Idempotency (#30, #240): re-running unchanged inputs leaves the evidence and
mapping graph untouched (the CLI still refreshes build-run metadata).
Reference evidence is append-only and content-addressed; when its current set
changes, every current observation-to-work mapping and the legacy ``citations``
materialization are deterministically rebuilt. ``--rebuild`` drops every
table except the rate-limited ``bhl_lookups`` cache.

Usage:
    python build_biblio_authority.py /path/to/output
    python build_biblio_authority.py /path/to/output --enrich-bhl --bhl-api-key YOUR_KEY
    python build_biblio_authority.py /path/to/output --rebuild
"""

from __future__ import annotations

import argparse
import hashlib
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
from urllib.parse import unquote

try:
    import requests
    from rapidfuzz import fuzz
    _HAS_BHL_DEPS = True
except ImportError:
    _HAS_BHL_DEPS = False

from pipeline.external import (
    CircuitBreaker,
    CircuitOpenError,
    is_transient,
    retry_with_backoff,
)

logger = logging.getLogger("corpus.biblio")

# Changes only when the deterministic observation -> work rules change. It is
# persisted beside every verdict so an operator can explain why a mapping was
# reconsidered independently of the package release number (#240).
REFERENCE_MAPPING_PRODUCER = "reference-mapping-v2"

# Cross-block escape hatch measured by the #155 audit. Short/generic titles
# are excluded; the threshold is public so the read-only QC tool uses the same
# title-length boundary when it reports the broader review population.
IDENTITY_TITLE_MIN_ALPHA = 25

# Module-level breaker — BHL is hit from a single CLI run, so per-process
# state is the right scope. Threshold is generous: BHL routinely 500s on
# papers with funky title characters, and we don't want one bad paper to
# shut the whole BHL pass down.
_BHL_BREAKER = CircuitBreaker("bhl", threshold=20, cooldown_s=120.0)

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
    """Normalize URL/prefix spelling to a bare lowercase DOI (#155)."""
    doi = unquote(doi.strip()).lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "http://dx.doi.org/",
                   "https://dx.doi.org/", "info:doi/", "doi:"):
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
            -- #51 — figure licensing + publishable gate.
            -- license: SPDX short identifier or one of the small custom
            --          vocabulary (public-domain, all-rights-reserved,
            --          publisher-permission, unknown).
            -- license_source: 'bibtex' | 'age_based_pd' | 'unknown'.
            -- publishable: derived at build time from license + age cutoff
            --              (default `licensing.pd_cutoff_years: 95`).
            license        TEXT,
            license_url    TEXT,
            license_source TEXT,
            publishable    INTEGER,
            -- #54 — PDF QC skip flag. serve=0 means package_for_serve
            -- excludes the paper from the served bundle; serve_reason
            -- is a short tag from the closed vocab in dev_docs/QC.md.
            serve          INTEGER NOT NULL DEFAULT 1,
            serve_reason   TEXT,
            -- #79 — set by bib/importer.py whenever a row is touched by
            -- a `corpus bib import`. NULL means the row's bib fields
            -- have never been blessed by a human-edited .bib; the
            -- format_citation MCP tool uses this to gate "from user .bib"
            -- (no warning) vs "grobid-reconciled" (warning) provenance.
            bib_imported_at REAL,
            -- #176 — operator override for OCR language packs, as a
            -- `+`-joined Tesseract pack list (e.g. `pol+eng`). Not a
            -- bibliographic fact and deliberately not the standard BibTeX
            -- `language` field, which means "language of the work" and
            -- which Zotero populates by default — an imported .bib must
            -- not silently start steering OCR. NULL for almost every row.
            ocrlang        TEXT,
            -- #186 — operator override for whether/how OCR runs on this
            -- document: force, redo, or skip-text. Like ocrlang, this is an
            -- instruction rather than bibliographic metadata.
            ocrmode        TEXT,
            -- #214 — curation fields. `doclang` is a BCP-47 tag recording
            -- what the paper *is* ("de-Latf": German set in Fraktur) where
            -- `ocrlang` above records what to *do* about it; `pagemap` is
            -- free text describing the scan's physical structure, so a
            -- `keeppages` range is auditable by eye. Both are round-tripped
            -- and read by nothing — no pipeline stage and no fingerprint.
            doclang        TEXT,
            pagemap        TEXT,
            -- #188 — physical page selection, e.g. `3--20`. Unlike the two
            -- above this one changes the document: every later stage sees
            -- only the selected pages, so it is fingerprinted.
            keeppages      TEXT,
            created_at     REAL NOT NULL,
            updated_at     REAL NOT NULL
        );

        -- v0.3 schema migration helpers (#51 + #54). The CREATE TABLE IF
        -- NOT EXISTS above is a no-op when the table already exists, so
        -- existing biblio_authority.sqlite files miss the new columns.
        -- The Python wrapper below adds them with ALTER TABLE.

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

        -- Immutable/re-derivable evidence from references.json (#240).
        -- Observations and memberships are append-only. Current-set selection
        -- lives separately so replacing/removing an artifact never mutates or
        -- deletes its raw or parsed evidence.
        CREATE TABLE IF NOT EXISTS reference_observations (
            observation_id     TEXT PRIMARY KEY,
            citing_corpus_hash TEXT NOT NULL,
            ordinal            INTEGER NOT NULL,
            grobid_xml_id      TEXT,
            raw_citation       TEXT,
            title              TEXT,
            year               INTEGER,
            journal            TEXT,
            doi                TEXT,
            authors_json       TEXT NOT NULL,
            first_seen_at      REAL NOT NULL
        );

        CREATE TABLE IF NOT EXISTS reference_observation_sets (
            corpus_hash        TEXT NOT NULL,
            source_fingerprint TEXT NOT NULL,
            observation_count INTEGER NOT NULL,
            first_seen_at      REAL NOT NULL,
            PRIMARY KEY (corpus_hash, source_fingerprint)
        );

        CREATE TABLE IF NOT EXISTS reference_observation_memberships (
            corpus_hash        TEXT NOT NULL,
            source_fingerprint TEXT NOT NULL,
            ordinal            INTEGER NOT NULL,
            observation_id     TEXT NOT NULL
                REFERENCES reference_observations(observation_id),
            PRIMARY KEY (corpus_hash, source_fingerprint, ordinal),
            FOREIGN KEY (corpus_hash, source_fingerprint)
                REFERENCES reference_observation_sets(
                    corpus_hash, source_fingerprint
                )
        );

        CREATE TABLE IF NOT EXISTS reference_current_sets (
            corpus_hash        TEXT PRIMARY KEY,
            source_fingerprint TEXT NOT NULL,
            selected_at        REAL NOT NULL,
            FOREIGN KEY (corpus_hash, source_fingerprint)
                REFERENCES reference_observation_sets(
                    corpus_hash, source_fingerprint
                )
        );

        CREATE TABLE IF NOT EXISTS observation_work (
            observation_id TEXT PRIMARY KEY
                REFERENCES reference_observations(observation_id),
            work_id         TEXT NOT NULL REFERENCES works(work_id),
            match_method    TEXT NOT NULL,
            match_score     REAL DEFAULT 1.0,
            producer_version TEXT NOT NULL,
            mapped_at       REAL NOT NULL
        );

        -- Corpus-paper -> canonical-work decisions made by bib.reconcile.
        -- This makes the legacy compatibility merge reviewable even though
        -- the superseded works row is removed.
        CREATE TABLE IF NOT EXISTS work_reconciliation_decisions (
            corpus_hash    TEXT NOT NULL,
            source_work_id TEXT NOT NULL,
            target_work_id TEXT NOT NULL,
            match_method   TEXT NOT NULL,
            match_score    REAL,
            producer_version TEXT NOT NULL,
            decided_at     REAL NOT NULL,
            PRIMARY KEY (
                corpus_hash, source_work_id, target_work_id, producer_version
            )
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

        -- Tracks the mtime of the per-paper artifacts we ingested
        -- (metadata.json for phase 1, references.json for phase 2) so
        -- subsequent runs can detect when a single paper was reprocessed
        -- and reconcile the cross-paper rows without a full --rebuild.
        -- artifact ∈ ('metadata', 'references').
        CREATE TABLE IF NOT EXISTS paper_artifacts_processed (
            corpus_hash    TEXT NOT NULL,
            artifact       TEXT NOT NULL,
            source_mtime   REAL NOT NULL,
            processed_at   REAL NOT NULL,
            PRIMARY KEY (corpus_hash, artifact)
        );

        CREATE INDEX IF NOT EXISTS idx_works_doi ON works(doi) WHERE doi IS NOT NULL;
        CREATE INDEX IF NOT EXISTS idx_works_corpus_hash ON works(corpus_hash) WHERE corpus_hash IS NOT NULL;
        CREATE INDEX IF NOT EXISTS idx_works_year ON works(year);
        CREATE INDEX IF NOT EXISTS idx_works_in_corpus ON works(in_corpus);
        CREATE INDEX IF NOT EXISTS idx_works_publishable ON works(publishable);
        CREATE INDEX IF NOT EXISTS idx_works_serve ON works(serve);
        CREATE INDEX IF NOT EXISTS idx_work_authors_surname ON work_authors(surname_normalized);
        CREATE INDEX IF NOT EXISTS idx_citations_cited ON citations(cited_work_id);
        CREATE INDEX IF NOT EXISTS idx_citations_citing ON citations(citing_work_id);
        CREATE INDEX IF NOT EXISTS idx_reference_observations_citing
            ON reference_observations(citing_corpus_hash);
        CREATE INDEX IF NOT EXISTS idx_reference_memberships_observation
            ON reference_observation_memberships(observation_id);
        CREATE INDEX IF NOT EXISTS idx_observation_work_work
            ON observation_work(work_id);
        CREATE INDEX IF NOT EXISTS idx_taxon_work_links_work ON taxon_work_links(work_id);
    """)
    _migrate_works_columns(conn)
    conn.commit()


# Column additions per release. ALTER TABLE ADD COLUMN is a no-op on
# new DBs (CREATE TABLE above already declared the columns) and a
# one-time migration on DBs built by older releases.
_V03_WORKS_COLUMNS = [
    ("license",        "TEXT"),
    ("license_url",    "TEXT"),
    ("license_source", "TEXT"),
    ("publishable",    "INTEGER"),
    ("serve",          "INTEGER NOT NULL DEFAULT 1"),
    ("serve_reason",   "TEXT"),
]
_V05_WORKS_COLUMNS = [
    ("bib_imported_at", "REAL"),  # #79
]
_V11_WORKS_COLUMNS = [
    ("ocrlang", "TEXT"),  # #176
]
_V12_WORKS_COLUMNS = [
    ("doclang", "TEXT"),    # #214
    ("pagemap", "TEXT"),    # #214
    ("keeppages", "TEXT"),  # #188
]
_V13_WORKS_COLUMNS = [
    ("ocrmode", "TEXT"),  # #186
]


def _migrate_works_columns(conn: sqlite3.Connection) -> None:
    """Idempotent ALTER TABLE for works.* additions from past releases."""
    have = {row[1] for row in conn.execute("PRAGMA table_info(works)")}
    for name, decl in (*_V03_WORKS_COLUMNS, *_V05_WORKS_COLUMNS,
                       *_V11_WORKS_COLUMNS, *_V12_WORKS_COLUMNS,
                       *_V13_WORKS_COLUMNS):
        if name not in have:
            conn.execute(f"ALTER TABLE works ADD COLUMN {name} {decl}")


# #51 — license vocab. Per dev_docs/LICENSING.md, the value is an SPDX
# short identifier extended with a small custom vocabulary. The
# `_PUBLISHABLE_LICENSES` set captures the licenses we treat as
# reusable in derived publications by default; everything else gets
# publishable=0 unless the operator overrides.
_PUBLISHABLE_LICENSES = frozenset({
    # SPDX (and our own "public-domain" tag)
    "public-domain", "cc0-1.0",
    "cc-by-1.0", "cc-by-2.0", "cc-by-2.5", "cc-by-3.0", "cc-by-4.0",
    "cc-by-sa-3.0", "cc-by-sa-4.0",
})
_NON_PUBLISHABLE_LICENSES = frozenset({
    "all-rights-reserved", "publisher-permission", "unknown",
})


def _seed_license_and_serve(conn: sqlite3.Connection, work_id: str, meta: dict) -> None:
    """Copy bib-derived license + serve + OCR + curation fields from
    metadata.json into works.*.

    ``ocrlang`` (#176) and ``ocrmode`` (#186) ride along here purely so
    ``corpus bib export`` can round-trip them. Nothing in the pipeline reads
    them back from this table — the scan stage reads them straight off the
    BibIndex, because it has to run before the authority DB exists.

    ``doclang`` and ``pagemap`` (#214) ride along for the same reason and
    go one step further: nothing reads them *at all*. They are a curator's
    record of what a document is and how the scan is put together, and
    keeping them inert is what makes them safe to correct — fixing a typo
    in a note must not reprocess a document.
    """
    license_v = meta.get("license")
    license_url = meta.get("license_url")
    serve_v = meta.get("serve")
    serve_reason = meta.get("serve_reason")
    ocrlang = meta.get("ocrlang")
    ocrmode = meta.get("ocrmode")
    doclang = meta.get("doclang")
    pagemap = meta.get("pagemap")
    keeppages = meta.get("keeppages")

    # Build the SET clause only for fields we have a value for.
    sets, params = [], []
    if license_v:
        sets.append("license = ?")
        params.append(license_v)
        sets.append("license_source = ?")
        params.append("bibtex")
    if license_url:
        sets.append("license_url = ?")
        params.append(license_url)
    if serve_v is not None:
        sets.append("serve = ?")
        params.append(int(serve_v))
    if serve_reason:
        sets.append("serve_reason = ?")
        params.append(serve_reason)
    if ocrlang:
        sets.append("ocrlang = ?")
        params.append(ocrlang)
    if ocrmode:
        sets.append("ocrmode = ?")
        params.append(ocrmode)
    if doclang:
        sets.append("doclang = ?")
        params.append(doclang)
    if pagemap:
        sets.append("pagemap = ?")
        params.append(pagemap)
    if keeppages:
        sets.append("keeppages = ?")
        params.append(keeppages)
    if sets:
        params.append(work_id)
        conn.execute(
            f"UPDATE works SET {', '.join(sets)} WHERE work_id = ?",
            params,
        )


def derive_publishable(license_v: Optional[str], year: Optional[int],
                       pd_cutoff_years: int = 95) -> tuple[Optional[int], str]:
    """Decide ``publishable`` (0/1/None) and ``license_source`` for a work.

    Logic (per #51):
      1. Explicit license string wins. SPDX/CC-BY family + ``public-domain``
         → publishable. ``all-rights-reserved`` / ``publisher-permission`` /
         ``unknown`` → not publishable.
      2. No license + a year + year is older than the configured PD cutoff
         → ``age_based_pd``; publishable.
      3. Otherwise: source ``unknown``; not publishable (conservative
         default per the issue body).

    Returns (publishable, license_source). publishable=None means
    "couldn't decide" (e.g. unrecognized license string).

    **This boolean is lossy and is not what should be shown to clients**
    (#154 §2): a ``0`` means "we could not establish public domain", which
    is *not* the same as "the rightsholder forbade it", yet both collapse
    to the same value. In the served siphonophore corpus 86% of works were
    ``publishable=0, license_source=unknown`` and **not one** was asserted
    ``all-rights-reserved``. Use :func:`clearance_state` for anything a
    client or a human reads; keep this boolean for the storage column and
    the strict-profile gate, where "could not establish" and "forbidden"
    do warrant the same conservative answer.
    """
    import datetime as _dt
    if license_v:
        norm = license_v.strip().lower()
        if norm in _PUBLISHABLE_LICENSES:
            return 1, "bibtex"
        if norm in _NON_PUBLISHABLE_LICENSES:
            return 0, "bibtex"
        # Unrecognized license string — leave publishable null but record
        # source. Warn (#154 §2): silently NULLing means a typo'd
        # `license = {CC-BY 4.0}` (space, not hyphen) blocks a figure under
        # strict profiles exactly like an asserted all-rights-reserved,
        # with nothing in the build log to explain it.
        logger.warning(
            "unrecognized license string %r — leaving publishable NULL. "
            "Under a strict output profile this blocks the work's figures "
            "as firmly as an explicit all-rights-reserved. Recognized "
            "values: %s",
            license_v,
            ", ".join(sorted(_PUBLISHABLE_LICENSES | _NON_PUBLISHABLE_LICENSES)),
        )
        return None, "bibtex"
    if year is not None:
        current = _dt.datetime.now().year
        if (current - int(year)) >= pd_cutoff_years:
            return 1, "age_based_pd"
    return 0, "unknown"


# The five states a work's publication clearance can actually be in
# (#154 §2). ``publishable`` collapses the last three into a single 0,
# which is why the flag reads as a prohibition when it usually means
# "unknown".
CLEARANCE_PUBLIC_DOMAIN = "public_domain"
CLEARANCE_LICENSED_OPEN = "licensed_open"
CLEARANCE_RESTRICTED = "restricted"
CLEARANCE_UNDETERMINED = "undetermined"
CLEARANCE_NO_RECORD = "no_record"


def clearance_state(work: Optional[Dict]) -> str:
    """Publication-clearance state for a works row (#154 §2).

    Distinguishes the states ``publishable`` flattens:

    * ``public_domain``  — asserted public-domain, or past the PD cutoff
      (``license_source = age_based_pd``).
    * ``licensed_open``  — an explicit open license (CC-BY family).
    * ``restricted``     — the rightsholder's terms were recorded and
      forbid republication (``all-rights-reserved``, etc.). **Positive
      evidence of restriction**, unlike the states below.
    * ``undetermined``   — a license string was recorded but not
      recognized (``publishable IS NULL``). Often a typo.
    * ``no_record``      — nothing on file: no license and not old enough
      to be age-based PD, or the work isn't in the authority DB at all.
      This is the overwhelmingly common case and carries **no** evidence
      of restriction either way.

    Derived from the stored columns rather than a new one, so no schema
    migration or authority-DB rebuild is required.
    """
    if not work:
        return CLEARANCE_NO_RECORD
    source = (work.get("license_source") or "").strip().lower()
    pub = work.get("publishable")
    license_v = (work.get("license") or "").strip().lower()

    if pub is None:
        # Recorded but unrecognized — the typo case.
        return CLEARANCE_UNDETERMINED if source else CLEARANCE_NO_RECORD
    if int(pub) == 1:
        if source == "age_based_pd":
            return CLEARANCE_PUBLIC_DOMAIN
        if license_v in {"public-domain", "cc0", "cc0-1.0"}:
            return CLEARANCE_PUBLIC_DOMAIN
        return CLEARANCE_LICENSED_OPEN
    # publishable == 0: restricted only when we actually recorded terms.
    if source == "bibtex" and license_v:
        return CLEARANCE_RESTRICTED
    return CLEARANCE_NO_RECORD


def _resolve_pd_cutoff_from_config(config_path: Optional[Path]) -> int:
    """Read ``licensing.pd_cutoff_years`` from the per-corpuscle config
    or fall back to the schema default (95).
    """
    default = 95
    if config_path is None:
        return default
    if not config_path.exists():
        return default
    try:
        import yaml
        raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        return int((raw.get("licensing") or {}).get("pd_cutoff_years", default))
    except Exception:
        return default


def apply_publishable_derivation(conn: sqlite3.Connection,
                                 pd_cutoff_years: int = 95) -> int:
    """Walk every work, set ``publishable`` + ``license_source`` based on
    ``license`` + ``year``. Returns count of updated rows.

    Idempotent: run after every build (or when ``licensing.pd_cutoff_years``
    in config.yaml changes).
    """
    rows = list(conn.execute(
        "SELECT work_id, license, year FROM works"
    ))
    n = 0
    for r in rows:
        publishable, source = derive_publishable(r[1], r[2], pd_cutoff_years)
        conn.execute(
            "UPDATE works SET publishable = ?, license_source = ? WHERE work_id = ?",
            (publishable, source, r[0]),
        )
        n += 1
    conn.commit()
    return n


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


def _artifact_state(
    conn: sqlite3.Connection, corpus_hash: str, artifact: str,
    current_mtime: float,
) -> Tuple[bool, bool]:
    """Return ``(seen, stale)`` for a per-paper artifact.

    ``seen`` is True when the row exists; ``stale`` is True when the
    on-disk file mtime is newer than the recorded mtime (or no mtime
    is recorded). Callers use ``not seen`` to take the first-ingest
    path and ``seen and stale`` to take the re-ingest path.
    """
    cur = conn.execute(
        "SELECT source_mtime FROM paper_artifacts_processed "
        "WHERE corpus_hash = ? AND artifact = ?",
        (corpus_hash, artifact),
    )
    row = cur.fetchone()
    if row is None:
        return False, True
    stored = row[0]
    if stored is None or current_mtime > stored:
        return True, True
    return True, False


def _record_artifact(
    conn: sqlite3.Connection, corpus_hash: str, artifact: str,
    current_mtime: float,
) -> None:
    """Stamp the artifact as processed at ``current_mtime``."""
    conn.execute(
        """INSERT INTO paper_artifacts_processed
           (corpus_hash, artifact, source_mtime, processed_at)
           VALUES (?, ?, ?, ?)
           ON CONFLICT(corpus_hash, artifact) DO UPDATE SET
               source_mtime = excluded.source_mtime,
               processed_at = excluded.processed_at""",
        (corpus_hash, artifact, current_mtime, time.time()),
    )


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
    cur = conn.execute(
        """SELECT work_id FROM works WHERE doi = ?
           ORDER BY in_corpus DESC,
                    (bib_imported_at IS NOT NULL) DESC,
                    (guid_type = 'doi') DESC,
                    work_id
           LIMIT 1""",
        (doi,),
    )
    row = cur.fetchone()
    return row[0] if row else None


def _doi_corruption_shape(left: str, right: str) -> Optional[str]:
    """Name a narrow OCR corruption relating two non-identical DOIs.

    This does not repair either DOI. It is only a candidate generator for a
    second, independent title check (#239). Hyphens may be either spurious or
    genuinely lost, so equality after removing them is symmetric. Glued text
    is a long alphabetic suffix on one otherwise-complete DOI.
    """
    left = normalize_doi(left or "")
    right = normalize_doi(right or "")
    if not left or not right or left == right:
        return None
    if left.replace("-", "") == right.replace("-", ""):
        return "doi_hyphenation"
    shorter, longer = sorted((left, right), key=len)
    if (
        longer.startswith(shorter)
        and shorter.count("(") > shorter.count(")")
        and longer[len(shorter):].startswith(")")
    ):
        return "doi_truncation"
    suffix = longer[len(shorter):] if longer.startswith(shorter) else ""
    if len(suffix) >= 4 and suffix.isalpha():
        return "doi_trailing_text"
    return None


def lookup_doi_variant_by_title(
    conn: sqlite3.Connection,
    doi: str,
    title: str,
) -> Optional[Tuple[str, str, float]]:
    """Resolve a corrupted DOI only when the title independently agrees.

    Returns ``(work_id, match_method, score)``. The established dual title
    threshold is reused: token-set ratio >= 85 *and* straight ratio >= 60.
    Candidate choice is deterministic and favors an in-corpus or more-cited
    canonical row before a lexical work-id tie-break. A DOI-shape resemblance
    alone can never merge works.
    """
    if not _HAS_BHL_DEPS or not doi or not title:
        return None
    normalized_title = normalize_for_key(title)
    if not normalized_title:
        return None
    rows = conn.execute(
        """SELECT w.work_id, w.doi, w.title, w.in_corpus,
                  COUNT(c.citing_work_id) AS cited_count
           FROM works w
           LEFT JOIN citations c ON c.cited_work_id = w.work_id
           WHERE w.doi IS NOT NULL AND w.title IS NOT NULL
           GROUP BY w.work_id, w.doi, w.title, w.in_corpus"""
    ).fetchall()
    candidates = []
    for work_id, candidate_doi, candidate_title, in_corpus, cited_count in rows:
        corruption = _doi_corruption_shape(doi, candidate_doi)
        if corruption is None:
            continue
        normalized_candidate = normalize_for_key(candidate_title or "")
        if not normalized_candidate:
            continue
        set_score = int(fuzz.token_set_ratio(
            normalized_title, normalized_candidate,
        ))
        ratio_score = int(fuzz.ratio(normalized_title, normalized_candidate))
        if set_score < 85 or ratio_score < 60:
            continue
        candidates.append((
            int(bool(in_corpus)), int(cited_count or 0), set_score, ratio_score,
            str(work_id), corruption,
        ))
    if not candidates:
        return None
    candidates.sort(key=lambda row: (-row[0], -row[1], -row[2], -row[3], row[4]))
    _in_corpus, _cited_count, set_score, _ratio_score, work_id, corruption = \
        candidates[0]
    return work_id, corruption + "_title", set_score / 100.0


def lookup_by_alias(conn: sqlite3.Connection, alias_key: str) -> Optional[str]:
    """Return work_id for a given alias key, or None."""
    cur = conn.execute(
        """SELECT wa.work_id
           FROM work_aliases wa JOIN works w ON w.work_id = wa.work_id
           WHERE wa.alias_key = ?
           ORDER BY w.in_corpus DESC,
                    (w.bib_imported_at IS NOT NULL) DESC,
                    wa.work_id
           LIMIT 1""",
        (alias_key,),
    )
    row = cur.fetchone()
    return row[0] if row else None


def _normalized_ref_author_set(authors: List[str]) -> frozenset[str]:
    """Return the non-empty normalized surname set from a reference."""
    return frozenset(
        normalized
        for author in authors
        if (normalized := normalize_for_key(extract_surname_from_ref_author(author)))
    )


def _in_corpus_identity_index(
    conn: sqlite3.Connection,
) -> Dict[int, List[Tuple[str, str, str, frozenset[str]]]]:
    """Index immutable corpus candidates once for a materialization pass."""
    authors_by_work: Dict[str, set[str]] = {}
    for work_id, surname in conn.execute(
        """SELECT wa.work_id, wa.surname_normalized
           FROM work_authors wa JOIN works w ON w.work_id = wa.work_id
           WHERE w.in_corpus = 1 AND wa.surname_normalized != ''"""
    ):
        authors_by_work.setdefault(work_id, set()).add(surname)
    index: Dict[int, List[Tuple[str, str, str, frozenset[str]]]] = {}
    for work_id, title, year, doi in conn.execute(
        """SELECT work_id, title, year, doi FROM works
           WHERE in_corpus = 1 AND title IS NOT NULL AND year IS NOT NULL
           ORDER BY work_id"""
    ):
        index.setdefault(year, []).append((
            work_id,
            title or "",
            doi or "",
            frozenset(authors_by_work.get(work_id, set())),
        ))
    return index


def lookup_in_corpus_by_identity(
    conn: sqlite3.Connection,
    title: str,
    year: Optional[int],
    authors: List[str],
    incoming_doi: str = "",
    candidate_index: Optional[
        Dict[int, List[Tuple[str, str, str, frozenset[str]]]]
    ] = None,
) -> Optional[Tuple[str, str, float]]:
    """Return a unique title/year/author-set corpus work (#155, #225).

    Author order can differ between a paper header and a reference, so this
    safely crosses the first-author block used by the fuzzy cascade. Every
    normalized surname must still agree, both titles must be substantive, and
    the match must be unique. Non-exact titles use the same dual fuzzy guards
    as the existing first-author cascade. A different DOI does not erase the
    match, but is named by the persisted method so the verdict remains honest.
    """
    if not title or year is None:
        return None
    normalized_title = normalize_for_key(title)
    if (
        sum(character.isalpha() for character in normalized_title)
        < IDENTITY_TITLE_MIN_ALPHA
    ):
        return None
    author_set = _normalized_ref_author_set(authors)
    if not author_set:
        return None
    doi = normalize_doi(incoming_doi) if incoming_doi else ""
    matches = []
    if candidate_index is None:
        candidate_rows = _in_corpus_identity_index(conn).get(year, [])
    else:
        candidate_rows = candidate_index.get(year, [])
    for work_id, candidate_title, candidate_doi, candidate_authors in candidate_rows:
        normalized_candidate_title = normalize_for_key(candidate_title or "")
        if (
            sum(character.isalpha() for character in normalized_candidate_title)
            < IDENTITY_TITLE_MIN_ALPHA
        ):
            continue
        if candidate_authors != author_set:
            continue
        if normalized_candidate_title == normalized_title:
            method = "title_year_authors_exact"
            score = 1.0
        elif _HAS_BHL_DEPS:
            set_score = int(fuzz.token_set_ratio(
                normalized_title, normalized_candidate_title,
            ))
            ratio_score = int(fuzz.ratio(
                normalized_title, normalized_candidate_title,
            ))
            if set_score < 85 or ratio_score < 60:
                continue
            method = "title_year_authors_fuzzy"
            score = set_score / 100.0
        else:
            continue
        normalized_candidate_doi = (
            normalize_doi(candidate_doi) if candidate_doi else ""
        )
        doi_conflict = bool(
            doi and normalized_candidate_doi and doi != normalized_candidate_doi
        )
        if doi_conflict:
            method += "_doi_conflict"
            score = min(score, 0.95)
        matches.append((work_id, method, score))
    return matches[0] if len(matches) == 1 else None


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
               WHERE wa.surname_normalized = ? AND w.year = ? AND wa.position = 0
               ORDER BY wa.work_id""",
            (norm_surname, year),
        )
    else:
        cur = conn.execute(
            """SELECT DISTINCT wa.work_id, w.title
               FROM work_authors wa JOIN works w ON wa.work_id = w.work_id
               WHERE wa.surname_normalized = ? AND wa.position = 0
               ORDER BY wa.work_id""",
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
    refreshed = 0
    batch = 0
    for hash_dir in sorted(docs_dir.iterdir()):
        if not hash_dir.is_dir():
            continue
        meta_path = hash_dir / "metadata.json"
        if not meta_path.exists():
            continue

        corpus_hash = hash_dir.name
        try:
            current_mtime = meta_path.stat().st_mtime
        except OSError as e:
            logger.warning("Skipping %s: %s", meta_path, e)
            continue

        # Skip when we've already seeded this corpus_hash AND the
        # source metadata.json hasn't been regenerated since. When
        # metadata is newer than the stored mtime, fall through to
        # refresh the row's title/year/journal/doi/license/serve
        # fields and rebuild its author list.
        seen, stale = _artifact_state(
            conn, corpus_hash, "metadata", current_mtime,
        )
        existing_work_id: Optional[str] = None
        if seen and not stale:
            continue
        if seen and stale:
            cur = conn.execute(
                "SELECT work_id FROM works WHERE corpus_hash = ?",
                (corpus_hash,),
            )
            row = cur.fetchone()
            existing_work_id = row[0] if row else None

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

        if existing_work_id is not None and existing_work_id != work_id:
            # Identity changed (e.g. metadata.json gained a DOI that
            # rewrites work_id from corpus_key → doi). Citations table
            # has FK references on the old work_id — re-seeding under
            # a new work_id would orphan them. Leave the prior row in
            # place; the operator can address with --rebuild if they
            # need the identity update.
            logger.warning(
                "metadata.json for %s now resolves to work_id %r, but "
                "this corpus_hash is already bound to %r. Keeping the "
                "existing row to preserve citation graph integrity; "
                "use --rebuild to re-key.",
                corpus_hash, work_id, existing_work_id,
            )
            _record_artifact(conn, corpus_hash, "metadata", current_mtime)
            continue

        inserted = insert_work(conn, work_id, guid_type, title, year, journal,
                               doi, corpus_hash, in_corpus=True,
                               source="corpus_paper")
        if inserted:
            insert_authors(conn, work_id, authors)
            # Register alias for dedup matching
            if first_surname:
                alias = make_alias_key(first_surname, year, title or meta.get("filename", ""))
                insert_alias(conn, alias, work_id)
            # #51 + #54 — propagate license / skip fields if the bib-derived
            # metadata.json carries them (bib.parser.bib_entry_to_metadata
            # surfaces license / licenseurl / serve / servereason).
            _seed_license_and_serve(conn, work_id, meta)
            count += 1
        elif existing_work_id is not None:
            # Same work_id — refresh fields and rebuild the author list.
            now = time.time()
            conn.execute(
                """UPDATE works SET title=?, year=?, journal=?,
                       doi=?, updated_at=? WHERE work_id=?""",
                (title, year, journal, doi or None, now, work_id),
            )
            conn.execute("DELETE FROM work_authors WHERE work_id=?", (work_id,))
            insert_authors(conn, work_id, authors)
            if first_surname:
                alias = make_alias_key(first_surname, year, title or meta.get("filename", ""))
                insert_alias(conn, alias, work_id)
            _seed_license_and_serve(conn, work_id, meta)
            refreshed += 1

        _record_artifact(conn, corpus_hash, "metadata", current_mtime)
        batch += 1
        if batch >= 100:
            conn.commit()
            batch = 0
            logger.info(
                "Phase 1 progress: %d corpus papers seeded, %d refreshed",
                count, refreshed,
            )

    conn.commit()
    logger.info(
        "Phase 1 complete: %d corpus papers seeded, %d refreshed",
        count, refreshed,
    )
    return count


# ── Phase 2: Ingest references + citation graph ─────────────────────

def _resolve_reference(conn: sqlite3.Connection, ref: dict,
                       enrich_bhl: bool = False,
                       bhl_api_key: str = "",
                       bhl_max_year: Optional[int] = None,
                       fallback_key: str = "",
                       identity_index: Optional[
                           Dict[int, List[
                               Tuple[str, str, str, frozenset[str]]
                           ]]
                       ] = None) -> Tuple[str, str, float]:
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

    # Old references.json artifacts may predate #226's TEI fix and carry the
    # journal in both fields. Preserve the journal but do not let it become a
    # work title or a fuzzy-reconciliation feature.
    if title and journal and normalize_for_key(title) == normalize_for_key(journal):
        title = ""

    # Parse first author surname from ref format ("F Johnson")
    first_surname = ""
    if authors_raw:
        first_surname = extract_surname_from_ref_author(authors_raw[0])

    # ── Cascade step 1: DOI exact ──────────────────────────────────
    doi = normalize_doi(doi_raw) if doi_raw else ""
    if doi:
        existing = lookup_by_doi(conn, doi)
        if existing:
            existing_in_corpus = conn.execute(
                "SELECT in_corpus FROM works WHERE work_id = ?", (existing,),
            ).fetchone()[0]
            if existing_in_corpus:
                return existing, "doi_exact", 1.0
        variant = lookup_doi_variant_by_title(conn, doi, title)
        if variant is not None:
            return variant
        identity_match = lookup_in_corpus_by_identity(
            conn, title, year, authors_raw, incoming_doi=doi,
            candidate_index=identity_index,
        )
        if identity_match is not None:
            return identity_match
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

    identity_match = lookup_in_corpus_by_identity(
        conn, title, year, authors_raw, candidate_index=identity_index,
    )
    if identity_match is not None:
        return identity_match

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
        # No author at all — use the raw evidence, or the observation's
        # content address when Grobid returned a completely empty record.
        # Phase 2 always supplies ``fallback_key``; the timestamp remains only
        # for direct legacy callers that have no occurrence identity.
        fallback = title or raw or fallback_key or f"unknown_{time.time_ns()}"
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
# Key is the query string; value is (results, error). Historical
# taxonomic literature is author-heavy (the same monograph cited
# under many slight title variants), so the same (author, year) query
# can recur dozens of times across cited references — without this
# cache, stage 2 would re-issue the same broad BHL query each time.
_BHL_QUERY_CACHE: Dict[str, Tuple[Optional[list], Optional[str]]] = {}


def _bhl_search(query: str, api_key: str) -> Tuple[Optional[list], Optional[str]]:
    """Call BHL PublicationSearch API. Returns (results_list, error_string).

    Retry + backoff via :mod:`external`. Per-process circuit breaker
    (:data:`_BHL_BREAKER`) trips after ~20 consecutive failures so a BHL
    outage doesn't burn the entire enrichment pass on every paper.
    """
    cached = _BHL_QUERY_CACHE.get(query)
    if cached is not None:
        return cached

    def do_request() -> requests.Response:
        _BHL_BREAKER.check()
        # Built-in inter-call rate-limit nudge — BHL throttles harshly
        # without it.
        time.sleep(1)
        try:
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
        except (requests.ConnectionError, requests.Timeout):
            _BHL_BREAKER.record_failure()
            raise
        try:
            r.raise_for_status()
        except requests.HTTPError as e:
            if is_transient(e):
                _BHL_BREAKER.record_failure()
            raise
        _BHL_BREAKER.record_success()
        return r

    try:
        r = retry_with_backoff(do_request, max_attempts=3, base_delay=2.0)
        result = (r.json().get("Result", []), None)
    except CircuitOpenError as e:
        result = (None, f"circuit open: {e}")
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


def _reference_observation_id(corpus_hash: str, ordinal: int, ref: dict) -> str:
    """Content-address one occurrence in a citing paper's bibliography."""
    canonical = json.dumps(
        ref, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
        default=str,
    )
    digest = hashlib.sha256(
        f"{corpus_hash}\0{ordinal}\0{canonical}".encode("utf-8")
    ).hexdigest()
    return f"refobs:{digest}"


def _ingest_reference_observations(
    conn: sqlite3.Connection,
    corpus_hash: str,
    references: List[dict],
    source_fingerprint: str,
) -> None:
    """Append one evidence set and select it as current for its paper."""
    now = time.time()
    conn.execute(
        """INSERT OR IGNORE INTO reference_observation_sets
           (corpus_hash, source_fingerprint, observation_count, first_seen_at)
           VALUES (?, ?, ?, ?)""",
        (corpus_hash, source_fingerprint, len(references), now),
    )
    for ordinal, ref in enumerate(references):
        observation_id = _reference_observation_id(corpus_hash, ordinal, ref)
        authors = ref.get("authors") or []
        authors_json = json.dumps(authors, ensure_ascii=False, separators=(",", ":"))
        conn.execute(
            """INSERT OR IGNORE INTO reference_observations
               (observation_id, citing_corpus_hash, ordinal, grobid_xml_id,
                raw_citation, title, year, journal, doi, authors_json,
                first_seen_at)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                observation_id, corpus_hash, ordinal, ref.get("xml_id", ""),
                ref.get("raw", ""), ref.get("title", ""), ref.get("year"),
                ref.get("journal", ""), ref.get("doi", ""), authors_json,
                now,
            ),
        )
        conn.execute(
            """INSERT OR IGNORE INTO reference_observation_memberships
               (corpus_hash, source_fingerprint, ordinal, observation_id)
               VALUES (?, ?, ?, ?)""",
            (corpus_hash, source_fingerprint, ordinal, observation_id),
        )
    conn.execute(
        """INSERT INTO reference_current_sets
           (corpus_hash, source_fingerprint, selected_at)
           VALUES (?, ?, ?)
           ON CONFLICT(corpus_hash) DO UPDATE SET
             source_fingerprint = excluded.source_fingerprint,
             selected_at = excluded.selected_at""",
        (corpus_hash, source_fingerprint, now),
    )


def _clear_derived_reference_materialization(conn: sqlite3.Connection) -> None:
    """Remove re-derivable mappings/ghosts, never raw observations.

    BHL identities survive because the rate-limited lookup cache records only
    the outcome, not enough response data to recreate the BHL work row. Their
    aliases make a subsequent mapping reuse the same externally established
    identity without another network call.
    """
    derived_ids = [
        row[0] for row in conn.execute(
            """SELECT work_id FROM works
               WHERE source IN ('cited_reference', 'taxon_authority')
                 AND guid_type != 'bhl' AND in_corpus = 0
                 AND bib_imported_at IS NULL"""
        )
    ]
    conn.execute("DELETE FROM citations")
    conn.execute("DELETE FROM observation_work")
    for work_id in derived_ids:
        conn.execute("DELETE FROM taxon_work_links WHERE work_id = ?", (work_id,))
        conn.execute("DELETE FROM work_aliases WHERE work_id = ?", (work_id,))
        conn.execute("DELETE FROM work_authors WHERE work_id = ?", (work_id,))
        conn.execute("DELETE FROM works WHERE work_id = ?", (work_id,))


def _active_reference_observations(conn: sqlite3.Connection):
    return conn.execute(
        """SELECT ro.observation_id, ro.citing_corpus_hash, ro.ordinal,
                  ro.grobid_xml_id, ro.raw_citation, ro.title, ro.year,
                  ro.journal, ro.doi, ro.authors_json
           FROM reference_current_sets current
           JOIN reference_observation_memberships member
             ON member.corpus_hash = current.corpus_hash
            AND member.source_fingerprint = current.source_fingerprint
           JOIN reference_observations ro
             ON ro.observation_id = member.observation_id
           ORDER BY ro.observation_id"""
    ).fetchall()


def _rebuild_reference_materialization(
    conn: sqlite3.Connection,
    *,
    enrich_bhl: bool,
    bhl_api_key: str,
    bhl_max_year: Optional[int],
) -> Tuple[int, int]:
    """Derive current mappings and the frozen ``citations`` view from evidence."""
    observations = _active_reference_observations(conn)
    _clear_derived_reference_materialization(conn)
    known_work_ids = {row[0] for row in conn.execute("SELECT work_id FROM works")}
    identity_index = _in_corpus_identity_index(conn)
    citing_ids = {
        row[0]: row[1] for row in conn.execute(
            "SELECT corpus_hash, work_id FROM works WHERE corpus_hash IS NOT NULL"
        )
    }
    n_mapped = 0
    n_new_works = 0
    now = time.time()
    for (
        observation_id, corpus_hash, _ordinal, xml_id, raw_citation,
        title, year, journal, doi, authors_json,
    ) in observations:
        citing_work_id = citing_ids.get(corpus_hash)
        if citing_work_id is None:
            logger.warning(
                "No citing work for active reference observation %s (%s)",
                observation_id, corpus_hash,
            )
            continue
        try:
            authors = json.loads(authors_json or "[]")
        except (json.JSONDecodeError, TypeError):
            authors = []
        ref = {
            "xml_id": xml_id or "",
            "raw": raw_citation or "",
            "title": title or "",
            "year": year,
            "journal": journal or "",
            "doi": doi or "",
            "authors": authors,
        }
        cited_work_id, match_method, match_score = _resolve_reference(
            conn, ref, enrich_bhl=enrich_bhl,
            bhl_api_key=bhl_api_key, bhl_max_year=bhl_max_year,
            fallback_key=observation_id,
            identity_index=identity_index,
        )
        if cited_work_id not in known_work_ids:
            known_work_ids.add(cited_work_id)
            n_new_works += 1
        if match_method == "new":
            match_method = "new_work"
        conn.execute(
            """INSERT INTO observation_work
               (observation_id, work_id, match_method, match_score,
                producer_version, mapped_at)
               VALUES (?, ?, ?, ?, ?, ?)""",
            (
                observation_id, cited_work_id, match_method, match_score,
                REFERENCE_MAPPING_PRODUCER, now,
            ),
        )
        insert_citation(
            conn, citing_work_id, cited_work_id, corpus_hash,
            grobid_xml_id=xml_id or "", raw_citation=raw_citation or "",
            match_method=match_method, match_score=match_score,
        )
        n_mapped += 1
    return n_mapped, n_new_works


def phase2_references(conn: sqlite3.Connection, output_dir: Path,
                      enrich_bhl: bool = False,
                      bhl_api_key: str = "",
                      bhl_max_year: Optional[int] = None) -> Tuple[int, int]:
    """Ingest immutable observations, then derive mappings and citations.

    The expensive resolver runs over the complete active observation set when
    any source set changes or the mapping producer version changes. An
    unchanged run is a no-op. Therefore adding one paper and building cleanly
    produce the same current mapping rather than preserving document-order
    decisions from an earlier database (#240).
    """
    docs_dir = output_dir / "documents"
    changed = False
    n_skipped = 0
    n_refreshed = 0
    present_reference_hashes = set()
    for hash_dir in sorted(docs_dir.iterdir()):
        if not hash_dir.is_dir():
            continue
        refs_path = hash_dir / "references.json"
        if not refs_path.exists():
            continue
        corpus_hash = hash_dir.name
        present_reference_hashes.add(corpus_hash)
        try:
            refs_data = json.loads(refs_path.read_bytes())
        except (json.JSONDecodeError, OSError) as e:
            logger.warning("Skipping %s: %s", refs_path, e)
            continue
        if not isinstance(refs_data, dict):
            logger.warning("Skipping %s: top-level JSON is not an object", refs_path)
            continue
        references = refs_data.get("references", []) or []
        if not isinstance(references, list) or not all(
            isinstance(ref, dict) for ref in references
        ):
            logger.warning("Skipping %s: references is not a list of objects", refs_path)
            continue
        canonical_set = json.dumps(
            references, ensure_ascii=False, sort_keys=True,
            separators=(",", ":"), default=str,
        )
        source_fingerprint = hashlib.sha256(
            canonical_set.encode("utf-8")
        ).hexdigest()
        prior = conn.execute(
            "SELECT source_fingerprint FROM reference_current_sets "
            "WHERE corpus_hash = ?",
            (corpus_hash,),
        ).fetchone()
        if prior is not None and prior[0] == source_fingerprint:
            n_skipped += 1
            continue
        _ingest_reference_observations(
            conn, corpus_hash, references, source_fingerprint,
        )
        try:
            _record_artifact(
                conn, corpus_hash, "references", refs_path.stat().st_mtime,
            )
        except OSError:
            pass
        changed = True
        n_refreshed += int(prior is not None)

    # A removed paper or references artifact leaves historical evidence in the
    # table but must not contribute to the current graph.
    for (corpus_hash,) in conn.execute(
        "SELECT corpus_hash FROM reference_current_sets"
    ).fetchall():
        if corpus_hash in present_reference_hashes:
            continue
        conn.execute(
            "DELETE FROM reference_current_sets WHERE corpus_hash = ?",
            (corpus_hash,),
        )
        changed = True

    active_count = len(_active_reference_observations(conn))
    # A legacy ``works.corpus_hash`` can represent only one member when
    # several documents share one canonical work identity (for example,
    # volumes carrying the same BHL DOI). Those observations cannot yet form
    # citation edges. Do not interpret the known structural gap as a stale
    # producer and rematerialize forever; the QC report exposes the unmapped
    # population for review.
    mappable_active_count = conn.execute(
        """SELECT COUNT(*)
           FROM reference_current_sets current
           JOIN reference_observation_memberships member
             ON member.corpus_hash = current.corpus_hash
            AND member.source_fingerprint = current.source_fingerprint
           JOIN reference_observations ro
             ON ro.observation_id = member.observation_id
           WHERE EXISTS (
             SELECT 1 FROM works w
             WHERE w.in_corpus = 1
               AND w.corpus_hash = ro.citing_corpus_hash
           )"""
    ).fetchone()[0]
    current_mapping_count = conn.execute(
        """SELECT COUNT(*)
           FROM observation_work ow
           JOIN reference_observation_memberships member
             ON member.observation_id = ow.observation_id
           JOIN reference_current_sets current
             ON current.corpus_hash = member.corpus_hash
            AND current.source_fingerprint = member.source_fingerprint
           WHERE ow.producer_version = ?""",
        (REFERENCE_MAPPING_PRODUCER,),
    ).fetchone()[0]
    if current_mapping_count != mappable_active_count:
        changed = True

    if changed:
        n_citations, n_new_works = _rebuild_reference_materialization(
            conn, enrich_bhl=enrich_bhl, bhl_api_key=bhl_api_key,
            bhl_max_year=bhl_max_year,
        )
    else:
        n_citations = n_new_works = 0
    conn.commit()
    logger.info(
        "Phase 2 complete: %d active observations, %d mapped this run, "
        "%d derived works created (%d source sets skipped, %d refreshed)",
        active_count, n_citations, n_new_works, n_skipped, n_refreshed,
    )
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
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Report what would be processed without opening the SQLite or "
             "calling external services.",
    )
    parser.add_argument(
        "--strict-network", action="store_true",
        help="Fail fast on the first transient external-service failure "
             "(BHL 5xx, connect error, timeout) instead of retrying.",
    )
    parser.add_argument(
        "--pd-cutoff-years", type=int, default=None,
        help="Public-domain cutoff for the publishable-derivation pass "
             "(#51). Default: read from config.yaml's `licensing."
             "pd_cutoff_years`, falling back to 95 (US copyright rule).",
    )
    parser.add_argument(
        "--config", type=Path, default=None,
        help="Path to per-corpuscle config.yaml. Read for "
             "`licensing.pd_cutoff_years` (#51). Defaults to "
             "<output_dir>/../config.yaml if found.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.strict_network:
        from pipeline.external import set_strict_network
        set_strict_network(True)

    if args.output is None:
        args.output = args.output_dir / "biblio_authority.sqlite"
    if args.taxonomy_db is None:
        args.taxonomy_db = args.output_dir / "taxonomy.sqlite"

    if args.dry_run:
        documents_dir = args.output_dir / "documents"
        # No documents/ means the corpuscle has not been extracted yet;
        # "0 hash dirs" is the honest plan. See pipeline/taxon_mentions.py.
        hash_dirs = (
            [d for d in documents_dir.iterdir() if d.is_dir()]
            if documents_dir.is_dir() else []
        )
        n_papers = len(hash_dirs)
        n_metadata = sum(1 for d in hash_dirs if (d / "metadata.json").exists())
        n_refs = sum(1 for d in hash_dirs if (d / "references.json").exists())
        action = "rebuild from scratch" if args.rebuild else "incrementally update"
        logger.info(
            "Dry-run: would %s %s. Source: %d hash dirs (%d with metadata.json, "
            "%d with references.json). Taxonomy: %s. BHL enrichment: %s. "
            "No SQLite or network writes.",
            action, args.output, n_papers, n_metadata, n_refs,
            args.taxonomy_db,
            "ON" if args.enrich_bhl else "off",
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
            # Preserve bhl_lookups across rebuilds — it's the cache of
            # rate-limited BHL API calls. Dropping it forces every
            # lookup to be re-paid the next time ``--enrich-bhl`` is
            # used, which takes hours and is unrelated to the schema
            # churn --rebuild is meant to address.
            logger.info("Rebuilding: dropping all tables except bhl_lookups")
            conn.executescript("""
                DROP TABLE IF EXISTS taxon_work_links;
                DROP TABLE IF EXISTS citations;
                DROP TABLE IF EXISTS observation_work;
                DROP TABLE IF EXISTS work_reconciliation_decisions;
                DROP TABLE IF EXISTS reference_current_sets;
                DROP TABLE IF EXISTS reference_observation_memberships;
                DROP TABLE IF EXISTS reference_observation_sets;
                DROP TABLE IF EXISTS reference_observations;
                DROP TABLE IF EXISTS work_aliases;
                DROP TABLE IF EXISTS work_authors;
                DROP TABLE IF EXISTS works;
                DROP TABLE IF EXISTS build_meta;
                DROP TABLE IF EXISTS paper_artifacts_processed;
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

        # #51 — derive publishable from license + age cutoff. Runs at
        # the end of every build so the column is always consistent
        # with the current license values + the configured cutoff.
        pd_cutoff = args.pd_cutoff_years
        if pd_cutoff is None:
            pd_cutoff = _resolve_pd_cutoff_from_config(args.config)
        n_derived = apply_publishable_derivation(conn, pd_cutoff_years=pd_cutoff)
        n_publishable = conn.execute(
            "SELECT COUNT(*) FROM works WHERE publishable = 1"
        ).fetchone()[0]
        logger.info(
            "  Publishable:  %d / %d works (pd_cutoff_years=%d)",
            n_publishable, n_derived, pd_cutoff,
        )

        conn.commit()
        rc = 0
    finally:
        conn.close()

    # #67: apply any sidecar staged by `corpus bib import` before the
    # first `corpus run`. Closing + reopening the connection keeps the
    # transaction boundaries clean.
    try:
        from .importer import apply_pending_bib
        applied = apply_pending_bib(args.output_dir, args.output)
        if applied is not None:
            logger.info(
                "Applied staged bib overrides: %d updated, %d no-change, %d no-match",
                applied.get("changed", 0),
                applied.get("no_changes", 0),
                applied.get("no_match", 0),
            )
    except Exception as e:
        logger.warning("Could not apply staged bib overrides: %s", e)

    return rc


if __name__ == "__main__":
    sys.exit(main())
