#!/usr/bin/env python3
"""Build a Darwin Core taxonomy SQLite snapshot from one of several sources.

The corpus pipeline consumes a single SQLite file (default
``resources/taxonomy.sqlite``) at run time. This tool produces it from:

* ``--source dwc PATH``   — a DwC Taxon CSV/TSV (recommended; download
                            from any DwC publisher, e.g. WoRMS, GBIF,
                            iNaturalist, ITIS, or hand-curated).
* ``--source dwca PATH``  — a Darwin Core Archive (.zip with meta.xml +
                            Taxon.tsv), or an extracted directory, or a
                            bare ``Taxon.tsv``.
* ``--source worms``      — walks the WoRMS REST API from a root AphiaID
                            and writes the same DwC schema. Network +
                            rate-limited (~0.3 s per call).

Optional ``--root-id <taxonID>`` prunes the snapshot to the subgraph
rooted at that taxon: useful for ingesting only Siphonophorae out of the
full WoRMS export, for example. Synonyms whose ``acceptedNameUsageID``
points into the kept set are preserved; everything else is dropped.

Schema (Darwin Core Taxon class, plus an indexed ``names`` table for
lookup):

    taxa
      taxon_id (PK), scientific_name, scientific_name_authorship,
      taxon_rank, taxonomic_status,
      parent_name_usage_id, accepted_name_usage_id, accepted_name,
      kingdom, phylum, class, "order", family, genus,
      source, source_record_id, citation, fetched_at
    names
      name, name_lowercase, taxon_id, name_type
        name_type ∈ {accepted, unaccepted, synonym}
    meta
      key, value     -- root_id, source, etc.

Usage::

    # From a downloaded DwC CSV/TSV
    python ingest_taxonomy.py --source dwc /path/to/Taxon.tsv \
        --root-id urn:lsid:marinespecies.org:taxname:1371
    # Live WoRMS walk
    python ingest_taxonomy.py --source worms --root-id 1371
    # Whole archive
    python ingest_taxonomy.py --source dwca /path/to/dwca.zip
"""

from __future__ import annotations

import argparse
import csv
import io
import logging
import signal
import sqlite3
import sys
import time
import zipfile
from collections import defaultdict, deque
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional


logger = logging.getLogger("ingest_taxonomy")


DEFAULT_OUTPUT = Path(__file__).resolve().parent / "resources" / "taxonomy.sqlite"


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------


# Optional rank columns we denormalize when the source provides them.
# Quoted because "order" is a SQL reserved word.
_DENORM_RANK_COLUMNS = ["kingdom", "phylum", "class", "order", "family", "genus"]


def create_schema(conn: sqlite3.Connection) -> None:
    """Create the DwC-aligned schema if it doesn't exist.

    Field names match Darwin Core Taxon-class terms verbatim (snake_case
    in SQL, camelCase in DwC). ``taxon_id`` is TEXT because DwC
    ``taxonID`` is a string identifier — for WoRMS we render the AphiaID
    as text, for GBIF it's already a string, and for hand-curated files
    the user can supply any unique key.
    """
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS taxa (
            taxon_id                       TEXT PRIMARY KEY,
            scientific_name                TEXT NOT NULL,
            scientific_name_authorship     TEXT,
            taxon_rank                     TEXT,
            taxonomic_status               TEXT,
            parent_name_usage_id           TEXT,
            accepted_name_usage_id         TEXT,
            accepted_name                  TEXT,
            kingdom                        TEXT,
            phylum                         TEXT,
            class                          TEXT,
            "order"                        TEXT,
            family                         TEXT,
            genus                          TEXT,
            source                         TEXT,
            source_record_id               TEXT,
            citation                       TEXT,
            fetched_at                     REAL
        );

        CREATE TABLE IF NOT EXISTS names (
            name            TEXT NOT NULL,
            name_lowercase  TEXT NOT NULL,
            taxon_id        TEXT NOT NULL,
            name_type       TEXT NOT NULL,
            FOREIGN KEY (taxon_id) REFERENCES taxa(taxon_id)
        );

        CREATE INDEX IF NOT EXISTS idx_names_lowercase ON names(name_lowercase);
        CREATE INDEX IF NOT EXISTS idx_names_taxon_id ON names(taxon_id);
        CREATE INDEX IF NOT EXISTS idx_taxa_parent ON taxa(parent_name_usage_id);
        CREATE INDEX IF NOT EXISTS idx_taxa_accepted ON taxa(accepted_name_usage_id);

        CREATE TABLE IF NOT EXISTS meta (
            key   TEXT PRIMARY KEY,
            value TEXT
        );
        """
    )
    conn.commit()


# ---------------------------------------------------------------------------
# Canonical record (the in-memory shape we accumulate before writing).
# Every source backend yields these; pruning + insert work on them.
# ---------------------------------------------------------------------------


def make_record(
    *,
    taxon_id: str,
    scientific_name: str,
    scientific_name_authorship: Optional[str] = None,
    taxon_rank: Optional[str] = None,
    taxonomic_status: Optional[str] = None,
    parent_name_usage_id: Optional[str] = None,
    accepted_name_usage_id: Optional[str] = None,
    accepted_name: Optional[str] = None,
    higher: Optional[Dict[str, Optional[str]]] = None,
    source: Optional[str] = None,
    source_record_id: Optional[str] = None,
    citation: Optional[str] = None,
    extra_names: Optional[List[str]] = None,
) -> Dict:
    """Build a normalized in-memory record. ``extra_names`` is for
    synonym names attached to an accepted record (typical of the WoRMS
    REST API, which returns synonyms as a separate /synonyms endpoint
    rather than as their own rows)."""
    rec = {
        "taxon_id": str(taxon_id),
        "scientific_name": scientific_name or "",
        "scientific_name_authorship": scientific_name_authorship,
        "taxon_rank": taxon_rank,
        "taxonomic_status": taxonomic_status,
        "parent_name_usage_id": str(parent_name_usage_id) if parent_name_usage_id is not None else None,
        "accepted_name_usage_id": str(accepted_name_usage_id) if accepted_name_usage_id is not None else None,
        "accepted_name": accepted_name,
        "source": source,
        "source_record_id": source_record_id,
        "citation": citation,
        "extra_names": list(extra_names or []),
    }
    higher = higher or {}
    for k in _DENORM_RANK_COLUMNS:
        rec[k] = higher.get(k)
    return rec


# ---------------------------------------------------------------------------
# Source: DwC CSV / TSV
# ---------------------------------------------------------------------------


# Common header aliases. DwC files in the wild come in three flavors:
# bare camelCase ("taxonID"), namespaced ("dwc:taxonID"), or with
# trailing "_id" Linnean-rank columns. We accept all of these.
_DWC_ALIASES = {
    "taxon_id":                   ["taxonID", "dwc:taxonID", "taxon_id", "id"],
    "scientific_name":            ["scientificName", "dwc:scientificName", "scientific_name"],
    "scientific_name_authorship": ["scientificNameAuthorship", "dwc:scientificNameAuthorship",
                                   "namePublishedInYear", "authorship"],
    "taxon_rank":                 ["taxonRank", "dwc:taxonRank", "rank"],
    "taxonomic_status":           ["taxonomicStatus", "dwc:taxonomicStatus", "status"],
    "parent_name_usage_id":       ["parentNameUsageID", "dwc:parentNameUsageID", "parent_id"],
    "accepted_name_usage_id":     ["acceptedNameUsageID", "dwc:acceptedNameUsageID",
                                   "valid_AphiaID", "accepted_id"],
    "accepted_name":              ["acceptedNameUsage", "dwc:acceptedNameUsage", "valid_name"],
    "kingdom":                    ["kingdom", "dwc:kingdom"],
    "phylum":                     ["phylum", "dwc:phylum"],
    "class":                      ["class", "dwc:class"],
    "order":                      ["order", "dwc:order"],
    "family":                     ["family", "dwc:family"],
    "genus":                      ["genus", "dwc:genus"],
    "citation":                   ["bibliographicCitation", "dwc:bibliographicCitation", "citation"],
}


def _build_header_index(fieldnames: Iterable[str]) -> Dict[str, str]:
    """Map our internal field name → the actual CSV header that supplies it.

    Picks the first alias hit per field; fields with no header in the
    file are simply absent from the returned dict.
    """
    available = {name: name for name in fieldnames if name}
    lower_lookup = {name.lower(): name for name in available}
    out: Dict[str, str] = {}
    for canonical, aliases in _DWC_ALIASES.items():
        for a in aliases:
            if a in available:
                out[canonical] = a
                break
            if a.lower() in lower_lookup:
                out[canonical] = lower_lookup[a.lower()]
                break
    return out


def _empty_to_none(v: Any) -> Optional[str]:
    if v is None:
        return None
    s = str(v).strip()
    return s or None


def iter_dwc_csv(path: Path) -> Iterator[Dict]:
    """Yield canonical records from a DwC-flavored CSV/TSV file.

    Sniffs the delimiter from the first 4 KiB; defaults to tab on
    failure (the most common DwC export). Required columns are
    ``taxonID`` and ``scientificName``; everything else is optional.
    """
    with open(path, "r", encoding="utf-8", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
            delimiter = dialect.delimiter
        except csv.Error:
            delimiter = "\t"
        reader = csv.DictReader(f, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"{path}: no header row found")
        idx = _build_header_index(reader.fieldnames)
        if "taxon_id" not in idx or "scientific_name" not in idx:
            raise ValueError(
                f"{path}: missing required DwC columns. Need at least "
                f"taxonID and scientificName; saw headers {reader.fieldnames!r}"
            )
        for row in reader:
            taxon_id = _empty_to_none(row.get(idx["taxon_id"]))
            sci = _empty_to_none(row.get(idx["scientific_name"]))
            if not taxon_id or not sci:
                continue
            higher = {
                k: _empty_to_none(row.get(idx[k])) if k in idx else None
                for k in _DENORM_RANK_COLUMNS
            }
            yield make_record(
                taxon_id=taxon_id,
                scientific_name=sci,
                scientific_name_authorship=_empty_to_none(row.get(idx["scientific_name_authorship"])) if "scientific_name_authorship" in idx else None,
                taxon_rank=_empty_to_none(row.get(idx["taxon_rank"])) if "taxon_rank" in idx else None,
                taxonomic_status=_empty_to_none(row.get(idx["taxonomic_status"])) if "taxonomic_status" in idx else None,
                parent_name_usage_id=_empty_to_none(row.get(idx["parent_name_usage_id"])) if "parent_name_usage_id" in idx else None,
                accepted_name_usage_id=_empty_to_none(row.get(idx["accepted_name_usage_id"])) if "accepted_name_usage_id" in idx else None,
                accepted_name=_empty_to_none(row.get(idx["accepted_name"])) if "accepted_name" in idx else None,
                higher=higher,
                source="dwc",
                source_record_id=taxon_id,
                citation=_empty_to_none(row.get(idx["citation"])) if "citation" in idx else None,
            )


# ---------------------------------------------------------------------------
# Source: Darwin Core Archive (.zip or extracted dir or bare Taxon.tsv)
# ---------------------------------------------------------------------------


def iter_dwca(path: Path) -> Iterator[Dict]:
    """Yield records from a DwC-A. Accepts a .zip, an extracted directory,
    or a bare ``Taxon.tsv``. Looks for any file whose lowercased name
    starts with ``taxon`` and ends in .tsv/.csv/.txt — DwC archives use
    ``taxon.txt`` by convention but we don't insist on it.

    For the .zip case we extract the matching member to a temp buffer.
    Avoids requiring a full DwC-A library (python-dwca-reader); we only
    need the Taxon core for taxonomy work.
    """
    p = Path(path)
    if p.is_file() and p.suffix.lower() == ".zip":
        with zipfile.ZipFile(p) as zf:
            taxon_member = next(
                (n for n in zf.namelist()
                 if Path(n).name.lower().startswith("taxon")
                 and Path(n).suffix.lower() in (".tsv", ".csv", ".txt")),
                None,
            )
            if taxon_member is None:
                raise ValueError(f"{p}: no Taxon.tsv-like file in archive")
            logger.info("Reading %s from %s", taxon_member, p)
            with zf.open(taxon_member) as raw:
                # csv.DictReader needs text; wrap.
                text = io.TextIOWrapper(raw, encoding="utf-8", newline="")
                yield from _iter_csv_handle(text, source_label=f"dwca:{p.name}")
        return
    if p.is_dir():
        cand = next(
            (q for q in p.iterdir()
             if q.is_file()
             and q.name.lower().startswith("taxon")
             and q.suffix.lower() in (".tsv", ".csv", ".txt")),
            None,
        )
        if cand is None:
            raise ValueError(f"{p}: no Taxon.tsv-like file in directory")
        logger.info("Reading %s", cand)
        yield from iter_dwc_csv(cand)
        return
    # Treat as a bare file.
    yield from iter_dwc_csv(p)


def _iter_csv_handle(handle, *, source_label: str) -> Iterator[Dict]:
    """Like ``iter_dwc_csv`` but operates on an already-open text handle."""
    sample = handle.read(4096)
    handle.seek(0)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        delimiter = dialect.delimiter
    except csv.Error:
        delimiter = "\t"
    reader = csv.DictReader(handle, delimiter=delimiter)
    if not reader.fieldnames:
        raise ValueError(f"{source_label}: no header row")
    idx = _build_header_index(reader.fieldnames)
    if "taxon_id" not in idx or "scientific_name" not in idx:
        raise ValueError(
            f"{source_label}: missing required DwC columns "
            f"(taxonID + scientificName); saw {reader.fieldnames!r}"
        )
    for row in reader:
        taxon_id = _empty_to_none(row.get(idx["taxon_id"]))
        sci = _empty_to_none(row.get(idx["scientific_name"]))
        if not taxon_id or not sci:
            continue
        higher = {
            k: _empty_to_none(row.get(idx[k])) if k in idx else None
            for k in _DENORM_RANK_COLUMNS
        }
        yield make_record(
            taxon_id=taxon_id,
            scientific_name=sci,
            scientific_name_authorship=_empty_to_none(row.get(idx["scientific_name_authorship"])) if "scientific_name_authorship" in idx else None,
            taxon_rank=_empty_to_none(row.get(idx["taxon_rank"])) if "taxon_rank" in idx else None,
            taxonomic_status=_empty_to_none(row.get(idx["taxonomic_status"])) if "taxonomic_status" in idx else None,
            parent_name_usage_id=_empty_to_none(row.get(idx["parent_name_usage_id"])) if "parent_name_usage_id" in idx else None,
            accepted_name_usage_id=_empty_to_none(row.get(idx["accepted_name_usage_id"])) if "accepted_name_usage_id" in idx else None,
            accepted_name=_empty_to_none(row.get(idx["accepted_name"])) if "accepted_name" in idx else None,
            higher=higher,
            source=source_label,
            source_record_id=taxon_id,
            citation=_empty_to_none(row.get(idx["citation"])) if "citation" in idx else None,
        )


# ---------------------------------------------------------------------------
# Source: WoRMS REST API
# ---------------------------------------------------------------------------


WORMS_API = "https://www.marinespecies.org/rest"


class _RateLimitedSession:
    """Wrap requests so we don't hammer the WoRMS API faster than
    ``min_interval`` seconds per call. The API docs ask for courtesy."""

    def __init__(self, min_interval: float = 0.3):
        import requests
        self.min_interval = min_interval
        self._last_call = 0.0
        self.session = requests.Session()

    def get(self, url: str, **kw):
        delay = self.min_interval - (time.monotonic() - self._last_call)
        if delay > 0:
            time.sleep(delay)
        self._last_call = time.monotonic()
        return self.session.get(url, **kw)


def _worms_rec_to_canonical(rec: Dict, *, name_type: str) -> Dict:
    """Map a WoRMS AphiaRecord to our canonical in-memory record."""
    # WoRMS uses ``valid_AphiaID``/``valid_name`` for accepted-name pointers
    # (its older API predates the DwC ``acceptedNameUsageID`` term).
    higher = {
        "kingdom": rec.get("kingdom"),
        "phylum": rec.get("phylum"),
        "class": rec.get("class"),
        "order": rec.get("order"),
        "family": rec.get("family"),
        "genus": rec.get("genus"),
    }
    accepted_id = rec.get("valid_AphiaID")
    return make_record(
        taxon_id=str(rec.get("AphiaID")),
        scientific_name=rec.get("scientificname") or "",
        scientific_name_authorship=rec.get("authority"),
        taxon_rank=rec.get("rank"),
        taxonomic_status=rec.get("status"),
        parent_name_usage_id=str(rec["parentNameUsageID"]) if rec.get("parentNameUsageID") is not None else None,
        accepted_name_usage_id=str(accepted_id) if accepted_id is not None else None,
        accepted_name=rec.get("valid_name"),
        higher=higher,
        source="worms",
        source_record_id=str(rec.get("AphiaID")),
        citation=rec.get("citation"),
    )


def iter_worms_walk(root_id: str, *, min_interval: float = 0.3) -> Iterator[Dict]:
    """Breadth-first walk from a WoRMS AphiaID. Yields canonical records.

    For accepted taxa we also fetch synonyms; each synonym is yielded as
    its own (unaccepted) record AND attached as an ``extra_names`` entry
    on the accepted record so name-based lookup hits the accepted taxon
    directly.
    """
    sess = _RateLimitedSession(min_interval=min_interval)
    visited: set = set()
    queue: deque = deque([root_id])
    while queue:
        aid = queue.popleft()
        if aid in visited:
            continue
        visited.add(aid)
        try:
            r = sess.get(f"{WORMS_API}/AphiaRecordByAphiaID/{aid}")
            if r.status_code == 204:
                continue
            r.raise_for_status()
            rec = r.json()
        except Exception as e:
            logger.warning("worms: record fetch failed for %s: %s", aid, e)
            continue
        if not rec:
            continue
        name_type = "accepted" if rec.get("status") == "accepted" else "unaccepted"
        canonical = _worms_rec_to_canonical(rec, name_type=name_type)

        # Synonyms — attach their names to the accepted record AND emit
        # them as their own unaccepted records (so cross-references via
        # accepted_name_usage_id resolve).
        if rec.get("status") == "accepted":
            try:
                syn_resp = _worms_paged(sess, f"{WORMS_API}/AphiaSynonymsByAphiaID/{aid}")
            except Exception as e:
                logger.warning("worms: synonyms fetch failed for %s: %s", aid, e)
                syn_resp = []
            for syn in syn_resp:
                syn_id = syn.get("AphiaID")
                if not syn_id:
                    continue
                syn_rec = _worms_rec_to_canonical(syn, name_type="synonym")
                # Force the link: the synonym's accepted_name points at
                # the parent record so the lookup resolves correctly.
                syn_rec["accepted_name_usage_id"] = canonical["taxon_id"]
                syn_rec["accepted_name"] = canonical["scientific_name"]
                syn_rec["taxonomic_status"] = "synonym"
                yield syn_rec
                if syn_rec["scientific_name"]:
                    canonical["extra_names"].append(syn_rec["scientific_name"])

        yield canonical

        # Enqueue children for recursion.
        try:
            kids = _worms_paged(sess, f"{WORMS_API}/AphiaChildrenByAphiaID/{aid}",
                                params={"marine_only": "false"})
        except Exception as e:
            logger.warning("worms: children fetch failed for %s: %s", aid, e)
            kids = []
        for child in kids:
            cid = str(child.get("AphiaID")) if child.get("AphiaID") is not None else None
            if cid and cid not in visited:
                queue.append(cid)


def _worms_paged(sess: _RateLimitedSession, url: str, params: Optional[Dict] = None) -> List[Dict]:
    """Page through a WoRMS endpoint. Stops on 204 or short page."""
    out: List[Dict] = []
    offset = 1
    while True:
        full_params = dict(params or {})
        full_params["offset"] = offset
        r = sess.get(url, params=full_params)
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


# ---------------------------------------------------------------------------
# Pruning to a subgraph rooted at a taxon_id
# ---------------------------------------------------------------------------


def prune_to_subgraph(records: List[Dict], root_id: str) -> List[Dict]:
    """Keep only records descended from ``root_id`` (plus their synonyms).

    Walks down ``parent_name_usage_id`` from the root using a children
    index. Synonyms are kept iff their ``accepted_name_usage_id`` falls
    in the descendant set. The root itself is included.
    """
    by_id: Dict[str, Dict] = {r["taxon_id"]: r for r in records}
    if root_id not in by_id:
        logger.warning(
            "prune: root_id %r not found among %d records — keeping nothing",
            root_id, len(records),
        )
        return []
    children: Dict[str, List[str]] = defaultdict(list)
    for r in records:
        p = r.get("parent_name_usage_id")
        if p:
            children[p].append(r["taxon_id"])

    keep: set = set()
    queue: deque = deque([root_id])
    while queue:
        cur = queue.popleft()
        if cur in keep:
            continue
        keep.add(cur)
        for c in children.get(cur, ()):
            if c not in keep:
                queue.append(c)

    # Pull in synonyms whose accepted_name_usage_id resolves into the kept set.
    for r in records:
        acc = r.get("accepted_name_usage_id")
        if acc and acc in keep:
            keep.add(r["taxon_id"])

    pruned = [by_id[t] for t in keep if t in by_id]
    logger.info(
        "prune: kept %d / %d records under root %r", len(pruned), len(records), root_id,
    )
    return pruned


# ---------------------------------------------------------------------------
# Insert
# ---------------------------------------------------------------------------


def insert_records(conn: sqlite3.Connection, records: Iterable[Dict]) -> Dict[str, int]:
    """Write canonical records to taxa + names. Idempotent (REPLACE).

    Synonym names attached via ``extra_names`` are inserted into
    ``names`` against the *accepted* taxon row, so a single name lookup
    resolves directly to the accepted record.
    """
    n_taxa = 0
    n_names = 0
    now = time.time()
    seen_name_keys: set = set()  # (name_lowercase, taxon_id, name_type) — dedup

    for rec in records:
        conn.execute(
            """
            INSERT OR REPLACE INTO taxa (
                taxon_id, scientific_name, scientific_name_authorship,
                taxon_rank, taxonomic_status,
                parent_name_usage_id, accepted_name_usage_id, accepted_name,
                kingdom, phylum, class, "order", family, genus,
                source, source_record_id, citation, fetched_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                rec["taxon_id"],
                rec["scientific_name"],
                rec.get("scientific_name_authorship"),
                rec.get("taxon_rank"),
                rec.get("taxonomic_status"),
                rec.get("parent_name_usage_id"),
                rec.get("accepted_name_usage_id"),
                rec.get("accepted_name"),
                rec.get("kingdom"),
                rec.get("phylum"),
                rec.get("class"),
                rec.get("order"),
                rec.get("family"),
                rec.get("genus"),
                rec.get("source"),
                rec.get("source_record_id"),
                rec.get("citation"),
                now,
            ),
        )
        n_taxa += 1

        # Primary name. ``name_type`` is derived from taxonomic_status:
        # we treat anything that isn't "accepted" as either "synonym"
        # (when there's an accepted_name_usage_id) or "unaccepted".
        status = (rec.get("taxonomic_status") or "").lower()
        if status == "accepted":
            primary_type = "accepted"
        elif rec.get("accepted_name_usage_id") and rec["accepted_name_usage_id"] != rec["taxon_id"]:
            primary_type = "synonym"
        else:
            primary_type = "unaccepted"

        primary_name = rec.get("scientific_name")
        if primary_name:
            # For synonyms, register the name against the *accepted* taxon
            # so name lookup follows synonymy in one query.
            target_id = (rec.get("accepted_name_usage_id")
                         if primary_type == "synonym"
                         else rec["taxon_id"])
            key = (primary_name.lower(), target_id, primary_type)
            if key not in seen_name_keys:
                conn.execute(
                    "INSERT INTO names (name, name_lowercase, taxon_id, name_type) "
                    "VALUES (?, ?, ?, ?)",
                    (primary_name, primary_name.lower(), target_id, primary_type),
                )
                seen_name_keys.add(key)
                n_names += 1

        # Extra names (e.g. synonyms attached to an accepted record by
        # the WoRMS walker). All filed against this record's taxon_id as
        # synonym, since this record is accepted.
        for extra in rec.get("extra_names") or []:
            key = (extra.lower(), rec["taxon_id"], "synonym")
            if key not in seen_name_keys:
                conn.execute(
                    "INSERT INTO names (name, name_lowercase, taxon_id, name_type) "
                    "VALUES (?, ?, ?, ?)",
                    (extra, extra.lower(), rec["taxon_id"], "synonym"),
                )
                seen_name_keys.add(key)
                n_names += 1

    conn.commit()
    return {"taxa": n_taxa, "names": n_names}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--source", required=True, choices=["dwc", "dwca", "worms"],
        help="Where to load taxonomy data from. 'dwc' is a single CSV/TSV "
             "with DwC Taxon headers; 'dwca' is a Darwin Core Archive "
             "(.zip / dir / Taxon.tsv); 'worms' walks the WoRMS REST API.",
    )
    parser.add_argument(
        "--input", type=Path, default=None,
        help="Path to the input file (required for --source dwc and dwca).",
    )
    parser.add_argument(
        "--root-id", type=str, default=None,
        help="Optional taxonID. For --source worms this is the AphiaID to "
             "start walking from (required). For dwc/dwca this prunes the "
             "ingest to descendants of this taxon (plus their synonyms).",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=DEFAULT_OUTPUT,
        help=f"Output SQLite path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--rebuild", action="store_true",
        help="Drop and recreate tables before inserting (otherwise REPLACE-merges).",
    )
    parser.add_argument(
        "--min-interval", type=float, default=0.3,
        help="WoRMS API: minimum seconds between calls (default: 0.3).",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.source in ("dwc", "dwca") and args.input is None:
        parser.error(f"--source {args.source} requires --input PATH")
    if args.source == "worms" and args.root_id is None:
        parser.error("--source worms requires --root-id <AphiaID>")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(args.output)

    def _sigint(signum, frame):  # noqa: ARG001
        logger.info("SIGINT received; committing and exiting")
        conn.commit()
        conn.close()
        sys.exit(130)
    signal.signal(signal.SIGINT, _sigint)

    try:
        if args.rebuild:
            logger.info("Rebuilding: dropping taxa + names tables")
            conn.executescript(
                "DROP TABLE IF EXISTS taxa; "
                "DROP TABLE IF EXISTS names; "
                "DROP TABLE IF EXISTS meta;"
            )
        create_schema(conn)

        # Load source records. WoRMS streams (the walk is itself a BFS
        # from the root, so pruning is implicit). DwC/DwCA loads fully
        # into memory because we need the full record set to compute
        # descendant pruning.
        if args.source == "worms":
            logger.info("Walking WoRMS from AphiaID %s ...", args.root_id)
            records = list(iter_worms_walk(args.root_id, min_interval=args.min_interval))
        elif args.source == "dwc":
            logger.info("Reading DwC CSV/TSV %s ...", args.input)
            records = list(iter_dwc_csv(args.input))
        else:  # dwca
            logger.info("Reading DwC-A %s ...", args.input)
            records = list(iter_dwca(args.input))
        logger.info("Loaded %d candidate records", len(records))

        # For dwc/dwca, prune to subgraph if requested. (For worms the
        # walk already produced only descendants of root_id.)
        if args.root_id and args.source != "worms":
            records = prune_to_subgraph(records, args.root_id)

        if not records:
            logger.error("No records to insert. Nothing written.")
            return 2

        stats = insert_records(conn, records)
        logger.info("Inserted %d taxa, %d names", stats["taxa"], stats["names"])

        # Stamp meta.
        conn.execute(
            "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
            ("source", args.source),
        )
        if args.root_id:
            conn.execute(
                "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
                ("root_id", args.root_id),
            )
        if args.input:
            conn.execute(
                "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
                ("source_path", str(args.input)),
            )
        conn.execute(
            "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
            ("last_ingest_ts", str(time.time())),
        )
        conn.commit()

        # Quick rank distribution for the log.
        cur = conn.execute(
            "SELECT taxon_rank, COUNT(*) FROM taxa GROUP BY taxon_rank "
            "ORDER BY COUNT(*) DESC LIMIT 10"
        )
        ranks = list(cur)
        logger.info("Done. Output: %s", args.output)
        logger.info("Rank distribution (top 10): %s", ranks)
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
