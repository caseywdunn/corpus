"""Canonical-work membership and document-local curation evidence.

``works.corpus_hash`` remains a representative for the frozen wire surface;
it is not the inventory. A license or page directive belongs to a particular
PDF and must never leak to another document sharing that work's DOI.
"""
from __future__ import annotations

import json

DOCUMENT_FIELDS = ("license", "license_url", "license_source", "publishable",
                   "serve", "serve_reason", "ocrlang", "ocrmode", "doclang",
                   "pagemap", "keeppages")


def consumed_metadata(meta):
    """Retain only fields consumed here, not build paths or service receipts."""
    keys = ("title", "year", "journal", "doi", "filename", "authors", *DOCUMENT_FIELDS)
    result = {key: meta[key] for key in keys if key in meta}
    if "authors" in result:
        result["authors"] = [{key: author.get(key, "") for key in ("surname", "forename")}
                             for author in result["authors"]]
    return result


def has_memberships(conn):
    return conn.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='work_documents'").fetchone() is not None


def create_schema(conn):
    conn.execute("""CREATE TABLE IF NOT EXISTS work_documents (
        corpus_hash TEXT PRIMARY KEY,
        work_id TEXT NOT NULL REFERENCES works(work_id),
        metadata_json TEXT NOT NULL,
        source_sha256 TEXT NOT NULL
    )""")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_work_documents_work ON work_documents(work_id)")


def work_map(conn):
    """Include scalar legacy rows until their next normal authority ingest."""
    result = dict(conn.execute("SELECT corpus_hash, work_id FROM works WHERE corpus_hash IS NOT NULL AND in_corpus=1"))
    if has_memberships(conn):
        result.update(conn.execute("SELECT corpus_hash, work_id FROM work_documents"))
    return result


def document_metadata(conn, corpus_hash):
    if not has_memberships(conn):
        return None
    row = conn.execute("SELECT metadata_json FROM work_documents WHERE corpus_hash=?", (corpus_hash,)).fetchone()
    return json.loads(row[0]) if row else None


def find_work(conn, corpus_hash):
    if has_memberships(conn):
        row = conn.execute("SELECT work_id FROM work_documents WHERE corpus_hash=?", (corpus_hash,)).fetchone()
        if row:
            return row[0]
    row = conn.execute("SELECT work_id FROM works WHERE corpus_hash=?", (corpus_hash,)).fetchone()
    return row[0] if row else None


def document_fields(meta):
    from .authority import derive_publishable
    values = {key: meta.get(key, 1 if key == "serve" else None) for key in DOCUMENT_FIELDS}
    publishable, source = derive_publishable(meta.get("license"), meta.get("year"))
    if "publishable" not in meta:
        values["publishable"] = publishable
    if "license_source" not in meta:
        values["license_source"] = source
    return values


def refresh_representative(conn, work_id, *, refresh_header=False):
    """Select a stable, preferably served representative; preserve all members."""
    rows = conn.execute("SELECT corpus_hash, metadata_json FROM work_documents WHERE work_id=? ORDER BY corpus_hash", (work_id,)).fetchall()
    if not rows:
        conn.execute("UPDATE works SET corpus_hash=NULL, in_corpus=0 WHERE work_id=?", (work_id,))
        return
    members = [(sha, json.loads(raw)) for sha, raw in rows]
    sha, meta = min(members, key=lambda item: (item[1].get("serve", 1) == 0, item[0]))
    fields = document_fields(meta)
    conn.execute(f"UPDATE works SET corpus_hash=?, in_corpus=1, {', '.join(k+'=?' for k in fields)} WHERE work_id=?",
                 (sha, *fields.values(), work_id))
    curated = conn.execute("SELECT bib_imported_at FROM works WHERE work_id=?", (work_id,)).fetchone()
    if refresh_header and not (curated and curated[0] is not None):
        from .authority import insert_authors, normalize_doi
        conn.execute("UPDATE works SET title=?, year=?, journal=?, doi=?, source='corpus_paper' WHERE work_id=?",
                     (meta.get("title") or "", meta.get("year"), meta.get("journal") or "", normalize_doi(meta.get("doi") or "") or None, work_id))
        conn.execute("DELETE FROM work_authors WHERE work_id=?", (work_id,))
        insert_authors(conn, work_id, [(a.get("surname", ""), a.get("forename", ""))
                                     for a in meta.get("authors", []) if a.get("surname")])


def move_memberships(conn, old_work_id, new_work_id):
    if has_memberships(conn):
        conn.execute("UPDATE work_documents SET work_id=? WHERE work_id=?", (new_work_id, old_work_id))
        if conn.execute("SELECT 1 FROM work_documents WHERE work_id=?", (new_work_id,)).fetchone():
            refresh_representative(conn, new_work_id)


def update_document_fields(conn, corpus_hash, changes):
    meta = document_metadata(conn, corpus_hash)
    if meta is None:
        return False
    meta.update({key: value for key, value in changes.items() if key in (*DOCUMENT_FIELDS, "year")})
    conn.execute("UPDATE work_documents SET metadata_json=? WHERE corpus_hash=?",
                 (json.dumps(meta, sort_keys=True, ensure_ascii=False), corpus_hash))
    work_id = conn.execute("SELECT work_id FROM work_documents WHERE corpus_hash=?", (corpus_hash,)).fetchone()[0]
    refresh_representative(conn, work_id)
    return True
