"""Bibliographic fields and deterministic source authority (#296, #301).

Document directives (including physical ``keeppages``) are deliberately
separate from publication locators. Source rows retain the supplied values;
canonical works are a materialized projection, never the source of evidence.
"""
from __future__ import annotations

import json
import time

LOCATOR_FIELDS = ("volume", "number", "pages", "eid", "articleno", "chapter", "booktitle", "publisher", "edition", "series")
BIBLIOGRAPHIC_FIELDS = ("title", "year", "journal", "doi", *LOCATOR_FIELDS)
SOURCE_FIELDS = (*BIBLIOGRAPHIC_FIELDS, "authors", "bib_key")


def create_schema(conn):
    conn.execute("""CREATE TABLE IF NOT EXISTS work_bib_sources (
        work_id TEXT NOT NULL REFERENCES works(work_id),
        source_id TEXT NOT NULL,
        origin TEXT NOT NULL,
        metadata_json TEXT NOT NULL,
        PRIMARY KEY (work_id, source_id)
    )""")


def record_source(conn, work_id, source_id, metadata, *, origin):
    """Record a current authoritative input; extracted observations live elsewhere."""
    payload = {key: metadata[key] for key in SOURCE_FIELDS if key in metadata}
    conn.execute("""INSERT INTO work_bib_sources VALUES (?,?,?,?)
        ON CONFLICT(work_id, source_id) DO UPDATE SET
        origin=excluded.origin, metadata_json=excluded.metadata_json""",
        (work_id, source_id, origin, json.dumps(payload, sort_keys=True, ensure_ascii=False)))


def sources(conn, work_id):
    rows = [(source_id, origin, json.loads(raw)) for source_id, origin, raw in conn.execute(
        "SELECT source_id,origin,metadata_json FROM work_bib_sources WHERE work_id=? ORDER BY source_id",
        (work_id,))]
    # An explicit re-import of an entry supersedes the build-time copy of
    # that same entry. Independent entries remain visible as conflicts.
    imported = {meta.get("bib_key") for _, origin, meta in rows if origin == "import" and meta.get("bib_key")}
    return sorted((row for row in rows if row[1] == "import" or row[2].get("bib_key") not in imported),
                  key=lambda row: (row[1] != "import", row[0]))


def conflicts(conn, work_id):
    rows = sources(conn, work_id)
    result = []
    for field in (*BIBLIOGRAPHIC_FIELDS, "authors"):
        values = [(source_id, meta[field]) for source_id, _, meta in rows if meta.get(field) not in (None, "", [])]
        if len({json.dumps(value, sort_keys=True, ensure_ascii=False) for _, value in values}) > 1:
            result.append({"field": field, "policy": "explicit_import_then_source_id",
                           "selected_source": values[0][0], "sources": [
                               {"source_id": source_id, "value": value} for source_id, value in values]})
    return result


def materialize(conn, work_id):
    """Apply supplied curated fields and the complete author list, in stable order."""
    rows = sources(conn, work_id)
    if not rows:
        conn.execute("UPDATE works SET bib_imported_at=NULL,bib_key=NULL,bib_source=NULL WHERE work_id=? AND bib_source='metadata'", (work_id,))
        return False
    from .authority import insert_authors, normalize_doi
    chosen = {}
    for _, _, meta in rows:
        for key, value in meta.items():
            if key not in chosen and value not in (None, "", []):
                chosen[key] = value
    values = {key: chosen.get(key) for key in BIBLIOGRAPHIC_FIELDS
              if any(key in meta for _, _, meta in rows)}
    if "doi" in values:
        values["doi"] = normalize_doi(values["doi"]) if values["doi"] else None
    values.update(bib_key=chosen.get("bib_key"), bib_source=rows[0][1])
    conn.execute(f"UPDATE works SET {', '.join(key+'=?' for key in values)}, bib_imported_at=COALESCE(bib_imported_at,?), updated_at=? WHERE work_id=?",
                 (*values.values(), time.time(), time.time(), work_id))
    if "authors" in chosen:
        conn.execute("DELETE FROM work_authors WHERE work_id=?", (work_id,))
        insert_authors(conn, work_id, [(a.get("surname", ""), a.get("forename", "")) for a in chosen["authors"]])
    return True


def capture_legacy(conn, work_id):
    """Carry pre-source-table explicit imports through reconciliation."""
    row = conn.execute("SELECT * FROM works WHERE work_id=?", (work_id,)).fetchone()
    names = [column[1] for column in conn.execute("PRAGMA table_info(works)")]
    meta = dict(zip(names, row)) if row else {}
    if meta.get("bib_imported_at") is None or sources(conn, work_id):
        return
    meta["authors"] = [{"surname": s, "forename": f or ""} for s, f in conn.execute(
        "SELECT surname,forename FROM work_authors WHERE work_id=? ORDER BY position", (work_id,))]
    record_source(conn, work_id, "legacy:" + work_id, meta, origin="import")


def move_sources(conn, old_id, new_id):
    for source_id, origin, raw in conn.execute(
            "SELECT source_id,origin,metadata_json FROM work_bib_sources WHERE work_id=?", (old_id,)).fetchall():
        record_source(conn, new_id, source_id, json.loads(raw), origin=origin)
    conn.execute("DELETE FROM work_bib_sources WHERE work_id=?", (old_id,))
    materialize(conn, new_id)
