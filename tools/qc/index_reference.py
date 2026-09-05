"""Logical index evidence for release comparison; no models or build writes."""
from collections import defaultdict
import hashlib
import json
import sqlite3


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                    separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _rows(conn, table, *, omit=()):
    # Table names come only from the static catalog below, not operator input.
    columns = [r[1] for r in conn.execute(f'PRAGMA table_info("{table}")') if r[1] not in omit]
    if not columns:
        return None
    rows = []
    selected = ", ".join('"' + name.replace('"', '""') + '"' for name in columns)
    for row in conn.execute(f'SELECT {selected} FROM "{table}"'):
        value = dict(zip(columns, row))
        if "bib_imported_at" in value:
            value["bib_imported_at"] = value["bib_imported_at"] is not None
        for key in columns:
            if key.endswith("_json") and value[key] is not None:
                value[key] = json.loads(value[key])
        rows.append(value)
    fingerprints = sorted(digest(row) for row in rows)
    return {"columns": columns, "row_count": len(rows), "sha256": digest(fingerprints)}


def sqlite_snapshot(path, tables):
    if not path.is_file():
        return None
    # mode=ro avoids creating missing databases and prevents accidental writes.
    conn = sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)
    try:
        return {table: _rows(conn, table, omit=omit) for table, omit in tables.items()}
    finally:
        conn.close()


def vector_snapshot(path):
    if not path.is_dir():
        return None
    import lancedb
    from pipeline.embeddings import lancedb_table_names
    db = lancedb.connect(str(path))
    tables = {}
    for name in sorted(lancedb_table_names(db)):
        table = db.open_table(name)
        by_document = defaultdict(list)
        # Streaming: don't materialize the corpus's full embedding matrix.
        for batch in table.search().limit(None).to_batches():
            for row in batch.to_pylist():
                # Transaction UUIDs differ across equivalent clean/incremental
                # builds. Every actual text, metadata and float value matters.
                row.pop("embedding_generation", None)
                sha = (row.get("metadata") or {}).get("pdf_hash", "<missing>")
                by_document[sha].append(digest(row))
        tables[name] = {
            "schema": str(table.schema.remove_metadata()),
            "documents": {sha: {"row_count": len(rows), "sha256": digest(sorted(rows))}
                          for sha, rows in sorted(by_document.items())},
            "row_count": sum(len(rows) for rows in by_document.values()),
        }
    return tables


def index_snapshot(build):
    from pipeline.embedding_state import embedding_identity
    return {
        "embedding_identity": embedding_identity(build) if (build / "documents").is_dir() else None,
        "bibliography": sqlite_snapshot(build / "biblio_authority.sqlite", {
            "works": ("created_at", "updated_at"),
            "work_documents": ("source_sha256",),
            "work_authors": (), "work_aliases": (), "citations": (),
            "reference_current_sets": ("selected_at",),
            "observation_work": ("mapped_at",), "taxon_work_links": (),
        }),
        "taxon_mentions": sqlite_snapshot(build / "taxon_mentions.sqlite", {
            "taxon_mentions": ("mention_id",),
        }),
        "taxonomy": sqlite_snapshot(build / "taxonomy.sqlite", {"taxa": ("fetched_at",), "names": ()}),
        "vectors": vector_snapshot(build / "vector_db/lancedb"),
    }


def figure_pixels(directory):
    """Hash decoded RGBA pixels, not PNG encoder metadata/compression."""
    from PIL import Image
    result = {}
    if not directory.is_dir():
        return result
    for path in sorted(directory.rglob("*.png")):
        # Legacy query-time crops are not build evidence.
        if "crops" in path.relative_to(directory).parts:
            continue
        with Image.open(path) as image:
            rgba = image.convert("RGBA")
            result[path.relative_to(directory).as_posix()] = {
                "size": list(rgba.size), "rgba_sha256": hashlib.sha256(rgba.tobytes()).hexdigest(),
            }
    return result
