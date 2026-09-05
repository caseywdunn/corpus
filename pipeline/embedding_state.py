"""Build-side embedding inputs and completion evidence (#271).

No model is loaded here. A marker is evidence only when its input fingerprint
and generation agree with the current artifacts and committed table rows.
See dev_docs/OVERVIEW.md, "Stage 2: Embedding".
"""

from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from .version import __version__

MARKER_VERSION = 3


def load_document_data(hash_dir: Path) -> dict:
    """Require valid chunks; missing optional metadata uses schema defaults.

    A malformed or missing chunk artifact is never an empty document: treating
    it as one would delete the last good rows after an interrupted extraction.
    """
    data = {}
    for name in ("summary", "metadata", "chunks"):
        path = hash_dir / f"{name}.json"
        if name != "chunks" and not path.exists():
            data[name] = {}
            continue
        value = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(value, dict):
            raise ValueError(f"{path}: expected an object")
        data[name] = value
    chunks = data["chunks"].get("chunks")
    if not isinstance(chunks, list) or any(not isinstance(c, dict) for c in chunks):
        raise ValueError(f"{hash_dir}/chunks.json: expected a chunks list")
    return data


def record_payloads(pdf_hash: str, doc: dict) -> list[dict]:
    """The exact text/metadata fields Stage 2 stores, without vectors."""
    metadata = doc["metadata"]
    paths = doc["summary"].get("relative_paths") or []
    authors = []
    for author in metadata.get("authors") or []:
        name = f"{(author.get('forename') or '').strip()} {(author.get('surname') or '').strip()}".strip()
        if name:
            authors.append(name)
    common = {
        "pdf_hash": pdf_hash,
        "filename": paths[0] if paths else "",
        "relative_paths": paths,
        "title": metadata.get("title") or "",
        "authors": authors,
        "year": metadata.get("year"),
        "journal": metadata.get("journal") or "",
        "doi": metadata.get("doi") or "",
        "total_pages": (doc["chunks"].get("metadata") or {}).get("total_pages"),
    }
    return [{
        "text": chunk.get("text") or "",
        "metadata": {
            **common,
            "chunk_id": chunk.get("chunk_id", ""),
            "section_class": chunk.get("section_class"),
            "headings": chunk.get("headings") or [],
        },
    } for chunk in doc["chunks"]["chunks"]]


def input_fingerprint(pdf_hash: str, doc: dict) -> str:
    payload = json.dumps(record_payloads(pdf_hash, doc), sort_keys=True,
                         ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def row_census(table) -> dict[str, Counter]:
    """One projected scan; never load vectors or text to verify completion."""
    columns = ["metadata.pdf_hash"]
    has_generation = "embedding_generation" in table.schema.names
    if has_generation:
        columns.append("embedding_generation")
    census = defaultdict(Counter)
    for batch in table.search().select(columns).limit(None).to_batches():
        for row in batch.to_pylist():
            census[row["metadata.pdf_hash"]][row.get("embedding_generation")] += 1
    return dict(census)


def read_marker(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def marker_problem(hash_dir: Path, marker: dict, census: dict, *,
                   model: str | None = None, dim: int | None = None,
                   table_name: str = "document_chunks",
                   fingerprint: str | None = None,
                   producer: dict | None = None) -> str | None:
    """Return a reason to rebuild, or None for verified current rows."""
    if (marker.get("marker_version") != MARKER_VERSION
            or marker.get("pipeline_version") != __version__
            or marker.get("status") != "completed"
            or marker.get("embedding_table") != table_name
            or marker.get("pdf_hash") != hash_dir.name):
        return "missing, legacy or incompatible completion marker"
    count = marker.get("chunks_count")
    generation = marker.get("embedding_generation")
    if (type(count) is not int or count < 0
            or not isinstance(generation, str) or not generation
            or not isinstance(marker.get("embedding_model"), str)
            or not marker["embedding_model"]
            or type(marker.get("embedding_dim")) is not int
            or marker["embedding_dim"] <= 0
            or not isinstance(marker.get("embedding_producer"), dict)
            or not marker["embedding_producer"].get("verification")):
        return "invalid completion marker"
    from .model_provenance import same_embedding_space
    if producer is not None and not same_embedding_space(marker["embedding_producer"], producer):
        return "embedding producer changed"
    if ((model is not None and marker["embedding_model"] != model)
            or (dim is not None and marker["embedding_dim"] != dim)):
        return "embedding model or dimension changed"
    if fingerprint is None:
        fingerprint = input_fingerprint(hash_dir.name, load_document_data(hash_dir))
    if marker.get("input_fingerprint") != fingerprint:
        return "embedding inputs changed"
    expected = Counter({generation: count}) if count else Counter()
    if census.get(hash_dir.name, Counter()) != expected:
        return "committed rows do not match the completion marker"
    return None


def validate_embedding_index(output_dir: Path) -> tuple[str | None, int | None]:
    """Validate the whole build index before bundling, not an arbitrary marker.

    A build with neither an index nor markers may be served without semantic
    search. A partially embedded build must not masquerade as that mode.
    Orphan markers alone are harmless; orphan rows are a failed prune.
    """
    import lancedb
    from .embeddings import lancedb_table_names

    vdb = output_dir / "vector_db"
    dirs = sorted(p for p in (output_dir / "documents").iterdir() if p.is_dir())
    index_path = vdb / "lancedb"
    if not index_path.exists():
        if any((vdb / f"{p.name}_embedded.done").exists() for p in dirs):
            raise ValueError("Embedding markers exist but the vector index is missing")
        return None, None
    db = lancedb.connect(str(index_path))
    if "document_chunks" not in lancedb_table_names(db):
        raise ValueError("Vector index is missing its document_chunks table")
    table = db.open_table("document_chunks")
    census = row_census(table)
    orphan_rows = set(census) - {p.name for p in dirs}
    if orphan_rows:
        raise ValueError(f"Vector index has orphan document rows: {sorted(orphan_rows)}")
    dim = table.schema.field("vector").type.list_size
    model = None
    producer = None
    for hash_dir in dirs:
        marker = read_marker(vdb / f"{hash_dir.name}_embedded.done")
        problem = marker_problem(hash_dir, marker, census, model=model, dim=dim, producer=producer)
        if problem:
            raise ValueError(f"Incomplete vector index: {hash_dir.name}: {problem}. "
                             "Run embedding before bundling.")
        model = marker["embedding_model"]
        producer = marker["embedding_producer"]
    return (model, dim) if model is not None else (None, None)


def embedding_identity(output_dir: Path, *, require_all=False):
    """Read a consistent identity across current document receipts, not a sample.

    Bundle callers additionally validate committed rows and all inputs with
    validate_embedding_index. This helper alone is not a completeness claim.
    """
    from .model_provenance import same_embedding_space
    result = None
    for hd in sorted((output_dir / "documents").iterdir()):
        if not hd.is_dir():
            continue
        marker = read_marker(output_dir / "vector_db" / f"{hd.name}_embedded.done")
        model, producer, dim = (marker.get(key) for key in ("embedding_model", "embedding_producer", "embedding_dim"))
        if not model:
            if require_all and (output_dir / "vector_db/lancedb").is_dir():
                raise ValueError(f"Missing embedding identity receipt for {hd.name}")
            continue  # Legacy build without identity; caller reports weaker proof.
        current = {"schema_version": 1, "model": model, "dimension": dim, "producer": producer}
        if result is not None and (model != result["model"] or dim != result["dimension"]
                                  or (producer != result["producer"] and not same_embedding_space(producer, result["producer"]))):
            raise ValueError("Embedding index contains inconsistent producer receipts")
        result = current
    if result and result["producer"] and result["producer"].get("verification") == "local-file-content":
        result["model"] = result["producer"]["model"]
    return result
