"""Embedding update contract on real LanceDB, with a deterministic tiny model."""

import json
from pathlib import Path
from types import SimpleNamespace

import lancedb
import pytest

from mcpsrv.bundle import package
from pipeline import embed
from pipeline.embedding_state import (
    input_fingerprint, load_document_data, marker_problem, read_marker,
    row_census, validate_embedding_index,
)
from pipeline.status import render_artifacts


def write_json(path, value):
    path.write_text(json.dumps(value))


@pytest.fixture
def build(tmp_path):
    root = tmp_path / "build"
    (root / "vector_db").mkdir(parents=True)
    model = embed.make_chunk_model(2)
    db = lancedb.connect(str(root / "vector_db" / "lancedb"))
    table = db.create_table("document_chunks", schema=model.to_arrow_schema())
    backend = SimpleNamespace(model_name="test-model", dim=2,
                              embed=lambda texts: [[float(len(t)), 1.0] for t in texts])

    def document(pdf_hash="abc", texts=("first", "second")):
        hd = root / "documents" / pdf_hash
        hd.mkdir(parents=True, exist_ok=True)
        write_json(hd / "summary.json", {"relative_paths": [f"{pdf_hash}.pdf"]})
        write_json(hd / "metadata.json", {"title": f"Paper {pdf_hash}"})
        write_json(hd / "chunks.json", {"chunks": [
            {"chunk_id": f"c{i}", "text": text} for i, text in enumerate(texts)
        ]})
        return hd

    def run(hd):
        return embed.embed_document(hd, table, model, backend)

    def problem(hd):
        marker = read_marker(root / "vector_db" / f"{hd.name}_embedded.done")
        return marker_problem(hd, marker, row_census(table), model=backend.model_name, dim=2)

    return SimpleNamespace(root=root, table=table, model=model, backend=backend,
                           document=document, run=run, problem=problem)


def logical_rows(table):
    rows = table.to_arrow().to_pylist()
    for row in rows:
        row.pop("embedding_generation", None)
    return sorted(rows, key=lambda r: (r["metadata"]["pdf_hash"], r["metadata"]["chunk_id"]))


def cli(monkeypatch, build, *args):
    monkeypatch.setattr(embed, "get_embedder", lambda *a, **kw: build.backend)
    monkeypatch.setattr("sys.argv", ["embed", str(build.root), *args])
    result = embed.main()
    # The CLI uses another handle; LanceDB defaults to snapshot reads.
    build.table.checkout_latest()
    return result


def test_replacement_is_one_commit_and_equals_clean_build(build, tmp_path):
    a, b = build.document(), build.document("def")
    build.run(a)
    build.run(b)
    others = [r for r in build.table.to_arrow().to_pylist() if r["metadata"]["pdf_hash"] == "def"]
    build.document(texts=("revised and shorter",))
    assert build.problem(a) == "embedding inputs changed"
    before = build.table.version
    assert build.run(a) == 1
    assert build.table.version == before + 1
    assert build.problem(a) is None
    assert [r for r in build.table.to_arrow().to_pylist() if r["metadata"]["pdf_hash"] == "def"] == others

    clean = lancedb.connect(str(tmp_path / "clean")).create_table(
        "document_chunks", schema=build.model.to_arrow_schema())
    for hd in (a, b):
        embed.embed_document(hd, clean, build.model, build.backend,
                             marker_file=tmp_path / f"{hd.name}.done")
    assert logical_rows(build.table) == logical_rows(clean)


@pytest.mark.parametrize("file,key,value", [
    ("metadata", "title", "Changed title"),
    ("metadata", "authors", [{"forename": "Ada", "surname": "Lovelace"}]),
    ("metadata", "year", 2020),
    ("metadata", "doi", "10.1234/test"),
    ("summary", "relative_paths", ["renamed.pdf", "alias.pdf"]),
    ("chunks", "metadata", {"total_pages": 20}),
])
def test_consumed_metadata_changes_invalidate_resume(build, monkeypatch, file, key, value):
    hd = build.document()
    build.run(hd)
    path = hd / f"{file}.json"
    data = json.loads(path.read_text())
    data[key] = value
    write_json(path, data)
    assert build.problem(hd) == "embedding inputs changed"
    assert cli(monkeypatch, build, "--resume") == 0
    assert build.problem(hd) is None
    assert build.table.count_rows() == 2


def test_noop_and_unconsumed_fields_do_not_reembed(build, monkeypatch):
    hd = build.document()
    build.run(hd)
    fp = input_fingerprint(hd.name, load_document_data(hd))
    path = hd / "summary.json"
    data = json.loads(path.read_text())
    data["processing_summary"] = {"timestamp": "later"}
    write_json(path, data)
    assert input_fingerprint(hd.name, load_document_data(hd)) == fp
    before = build.table.version
    build.backend.embed = lambda _: pytest.fail("unchanged document was re-embedded")
    assert cli(monkeypatch, build, "--resume") == 0
    assert build.table.version == before


def test_empty_chunks_remove_old_rows_and_get_real_marker(build):
    hd = build.document()
    build.run(hd)
    build.document(texts=())
    build.backend.embed = lambda _: pytest.fail("empty input called model")
    assert build.run(hd) == 0
    assert build.table.count_rows() == 0
    assert build.problem(hd) is None
    assert validate_embedding_index(build.root) == ("test-model", 2)


@pytest.mark.parametrize("content", [None, "{", "{}", '[]', '{"chunks": null}', '{"chunks": [null]}'])
def test_missing_or_malformed_chunks_preserve_old_rows(build, content):
    hd = build.document()
    build.run(hd)
    before = build.table.version
    path = hd / "chunks.json"
    if content is None:
        path.unlink()
    else:
        path.write_text(content)
    with pytest.raises((ValueError, FileNotFoundError)):
        build.run(hd)
    assert build.table.version == before
    assert build.table.count_rows() == 2


def test_failed_backend_and_short_vector_batch_preserve_old_rows(build):
    hd = build.document()
    build.run(hd)
    before = build.table.version
    build.document(texts=("changed", "again"))
    build.backend.embed = lambda _: []
    with pytest.raises(embed.EmbeddingError):
        build.run(hd)
    assert build.table.version == before
    assert build.problem(hd) == "embedding inputs changed"


def test_input_changes_during_embedding_do_not_commit(build):
    hd = build.document()
    build.run(hd)
    before = build.table.version

    def concurrent_edit(texts):
        build.document(texts=("changed under the embedder",))
        return [[1.0, 0.0] for _ in texts]

    build.backend.embed = concurrent_edit
    with pytest.raises(embed.EmbeddingError, match="inputs changed"):
        build.run(hd)
    assert build.table.version == before


def test_failed_transaction_keeps_old_rows_and_marker(build, monkeypatch):
    hd = build.document()
    build.run(hd)
    before = build.table.version
    marker = build.root / "vector_db" / "abc_embedded.done"
    previous = marker.read_bytes()
    build.document(texts=("replacement",))

    def fail(*args, **kwargs):
        raise RuntimeError("transaction failed")

    monkeypatch.setattr(build.table, "merge_insert", fail)
    with pytest.raises(RuntimeError, match="transaction failed"):
        build.run(hd)
    assert build.table.version == before
    assert marker.read_bytes() == previous


def test_failure_after_commit_cannot_pass_resume(build, monkeypatch):
    hd = build.document()
    build.run(hd)
    marker = build.root / "vector_db" / "abc_embedded.done"
    previous = marker.read_bytes()
    # Same inputs AND count: only the committed generation detects this gap.
    with monkeypatch.context() as patch:
        patch.setattr(embed.os, "replace", lambda *a: (_ for _ in ()).throw(OSError("interrupted")))
        with pytest.raises(OSError, match="interrupted"):
            build.run(hd)
    assert marker.read_bytes() == previous
    assert build.problem(hd) == "committed rows do not match the completion marker"
    assert sorted(p.name for p in marker.parent.iterdir()) == ["abc_embedded.done", "lancedb"]
    assert cli(monkeypatch, build, "--resume") == 0
    assert build.problem(hd) is None
    assert build.table.count_rows() == 2


def test_legacy_duplicate_rows_replaced_without_losing_other_papers(build, tmp_path):
    hd = build.document()
    build.run(hd)
    old_rows = logical_rows(build.table)
    old_rows.append({**old_rows[0], "metadata": {**old_rows[0]["metadata"], "pdf_hash": "def"}})
    legacy = lancedb.connect(str(tmp_path / "legacy")).create_table("document_chunks", old_rows + old_rows[:1])
    embed.embed_document(hd, legacy, build.model, build.backend)
    assert legacy.count_rows() == 3
    assert legacy.count_rows("metadata.pdf_hash = 'def'") == 1
    assert marker_problem(hd, read_marker(build.root / "vector_db" / "abc_embedded.done"), row_census(legacy)) is None


@pytest.mark.parametrize("damage", ["legacy", "missing", "rows", "version", "model"])
def test_bad_completion_evidence_is_rejected(build, damage):
    hd = build.document()
    build.run(hd)
    path = build.root / "vector_db" / "abc_embedded.done"
    marker = read_marker(path)
    if damage == "legacy":
        write_json(path, {"status": "completed"})
    elif damage == "missing":
        path.unlink()
    elif damage == "rows":
        build.table.delete("metadata.pdf_hash = 'abc'")
    elif damage == "version":
        write_json(path, {**marker, "pipeline_version": "old"})
    else:
        write_json(path, {**marker, "embedding_model": "different"})
    assert build.problem(hd) is not None


def test_same_dim_model_switch_needs_rebuild(build, monkeypatch):
    hd = build.document()
    build.run(hd)
    before = build.table.version
    build.backend.model_name = "different-same-dim"
    assert cli(monkeypatch, build, "--resume") == 2
    assert build.table.version == before
    assert cli(monkeypatch, build, "--rebuild", "--resume") == 0
    assert validate_embedding_index(build.root) == ("different-same-dim", 2)


def test_selective_rebuild_rejected_before_drop(build, monkeypatch):
    build.run(build.document())
    before = build.table.version
    with pytest.raises(SystemExit):
        cli(monkeypatch, build, "--rebuild", "--pdf-hash", "abc")
    assert build.table.version == before


def test_dry_run_uses_real_evidence_without_loading_model(build, monkeypatch, caplog):
    hd = build.document()
    build.run(hd)
    before = build.table.version
    monkeypatch.setattr(embed, "get_embedder", lambda *a, **k: pytest.fail("dry-run loaded model"))
    monkeypatch.setattr("sys.argv", ["embed", str(build.root), "--dry-run", "--resume", "--model", "test-model"])
    with caplog.at_level("INFO"):
        assert embed.main() == 0
        assert "would embed 0 hash(es); 1 verified current" in caplog.text
        caplog.clear()
        build.document(texts=("changed",))
        assert embed.main() == 0
        assert "would embed 1 hash(es); 0 verified current" in caplog.text
    assert build.table.version == before


@pytest.mark.parametrize("damage", ["second_marker", "second_inputs", "mixed_model", "orphan_rows"])
def test_bundle_checks_whole_index_before_touching_destination(build, tmp_path, damage):
    a, b = build.document(), build.document("def")
    build.run(a)
    build.run(b)
    marker = build.root / "vector_db" / "def_embedded.done"
    if damage == "second_marker":
        marker.unlink()
    elif damage == "second_inputs":
        build.document("def", texts=("changed",))
    elif damage == "mixed_model":
        write_json(marker, {**read_marker(marker), "embedding_model": "other"})
    else:
        row = build.table.to_arrow().to_pylist()[0]
        row["metadata"]["pdf_hash"] = "orphan"
        build.table.add([row])
    dst = tmp_path / "serve"
    dst.mkdir()
    sentinel = dst / "bundle_manifest.json"
    sentinel.write_text("previous valid bundle")
    with pytest.raises(ValueError):
        package(build.root, dst, "test", False, False)
    assert sentinel.read_text() == "previous valid bundle"
    assert list(dst.iterdir()) == [sentinel]
    assert "not verified current" in render_artifacts(build.root)


def test_orphan_marker_does_not_set_bundle_model(build):
    build.run(build.document())
    write_json(build.root / "vector_db" / "000_embedded.done", {"embedding_model": "wrong"})
    assert validate_embedding_index(build.root) == ("test-model", 2)
    assert "verified current: test-model" in render_artifacts(build.root)


def test_stage1_does_not_publish_embedding_markers():
    root = Path(embed.__file__).parent
    for name in ("main.py", "chunking.py"):
        source = (root / name).read_text()
        assert "_embedded.done" not in source
        assert "ingest_to_vector_db" not in source


def test_missing_table_rejects_markers_and_resume_recreates_rows(build, monkeypatch):
    hd = build.document()
    build.run(hd)
    db = lancedb.connect(str(build.root / "vector_db" / "lancedb"))
    db.drop_table("document_chunks")
    with pytest.raises(ValueError, match="missing its document_chunks table"):
        validate_embedding_index(build.root)
    assert cli(monkeypatch, build, "--resume") == 0
    assert validate_embedding_index(build.root) == ("test-model", 2)


def test_missing_index_rejects_existing_receipts(build):
    build.run(build.document())
    index_path = build.root / "vector_db" / "lancedb"
    index_path.rename(build.root / "saved-index")
    with pytest.raises(ValueError, match="vector index is missing"):
        validate_embedding_index(build.root)


def test_index_free_bundle_supported_but_cannot_keep_old_destination_index(tmp_path):
    root = tmp_path / "build"
    hd = root / "documents" / "abc"
    hd.mkdir(parents=True)
    write_json(hd / "chunks.json", {"chunks": [{"text": "keyword-only"}]})
    assert validate_embedding_index(root) == (None, None)
    dst = tmp_path / "serve"
    manifest = package(root, dst, "test", False, False)
    assert manifest["embedding_model"] is None
    (dst / "vector_db" / "lancedb").mkdir(parents=True)
    sentinel = dst / "bundle_manifest.json"
    previous = sentinel.read_bytes()
    with pytest.raises(ValueError, match="fresh served directory"):
        package(root, dst, "test", False, False)
    assert sentinel.read_bytes() == previous


def test_custom_table_receipt_does_not_prove_default_table_complete(build, tmp_path):
    hd = build.document(texts=())
    custom = lancedb.connect(str(tmp_path / "custom")).create_table(
        "another_table", schema=build.model.to_arrow_schema())
    embed.embed_document(hd, custom, build.model, build.backend)
    assert build.problem(hd) is not None
