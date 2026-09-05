"""Build and query vectors must belong to the same embedding space."""
import json
from types import SimpleNamespace

import pytest

from mcpsrv.bundle import package
from mcpsrv.indexes import CorpusIndex
from pipeline import embed
from pipeline.embedding_state import read_marker, validate_embedding_index
from pipeline.model_provenance import embedding_producer
from tests.test_embedding_updates import build, cli  # noqa: F401 — pytest fixture


def identity(revision="r1"):
    return {"model": "test-model", "verification": "repository-revision", "revision": revision,
            "recipe": "sentence-transformers-unit-l2-v1"}


def bundle(build):
    build.backend.producer = identity()
    build.run(build.document())
    served = build.root / "_serve"
    package(build.root, served, "test", False, False)
    return served


def test_revision_change_requires_whole_index_rebuild(build, monkeypatch):
    build.backend.producer = identity()
    build.run(build.document())
    build.run(build.document("def"))
    version = build.table.version
    build.backend.producer = identity("r2")
    assert cli(monkeypatch, build, "--resume") == 2
    assert build.table.version == version
    assert cli(monkeypatch, build, "--rebuild", "--resume") == 0
    assert validate_embedding_index(build.root) == ("test-model", 2)


def test_mixed_producers_cannot_publish_over_previous_bundle(build):
    served = bundle(build)
    before = (served / "embedding_producer.json").read_bytes()
    build.backend.producer = identity("r2")
    build.run(build.document("def"))
    with pytest.raises(ValueError, match="producer changed"):
        package(build.root, served, "test", False, False)
    assert (served / "embedding_producer.json").read_bytes() == before


def test_legacy_migration_cannot_leave_unselected_old_rows(build, monkeypatch):
    for sha in ("abc", "def"):
        build.run(build.document(sha))
        path = build.root / "vector_db" / f"{sha}_embedded.done"
        marker = read_marker(path)
        marker.pop("embedding_producer")
        marker["marker_version"] = 2
        path.write_text(json.dumps(marker))
    version = build.table.version
    assert cli(monkeypatch, build, "--resume", "--pdf-hash", "abc") == 2
    assert build.table.version == version
    assert cli(monkeypatch, build, "--resume") == 0
    assert validate_embedding_index(build.root) == ("test-model", 2)


def test_query_pins_build_model_and_revision_without_bundle_writes(build, monkeypatch):
    served = bundle(build)
    before = {p.relative_to(served): (p.stat().st_size, p.stat().st_mtime_ns) for p in served.rglob("*") if p.is_file()}
    calls = []
    def loader(model, **kwargs):
        calls.append((model, kwargs))
        return build.backend
    monkeypatch.setattr("mcpsrv.indexes.get_embedder", loader)
    idx = CorpusIndex(served)
    idx.load()
    assert not calls  # still lazy
    assert idx.get_topic_searcher()[0] is build.backend
    assert calls == [("test-model", {"revision": "r1"})]
    after = {p.relative_to(served): (p.stat().st_size, p.stat().st_mtime_ns) for p in served.rglob("*") if p.is_file()}
    assert before == after


def test_same_dimension_wrong_model_or_loaded_revision_is_refused(build, monkeypatch):
    served = bundle(build)
    idx = CorpusIndex(served, embedding_model="different-same-dimension")
    idx.load()
    monkeypatch.setattr("mcpsrv.indexes.get_embedder", lambda *a, **k: pytest.fail("wrong model loaded"))
    assert idx.get_topic_searcher() == (None, None)
    assert "identity mismatch" in idx._topic_search_status
    idx = CorpusIndex(served)
    idx.load()
    monkeypatch.setattr("mcpsrv.indexes.get_embedder", lambda *a, **k: SimpleNamespace(dim=2, model_name="test-model", producer=identity("wrong")))
    assert idx.get_topic_searcher() == (None, None)
    assert "producer mismatch" in idx._topic_search_status


def test_legacy_bundle_uses_manifest_model_and_warns_about_weaker_proof(build, monkeypatch, caplog):
    served = bundle(build)
    (served / "embedding_producer.json").unlink()
    idx = CorpusIndex(served)
    idx.load()
    calls = []
    def loader(model, **kwargs):
        calls.append((model, kwargs))
        return build.backend
    monkeypatch.setattr("mcpsrv.indexes.get_embedder", loader)
    assert idx.get_topic_searcher()[0] is build.backend
    assert calls == [("test-model", {})]
    assert "Legacy embedding identity" in caplog.text


def test_corrupt_sidecar_does_not_fall_back_to_an_unverified_model(build, monkeypatch):
    served = bundle(build)
    (served / "embedding_producer.json").write_text("{}")
    idx = CorpusIndex(served)
    assert idx.load() == 1
    monkeypatch.setattr("mcpsrv.indexes.get_embedder", lambda *a, **k: pytest.fail("bad identity loaded"))
    assert idx.get_topic_searcher() == (None, None)
    assert "Malformed" in idx._topic_search_status


def test_build_without_all_producer_receipts_refuses_query_embedding(build, monkeypatch):
    build.run(build.document())
    (build.root / "vector_db/abc_embedded.done").unlink()
    idx = CorpusIndex(build.root)
    idx.load()
    monkeypatch.setattr("mcpsrv.indexes.get_embedder", lambda *a, **k: pytest.fail("unverified build loaded"))
    assert idx.get_topic_searcher() == (None, None)
    assert "Missing embedding identity" in idx._topic_search_status


def test_local_model_relocation_checks_content_and_keeps_paths_private(build, monkeypatch, tmp_path):
    original, replica = tmp_path / "original-model", tmp_path / "replica-model"
    for directory in (original, replica):
        directory.mkdir()
        (directory / "weights.bin").write_bytes(b"same model weights")
    build.backend.model_name = str(original)
    build.backend.producer = embedding_producer(str(original))
    build.run(build.document())
    served = build.root / "_serve"
    package(build.root, served, "test", False, False)
    assert str(tmp_path) not in (served / "embedding_producer.json").read_text()
    assert str(tmp_path) not in (served / "bundle_manifest.json").read_text()
    def loader(model, **kwargs):
        assert model == str(replica)
        return SimpleNamespace(dim=2, model_name=model, producer=embedding_producer(model))
    monkeypatch.setattr("mcpsrv.indexes.get_embedder", loader)
    idx = CorpusIndex(served, embedding_model=str(replica))
    idx.load()
    assert idx.get_topic_searcher()[0] is not None
    (replica / "weights.bin").write_bytes(b"different weights")
    idx = CorpusIndex(served, embedding_model=str(replica))
    idx.load()
    assert idx.get_topic_searcher() == (None, None)
