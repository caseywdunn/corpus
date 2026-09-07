"""Vision identity changes invalidate only their actual consumers (#174)."""
import json
import os
from types import SimpleNamespace

from pipeline import model_provenance, stages
from pipeline.build_inputs import config_fingerprints
from pipeline.vision import LocalVLMBackend
from tests.test_config_updates import _install_figure_stubs
from tests.test_metadata_resume import corpus  # noqa: F401 — pytest fixture


def test_default_model_is_explicit_and_cached_revision_changes_only_figure_consumers(monkeypatch):
    revision = ["old"]
    monkeypatch.setattr(model_provenance, "_cached_revision", lambda _: revision[0])
    before = config_fingerprints({}, panel_mode="vision-local")
    assert before["figure_materialization"]["figures.model"] == model_provenance.DEFAULT_VISION_MODELS["vision-local"]
    revision[0] = "new"
    after = config_fingerprints({}, panel_mode="vision-local")
    assert {key for key in before if before[key] != after[key]} == {"figure_materialization", "figure_crossref"}


def test_local_weights_use_bytes_not_size_or_mtime(tmp_path):
    model = tmp_path / "model"
    model.mkdir()
    weights = model / "model.safetensors"
    weights.write_bytes(b"first weights")
    before = model_provenance.vision_producer("vision-local", str(model))
    stat = weights.stat()
    weights.write_bytes(b"other weights")
    os.utime(weights, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    after = model_provenance.vision_producer("vision-local", str(model))
    assert before["verification"] == "local-file-content"
    assert before["model"] == "model"
    assert str(tmp_path) not in json.dumps(before)
    assert before["sha256"] != after["sha256"]


def test_loaded_revision_wins_over_stale_cache_without_loading_models(monkeypatch):
    monkeypatch.setattr(model_provenance, "_cached_revision", lambda _: "cached-old")
    backend = LocalVLMBackend.__new__(LocalVLMBackend)
    backend._model_id = "example/model"
    backend._model = SimpleNamespace(config=SimpleNamespace(_commit_hash="actually-loaded"))
    backend._max_new_tokens, backend._max_pixels, backend._min_pixels = 1024, 1003520, 3136
    resolved = backend.producer
    assert resolved["revision"] == "actually-loaded"
    fps = config_fingerprints({}, panel_mode="vision-local", vision_model="example/model",
                              resolved_vision_producer=resolved)
    assert fps["figure_materialization"]["figures.producer"] == resolved


def test_remote_identity_and_cache_miss_are_not_claimed_weight_verification(monkeypatch):
    def no_lookup(_):
        raise AssertionError("remote producer must not inspect HuggingFace")
    monkeypatch.setattr(model_provenance, "_cached_revision", no_lookup)
    assert model_provenance.vision_producer("vision-claude")["verification"] == "remote-model-id-only"
    monkeypatch.setattr(model_provenance, "_cached_revision", lambda _: None)
    assert model_provenance.vision_producer("vision-local")["verification"] == "unverified-cache-miss"
    old = config_fingerprints({}, panel_mode="vision-claude")
    new = config_fingerprints({"figures": {"producer_id": "new deployment"}}, panel_mode="vision-claude")
    assert old["figure_materialization"] != new["figure_materialization"]
    assert old["docling_extraction"] == new["docling_extraction"]


def test_revision_change_runs_through_both_resume_gates(corpus, monkeypatch):
    _install_figure_stubs(monkeypatch)
    revision = ["old"]
    monkeypatch.setattr(model_provenance, "_cached_revision", lambda _: revision[0])
    corpus.run("--figure-panels", "vision-local")
    hd = corpus.hd()
    old = stages._load_pipeline_state(hd)["stages"]
    corpus.run("--figure-panels", "vision-local")
    assert stages._load_pipeline_state(hd)["stages"] == old
    revision[0] = "new"
    corpus.run("--figure-panels", "vision-local")
    after = stages._load_pipeline_state(hd)["stages"]
    assert after["figure_materialization"] != old["figure_materialization"]
    for stage in ("scan_detection", "pdf_preparation", "docling_extraction", "metadata_extraction", "text_chunking"):
        assert old[stage] == after[stage]
    receipt = after["figure_materialization"]["input_fingerprint"]["config"]["figures.producer"]
    assert receipt["revision"] == "new"
    assert json.loads((hd / "figures.json").read_text())["figures"][0]["filename"] == "split.png"
