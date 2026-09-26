"""The fragment policy reaches both Stage 1 resume gates without PDF work."""

import json

import pytest

from pipeline import chunking, key_context, runner, stages, treatment_context
from tests import test_metadata_resume

corpus = test_metadata_resume.corpus


@pytest.mark.parametrize("module,constant,field", [
    (key_context, "KEY_BRANCH_CONTEXT_POLICY", "key_branch_context_policy"),
    (treatment_context, "TREATMENT_CONTEXT_POLICY", "treatment_context_policy"),
])
def test_context_policy_rechunks_through_discovery_and_stage_resume(corpus, monkeypatch, module, constant, field):
    current = getattr(module, constant)
    with monkeypatch.context() as old:
        old.setattr(module, constant, "legacy-context-policy")
        old.setattr(chunking, constant, "legacy-context-policy")
        corpus.run()
    directory = corpus.hd()
    before = stages._load_pipeline_state(directory)["stages"]
    assert json.loads((directory / "chunks.json").read_text())[field] != current
    with monkeypatch.context() as patch:
        for name in ("detect_scan_type", "prepare_pdf", "extract_docling_content",
                     "extract_metadata", "_pass3a_annotate_rois"):
            patch.setattr(runner, name, lambda *a, **k: pytest.fail("unrelated work reran"))
        corpus.run()
    after = stages._load_pipeline_state(directory)["stages"]
    for stage in ("scan_detection", "pdf_preparation", "metadata_extraction",
                  "docling_extraction", "figure_materialization"):
        assert before[stage] == after[stage]
    assert before["text_chunking"] != after["text_chunking"]
    actual = json.loads((directory / "chunks.json").read_text())
    assert actual[field] == current
    reasons = json.loads((directory / "summary.json").read_text())["processing_summary"]["rerun_reasons"]
    assert f"config.chunking.{field}" in reasons["text_chunking"]
    clean = corpus.run(destination=corpus.output.parent / "clean")
    assert actual == json.loads((corpus.hd(destination=clean) / "chunks.json").read_text())
    receipts = (directory / "pipeline_state.json").read_bytes()
    corpus.run()
    assert (directory / "pipeline_state.json").read_bytes() == receipts
