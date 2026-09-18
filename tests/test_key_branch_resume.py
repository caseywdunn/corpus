"""The fragment policy reaches both Stage 1 resume gates without PDF work."""

import json

import pytest

from pipeline import chunking, key_context, runner, stages
from tests import test_metadata_resume

corpus = test_metadata_resume.corpus


def test_branch_policy_rechunks_through_discovery_and_stage_resume(corpus, monkeypatch):
    current = key_context.KEY_BRANCH_CONTEXT_POLICY
    with monkeypatch.context() as old:
        old.setattr(key_context, "KEY_BRANCH_CONTEXT_POLICY", "legacy-key-policy")
        old.setattr(chunking, "KEY_BRANCH_CONTEXT_POLICY", "legacy-key-policy")
        corpus.run()
    directory = corpus.hd()
    before = stages._load_pipeline_state(directory)["stages"]
    assert json.loads((directory / "chunks.json").read_text())["key_branch_context_policy"] != current
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
    assert actual["key_branch_context_policy"] == current
    reasons = json.loads((directory / "summary.json").read_text())["processing_summary"]["rerun_reasons"]
    assert "config.chunking.key_branch_context_policy" in reasons["text_chunking"]
    clean = corpus.run(destination=corpus.output.parent / "clean")
    assert actual == json.loads((corpus.hd(destination=clean) / "chunks.json").read_text())
    receipts = (directory / "pipeline_state.json").read_bytes()
    corpus.run()
    assert (directory / "pipeline_state.json").read_bytes() == receipts
