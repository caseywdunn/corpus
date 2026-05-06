"""Unit tests for per-stage resume helpers backed by pipeline_state.json.

The resume signal is a structured completion record (pipeline_version +
input_fingerprint + completed_at) — not file existence. These tests
exercise that contract directly.
"""
from __future__ import annotations

import json
from pathlib import Path

from pipeline import PIPELINE_VERSION
from process_corpus import (
    _all_stage_artifacts_complete,
    _load_pipeline_state,
    _record_stage_completion,
    _should_run_stage,
    _stage,
    _stage_recorded_complete,
)


# ---------------------------------------------------------------------------
# _record_stage_completion / _load_pipeline_state
# ---------------------------------------------------------------------------


def test_record_creates_pipeline_state_file(tmp_path):
    _record_stage_completion(tmp_path, "scan_detection")
    assert (tmp_path / "pipeline_state.json").exists()


def test_record_captures_pipeline_version_and_fingerprint(tmp_path):
    fp = {"taxonomy": {"sha256": "abc"}}
    _record_stage_completion(tmp_path, "taxa_and_lexicon_extraction",
                             input_fingerprint=fp)
    state = json.loads((tmp_path / "pipeline_state.json").read_text())
    rec = state["stages"]["taxa_and_lexicon_extraction"]
    assert rec["pipeline_version"] == PIPELINE_VERSION
    assert rec["input_fingerprint"] == fp
    assert rec["completed_at"]
    assert state["pipeline_version_latest"] == PIPELINE_VERSION


def test_record_overwrites_prior_completion(tmp_path):
    _record_stage_completion(tmp_path, "docling_extraction")
    first = _load_pipeline_state(tmp_path)["stages"]["docling_extraction"]
    _record_stage_completion(tmp_path, "docling_extraction")
    second = _load_pipeline_state(tmp_path)["stages"]["docling_extraction"]
    # Same stage name; the second write replaces the first record.
    assert first["pipeline_version"] == second["pipeline_version"]


def test_load_returns_empty_skeleton_when_missing(tmp_path):
    state = _load_pipeline_state(tmp_path)
    assert state == {"stages": {}}


def test_load_returns_empty_skeleton_when_malformed(tmp_path):
    (tmp_path / "pipeline_state.json").write_text("not json")
    state = _load_pipeline_state(tmp_path)
    assert state == {"stages": {}}


# ---------------------------------------------------------------------------
# _stage_recorded_complete
# ---------------------------------------------------------------------------


def test_recorded_complete_true_after_record(tmp_path):
    _record_stage_completion(tmp_path, "scan_detection")
    assert _stage_recorded_complete(tmp_path, "scan_detection")


def test_recorded_complete_false_when_no_record(tmp_path):
    assert not _stage_recorded_complete(tmp_path, "scan_detection")


def test_recorded_complete_false_when_pipeline_version_mismatch(tmp_path):
    state = {
        "stages": {
            "scan_detection": {
                "completed_at": "2026-01-01T00:00:00Z",
                "pipeline_version": "0.0.0-old",
                "input_fingerprint": {},
            }
        }
    }
    (tmp_path / "pipeline_state.json").write_text(json.dumps(state))
    assert not _stage_recorded_complete(tmp_path, "scan_detection")


def test_recorded_complete_false_on_fingerprint_mismatch(tmp_path):
    _record_stage_completion(tmp_path, "taxa_and_lexicon_extraction",
                             input_fingerprint={"taxonomy": {"sha256": "old"}})
    assert not _stage_recorded_complete(
        tmp_path, "taxa_and_lexicon_extraction",
        expected_fingerprint={"taxonomy": {"sha256": "new"}},
    )


def test_recorded_complete_true_when_fingerprint_matches(tmp_path):
    fp = {"anatomy": {"sha256": "abc"}}
    _record_stage_completion(tmp_path, "taxa_and_lexicon_extraction",
                             input_fingerprint=fp)
    assert _stage_recorded_complete(
        tmp_path, "taxa_and_lexicon_extraction", expected_fingerprint=fp,
    )


# ---------------------------------------------------------------------------
# _should_run_stage
# ---------------------------------------------------------------------------


def test_should_run_when_resume_off(tmp_path):
    """Without --resume, every stage runs even when recorded complete."""
    _record_stage_completion(tmp_path, "scan_detection")
    ps = {}
    assert _should_run_stage(
        "scan_detection", hash_dir=tmp_path,
        resume=False, processing_summary=ps,
    )
    assert ps.get("skipped_stages") in (None, [])


def test_should_run_when_no_record(tmp_path):
    ps = {}
    assert _should_run_stage(
        "scan_detection", hash_dir=tmp_path,
        resume=True, processing_summary=ps,
    )
    assert ps.get("skipped_stages") in (None, [])


def test_should_skip_when_resume_on_and_recorded(tmp_path):
    _record_stage_completion(tmp_path, "scan_detection")
    ps = {}
    assert not _should_run_stage(
        "scan_detection", hash_dir=tmp_path,
        resume=True, processing_summary=ps,
    )
    assert ps["skipped_stages"] == ["scan_detection"]


def test_should_run_on_fingerprint_drift(tmp_path):
    """Lexicon swap forces re-run even with --resume + matching version."""
    _record_stage_completion(tmp_path, "taxa_and_lexicon_extraction",
                             input_fingerprint={"anatomy": {"sha256": "old"}})
    ps = {}
    assert _should_run_stage(
        "taxa_and_lexicon_extraction", hash_dir=tmp_path,
        resume=True, processing_summary=ps,
        expected_fingerprint={"anatomy": {"sha256": "new"}},
    )


def test_should_skip_records_each_stage(tmp_path):
    _record_stage_completion(tmp_path, "foo")
    _record_stage_completion(tmp_path, "bar")
    ps = {}
    _should_run_stage("foo", hash_dir=tmp_path, resume=True, processing_summary=ps)
    _should_run_stage("bar", hash_dir=tmp_path, resume=True, processing_summary=ps)
    assert ps["skipped_stages"] == ["foo", "bar"]


# ---------------------------------------------------------------------------
# _stage context manager records completion
# ---------------------------------------------------------------------------


def test_stage_records_completion_on_success(tmp_path):
    ps = {"stage_timings": [], "stage_failures": []}
    with _stage(ps, "scan_detection", hash_dir=tmp_path):
        pass
    assert _stage_recorded_complete(tmp_path, "scan_detection")


def test_stage_does_not_record_on_failure(tmp_path):
    ps = {"stage_timings": [], "stage_failures": []}
    try:
        with _stage(ps, "scan_detection", hash_dir=tmp_path):
            raise RuntimeError("boom")
    except RuntimeError:
        pass
    assert not _stage_recorded_complete(tmp_path, "scan_detection")
    # Failure is structured rather than silent.
    assert ps["stage_failures"][0]["stage"] == "scan_detection"


def test_stage_records_input_fingerprint(tmp_path):
    ps = {"stage_timings": [], "stage_failures": []}
    fp = {"taxonomy": {"sha256": "xyz"}}
    with _stage(ps, "taxa_and_lexicon_extraction", hash_dir=tmp_path,
                input_fingerprint=fp):
        pass
    rec = _load_pipeline_state(tmp_path)["stages"]["taxa_and_lexicon_extraction"]
    assert rec["input_fingerprint"] == fp


# ---------------------------------------------------------------------------
# _all_stage_artifacts_complete
# ---------------------------------------------------------------------------


_CORE = ("scan_detection", "pdf_preparation", "docling_extraction",
         "metadata_extraction", "text_chunking")


def _record_all_core(hd: Path) -> None:
    for name in _CORE:
        _record_stage_completion(hd, name)


def test_all_complete_when_core_recorded(tmp_path):
    _record_all_core(tmp_path)
    assert _all_stage_artifacts_complete(tmp_path)


def test_all_complete_taxa_required_recorded(tmp_path):
    _record_all_core(tmp_path)
    _record_stage_completion(tmp_path, "taxa_and_lexicon_extraction")
    assert _all_stage_artifacts_complete(
        tmp_path,
        expected_stages=list(_CORE) + ["taxa_and_lexicon_extraction"],
    )


def test_all_complete_taxa_required_missing(tmp_path):
    _record_all_core(tmp_path)
    # taxa_and_lexicon_extraction not recorded
    assert not _all_stage_artifacts_complete(
        tmp_path,
        expected_stages=list(_CORE) + ["taxa_and_lexicon_extraction"],
    )


def test_all_complete_chunks_missing(tmp_path):
    _record_all_core(tmp_path)
    # Manually drop the chunking record to simulate partial state.
    state = _load_pipeline_state(tmp_path)
    del state["stages"]["text_chunking"]
    (tmp_path / "pipeline_state.json").write_text(json.dumps(state))
    assert not _all_stage_artifacts_complete(tmp_path)


def test_all_complete_false_when_pipeline_version_drifts(tmp_path):
    """A corpus produced under an older PIPELINE_VERSION must not be
    treated as complete after a breaking release."""
    state = {
        "stages": {
            name: {
                "completed_at": "2026-01-01T00:00:00Z",
                "pipeline_version": "0.0.0-old",
                "input_fingerprint": {},
            }
            for name in _CORE
        }
    }
    (tmp_path / "pipeline_state.json").write_text(json.dumps(state))
    assert not _all_stage_artifacts_complete(tmp_path)
