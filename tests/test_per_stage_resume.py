"""Unit tests for per-stage resume helpers (#28)."""
from __future__ import annotations

from pathlib import Path

import pytest

from process_corpus import (
    _all_stage_artifacts_complete,
    _should_run_stage,
    _stage_artifacts_present,
)


# ---------------------------------------------------------------------------
# _stage_artifacts_present
# ---------------------------------------------------------------------------


def test_artifacts_present_empty_returns_true():
    assert _stage_artifacts_present([])


def test_artifacts_present_all_exist_nonzero(tmp_path):
    p1, p2 = tmp_path / "a.json", tmp_path / "b.json"
    p1.write_text("{}")
    p2.write_text("{}")
    assert _stage_artifacts_present([p1, p2])


def test_artifacts_present_missing_returns_false(tmp_path):
    p = tmp_path / "missing.json"
    assert not _stage_artifacts_present([p])


def test_artifacts_present_zero_byte_returns_false(tmp_path):
    p = tmp_path / "empty.json"
    p.write_text("")
    assert not _stage_artifacts_present([p])


# ---------------------------------------------------------------------------
# _should_run_stage
# ---------------------------------------------------------------------------


def test_should_run_when_resume_off(tmp_path):
    """Without --resume, every stage runs even if artifact present."""
    artifact = tmp_path / "x.json"
    artifact.write_text("{}")
    ps = {}
    assert _should_run_stage("foo", [artifact], resume=False, processing_summary=ps)
    assert ps.get("skipped_stages") in (None, [])


def test_should_run_when_artifact_missing(tmp_path):
    artifact = tmp_path / "missing.json"
    ps = {}
    assert _should_run_stage("foo", [artifact], resume=True, processing_summary=ps)
    assert ps.get("skipped_stages") in (None, [])


def test_should_skip_when_resume_on_and_artifact_present(tmp_path):
    artifact = tmp_path / "x.json"
    artifact.write_text("{}")
    ps = {}
    assert not _should_run_stage("foo", [artifact], resume=True, processing_summary=ps)
    assert ps["skipped_stages"] == ["foo"]


def test_should_skip_records_each_stage(tmp_path):
    p1 = tmp_path / "a.json"
    p2 = tmp_path / "b.json"
    p1.write_text("{}")
    p2.write_text("{}")
    ps = {}
    _should_run_stage("foo", [p1], resume=True, processing_summary=ps)
    _should_run_stage("bar", [p2], resume=True, processing_summary=ps)
    assert ps["skipped_stages"] == ["foo", "bar"]


# ---------------------------------------------------------------------------
# _all_stage_artifacts_complete
# ---------------------------------------------------------------------------


def _write_all_required(hd: Path) -> None:
    """Populate the unconditional required artifacts."""
    for name in ("scan_detection.json", "processed.pdf", "text.json",
                 "figures.json", "metadata.json", "chunks.json"):
        (hd / name).write_text("{}" if name.endswith(".json") else "x")


def test_all_complete_when_no_taxa_no_anatomy(tmp_path):
    _write_all_required(tmp_path)
    assert _all_stage_artifacts_complete(tmp_path)


def test_all_complete_taxa_required_present(tmp_path):
    _write_all_required(tmp_path)
    (tmp_path / "taxa.json").write_text("{}")
    assert _all_stage_artifacts_complete(tmp_path, expect_taxa=True)


def test_all_complete_taxa_required_missing(tmp_path):
    _write_all_required(tmp_path)
    # taxa.json deliberately not written
    assert not _all_stage_artifacts_complete(tmp_path, expect_taxa=True)


def test_all_complete_anatomy_required_missing(tmp_path):
    _write_all_required(tmp_path)
    (tmp_path / "taxa.json").write_text("{}")
    # anatomy.json missing
    assert not _all_stage_artifacts_complete(
        tmp_path, expect_taxa=True, expect_anatomy=True,
    )


def test_all_complete_chunks_missing(tmp_path):
    _write_all_required(tmp_path)
    (tmp_path / "chunks.json").unlink()
    assert not _all_stage_artifacts_complete(tmp_path)


def test_all_complete_zero_byte_metadata_treated_as_missing(tmp_path):
    _write_all_required(tmp_path)
    (tmp_path / "metadata.json").write_text("")  # zero-byte
    assert not _all_stage_artifacts_complete(tmp_path)
