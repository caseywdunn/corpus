"""Extraction failure records must not outlive the failure (#—, v1.4 rebuild).

The 2026-09-08 siphonophore rebuild finished with all 1775 documents complete
and all 28 array tasks exiting 0, while `stage1_failures.json` still named one
document as crashed by signal 9. It was written by an earlier run and nothing
ever removed it, so the finished corpuscle carried a failure record that reads
exactly like a live one.

The obvious fix — clear the file when a run has no failures — is wrong, and
these tests pin why: 28 array tasks share one output_dir, so a clean shard
would delete the record a still-failing shard just wrote.
"""
from __future__ import annotations

import json

import pytest

from pipeline.main import (
    _FAILURES_DIR,
    _LEGACY_FAILURES_FILE,
    _failure_record_name,
    _record_stage1_failures,
)

FAIL = [{"pdf_hash": "abc123", "pdf_path": "/lib/x.pdf",
         "exitcode": -9, "signal": 9}]


def _records(root):
    d = root / _FAILURES_DIR
    return sorted(p.name for p in d.glob("*.json")) if d.is_dir() else []


# ── the bug ────────────────────────────────────────────────────────────

def test_a_shard_that_now_succeeds_retracts_its_own_record(tmp_path):
    """The actual observed bug: task 8 OOMed, was re-run, succeeded — and
    the corpuscle still claimed a failure."""
    _record_stage1_failures(tmp_path, 8, FAIL)
    assert _records(tmp_path) == ["batch-00008.json"]
    _record_stage1_failures(tmp_path, 8, [])          # the retry
    assert _records(tmp_path) == []


def test_a_failure_is_recorded_with_its_shard(tmp_path):
    path = _record_stage1_failures(tmp_path, 3, FAIL)
    saved = json.loads(path.read_text())
    assert saved["n_failed"] == 1
    assert saved["batch_index"] == 3
    assert saved["failures"][0]["pdf_hash"] == "abc123"
    assert saved["pipeline_version"]


# ── why the obvious fix is wrong ───────────────────────────────────────

def test_a_clean_shard_never_erases_another_shards_failure(tmp_path):
    """28 array tasks share one output_dir. If clearing were global, the
    27 that finished cleanly would race to delete task 8's record and the
    build would report success with a crashed document in it."""
    _record_stage1_failures(tmp_path, 8, FAIL)
    for other in (0, 1, 2, 27):
        _record_stage1_failures(tmp_path, other, [])
    assert _records(tmp_path) == ["batch-00008.json"]
    surviving = json.loads(
        (tmp_path / _FAILURES_DIR / "batch-00008.json").read_text())
    assert surviving["failures"][0]["pdf_hash"] == "abc123"


def test_two_failing_shards_both_keep_their_records(tmp_path):
    _record_stage1_failures(tmp_path, 4, FAIL)
    _record_stage1_failures(tmp_path, 9, [{"pdf_hash": "def456",
                                           "exitcode": 1, "signal": None}])
    assert _records(tmp_path) == ["batch-00004.json", "batch-00009.json"]


def test_shard_names_sort_naturally(tmp_path):
    """A listing of 28 shards should not put batch-10 before batch-2."""
    for i in (2, 10, 27):
        _record_stage1_failures(tmp_path, i, FAIL)
    assert _records(tmp_path) == [
        "batch-00002.json", "batch-00010.json", "batch-00027.json"]


# ── an unsharded run speaks for the whole corpus ───────────────────────

def test_an_unsharded_success_clears_every_stale_record(tmp_path):
    """It processed every document, so its verdict supersedes whatever an
    earlier sharded run left behind."""
    _record_stage1_failures(tmp_path, 8, FAIL)
    _record_stage1_failures(tmp_path, 12, FAIL)
    (tmp_path / _LEGACY_FAILURES_FILE).write_text('{"n_failed": 1}')
    _record_stage1_failures(tmp_path, None, [])
    assert _records(tmp_path) == []
    assert not (tmp_path / _LEGACY_FAILURES_FILE).exists()


def test_a_sharded_success_leaves_the_legacy_file_alone(tmp_path):
    """A shard cannot speak for slices it did not read, and the legacy
    file covers all of them."""
    (tmp_path / _LEGACY_FAILURES_FILE).write_text('{"n_failed": 1}')
    _record_stage1_failures(tmp_path, 8, [])
    assert (tmp_path / _LEGACY_FAILURES_FILE).exists()


def test_an_unsharded_failure_records_under_its_own_name(tmp_path):
    path = _record_stage1_failures(tmp_path, None, FAIL)
    assert path.name == "all.json"
    assert json.loads(path.read_text())["batch_index"] is None


# ── it must not be able to break a build ───────────────────────────────

def test_an_unwritable_record_is_reported_not_raised(tmp_path, monkeypatch):
    """This runs at the end of a multi-hour extract. Losing the run
    because the record could not be written would be worse than the
    missing record."""
    def boom(*a, **k):
        raise OSError("read-only filesystem")
    monkeypatch.setattr("pathlib.Path.mkdir", boom)
    assert _record_stage1_failures(tmp_path, 1, FAIL) is None


def test_clearing_a_record_that_is_not_there_is_fine(tmp_path):
    assert _record_stage1_failures(tmp_path, 5, []) is None


@pytest.mark.parametrize("idx,expected", [
    (None, "all.json"), (0, "batch-00000.json"), (27, "batch-00027.json")])
def test_record_names(idx, expected):
    assert _failure_record_name(idx) == expected
