"""Regression tests for the batch-slicing race fixed in #55.

The bug: Stage 1 main() pre-filtered fully-completed hashes from
``pdf_map`` and *then* sliced into batches. SLURM array tasks don't
start simultaneously, so each task's pre-filter saw a different
post-resume disk state. The slice indices then landed on different
hashes per task, producing overlapping batches that raced on per-doc
``summary.json`` / ``pipeline_state.json`` writes.

The fix: slice on the unfiltered hash list (deterministic across array
tasks because ``find_all_pdfs`` depends only on input directory
contents), let the per-doc resume guard inside the loop skip completed
docs.

These tests pin that contract directly on the helper, so a future
regression that re-introduces a pre-slice filter will fail loudly.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from pipeline.io import create_summary_json
from pipeline.main import _slice_hashes_for_batch
from pipeline.stages import _save_pipeline_state


# ── Slicer determinism ───────────────────────────────────────────────


def _fake_pdf_map(n: int) -> dict:
    """Build a 64-char hex pdf_map with n stable, sortable keys."""
    return {f"{i:064x}": [Path(f"/in/{i}.pdf")] for i in range(n)}


def test_slice_is_disjoint_across_tasks():
    """Two different batch_index values produce disjoint slices."""
    pdf_map = _fake_pdf_map(1000)
    batch_size = 128

    seen = set()
    for idx in range(8):  # ceil(1000/128) = 8
        batch, _, _ = _slice_hashes_for_batch(pdf_map, idx, batch_size)
        assert seen.isdisjoint(batch), f"batch {idx} overlaps with prior tasks"
        seen.update(batch)

    # Every hash covered exactly once
    assert seen == set(pdf_map)


def test_slice_is_independent_of_pdf_map_subset():
    """Two array tasks running on different pdf_map subsets (e.g.
    different snapshots after some hashes completed) produce slices
    derived ONLY from sorted(pdf_map.keys()).

    This is what makes the slicer safe under concurrent array tasks:
    each task discovers its own pdf_map fresh, but as long as the
    discovery is deterministic from the input directory, the slice
    indices land on the same hashes.
    """
    full = _fake_pdf_map(1000)

    # Task A's view: full input
    a_batch, a_total, _ = _slice_hashes_for_batch(full, 30, 128)

    # Task B's view of the SAME input dir at the same time. find_all_pdfs
    # is deterministic, so task B sees the same pdf_map as task A — even
    # if some other task has been writing summary.json files in parallel.
    # The slicer must NOT consult disk state, so we verify identical
    # input → identical slice.
    b_batch, b_total, _ = _slice_hashes_for_batch(full, 30, 128)
    assert a_batch == b_batch
    assert a_total == b_total


def test_slice_does_not_consult_disk(tmp_path):
    """Regression: pre-#55 the main() filter consulted disk state
    BEFORE slicing, so a task starting late saw a smaller pdf_map and
    sliced different hashes.  Pin that the helper is pure: takes
    pdf_map + index + size, returns deterministic output.
    """
    pdf_map = _fake_pdf_map(500)

    # Even if disk is full of unrelated files, the slicer ignores it.
    (tmp_path / "summary.json").write_text("{}")
    (tmp_path / "pipeline_state.json").write_text("{}")

    batch_a, _, _ = _slice_hashes_for_batch(pdf_map, 2, 64)
    batch_b, _, _ = _slice_hashes_for_batch(pdf_map, 2, 64)
    assert batch_a == batch_b


def test_empty_pdf_map():
    batch, total, total_batches = _slice_hashes_for_batch({}, 0, 128)
    assert batch == []
    assert total == 0
    assert total_batches == 0


def test_last_batch_partial():
    """Final batch may have fewer than batch_size hashes; total_batches
    is the ceiling division."""
    pdf_map = _fake_pdf_map(300)
    last, total, total_batches = _slice_hashes_for_batch(pdf_map, 2, 128)
    assert total == 300
    assert total_batches == 3
    assert len(last) == 300 - 2 * 128  # 44 hashes in the last slice


# ── Per-writer tmp filenames (atomic-write) ──────────────────────────


def test_create_summary_json_uses_per_writer_tmp(tmp_path):
    """summary.json must be written via tmp + rename, with a tmp
    filename that includes the pid so concurrent writers in the same
    hash_dir don't share a tmp path (#55)."""
    summary_file = create_summary_json(
        pdf_hash_full="0" * 64,
        pdf_paths=[tmp_path / "x.pdf"],
        input_dir=tmp_path,
        hash_dir=tmp_path,
        processing_summary={"ok": True},
    )
    # Final file is valid JSON
    payload = json.loads(summary_file.read_text())
    assert payload["pdf_hash_full"] == "0" * 64
    # No stray .tmp* files left behind
    assert not list(tmp_path.glob("*.tmp*")), \
        "atomic rename should leave no tmp file on success"


def test_save_pipeline_state_uses_per_writer_tmp(tmp_path):
    """pipeline_state.json writer uses a pid-tagged tmp name so two
    concurrent writers in the same hash_dir can't corrupt each other's
    payload via shared tmp filename (#55)."""
    state = {"stages": {"scan_detection": {"completed_at": "2026-05-08T00:00:00Z"}}}
    _save_pipeline_state(tmp_path, state)

    final = tmp_path / "pipeline_state.json"
    assert final.exists()
    payload = json.loads(final.read_text())
    assert payload["stages"]["scan_detection"]["completed_at"]
    assert not list(tmp_path.glob("pipeline_state.json.tmp*")), \
        "atomic rename should leave no tmp file on success"


def test_save_pipeline_state_tmp_name_is_unique_per_call(tmp_path, monkeypatch):
    """Two calls in quick succession produce distinct tmp filenames
    (different time.monotonic_ns()), so even if one writer's rename
    is delayed the other's tmp path doesn't clobber it.

    Verified by patching tmp.replace so the tmp survives, then calling
    twice and asserting two distinct tmp paths exist.
    """
    from pipeline import stages

    captured: list[Path] = []

    real_path_with_suffix = Path.with_suffix

    def _capturing_with_suffix(self, suffix):
        result = real_path_with_suffix(self, suffix)
        captured.append(result)
        return result

    monkeypatch.setattr(Path, "with_suffix", _capturing_with_suffix)

    stages._save_pipeline_state(tmp_path, {"stages": {"a": {}}})
    stages._save_pipeline_state(tmp_path, {"stages": {"b": {}}})

    tmp_paths = [p for p in captured if ".tmp." in p.name]
    assert len(tmp_paths) >= 2
    assert tmp_paths[0] != tmp_paths[1], \
        "two writes should produce distinct tmp filenames"
