"""Tests for the end-of-run `run.log` writer in `corpus run` (#57).

The helper lives at ``pipeline.cli._write_run_log``. It captures the
same content as ``corpus status --report`` (rollup + cross-paper
artifact presence) and writes it to ``<output_dir>/run.log`` so the
operator has a self-contained record without re-running status.

These tests pin two contracts:
- Empty ``output_dir`` (no ``documents/`` yet) → no log written,
  returns ``None``. Caller falls back to the legacy success message.
- Populated ``output_dir`` → log file present with a header (version
  + timestamp + path) and a non-empty body.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline.cli import _write_run_log
from pipeline.version import __version__


def _write_summary(documents_dir: Path, h: str, processing_summary: dict) -> None:
    hd = documents_dir / h
    hd.mkdir(parents=True, exist_ok=True)
    summary = {"pdf_hash": h, "processing_summary": processing_summary}
    (hd / "summary.json").write_text(json.dumps(summary))


@pytest.fixture
def populated_output_dir(tmp_path: Path) -> Path:
    """A minimal output_dir with one clean paper, mirroring the
    fixture pattern in test_corpus_status.py."""
    output_dir = tmp_path / "output"
    docs = output_dir / "documents"
    docs.mkdir(parents=True)
    _write_summary(docs, "aaa", {
        "stage_timings": [
            {"stage": "scan_detection", "ok": True, "duration_s": 0.1,
             "started_at": "2026-05-05T00:00:00Z",
             "ended_at": "2026-05-05T00:00:00Z"},
        ],
        "stage_failures": [],
        "quality_flags": [],
    })
    return output_dir


def test_returns_none_when_documents_missing(tmp_path: Path):
    """If the orchestrator bailed before extract ran, there's no
    documents/ tree to report on — the helper returns None so the
    caller can fall back to the legacy success message."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    assert _write_run_log(output_dir) is None
    assert not (output_dir / "run.log").exists()


def test_writes_run_log_with_header_and_body(populated_output_dir: Path):
    log_path = _write_run_log(populated_output_dir)
    assert log_path is not None
    assert log_path == populated_output_dir / "run.log"
    assert log_path.exists()

    content = log_path.read_text(encoding="utf-8")
    # Header: timestamp + version + resolved path
    assert content.startswith("# corpus run report\n")
    assert f"# corpus version {__version__}\n" in content
    assert str(populated_output_dir.resolve()) in content
    # Body: rollup text. The minimal one-paper fixture above runs
    # scan_detection cleanly, so the per-stage pass-rate output
    # mentions it.
    assert "scan_detection" in content


def test_overwrites_prior_log(populated_output_dir: Path):
    """Each `corpus run` replaces the prior log with the current
    snapshot — users wanting history are expected to git the file
    or copy it aside. Verify the file is rewritten, not appended."""
    log_path = populated_output_dir / "run.log"
    log_path.write_text("stale content from a prior run\n", encoding="utf-8")
    _write_run_log(populated_output_dir)
    content = log_path.read_text(encoding="utf-8")
    assert "stale content" not in content
    assert content.startswith("# corpus run report\n")
