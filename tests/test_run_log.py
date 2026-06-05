"""Tests for the end-of-run `run.log` writer in `corpus run` (#57, #90).

The helper lives at ``pipeline.cli._write_run_log``. It captures the
same content as ``corpus status --report`` (rollup + cross-paper
artifact presence) and writes it to ``<output_dir>/run.log`` so the
operator has a self-contained record without re-running status.

#90 extends it into a central per-invocation log: a timestamped archive
copy under ``<output_dir>/runs/<ts>/run.log`` (the top-level copy is
kept as "latest" for #57 back-compat), plus argv, the resolved config,
dependency-stack versions, and aggregate stage success/failure counts.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline.cli import _write_run_log
from pipeline.config_schema import validate_config
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


def test_overwrites_prior_top_level_log(populated_output_dir: Path):
    """Each `corpus run` replaces the top-level "latest" copy with the
    current snapshot (the per-invocation archive under runs/ keeps the
    history). Verify the top-level file is rewritten, not appended."""
    log_path = populated_output_dir / "run.log"
    log_path.write_text("stale content from a prior run\n", encoding="utf-8")
    _write_run_log(populated_output_dir)
    content = log_path.read_text(encoding="utf-8")
    assert "stale content" not in content
    assert content.startswith("# corpus run report\n")


# ---------------------------------------------------------------------------
# #90 — central per-invocation log
# ---------------------------------------------------------------------------


def test_writes_timestamped_archive_copy(populated_output_dir: Path):
    """The per-invocation archive lands under runs/<ts>/run.log with
    content identical to the top-level "latest" copy."""
    top = _write_run_log(populated_output_dir)
    runs_dir = populated_output_dir / "runs"
    assert runs_dir.is_dir()
    archives = list(runs_dir.glob("*/run.log"))
    assert len(archives) == 1
    # The single timestamped dir name is a compact ISO stamp.
    assert archives[0].parent.name.isdigit() is False  # contains a 'T'
    assert "T" in archives[0].parent.name
    assert archives[0].read_text(encoding="utf-8") == top.read_text(encoding="utf-8")


def test_log_includes_argv_config_and_dep_versions(populated_output_dir: Path):
    cfg = validate_config({"output_dir": str(populated_output_dir)})
    argv = ["corpus", "run", "--config", "demo/config.yaml"]
    log_path = _write_run_log(populated_output_dir, cfg=cfg, argv=argv)
    content = log_path.read_text(encoding="utf-8")

    # argv echoed in the header
    assert "# argv corpus run --config demo/config.yaml" in content
    # resolved config block (a known default key)
    assert "Resolved config:" in content
    assert "output_dir:" in content
    # dependency versions — python is always present
    assert "Dependency versions:" in content
    assert "python" in content
    # aggregate stage counts
    assert "Run summary:" in content
    assert "stage runs ok:" in content


def test_run_summary_counts_match_rollup(populated_output_dir: Path):
    """The one clean scan_detection stage in the fixture → 1 ok, 0 fail."""
    log_path = _write_run_log(populated_output_dir)
    content = log_path.read_text(encoding="utf-8")
    assert "stage runs ok:      1" in content
    assert "stage runs failed:  0" in content
