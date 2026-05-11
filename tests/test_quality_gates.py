"""Unit tests for the quality-gate helper (#36)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline import stages as process_corpus
from pipeline.stages import _run_quality_gates


@pytest.fixture
def fake_hash_dir(tmp_path: Path) -> Path:
    """A hash dir scaffold with all artifacts written; tests overwrite
    individual files to trigger specific gates.
    """
    hd = tmp_path / "abc123"
    hd.mkdir()
    (hd / "figures").mkdir()
    # Default: clean paper that fails no gates.
    (hd / "text.json").write_text(json.dumps({
        "title": "T", "text": "x" * 5000, "pages": 10,
    }))
    (hd / "chunks.json").write_text(json.dumps({
        "total_chunks": 3,
        "chunks": [{"text": "x" * 200} for _ in range(3)],
    }))
    (hd / "figures.json").write_text(json.dumps({"figures": [], "total_figures": 0}))
    (hd / "references.json").write_text(json.dumps({
        "references": [{"raw": "ref"}], "total_references": 1,
    }))
    (hd / "scan_detection.json").write_text(json.dumps({
        "needs_ocr": False, "gibberish_score": 0.0,
    }))
    return hd


def _gate_names(flags):
    return {f["gate"] for f in flags}


def test_clean_paper_produces_no_flags(fake_hash_dir):
    assert _run_quality_gates(fake_hash_dir) == []


def test_empty_text_flagged(fake_hash_dir):
    (fake_hash_dir / "text.json").write_text(json.dumps({
        "title": "", "text": "tiny", "pages": 10,
    }))
    flags = _run_quality_gates(fake_hash_dir)
    assert "empty_text" in _gate_names(flags)


def test_low_text_density_flagged(fake_hash_dir):
    # 1000 chars / 50 pages = 20 chars/page < 200 default threshold
    (fake_hash_dir / "text.json").write_text(json.dumps({
        "title": "", "text": "x" * 1000, "pages": 50,
    }))
    flags = _run_quality_gates(fake_hash_dir)
    assert "low_text_density" in _gate_names(flags)


def test_zero_references_unexpected_flagged(fake_hash_dir):
    (fake_hash_dir / "references.json").write_text(json.dumps({
        "references": [], "total_references": 0,
    }))
    # text.json default has pages=10, which is ≥ 5
    flags = _run_quality_gates(fake_hash_dir)
    assert "zero_references_unexpected" in _gate_names(flags)


def test_zero_references_short_paper_not_flagged(fake_hash_dir):
    (fake_hash_dir / "text.json").write_text(json.dumps({
        "title": "", "text": "x" * 5000, "pages": 2,
    }))
    (fake_hash_dir / "references.json").write_text(json.dumps({
        "references": [], "total_references": 0,
    }))
    flags = _run_quality_gates(fake_hash_dir)
    assert "zero_references_unexpected" not in _gate_names(flags)


def test_single_token_chunks_flagged(fake_hash_dir):
    (fake_hash_dir / "chunks.json").write_text(json.dumps({
        "total_chunks": 5,
        "chunks": [{"text": "ab"} for _ in range(5)],
    }))
    flags = _run_quality_gates(fake_hash_dir)
    assert "single_token_chunks" in _gate_names(flags)


def test_gibberish_after_ocr_only_runs_when_ocr_used(fake_hash_dir):
    # _gibberish_score flags tokens of length ≤ 2 — feed all-2-letter
    # tokens to drive the score to ~1.0.
    gib = " ".join(["ab", "cd", "ef", "gh", "ij", "kl"]) * 200
    (fake_hash_dir / "scan_detection.json").write_text(json.dumps({
        "needs_ocr": True, "gibberish_score": 0.0,
    }))
    (fake_hash_dir / "text.json").write_text(json.dumps({
        "title": "", "text": gib, "pages": 5,
    }))
    flags = _run_quality_gates(fake_hash_dir)
    assert "gibberish_after_ocr" in _gate_names(flags)


def test_gibberish_gate_skipped_for_born_digital(fake_hash_dir):
    # Same gibberish text but needs_ocr=False — gate should not fire
    gib = " ".join(["ab", "cd", "ef", "gh", "ij", "kl"]) * 200
    (fake_hash_dir / "text.json").write_text(json.dumps({
        "title": "", "text": gib, "pages": 5,
    }))
    flags = _run_quality_gates(fake_hash_dir)
    assert "gibberish_after_ocr" not in _gate_names(flags)


def test_missing_artifacts_do_not_crash(tmp_path):
    # Empty hash dir — every gate that needs a file should silently no-op
    hd = tmp_path / "empty"
    hd.mkdir()
    flags = _run_quality_gates(hd)
    # Empty body → "empty_text" still fires, which is correct
    assert "empty_text" in _gate_names(flags)
