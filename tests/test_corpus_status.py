"""Unit tests for corpus_status aggregation + filtering (#40)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from corpus_status import aggregate, filtered_hashes, render_text


def _write_summary(documents_dir: Path, h: str, processing_summary: dict) -> None:
    hd = documents_dir / h
    hd.mkdir(parents=True, exist_ok=True)
    summary = {
        "pdf_hash": h,
        "processing_summary": processing_summary,
    }
    (hd / "summary.json").write_text(json.dumps(summary))


@pytest.fixture
def documents_dir(tmp_path: Path) -> Path:
    docs = tmp_path / "documents"
    docs.mkdir()

    # paper aaa: clean — all stages succeed, no flags
    _write_summary(docs, "aaa", {
        "stage_timings": [
            {"stage": "scan_detection", "ok": True, "duration_s": 0.1,
             "started_at": "2026-05-05T00:00:00Z", "ended_at": "2026-05-05T00:00:00Z"},
            {"stage": "docling_extraction", "ok": True, "duration_s": 5.0,
             "started_at": "2026-05-05T00:00:01Z", "ended_at": "2026-05-05T00:00:06Z"},
        ],
        "stage_failures": [],
        "quality_flags": [],
    })

    # paper bbb: docling timeout
    _write_summary(docs, "bbb", {
        "stage_timings": [
            {"stage": "scan_detection", "ok": True, "duration_s": 0.1,
             "started_at": "2026-05-05T00:00:00Z", "ended_at": "2026-05-05T00:00:00Z"},
            {"stage": "docling_extraction", "ok": False, "duration_s": 600.0,
             "started_at": "2026-05-05T00:00:01Z", "ended_at": "2026-05-05T00:10:01Z"},
        ],
        "stage_failures": [
            {"stage": "docling_extraction", "reason_code": "timeout",
             "reason_detail": "TimeoutExpired: …",
             "started_at": "2026-05-05T00:00:01Z", "ended_at": "2026-05-05T00:10:01Z",
             "duration_s": 600.0},
        ],
        "quality_flags": [],
    })

    # paper ccc: Grobid external_unavailable + gibberish_after_ocr flag
    _write_summary(docs, "ccc", {
        "stage_timings": [
            {"stage": "metadata_extraction", "ok": False, "duration_s": 5.0,
             "started_at": "2026-05-05T00:00:00Z", "ended_at": "2026-05-05T00:00:05Z"},
        ],
        "stage_failures": [
            {"stage": "metadata_extraction", "reason_code": "external_unavailable",
             "reason_detail": "GrobidUnavailableError: …",
             "started_at": "2026-05-05T00:00:00Z", "ended_at": "2026-05-05T00:00:05Z",
             "duration_s": 5.0},
        ],
        "quality_flags": [
            {"gate": "gibberish_after_ocr", "severity": "error",
             "detail": "score=0.6", "metric": 0.6},
        ],
    })

    # paper ddd: same Grobid failure as ccc (so reason×stage count = 2)
    _write_summary(docs, "ddd", {
        "stage_timings": [],
        "stage_failures": [
            {"stage": "metadata_extraction", "reason_code": "external_unavailable",
             "reason_detail": "GrobidUnavailableError: …",
             "started_at": "2026-05-05T00:00:00Z", "ended_at": "2026-05-05T00:00:05Z",
             "duration_s": 5.0},
        ],
        "quality_flags": [
            {"gate": "empty_text", "severity": "error", "detail": "0 chars", "metric": 0},
        ],
    })

    # paper eee: legacy errors[] only — pre-#34 paper
    _write_summary(docs, "eee", {
        "errors": ["something went wrong"],
        "processing_steps": ["scan_detection"],
    })

    return docs


def test_aggregate_counts(documents_dir):
    r = aggregate(documents_dir)
    assert r["total_documents"] == 5
    assert r["documents_with_summary"] == 5

    # Stages: scan_detection appears in 2 papers (aaa, bbb), both ok
    assert r["stages"]["scan_detection"] == {"ok": 2, "fail": 0, "total": 2}
    # docling: aaa ok, bbb fail
    assert r["stages"]["docling_extraction"] == {"ok": 1, "fail": 1, "total": 2}
    # metadata_extraction: ccc fail (ddd has no timing recorded)
    assert r["stages"]["metadata_extraction"] == {"ok": 0, "fail": 1, "total": 1}


def test_aggregate_failures_grouped(documents_dir):
    r = aggregate(documents_dir)
    # external_unavailable × metadata_extraction appears in ccc + ddd → 2
    assert r["failures_by_reason_stage"][("external_unavailable", "metadata_extraction")] == 2
    # timeout × docling_extraction appears in bbb → 1
    assert r["failures_by_reason_stage"][("timeout", "docling_extraction")] == 1


def test_aggregate_quality_flags(documents_dir):
    r = aggregate(documents_dir)
    assert r["quality_flags"]["gibberish_after_ocr"] == 1
    assert r["quality_flags"]["empty_text"] == 1
    assert r["papers_by_gate"]["gibberish_after_ocr"] == {"ccc"}


def test_aggregate_legacy_errors_only(documents_dir):
    r = aggregate(documents_dir)
    assert r["papers_with_legacy_errors_only"] == ["eee"]


def test_filter_by_reason(documents_dir):
    r = aggregate(documents_dir)
    assert filtered_hashes(r, filter_reason="external_unavailable") == ["ccc", "ddd"]
    assert filtered_hashes(r, filter_reason="timeout") == ["bbb"]


def test_filter_by_gate(documents_dir):
    r = aggregate(documents_dir)
    assert filtered_hashes(r, filter_gate="empty_text") == ["ddd"]


def test_filter_intersection(documents_dir):
    r = aggregate(documents_dir)
    # external_unavailable AND gibberish_after_ocr → only ccc
    assert filtered_hashes(
        r,
        filter_reason="external_unavailable",
        filter_gate="gibberish_after_ocr",
    ) == ["ccc"]


def test_no_filter_returns_failure_or_flag_papers(documents_dir):
    r = aggregate(documents_dir)
    # bbb (failure), ccc (both), ddd (both). aaa clean, eee legacy-only.
    assert filtered_hashes(r) == ["bbb", "ccc", "ddd"]


def test_render_text_does_not_crash(documents_dir):
    r = aggregate(documents_dir)
    out = render_text(r)
    assert "Stages" in out
    assert "Failures" in out
    assert "Quality flags" in out
    # Legacy note should mention eee
    assert "legacy errors" in out


# ---------------------------------------------------------------------------
# Stale-fingerprint detection (#29)
# ---------------------------------------------------------------------------


def _write_taxa(documents_dir: Path, h: str, sha: str) -> None:
    hd = documents_dir / h
    hd.mkdir(parents=True, exist_ok=True)
    (hd / "taxa.json").write_text(json.dumps({
        "total_mentions": 0, "unique_taxa": 0, "mentions": [],
        "input_fingerprint": {"path": "/tmp/taxonomy.sqlite", "sha256": sha, "size": 100},
    }))


def _write_anatomy(documents_dir: Path, h: str, sha: str) -> None:
    hd = documents_dir / h
    hd.mkdir(parents=True, exist_ok=True)
    (hd / "anatomy.json").write_text(json.dumps({
        "total_mentions": 0, "unique_terms": 0, "mentions": [],
        "input_fingerprint": {"path": "/tmp/anatomy.yaml", "sha256": sha, "size": 100},
    }))


def test_detect_stale_taxonomy(tmp_path):
    from corpus_status import detect_stale_fingerprints

    docs = tmp_path / "documents"
    docs.mkdir()
    # Current taxonomy has a known SHA-256
    cur_taxonomy = tmp_path / "taxonomy.sqlite"
    cur_taxonomy.write_bytes(b"current")
    import hashlib
    cur_sha = hashlib.sha256(b"current").hexdigest()
    old_sha = hashlib.sha256(b"old").hexdigest()

    _write_taxa(docs, "fresh", cur_sha)   # matches current — not stale
    _write_taxa(docs, "stale", old_sha)   # different sha — stale

    stale = detect_stale_fingerprints(docs, taxonomy_path=cur_taxonomy)
    assert stale["taxonomy_stale"] == ["stale"]
    assert stale["anatomy_stale"] == []


def test_detect_stale_anatomy(tmp_path):
    from corpus_status import detect_stale_fingerprints

    docs = tmp_path / "documents"
    docs.mkdir()
    cur_anatomy = tmp_path / "anatomy.yaml"
    cur_anatomy.write_bytes(b"current")
    import hashlib
    cur_sha = hashlib.sha256(b"current").hexdigest()
    old_sha = hashlib.sha256(b"old").hexdigest()

    _write_anatomy(docs, "fresh", cur_sha)
    _write_anatomy(docs, "stale1", old_sha)
    _write_anatomy(docs, "stale2", old_sha)

    stale = detect_stale_fingerprints(docs, anatomy_lexicon_path=cur_anatomy)
    assert sorted(stale["anatomy_stale"]) == ["stale1", "stale2"]


def test_pre_29_papers_without_fingerprint_not_flagged(tmp_path):
    from corpus_status import detect_stale_fingerprints

    docs = tmp_path / "documents"
    docs.mkdir()
    hd = docs / "old_paper"
    hd.mkdir()
    # taxa.json without input_fingerprint (pre-#29 layout)
    (hd / "taxa.json").write_text(json.dumps({
        "total_mentions": 0, "unique_taxa": 0, "mentions": [],
    }))

    cur_taxonomy = tmp_path / "taxonomy.sqlite"
    cur_taxonomy.write_bytes(b"x")

    stale = detect_stale_fingerprints(docs, taxonomy_path=cur_taxonomy)
    # No fingerprint to compare against → skipped, not flagged
    assert stale["taxonomy_stale"] == []


def test_filter_stale_taxonomy(tmp_path):
    from corpus_status import filtered_hashes

    rollup = aggregate(tmp_path / "nonexistent_documents") if False else {
        "papers_by_reason_stage": {},
        "papers_by_gate": {},
        "papers_with_failures": set(),
        "papers_with_quality_flags": set(),
    }
    stale = {"taxonomy_stale": ["aaa", "bbb"], "anatomy_stale": ["ccc"]}

    assert filtered_hashes(rollup, filter_stale="taxonomy", stale=stale) == ["aaa", "bbb"]
    assert filtered_hashes(rollup, filter_stale="anatomy", stale=stale) == ["ccc"]
    assert filtered_hashes(rollup, filter_stale="any", stale=stale) == ["aaa", "bbb", "ccc"]
