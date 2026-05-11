"""Unit tests for the per-stage timing + structured failure recorder (#34)."""
from __future__ import annotations

import subprocess

import pytest
import requests

from pipeline.stages import (
    _HugeDocumentError,
    _classify_exception,
    _pdf_page_count,
    _stage,
    _utcnow_iso,
)


def test_utcnow_iso_format():
    ts = _utcnow_iso()
    # YYYY-MM-DDTHH:MM:SSZ — second precision, no microseconds, Z suffix
    assert len(ts) == 20
    assert ts.endswith("Z")
    assert "T" in ts


@pytest.mark.parametrize("exc, expected_code", [
    (subprocess.TimeoutExpired("cmd", 5), "timeout"),
    (requests.ConnectionError("conn"), "external_unavailable"),
    (requests.Timeout("slow"), "external_unavailable"),
    (Exception("PDF is encrypted"), "unsupported_format"),
    (Exception("password-protected"), "unsupported_format"),
    (Exception("Invalid PDF: corrupt header"), "corrupted"),
    (Exception("syntax error in xref table"), "corrupted"),
    (_HugeDocumentError("PDF has 600 pages"), "too_large"),
    (ValueError("anything else"), "crash"),
    (RuntimeError("kaboom"), "crash"),
])
def test_classify_exception(exc, expected_code):
    code, detail = _classify_exception(exc)
    assert code == expected_code
    assert detail  # non-empty


def test_pdf_page_count_returns_int_for_real_pdf():
    """Demo PDFs are real born-digital papers; page_count should resolve."""
    from pathlib import Path
    demo = Path(__file__).resolve().parent.parent / "demo"
    if not demo.is_dir():
        pytest.skip("no demo/ corpus available")
    sample = next(demo.glob("*.pdf"), None)
    if sample is None:
        pytest.skip("no PDFs under demo/")
    n = _pdf_page_count(sample)
    assert isinstance(n, int)
    assert n > 0


def test_pdf_page_count_returns_none_for_missing_file(tmp_path):
    assert _pdf_page_count(tmp_path / "nope.pdf") is None


def test_stage_records_success_timing_only():
    ps = {}
    with _stage(ps, "demo_ok"):
        pass
    assert len(ps["stage_timings"]) == 1
    t = ps["stage_timings"][0]
    assert t["stage"] == "demo_ok"
    assert t["ok"] is True
    assert t["duration_s"] >= 0
    assert "stage_failures" not in ps


def test_stage_records_failure_in_both_lists_and_reraises():
    ps = {}
    with pytest.raises(RuntimeError, match="boom"):
        with _stage(ps, "demo_fail"):
            raise RuntimeError("boom")
    assert len(ps["stage_timings"]) == 1
    assert ps["stage_timings"][0]["ok"] is False
    assert len(ps["stage_failures"]) == 1
    f = ps["stage_failures"][0]
    assert f["stage"] == "demo_fail"
    assert f["reason_code"] == "crash"
    assert "RuntimeError" in f["reason_detail"]
    assert f["duration_s"] >= 0


def test_stage_failure_uses_classified_reason_code():
    ps = {}
    with pytest.raises(subprocess.TimeoutExpired):
        with _stage(ps, "demo_timeout"):
            raise subprocess.TimeoutExpired("cmd", 5)
    assert ps["stage_failures"][0]["reason_code"] == "timeout"


def test_huge_document_gate_emits_too_large():
    """End-to-end: simulate the huge-document gate by raising
    _HugeDocumentError inside _stage; confirm it surfaces as
    reason_code=too_large with the page count in the detail.
    """
    ps = {}
    with pytest.raises(_HugeDocumentError):
        with _stage(ps, "huge_document_check"):
            raise _HugeDocumentError("PDF has 600 pages (max: 500)")
    failure = ps["stage_failures"][0]
    assert failure["stage"] == "huge_document_check"
    assert failure["reason_code"] == "too_large"
    assert "600 pages" in failure["reason_detail"]


def test_stage_appends_in_order_for_multiple_stages():
    ps = {}
    with _stage(ps, "first"):
        pass
    with _stage(ps, "second"):
        pass
    with pytest.raises(RuntimeError):
        with _stage(ps, "third"):
            raise RuntimeError("x")
    assert [t["stage"] for t in ps["stage_timings"]] == ["first", "second", "third"]
    assert [t["ok"] for t in ps["stage_timings"]] == [True, True, False]
    assert [f["stage"] for f in ps["stage_failures"]] == ["third"]
