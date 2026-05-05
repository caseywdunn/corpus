"""Unit tests for the per-stage timing + structured failure recorder (#34)."""
from __future__ import annotations

import subprocess

import pytest
import requests

from process_corpus import _classify_exception, _stage, _utcnow_iso


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
    (ValueError("anything else"), "crash"),
    (RuntimeError("kaboom"), "crash"),
])
def test_classify_exception(exc, expected_code):
    code, detail = _classify_exception(exc)
    assert code == expected_code
    assert detail  # non-empty


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
