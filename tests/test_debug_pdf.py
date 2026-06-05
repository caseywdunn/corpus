"""Tests for `corpus debug-pdf` — the single-PDF debug runner (#92).

The command drives one PDF through ``run_pdf_processing_pipeline`` with
stage tracing and no bundle. The actual extraction pipeline needs
docling + models, so these tests stub ``run_pdf_processing_pipeline``
and pin the command's *wiring*: hash-dir layout, summary.json write,
exit code reflecting stage status, and the not-a-file guard.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pytest

import pipeline.cli as cli
from pipeline.io import calculate_pdf_hash, short_hash


def _args(pdf: Path, output_dir: Path, grobid_url=None, verbose=0):
    return argparse.Namespace(
        pdf=pdf, output_dir=output_dir, grobid_url=grobid_url, verbose=verbose,
    )


@pytest.fixture
def fake_pdf(tmp_path: Path) -> Path:
    # calculate_pdf_hash only reads bytes; content needn't be real PDF.
    p = tmp_path / "hard_document.pdf"
    p.write_bytes(b"%PDF-1.4 fake bytes for hashing\n")
    return p


def _stub_pipeline(summary: dict):
    """Return a stub for run_pdf_processing_pipeline that records its
    call and returns the given summary."""
    calls = {}

    def _stub(pdf_path, hash_dir, temp_dir, **kwargs):
        calls["pdf_path"] = pdf_path
        calls["hash_dir"] = hash_dir
        calls["kwargs"] = kwargs
        return summary

    _stub.calls = calls
    return _stub


def test_debug_pdf_success(fake_pdf, tmp_path, monkeypatch, capsys):
    out = tmp_path / "debug-out"
    summary = {
        "status": "success",
        "stage_timings": [
            {"stage": "scan_detection", "ok": True, "duration_s": 0.1},
            {"stage": "docling_extraction", "ok": True, "duration_s": 2.5},
        ],
        "stage_failures": [],
        "skipped_stages": [],
        "files_created": ["text.json", "chunks.json"],
    }
    stub = _stub_pipeline(summary)
    monkeypatch.setattr("pipeline.runner.run_pdf_processing_pipeline", stub)

    rc = cli._cmd_debug_pdf(_args(fake_pdf, out))
    assert rc == cli.EXIT_OK

    # Hash-dir layout mirrors documents/<short-hash>/.
    sh = short_hash(calculate_pdf_hash(fake_pdf))
    hash_dir = out / sh
    assert hash_dir.is_dir()
    assert (hash_dir / "summary.json").is_file()
    assert stub.calls["hash_dir"] == hash_dir
    # No grobid by default → bibliography stage gets no client.
    assert stub.calls["kwargs"]["grobid_client"] is None
    assert stub.calls["kwargs"]["resume"] is False

    captured = capsys.readouterr().out
    assert "Stage trace" in captured
    assert "scan_detection" in captured
    assert "docling_extraction" in captured


def test_debug_pdf_failure_exit_code(fake_pdf, tmp_path, monkeypatch):
    summary = {
        "status": "failed",
        "stage_timings": [
            {"stage": "docling_extraction", "ok": False, "duration_s": 1.0},
        ],
        "stage_failures": [
            {"stage": "docling_extraction", "error": "docling crashed"},
        ],
        "skipped_stages": [],
        "files_created": [],
    }
    monkeypatch.setattr(
        "pipeline.runner.run_pdf_processing_pipeline", _stub_pipeline(summary),
    )
    rc = cli._cmd_debug_pdf(_args(fake_pdf, tmp_path / "out"))
    assert rc == cli.EXIT_GENERIC


def test_debug_pdf_grobid_url_wires_client(fake_pdf, tmp_path, monkeypatch):
    summary = {"status": "success", "stage_timings": [], "stage_failures": [],
               "skipped_stages": [], "files_created": []}
    stub = _stub_pipeline(summary)
    monkeypatch.setattr("pipeline.runner.run_pdf_processing_pipeline", stub)

    sentinel = object()
    monkeypatch.setattr(
        "pipeline.grobid_client.GrobidClient", lambda base_url: sentinel,
    )
    cli._cmd_debug_pdf(_args(fake_pdf, tmp_path / "out",
                             grobid_url="http://localhost:8070"))
    assert stub.calls["kwargs"]["grobid_client"] is sentinel


def test_debug_pdf_rejects_non_file(tmp_path):
    missing = tmp_path / "nope.pdf"
    rc = cli._cmd_debug_pdf(_args(missing, tmp_path / "out"))
    assert rc == cli.EXIT_CONFIG_ERROR


def test_debug_pdf_registered_in_parser():
    """The subparser is wired into the router with the right handler."""
    parser = cli._build_parser()
    args = parser.parse_args(["debug-pdf", "some.pdf", "--grobid-url", "http://x"])
    assert args.func is cli._cmd_debug_pdf
    assert args.pdf == Path("some.pdf")
    assert args.grobid_url == "http://x"
