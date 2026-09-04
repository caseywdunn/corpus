"""Pages the per-page OCR timeout silently blanks (#254).

``--tesseract-timeout`` does not fail a page. Its documented behaviour is to
give up on OCR, copy the un-OCR'd image into the output and exit 0 — so the
page survives visually, carries an empty text layer, and the document is
recorded ``status=success`` with an empty ``stage_failures``. Two full builds
of the same growing library lost text on ~9.5% of documents this way,
and the only trace was a warning in one array task's log.

The join that makes it detectable without a heuristic: ocrmypdf *names* the
pages it abandoned, one line each, page number first. Intersect those with the
pages that ended up with no text and what is left is exactly the lost pages —
a plate or a blank verso is empty but was never named, and a named page that
kept a pre-existing text layer is named but not empty.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pipeline.scan import _ocr_timeout_pages, _report_ocr_page_loss
from pipeline.stages import _run_quality_gates


# The real stderr shape, from Johnson_Widder2001.pdf on Bouchet.
STDERR = """\
    1 [tesseract] took too long to OCR - skipping
    2 [tesseract] took too long to OCR - skipping
   13 [tesseract] took too long to OCR - skipping
   14 page is facing up, no rotation needed
"""


# --- reading the page numbers back out ---------------------------------------


def test_the_abandoned_pages_are_named_not_merely_counted():
    assert _ocr_timeout_pages(STDERR) == [1, 2, 13]


def test_unrelated_stderr_names_nothing():
    assert _ocr_timeout_pages("   4 [tesseract] possibly poor OCR\n") == []
    assert _ocr_timeout_pages("") == []
    assert _ocr_timeout_pages(None) == []


def test_the_tag_may_change_spelling_between_ocrmypdf_releases():
    """What identifies the line is the page number and the message; the
    bracketed tag has moved before and must not be load-bearing."""
    assert _ocr_timeout_pages("  7 [tesseract eng] took too long to OCR\n") == [7]
    assert _ocr_timeout_pages("  7 took too long to OCR - skipping\n") == [7]


def test_the_verbose_log_format_is_read_too():
    """ocrmypdf puts a level and logger name ahead of the page number under
    -v, so the number cannot be anchored to the start of the line."""
    assert _ocr_timeout_pages(
        "[2026-08-28 21:44:40] - ocrmypdf._exec.tesseract - WARNING - "
        "    7 [tesseract] took too long to OCR - skipping\n") == [7]


# --- the intersection ---------------------------------------------------------


def _pdf(tmp_path: Path, texts) -> Path:
    """A PDF whose i-th page carries ``texts[i]`` (empty string = blank)."""
    fitz = pytest.importorskip("fitz")
    out = tmp_path / "processed.pdf"
    doc = fitz.open()
    for text in texts:
        page = doc.new_page()
        if text:
            page.insert_text((72, 72), text)
    doc.save(out)
    doc.close()
    return out


def test_a_blank_page_nobody_complained_about_is_not_data_loss(tmp_path):
    """Plates and blank versos are empty by nature. Flagging those at `error`
    would bury the real signal, which is the entire point of the join."""
    pdf = _pdf(tmp_path, ["text", "", "text"])
    out = _report_ocr_page_loss(pdf, "plate.pdf", stderr="")
    assert out["pages_without_text"] == [2]
    assert out["pages_blanked"] == []
    assert out["ocr_no_text_recovered"] is False


def test_a_named_page_that_kept_its_text_is_not_data_loss(tmp_path):
    """The timeout copies the page through, so a page that already had a text
    layer survives it — named, but nothing was lost."""
    pdf = _pdf(tmp_path, ["text", "text"])
    out = _report_ocr_page_loss(
        pdf, "mixed.pdf", stderr="  2 [tesseract] took too long to OCR\n")
    assert out["pages_ocr_timed_out"] == [2]
    assert out["pages_blanked"] == []


def test_empty_and_named_is_a_page_we_lost(tmp_path):
    pdf = _pdf(tmp_path, ["", "", "text", ""])
    out = _report_ocr_page_loss(
        pdf, "lost.pdf",
        stderr="  1 [tesseract] took too long to OCR\n"
               "  2 [tesseract] took too long to OCR\n")
    assert out["pages_blanked"] == [1, 2]
    assert out["pages_without_text"] == [1, 2, 4]   # page 4 is a plate
    assert out["page_count"] == 4


def test_success_with_every_page_textless_is_an_explicit_ocr_failure(tmp_path):
    pdf = _pdf(tmp_path, ["", "", ""])
    out = _report_ocr_page_loss(pdf, "textless.pdf", stderr="")
    assert out["pages_without_text"] == [1, 2, 3]
    assert out["ocr_no_text_recovered"] is True


def test_an_unreadable_pdf_still_reports_what_stderr_said(tmp_path):
    broken = tmp_path / "broken.pdf"
    broken.write_bytes(b"not a pdf")
    out = _report_ocr_page_loss(broken, "broken.pdf", stderr=STDERR)
    assert out["pages_blanked"] == [1, 2, 13]


# --- and out through the quality gate -----------------------------------------


@pytest.fixture
def hash_dir(tmp_path: Path) -> Path:
    """A document that passes every other gate, so a firing gate is this one."""
    hd = tmp_path / "3641b6b639c7"
    hd.mkdir()
    (hd / "figures").mkdir()
    (hd / "text.json").write_text(json.dumps({"text": "x" * 5000, "pages": 14}))
    (hd / "chunks.json").write_text(json.dumps(
        {"chunks": [{"text": "x" * 200} for _ in range(3)]}))
    (hd / "figures.json").write_text(json.dumps({"figures": []}))
    (hd / "references.json").write_text(json.dumps({"total_references": 4}))
    (hd / "scan_detection.json").write_text(json.dumps(
        {"needs_ocr": True, "gibberish_score": 0.0}))
    return hd


def _gates(hd: Path):
    return {f["gate"]: f for f in _run_quality_gates(hd)}


def test_the_gate_is_silent_when_nothing_was_blanked(hash_dir):
    assert "ocr_pages_blanked" not in _gates(hash_dir)
    assert "ocr_no_text_recovered" not in _gates(hash_dir)


def test_all_textless_pages_reach_an_error_gate_even_without_stderr(hash_dir):
    """A zero exit code and no timeout warning are not evidence of text."""
    scan = json.loads((hash_dir / "scan_detection.json").read_text())
    scan.update({
        "pages_without_text": list(range(1, 15)),
        "page_count": 14,
        "ocr_no_text_recovered": True,
        "pages_blanked": [],
    })
    (hash_dir / "scan_detection.json").write_text(json.dumps(scan))

    flag = _gates(hash_dir)["ocr_no_text_recovered"]
    assert flag["severity"] == "error"
    assert flag["metric"] == 14


def test_old_artifact_shape_also_derives_the_all_textless_gate(hash_dir):
    scan = json.loads((hash_dir / "scan_detection.json").read_text())
    scan.update({"pages_without_text": [1, 2], "page_count": 2})
    (hash_dir / "scan_detection.json").write_text(json.dumps(scan))
    assert "ocr_no_text_recovered" in _gates(hash_dir)


def test_blanked_pages_reach_the_quality_gates_as_an_error(hash_dir):
    """Johnson_Widder2001: 13 of 14 pages gone, 83,293 characters down to 671,
    and the two gates that did fire — low_text_density and
    zero_references_unexpected — are symptoms that name no cause."""
    scan = json.loads((hash_dir / "scan_detection.json").read_text())
    scan["pages_blanked"] = list(range(1, 14))
    (hash_dir / "scan_detection.json").write_text(json.dumps(scan))

    flag = _gates(hash_dir)["ocr_pages_blanked"]
    assert flag["severity"] == "error"
    assert flag["metric"] == 13
    assert "13/14" in flag["detail"]


def test_the_detail_does_not_print_a_thousand_page_numbers(hash_dir):
    scan = json.loads((hash_dir / "scan_detection.json").read_text())
    scan["pages_blanked"] = list(range(1, 200))
    (hash_dir / "scan_detection.json").write_text(json.dumps(scan))
    assert len(_gates(hash_dir)["ocr_pages_blanked"]["detail"]) < 200
