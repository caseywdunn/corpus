"""Keep readable rotated source pages during automatic force OCR (#346)."""

from pathlib import Path
import shutil
from types import SimpleNamespace

import fitz

from pipeline import scan
from pipeline import native_text_recovery, source_spacing_recovery


def _pdf_with_rotated_plate(path: Path) -> None:
    with fitz.open() as pdf:
        before = pdf.new_page()
        before.insert_text((25, 40), "Prose before the landscape plate.")
        plate = pdf.new_page()
        plate.insert_text((25, 40),
                          "Fig. 210. Polymorphism in a bryozoan colony.\n"
                          "The original source identifies the zooid and its appendages.\n"
                          "The plate contains labels for the structure of the colony.\n"
                          "These captions have readable English words and figure numbers.")
        plate.set_rotation(270)
        prose = pdf.new_page()
        prose.insert_text((25, 40), "This adjacent prose page should still be OCRed.")
        other_orientation = pdf.new_page()
        other_orientation.insert_text((25, 40),
                                      "Fig. 211. A readable figure caption with many English words.\n"
                                      "This plate is rotated ninety degrees, not the source failure.\n"
                                      "Its existing text must not trigger the narrow exception.\n"
                                      "The surrounding prose is present and readable.")
        other_orientation.set_rotation(90)
        after = pdf.new_page()
        after.insert_text((25, 40), "Prose after the ninety degree plate.")
        damaged = pdf.new_page()
        damaged.insert_text((25, 40), "x y z 1 2 3")
        damaged.set_rotation(270)
        pdf.save(path)


def test_automatic_ocr_preserves_only_readable_rotated_pages(monkeypatch, tmp_path):
    source, output = tmp_path / "source.pdf", tmp_path / "prepared.pdf"
    _pdf_with_rotated_plate(source)
    detection = {
        "needs_ocr": True, "has_text": True, "ocr_mode": "force_ocr",
        "ocrlang_honored": True, "tesseract_packs": ["eng"],
        "keeppages_selected": [5, 6, 7, 8, 9, 10],
    }
    assert scan._readable_rotated_pages(source, detection) == (6, [2])

    captured = {}

    def fake_ocr(command, _timeout):
        captured["command"] = command
        shutil.copyfile(source, output)
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(scan.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(scan, "_run_ocr", fake_ocr)
    monkeypatch.setattr(scan, "_report_ocr_page_loss", lambda *_args: {})
    monkeypatch.setattr(native_text_recovery, "inspect_native_text_regions",
                        lambda *_args: {"candidate_count": 0, "source_exponents": {"candidate_count": 0}})
    monkeypatch.setattr(source_spacing_recovery, "inspect_source_spacing",
                        lambda *_args: {"candidate_count": 0, "unresolved": []})

    outcome = scan.prepare_pdf(source, detection, output)
    command = captured["command"]
    assert command[command.index("--pages") + 1] == "1,3-6"
    assert outcome["rotated_native_text_preservation"] == {
        "policy": scan.ROTATED_NATIVE_TEXT_POLICY,
        "pages": [2],
        "original_pages": [6],
    }
    assert command[-2:] == [str(source), str(output)]


def test_explicit_force_keeps_all_pages_selected(tmp_path):
    source = tmp_path / "source.pdf"
    _pdf_with_rotated_plate(source)
    detection = {
        "needs_ocr": True, "has_text": True, "ocr_mode": "force_ocr",
        "ocrmode_honored": True,
    }
    assert scan._readable_rotated_pages(source, detection) == (0, [])
