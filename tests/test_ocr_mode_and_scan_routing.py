"""OCR mode overrides and the scan-routing regressions in #186/#264/#266."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from bib import (
    BibIndex,
    bib_entry_to_metadata,
    entry_ocrmode,
    ocrmode_for_pdf,
    parse_bibtex,
)
from bib.authority import create_schema
from bib.export import export_bibtex
from bib.importer import import_bibtex
from pipeline import scan
from pipeline.stages import _OCR_DEPENDENT_STAGES, _expected_fingerprints_for_run


INSTALLED = frozenset({"eng", "chi_sim", "ell", "rus"})


class _Page:
    def __init__(self, text: str):
        self._text = text

    def get_text(self, *args, **kwargs):
        return self._text


class _Document:
    def __init__(self, text: str):
        self._pages = [_Page(text)]

    def __len__(self):
        return len(self._pages)

    def __getitem__(self, index):
        return self._pages[index]

    def close(self):
        pass


def _stub_detection(monkeypatch, text: str) -> Path:
    import fitz

    monkeypatch.setattr(fitz, "open", lambda *_a, **_k: _Document(text))
    monkeypatch.setattr(scan, "_available_tesseract_langs", lambda: INSTALLED)
    monkeypatch.setattr(scan, "_scanned_page_fraction", lambda *_a, **_k: 0.0)
    monkeypatch.setattr(scan, "detect_vertical_cjk", lambda *_a, **_k: None)
    monkeypatch.setitem(scan.CONFIG["ocr"], "probe_language_by_ocr", False)
    return Path("fixture.pdf")


def test_no_text_layer_uses_force_not_skip_text(monkeypatch):
    pdf = _stub_detection(monkeypatch, "")
    out = scan.detect_scan_type(pdf)
    assert out["detection_reason"] == "no_text_layer"
    assert out["ocr_mode"] == "force_ocr"


def test_vendor_boilerplate_uses_force_not_skip_text(monkeypatch):
    marker = scan._VENDOR_BOILERPLATE[0]
    pdf = _stub_detection(monkeypatch, (marker + " ") * 100)
    out = scan.detect_scan_type(pdf)
    assert out["detection_reason"] == "vendor_boilerplate_only"
    assert out["ocr_mode"] == "force_ocr"


def test_non_latin_pin_triggers_visual_check_below_gibberish_floor(monkeypatch):
    pdf = _stub_detection(monkeypatch, "legacy ascii glyph stream " * 300)
    monkeypatch.setattr(scan, "_detect_language", lambda _text: ("en", 1.0))
    monkeypatch.setattr(scan, "_gibberish_score", lambda _text: 0.295)
    monkeypatch.setattr(scan, "_text_layer_scripts", lambda _text: {"Latin": 1.0})
    monkeypatch.setattr(scan, "_visual_page_script", lambda _pdf: "Han")

    out = scan.detect_scan_type(pdf, ocrlang="chi_sim+eng")
    assert out["detection_reason"] == "visual_script_mismatch"
    assert out["file_type"] == "broken_text_layer"
    assert out["ocr_mode"] == "force_ocr"
    assert out["visual_script"] == "Han"
    assert out["ocrlang_applied"] is True


def test_latin_pin_does_not_expand_the_visual_probe_population(monkeypatch):
    pdf = _stub_detection(monkeypatch, "ordinary latin prose " * 300)
    monkeypatch.setattr(scan, "_detect_language", lambda _text: ("en", 1.0))
    monkeypatch.setattr(scan, "_gibberish_score", lambda _text: 0.295)
    monkeypatch.setattr(scan, "_text_layer_scripts", lambda _text: {"Latin": 1.0})

    def should_not_run(_pdf):
        raise AssertionError("ordinary Latin pin should not pay for OSD")

    monkeypatch.setattr(scan, "_visual_page_script", should_not_run)
    out = scan.detect_scan_type(pdf, ocrlang="eng")
    assert out["detection_reason"] == "clean_text_layer"
    assert out["ocrlang_honored"] is True
    assert out["ocrlang_applied"] is False


@pytest.mark.parametrize("raw,internal", [
    ("force", "force_ocr"),
    ("redo", "redo_ocr"),
    ("skip-text", "skip_text"),
])
def test_ocrmode_forces_ocr_and_selects_the_requested_mode(raw, internal):
    detected = {
        "filename": "clean.pdf",
        "needs_ocr": False,
        "ocr_mode": None,
    }
    out = scan._apply_ocrmode_override(raw, detected)
    assert out["needs_ocr"] is True
    assert out["ocr_mode"] == internal
    assert out["ocrmode_honored"] is True
    assert out["ocrmode_detection_needs_ocr"] is False
    assert out["ocrmode_detection_mode"] is None


def test_unknown_ocrmode_is_recorded_but_ignored(caplog):
    detected = {"filename": "clean.pdf", "needs_ocr": False, "ocr_mode": None}
    out = scan._apply_ocrmode_override("aggressive", detected)
    assert out["needs_ocr"] is False
    assert out["ocr_mode"] is None
    assert out["ocrmode_honored"] is False
    assert "force, redo, skip-text" in caplog.text


def test_ocrmode_bib_lookup_and_metadata_round_trip():
    entry = parse_bibtex(
        "@article{K, file={clean.pdf}, ocrmode={force}, ocrlang={ell+eng}}"
    )[0]
    index = BibIndex([entry])
    assert entry_ocrmode(entry) == "force"
    assert ocrmode_for_pdf(index, "clean.pdf") == "force"
    assert ocrmode_for_pdf(index, "absent.pdf") is None
    assert bib_entry_to_metadata(entry, "clean.pdf")["ocrmode"] == "force"


def test_ocrmode_fingerprints_every_descendant_stage():
    fps = _expected_fingerprints_for_run(
        ocrlang=None, ocrmode="force", keeppages=None,
    )
    for stage in _OCR_DEPENDENT_STAGES:
        assert fps[stage]["ocrmode"] == "force", stage


def test_ocrmode_fingerprint_argument_is_required():
    with pytest.raises(TypeError, match="ocrmode"):
        _expected_fingerprints_for_run(ocrlang=None, keeppages=None)


def _authority_db(path: Path) -> Path:
    conn = sqlite3.connect(path)
    create_schema(conn)
    conn.execute(
        "INSERT INTO works (work_id, guid_type, title, year, corpus_hash, "
        "in_corpus, source, confidence, ocrmode, created_at, updated_at) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        ("w1", "corpus_key", "A paper", 2017, "abc123", 1,
         "corpus_paper", 1.0, "force", 0, 0),
    )
    conn.commit()
    conn.close()
    return path


def test_authority_export_and_import_round_trip_ocrmode(tmp_path):
    db = _authority_db(tmp_path / "authority.sqlite")
    exported = export_bibtex(db)
    assert "ocrmode = {force}" in exported

    edited = tmp_path / "edited.bib"
    edited.write_text(exported.replace("ocrmode = {force}", "ocrmode = {redo}"))
    import_bibtex(db, edited)
    conn = sqlite3.connect(db)
    assert conn.execute(
        "SELECT ocrmode FROM works WHERE work_id = 'w1'"
    ).fetchone()[0] == "redo"
    conn.close()


def test_authority_migration_adds_ocrmode_to_an_existing_db(tmp_path):
    db = _authority_db(tmp_path / "old.sqlite")
    conn = sqlite3.connect(db)
    conn.execute("ALTER TABLE works DROP COLUMN ocrmode")
    create_schema(conn)
    have = {row[1] for row in conn.execute("PRAGMA table_info(works)")}
    conn.close()
    assert "ocrmode" in have


def test_dry_run_does_not_persist_the_ocrmode_schema_migration(tmp_path):
    db = _authority_db(tmp_path / "old-dry-run.sqlite")
    conn = sqlite3.connect(db)
    conn.execute("ALTER TABLE works DROP COLUMN ocrmode")
    conn.commit()
    conn.close()
    edit = tmp_path / "force.bib"
    edit.write_text("@article{K, work_id={w1}, ocrmode={force}}")

    result = import_bibtex(db, edit, dry_run=True)
    assert result["changed"] == 1
    conn = sqlite3.connect(db)
    have = {row[1] for row in conn.execute("PRAGMA table_info(works)")}
    conn.close()
    assert "ocrmode" not in have
