"""Per-document OCR language override, and the OSD-precedence fix (#176).

Two behaviors, both about which Tesseract packs a paper gets OCR'd with:

1. `ocrlang` in the paper's bib entry pins the packs outright, overruling
   both langdetect and Tesseract OSD. It is the operator's only recourse
   against a detection that is confident and wrong.
2. A non-Latin OSD verdict no longer discards a high-confidence Latin
   language detection — the two are unioned.

Everything here is unit-level: `_available_tesseract_langs` is patched, so
no Tesseract install is required and the tests run in T0.
"""
from __future__ import annotations

import sqlite3
import time
from pathlib import Path

import pytest

from bib import entry_ocrlang, ocrlang_for_pdf, parse_bibtex
from bib.parser import BibIndex, bib_entry_to_metadata
from pipeline import scan


# Packs a well-provisioned host would have. Deliberately excludes `kor`,
# so the "requested but not installed" path has something to fail on.
INSTALLED = frozenset({
    "eng", "deu", "deu_latf", "fra", "spa", "ita", "pol", "rus", "ell", "lat",
})


@pytest.fixture(autouse=True)
def _installed_packs(monkeypatch):
    monkeypatch.setattr(scan, "_available_tesseract_langs", lambda: INSTALLED)


# ---------------------------------------------------------------------------
# Reading the tag off the bib
# ---------------------------------------------------------------------------


def test_entry_ocrlang_reads_the_field():
    entry = parse_bibtex("@article{K, title={T}, ocrlang={pol+eng}}")[0]
    assert entry_ocrlang(entry) == "pol+eng"


def test_entry_ocrlang_none_when_absent():
    entry = parse_bibtex("@article{K, title={T}}")[0]
    assert entry_ocrlang(entry) is None
    assert entry_ocrlang(None) is None


def test_standard_language_field_is_not_an_ocr_override():
    """`language` must never steer OCR — reference managers emit it by
    default, so honoring it would let an ordinary Zotero export silently
    change how papers are OCR'd."""
    entry = parse_bibtex("@article{K, title={T}, language={polish}}")[0]
    assert entry_ocrlang(entry) is None


def test_ocrlang_for_pdf_matches_on_the_file_field():
    idx = BibIndex(parse_bibtex(
        "@article{A, file={Laska2014.pdf}, ocrlang={pol+eng}}\n"
        "@article{B, file={Other.pdf}}"
    ))
    assert ocrlang_for_pdf(idx, "Laska2014.pdf") == "pol+eng"
    assert ocrlang_for_pdf(idx, "Other.pdf") is None
    assert ocrlang_for_pdf(idx, "Absent.pdf") is None
    assert ocrlang_for_pdf(None, "Laska2014.pdf") is None


def test_ocrlang_reaches_metadata_json():
    entry = parse_bibtex("@article{K, title={T}, ocrlang={ell+eng}}")[0]
    assert bib_entry_to_metadata(entry, "K.pdf")["ocrlang"] == "ell+eng"


# ---------------------------------------------------------------------------
# Parsing and validating the value
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("raw,expected", [
    ("pol+eng", ["pol", "eng"]),
    ("POL+ENG", ["pol", "eng"]),
    ("pol, eng", ["pol", "eng"]),
    ("pol eng", ["pol", "eng"]),
    ("pol+pol+eng", ["pol", "eng"]),   # deduped, order preserved
    ("  pol+eng  ", ["pol", "eng"]),
    ("", []),
    (None, []),
])
def test_parse_ocrlang(raw, expected):
    assert scan._parse_ocrlang(raw) == expected


def test_resolve_pin_keeps_installed_drops_the_rest():
    honored, dropped = scan._resolve_ocrlang_pin("pol+kor+eng")
    assert honored == ["pol", "eng"]
    assert dropped == ["kor"]


def test_resolve_pin_all_missing():
    honored, dropped = scan._resolve_ocrlang_pin("kor+jpn")
    assert honored == []
    assert dropped == ["kor", "jpn"]


# ---------------------------------------------------------------------------
# The OSD-precedence fix
# ---------------------------------------------------------------------------


def test_osd_no_longer_discards_a_confident_detection():
    """Laska2014: Polish/English paper, langdetect `en` at p=1.00, OSD says
    Cyrillic. It used to resolve to ['rus'] alone."""
    packs = scan._resolve_tesseract_packs("en", visual_script="Cyrillic")
    assert packs == ["rus", "eng"]


def test_the_non_english_case_that_was_silently_losing_its_pack():
    """Alvarino 1976b: Spanish, OSD says Greek. The English papers were
    rescued by _compose_ocr_langs appending `eng`; a Spanish one just lost
    `spa` and OCR'd without it."""
    assert scan._resolve_tesseract_packs("es", visual_script="Greek") == ["ell", "spa"]


def test_osd_still_leads():
    """OSD first in the list — it reads the page image, so it stays right
    on exactly the documents langdetect gets wrong."""
    packs = scan._resolve_tesseract_packs("fr", visual_script="Cyrillic")
    assert packs[0] == "rus"


def test_latin_osd_verdict_is_not_a_disagreement():
    """_SCRIPT_TO_TESSERACT has no "Latin" key, so an OSD verdict of Latin
    never forms a union — behavior here is unchanged."""
    assert scan._resolve_tesseract_packs("fr", visual_script="Latin") == ["fra"]


def test_german_still_gets_the_fraktur_companion():
    assert scan._resolve_tesseract_packs("de") == ["deu", "deu_latf"]


def test_no_duplicate_when_osd_and_detection_agree():
    packs = scan._resolve_tesseract_packs("el", visual_script="Greek")
    assert packs == ["ell"]


# ---------------------------------------------------------------------------
# Applying the pin during detection
# ---------------------------------------------------------------------------


def _detected(**over):
    base = {
        "filename": "Laska2014.pdf",
        "file_type": "broken_text_layer",
        "needs_ocr": True,
        "detected_language": "en",
        "language_confidence": 1.0,
        "visual_script": "Cyrillic",
    }
    base.update(over)
    return base


def test_pin_overrides_both_detection_and_osd():
    r = scan._annotate_pack_availability("pol+eng", _detected())
    assert r["tesseract_packs"] == ["pol", "eng"]
    assert r["ocrlang_honored"] is True
    assert r["ocrlang_requested"] == "pol+eng"


def test_pin_leaves_the_detection_record_intact():
    """The override must stay auditable: what the pipeline believed and
    what the operator overruled it with, side by side in one file."""
    r = scan._annotate_pack_availability("pol+eng", _detected())
    assert r["detected_language"] == "en"
    assert r["language_confidence"] == 1.0
    assert r["visual_script"] == "Cyrillic"


def test_pin_may_be_a_single_pack_with_no_eng_appended():
    """Pinning `ell` alone, with no Latin pack competing for the same
    glyphs, is most of the point of the override."""
    r = scan._annotate_pack_availability("ell", _detected())
    assert r["tesseract_packs"] == ["ell"]


def test_partially_installed_pin_is_honored_for_what_survives(caplog):
    r = scan._annotate_pack_availability("kor+eng", _detected())
    assert r["tesseract_packs"] == ["eng"]
    assert r["ocrlang_honored"] is True
    assert r["ocrlang_dropped"] == ["kor"]
    assert "kor" in caplog.text


def test_unusable_pin_falls_back_rather_than_degrading_ocr(caplog):
    """A typo must not produce worse OCR than no tag at all."""
    r = scan._annotate_pack_availability("koreann", _detected())
    assert r["ocrlang_honored"] is False
    assert r["ocrlang_dropped"] == ["koreann"]
    # Fell through to normal resolution — the union from the fix above.
    assert r["tesseract_packs"] == ["rus", "eng"]
    assert "ignoring the override" in caplog.text


def test_no_ocrlang_keys_when_no_tag():
    r = scan._annotate_pack_availability(None, _detected())
    assert "ocrlang_requested" not in r
    assert "ocrlang_honored" not in r


def test_pin_does_not_force_ocr_on_a_born_digital_paper():
    """ocrlang selects packs; it does not decide whether to OCR."""
    r = scan._annotate_pack_availability(
        "ell+eng",
        _detected(file_type="born_digital", needs_ocr=False, visual_script=None),
    )
    assert r["needs_ocr"] is False


# ---------------------------------------------------------------------------
# tesseract_packs must mirror the real ocrmypdf invocation (#176)
# ---------------------------------------------------------------------------


def test_recorded_packs_include_the_eng_suffix():
    """The field used to record targeted resolution only, while
    _compose_ocr_langs appended `eng` on top — so a paper OCR'd with
    `rus+eng` was filed as `['rus']`. That gap is what made the OSD bug
    above look worse than it was."""
    r = scan._annotate_pack_availability(None, _detected(visual_script=None,
                                                         detected_language="fr"))
    assert r["tesseract_packs"] == ["fra", "eng"]


def test_pack_available_still_tracks_targeted_resolution():
    """tesseract_packs is never empty now, so the availability flag has to
    be judged before the fallback union is folded in — otherwise every
    document looks covered."""
    r = scan._annotate_pack_availability(None, _detected(
        detected_language="ko", visual_script=None))
    assert r["tesseract_pack_available"] is False
    assert r["tesseract_packs"]           # fallback union, non-empty


# ---------------------------------------------------------------------------
# What actually reaches ocrmypdf
# ---------------------------------------------------------------------------


def _capture_ocr_cmd(monkeypatch, tmp_path, detection_result):
    """Run prepare_pdf against a stubbed ocrmypdf and return its argv."""
    captured = {}

    class _Result:
        returncode = 0
        stdout = ""
        stderr = ""

    out = tmp_path / "out.pdf"

    def fake_run(cmd, **kw):
        captured["cmd"] = cmd
        # Write the output path this test asked for. Deriving it from
        # `cmd[-1]` created a file named after whatever argument happened to
        # land last — `--jobs 6` put a stray `6` in the repo root, tracked,
        # for three days.
        out.write_bytes(b"%PDF-1.4\n")
        return _Result()

    monkeypatch.setattr(scan.shutil, "which", lambda n: f"/usr/bin/{n}")
    # prepare_pdf runs ocrmypdf through _run_ocr, not subprocess.run, so the
    # Tesseract grandchildren die with it on a timeout (#209, #254).
    monkeypatch.setattr(scan, "_run_ocr", lambda cmd, timeout: fake_run(cmd))
    monkeypatch.setattr(scan, "_report_ocr_page_loss", lambda *a, **k: {})

    src = tmp_path / "in.pdf"
    src.write_bytes(b"%PDF-1.4\n")
    scan.prepare_pdf(src, detection_result, out)
    return captured["cmd"]


def test_pinned_packs_reach_ocrmypdf_verbatim(monkeypatch, tmp_path):
    det = scan._annotate_pack_availability(
        "ell", _detected(ocr_mode="force_ocr"))
    cmd = _capture_ocr_cmd(monkeypatch, tmp_path, det)
    assert cmd[cmd.index("-l") + 1] == "ell"


def test_without_a_pin_the_union_reaches_ocrmypdf(monkeypatch, tmp_path):
    det = scan._annotate_pack_availability(
        None, _detected(ocr_mode="force_ocr"))
    cmd = _capture_ocr_cmd(monkeypatch, tmp_path, det)
    assert cmd[cmd.index("-l") + 1] == "rus+eng"


# ---------------------------------------------------------------------------
# Bib round-trip
# ---------------------------------------------------------------------------


def _works_db(tmp_path: Path) -> Path:
    """Minimal biblio_authority.sqlite built through the real migration."""
    from bib import authority

    db = tmp_path / "biblio_authority.sqlite"
    conn = sqlite3.connect(db)
    authority.create_schema(conn)
    now = time.time()
    conn.execute(
        "INSERT INTO works (work_id, guid_type, title, year, journal, doi, "
        "corpus_hash, in_corpus, source, confidence, ocrlang, created_at, "
        "updated_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
        ("10.1/aaa", "doi", "A paper", 2014, "J", "10.1/aaa", "aaa", 1,
         "corpus_paper", 1.0, "pol+eng", now, now),
    )
    conn.commit()
    conn.close()
    return db


def test_export_round_trips_the_override(tmp_path):
    """export → hand-edit → import must not drop the operator's tag."""
    from bib.export import export_bibtex

    db = _works_db(tmp_path)
    out = export_bibtex(db, documents_dir=None)
    assert "ocrlang = {pol+eng}" in out

    reparsed = parse_bibtex(out)[0]
    assert entry_ocrlang(reparsed) == "pol+eng"


def test_migration_adds_the_column_to_an_older_db(tmp_path):
    """Existing corpuscles must migrate rather than fail on export."""
    from bib import authority

    db = tmp_path / "old.sqlite"
    conn = sqlite3.connect(db)
    authority.create_schema(conn)
    conn.execute("ALTER TABLE works DROP COLUMN ocrlang")
    assert "ocrlang" not in {
        r[1] for r in conn.execute("PRAGMA table_info(works)")}

    authority._migrate_works_columns(conn)
    assert "ocrlang" in {
        r[1] for r in conn.execute("PRAGMA table_info(works)")}
    conn.close()
