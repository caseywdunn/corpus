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
    fps = _expected_fingerprints_for_run(config_fingerprints=None, metadata_fingerprint=None,
        ocrlang=None, ocrmode="force", keeppages=None,
    )
    for stage in _OCR_DEPENDENT_STAGES:
        assert fps[stage]["ocrmode"] == "force", stage


def test_ocrmode_fingerprint_argument_is_required():
    with pytest.raises(TypeError, match="ocrmode"):
        _expected_fingerprints_for_run(config_fingerprints=None, metadata_fingerprint=None, ocrlang=None, keeppages=None)


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


# --- #172 / #266: the OSD verdict reaches pack selection -----------------
#
# The probe renders pages and runs Tesseract OSD on each one to pick its
# probe packs. Those verdicts used to be dropped, so `visual_script` was
# recorded as null on every path and pack resolution fell back to a
# language read off the text layer — which on a legacy CJK font remapped
# into ASCII reports 100% Latin with total confidence.


def _probe(*scripts, votes=(("en", 1.0, 1),)):
    """A stub probe result: these corroborated verdicts, these votes."""
    return scan._ProbeResult(votes=tuple(votes), scripts=tuple(scripts),
                             osd_scripts=tuple(scripts))


@pytest.mark.parametrize("scripts,expected", [
    # A single non-Latin page decides it: a Latin title page in front of a
    # Chinese body is the normal shape of this material.
    (("Han", "Latin"), "Han"),
    (("Latin", "Han"), "Han"),
    # Most-seen non-Latin script wins over a rarer one.
    (("Cyrillic", "Cyrillic", "Japanese"), "Cyrillic"),
    # Fraktur collapses to Latin rather than steering packs to deu_latf,
    # because OSD's Latin/Fraktur call is the noisy one.
    (("Fraktur", "Latin"), "Latin"),
    (("Fraktur",), "Latin"),
    (("Latin",), "Latin"),
    # "OSD failed everywhere" is a different fact from "every page Latin".
    ((None, None), None),
    ((), None),
])
def test_dominant_visual_script(scripts, expected):
    assert scan._dominant_visual_script(list(scripts)) == expected


def test_raster_scan_selects_packs_from_page_images_not_the_text_layer(monkeypatch):
    """#266: Lin & Zhang 1991 — a Chinese scan whose text layer is 100%
    Latin mojibake, which resolved to `eng` and lost its whole body."""
    pdf = _stub_detection(monkeypatch, "legacy ascii glyph stream " * 300)
    monkeypatch.setattr(scan, "_scanned_page_fraction", lambda *_a, **_k: 1.0)
    monkeypatch.setattr(scan, "_detect_language", lambda _text: ("en", 1.0))
    monkeypatch.setattr(scan, "_gibberish_score", lambda _text: 0.295)
    monkeypatch.setattr(scan, "_text_layer_scripts", lambda _text: {"Latin": 1.0})
    monkeypatch.setitem(scan.CONFIG["ocr"], "probe_language_by_ocr", True)
    monkeypatch.setattr(scan, "_probe_language_by_ocr",
                        lambda *_a, **_k: _probe("Han", "Latin"))

    out = scan.detect_scan_type(pdf)

    assert out["detection_reason"] == "raster_page_images"
    assert out["needs_ocr"] is True
    assert out["visual_script"] == "Han"
    assert out["visual_scripts"] == ["Han", "Latin"]
    # The hint has to move with the verdict, or it narrows the fallback
    # union to exactly the family the pages are not written in.
    assert out["script_hint"] == "Han"
    # The point of the fix: Chinese packs, despite langdetect's confident
    # `en` off the mojibake layer.
    assert "chi_sim" in out["tesseract_packs"]


def test_a_latin_scan_keeps_its_language_derived_packs(monkeypatch):
    """The converse: OSD agreeing with the layer must not disturb a
    correct language guess."""
    pdf = _stub_detection(monkeypatch, "ordinary latin prose " * 300)
    monkeypatch.setattr(scan, "_scanned_page_fraction", lambda *_a, **_k: 1.0)
    monkeypatch.setattr(scan, "_detect_language", lambda _text: ("ru", 1.0))
    monkeypatch.setattr(scan, "_gibberish_score", lambda _text: 0.10)
    monkeypatch.setattr(scan, "_text_layer_scripts", lambda _text: {"Latin": 1.0})
    monkeypatch.setitem(scan.CONFIG["ocr"], "probe_language_by_ocr", True)
    monkeypatch.setattr(scan, "_probe_language_by_ocr",
                        lambda *_a, **_k: _probe("Latin", "Latin",
                                                 votes=(("ru", 1.0, 2),)))

    out = scan.detect_scan_type(pdf)

    assert out["visual_script"] == "Latin"
    assert out["script_hint"] == "Latin"
    assert out["tesseract_packs"] == ["rus", "eng"]


def test_no_text_layer_records_the_probes_script_verdict(monkeypatch):
    """With no text layer, OSD is the only script evidence there is."""
    pdf = _stub_detection(monkeypatch, "")
    monkeypatch.setitem(scan.CONFIG["ocr"], "probe_language_by_ocr", True)
    monkeypatch.setattr(scan, "_probe_language_by_ocr",
                        lambda *_a, **_k: _probe("Han", None))

    out = scan.detect_scan_type(pdf)

    assert out["detection_reason"] == "no_text_layer"
    assert out["visual_script"] == "Han"
    assert out["visual_scripts"] == ["Han", None]
    assert "chi_sim" in out["tesseract_packs"]


def test_gibberish_path_asks_the_page_images_for_a_script(monkeypatch):
    """A branch that has just declared the text layer garbage must not
    then select packs from a language read off that garbage alone."""
    pdf = _stub_detection(monkeypatch, "AKAllEMH5I HAYK " * 200)
    monkeypatch.setattr(scan, "_detect_language", lambda _text: ("en", 1.0))
    monkeypatch.setattr(scan, "_gibberish_score", lambda _text: 0.90)
    monkeypatch.setattr(scan, "_text_layer_scripts", lambda _text: {"Latin": 1.0})
    monkeypatch.setattr(scan, "_visual_page_script", lambda _pdf: "Han")

    out = scan.detect_scan_type(pdf)

    assert out["detection_reason"] == "gibberish_score_above_threshold"
    assert out["visual_script"] == "Han"
    assert out["script_hint"] == "Han"
    assert "chi_sim" in out["tesseract_packs"]


@pytest.mark.parametrize("verdict", ["Cyrillic", "Greek", "Thai", "Arabic",
                                     "Bengali", "Devanagari"])
def test_only_cjk_verdicts_may_override_the_text_layer(monkeypatch, verdict):
    """The exclusions are measured, not cautious. Cyrillic and Greek are
    circular — OCRing a Latin page under `rus` transcribes its letters as
    Cyrillic lookalikes, so the confirming characters are manufactured by
    the check. Thai and Arabic are simply wrong here: all 25 Thai and 28
    Arabic verdicts in the sweep were on Latin-script papers."""
    pdf = _stub_detection(monkeypatch, "AKAllEMH5I HAYK " * 200)
    monkeypatch.setattr(scan, "_detect_language", lambda _text: ("en", 1.0))
    monkeypatch.setattr(scan, "_gibberish_score", lambda _text: 0.90)
    monkeypatch.setattr(scan, "_text_layer_scripts", lambda _text: {"Latin": 1.0})
    monkeypatch.setattr(scan, "_visual_page_script", lambda _pdf: verdict)

    out = scan.detect_scan_type(pdf)

    assert out["visual_script"] == "Latin"
    assert out["tesseract_packs"] == ["eng"]
    # Recorded even so: the verdict was made and is part of the evidence.
    assert out["osd_page_scripts"] == [verdict]


def test_one_page_of_many_does_not_decide_a_script(monkeypatch):
    """Boone 1933 is English across 50 pages and yields a single Japanese
    verdict; a real mixed-script volume shows the script on more than one
    sampled page, or on half of a short one."""
    assert scan._dominant_visual_script(
        ["Japanese", "Latin", "Latin", "Latin", "Latin"]) == "Latin"
    # Two pages is evidence.
    assert scan._dominant_visual_script(
        ["Japanese", "Japanese", "Latin", "Latin", "Latin"]) == "Japanese"
    # So is one of two, which is Lin & Zhang 1991's whole paper.
    assert scan._dominant_visual_script(["Han", "Latin"]) == "Han"


@pytest.mark.parametrize("share,expected", [
    (0.577, "Han"),   # Lin & Zhang 1991's Chinese page
    (0.244, "Han"),   # Lindsay 2006's thinnest Japanese body page
    (0.143, "Latin"),  # Boone 1933 p17 — an English page misread
    (0.0, "Latin"),
])
def test_the_confirmation_floor_sits_between_the_measured_populations(
    share, expected, monkeypatch,
):
    monkeypatch.setattr(scan, "_script_char_share", lambda _t, _s: share)
    assert scan._confirm_page_script("Han", "text") == expected


def test_vendor_boilerplate_path_does_not_trust_the_banners_language(monkeypatch):
    """The banner is the scanning vendor's English, so it is the one
    signal on this path guaranteed not to describe the content."""
    marker = scan._VENDOR_BOILERPLATE[0]
    pdf = _stub_detection(monkeypatch, (marker + " ") * 100)
    monkeypatch.setattr(scan, "_visual_page_script", lambda _pdf: "Han")

    out = scan.detect_scan_type(pdf)

    assert out["detection_reason"] == "vendor_boilerplate_only"
    assert out["visual_script"] == "Han"
    assert "chi_sim" in out["tesseract_packs"]


# --- An OSD verdict is only acted on if the page's own OCR bears it out.
#
# OSD is the only signal that survives a corrupt text layer, which is why
# it is asked at all — but on this material it is wrong often and
# confidently. Run over the reference library, the bare check called 424
# of 1,580 Latin-text-layer documents non-Latin: Fewkes 1882a as Thai,
# Alvariño 1964 as Cyrillic, Bigelow & Sears 1939 as Bengali. Acting on
# the verdict alone is the regression `_resolve_tesseract_packs` records,
# where 68 papers lost their correct pack.


def test_script_char_share_measures_the_claimed_script():
    assert scan._script_char_share("第 16 卷 第 4 期", "Han") == 1.0
    assert scan._script_char_share("ordinary latin prose", "Han") == 0.0
    assert scan._script_char_share("Ta6auna cranguk", "Cyrillic") == 0.0
    assert scan._script_char_share("Таблица", "Cyrillic") == 1.0
    # Kana and kanji both count as Japanese: a page of kanji with no kana
    # still needs `jpn`.
    assert scan._script_char_share("相模湾に出現する", "Japanese") == 1.0
    # A script with no range table cannot be corroborated either way.
    assert scan._script_char_share("anything", "Latin") == 0.0
    assert scan._script_char_share("", "Han") == 0.0


def test_a_misfiring_osd_verdict_is_rejected():
    """`tha` transcribing a Latin page emits no Thai codepoints at all,
    so the misfire scores 0.0 rather than merely low."""
    assert scan._confirm_page_script(
        "Thai", "ON THE ACALEPHAE OF THE COAST OF NEW ENGLAND") == "Latin"
    assert scan._confirm_page_script("Cyrillic", "Estudio cuantitativo") == "Latin"


def test_a_corroborated_osd_verdict_survives():
    # A real content page mixes the script with Latin taxon names and
    # authorities, which is why the floor is 0.10 and not higher. Force-
    # OCR'd under their resolved packs, Lin & Zhang 1991's Chinese page
    # scores 0.577 Han and Lindsay 2006's six body pages 0.244-0.469
    # Japanese, against 0.018 and 0.032 for their Latin-only reference
    # pages — so the floor sits in a measured gap here too.
    mixed = "中 国 海域 管 水 母 新 纪录 " + "NEW RECORDS OF SIPHONOPHORES FROM CHINA SEA "
    assert 0.10 <= scan._script_char_share(mixed, "Han") < 0.5
    assert scan._confirm_page_script("Han", mixed) == "Han"


@pytest.mark.parametrize("verdict", ["Latin", "Fraktur", None])
def test_latin_family_verdicts_pass_through_uncorroborated(verdict):
    """Nothing downstream acts on them, and there is no range table to
    check them against."""
    assert scan._confirm_page_script(verdict, "any text at all") == verdict


def test_a_rejected_verdict_stays_visible_in_the_record(monkeypatch):
    """The corroborated script drives packs, but the raw OSD verdict is
    recorded too — a rejected one must not look as though OSD never ran."""
    pdf = _stub_detection(monkeypatch, "ordinary latin prose " * 300)
    monkeypatch.setattr(scan, "_scanned_page_fraction", lambda *_a, **_k: 1.0)
    monkeypatch.setattr(scan, "_detect_language", lambda _text: ("en", 1.0))
    monkeypatch.setattr(scan, "_gibberish_score", lambda _text: 0.10)
    monkeypatch.setattr(scan, "_text_layer_scripts", lambda _text: {"Latin": 1.0})
    monkeypatch.setitem(scan.CONFIG["ocr"], "probe_language_by_ocr", True)
    monkeypatch.setattr(scan, "_probe_language_by_ocr", lambda *_a, **_k: scan._ProbeResult(
        votes=(("en", 1.0, 3),), scripts=("Latin", "Latin", "Latin"),
        osd_scripts=("Thai", "Latin", "Thai"),
    ))

    out = scan.detect_scan_type(pdf)

    assert out["visual_script"] == "Latin"
    assert out["osd_page_scripts"] == ["Thai", "Latin", "Thai"]
    assert out["script_hint"] == "Latin"
    assert out["tesseract_packs"] == ["eng"]


def test_unmappable_text_layer_is_routed_to_ocr(monkeypatch):
    """#266: a font with no ToUnicode table extracts as glyph indices,
    which are not letters — so the gibberish score cannot see them and
    the paper ships as `clean_text_layer` with its body gone. Hunt et al.
    2001 was 84% unmappable and had 896 usable letters in 59,056
    characters."""
    pdf = _stub_detection(monkeypatch, ("\x01\x02\x03\x04 " * 400) + "Lensia")
    monkeypatch.setattr(scan, "_detect_language", lambda _text: ("hu", 1.0))
    monkeypatch.setitem(scan.CONFIG["ocr"], "probe_language_by_ocr", True)
    monkeypatch.setattr(scan, "_probe_language_by_ocr",
                        lambda *_a, **_k: _probe("Latin", votes=(("en", 1.0, 3),)))

    out = scan.detect_scan_type(pdf)

    assert out["detection_reason"] == "unmappable_text_layer"
    assert out["file_type"] == "broken_text_layer"
    assert out["needs_ocr"] is True
    assert out["ocr_mode"] == "force_ocr"
    assert out["unmappable_char_fraction"] > 0.5
    # `hu` came off characters that are definitionally not language.
    assert out["text_layer_language"] == "hu"
    assert out["detected_language"] == "en"
    assert out["tesseract_packs"] == ["eng"]


def test_a_stray_unmappable_glyph_does_not_condemn_a_clean_paper(monkeypatch):
    """Ishikawa et al. 2004 sits at 0.012: a µ and a ≥ that failed to
    map, in an otherwise clean English paper."""
    pdf = _stub_detection(monkeypatch, ("clean english prose here " * 200) + "\x01\x06")
    monkeypatch.setattr(scan, "_detect_language", lambda _text: ("en", 1.0))

    out = scan.detect_scan_type(pdf)

    assert out["detection_reason"] == "clean_text_layer"
    assert out["needs_ocr"] is False
    # Recorded on the clean path too, so the number that cleared the gate
    # is available when a paper looks present and reads empty.
    assert out["unmappable_char_fraction"] < 0.005


def test_the_unmappable_threshold_is_configurable(monkeypatch):
    pdf = _stub_detection(monkeypatch, "\x01\x02ab " * 300)
    monkeypatch.setattr(scan, "_detect_language", lambda _text: ("en", 1.0))
    monkeypatch.setattr(scan, "_gibberish_score", lambda _text: 0.10)
    monkeypatch.setitem(scan.CONFIG["ocr"], "probe_language_by_ocr", True)
    monkeypatch.setattr(scan, "_probe_language_by_ocr",
                        lambda *_a, **_k: _probe("Latin"))
    assert scan.detect_scan_type(pdf)["unmappable_char_fraction"] == 0.5

    monkeypatch.setitem(scan.CONFIG["ocr"], "unmappable_char_max", 0.05)
    assert scan.detect_scan_type(pdf)["detection_reason"] == "unmappable_text_layer"
    monkeypatch.setitem(scan.CONFIG["ocr"], "unmappable_char_max", 0.60)
    assert scan.detect_scan_type(pdf)["detection_reason"] != "unmappable_text_layer"
