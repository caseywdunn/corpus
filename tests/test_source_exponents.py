"""Original scan evidence survives OCR without guessing a unit's power."""
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import fitz
from PIL import Image
import pytest

from pipeline import source_exponents as exponents

FIXTURE = Path(__file__).parent / "fixtures/text_integrity/original_exponents"
CAPTURE = json.loads((FIXTURE / "capture.json").read_text())


def _prepared(tmp_path, monkeypatch):
    from docling_core.types.doc import DocItemLabel, DoclingDocument, ProvenanceItem, Size
    from pipeline import scientific_text
    path = tmp_path / "prepared.pdf"
    with fitz.open() as pdf:
        for _ in range(8):
            pdf.new_page(width=CAPTURE["page_size"][0], height=CAPTURE["page_size"][1])
        pdf.save(path)
    document = DoclingDocument(name="captured-prepared-item")
    document.add_page(page_no=8, size=Size(width=CAPTURE["page_size"][0], height=CAPTURE["page_size"][1]))
    entry = CAPTURE["prepared_item"]
    item = document.add_text(label=DocItemLabel.TEXT, text=entry["text"],
                            prov=ProvenanceItem.model_validate(entry["prov"][0]))
    item.orig = entry["orig"]
    monkeypatch.setattr(scientific_text, "_source_lines", lambda page: deepcopy(CAPTURE["prepared_fragments"]))
    return document, path, deepcopy(CAPTURE["receipt"])


def test_original_capture_pins_both_digit_readings_and_actual_prepared_corruption():
    candidates = CAPTURE["receipt"]["source_exponents"]["candidates"]
    assert [c["status"] for c in candidates] == ["confirmed", "unresolved"]
    assert [c["glyph"]["text"] for c in candidates] == ["3", "3"]
    assert all(not c["glyph"]["raised"] for c in candidates)
    assert [o["text"] for o in candidates[0]["digit_ocr"]["observations"]] == ["3\n", "3\n"]
    assert candidates[1]["digit_ocr"]["observations"][0]["text"] == ""
    assert "2000  мм?" in CAPTURE["prepared_item"]["text"]
    assert "0.3-10 мм'" in CAPTURE["prepared_item"]["text"]
    for crop in CAPTURE["crops"]:
        assert hashlib.sha256((FIXTURE / crop["file"]).read_bytes()).hexdigest() == crop["sha256"]


def test_recorded_raster_geometry_supports_raised_candidates_without_font_flags():
    chars = [c for row in CAPTURE["source_fragments"] for c in row]
    assert [c["prefix"] for c in exponents.source_candidates(chars)] == ["2000мм", "0.3-10мм"]
    for candidate in CAPTURE["receipt"]["source_exponents"]["candidates"]:
        assert exponents.raised_ink_agrees(candidate["source_ink"]["base"]["ink_bbox"],
                                          candidate["source_ink"]["digit"]["ink_bbox"])


@pytest.mark.parametrize("change", ["baseline", "full_size", "affiliation", "year", "far_away", "quantity_other_line"])
def test_candidate_negative_controls_do_not_request_ocr(change):
    chars = deepcopy(CAPTURE["source_fragments"][0])
    base, digit = chars[-2:]
    if change == "baseline":
        digit["origin"][1] = base["origin"][1]
        digit["raised"] = True  # A bogus flag is not source geometry.
    elif change == "full_size":
        digit["size"] = base["size"]
    elif change == "affiliation":
        for c in chars[:-1]:
            c["text"] = "A"
    elif change == "year":
        base["text"] = "9"
    elif change == "quantity_other_line":
        chars[1]["origin"][1] -= 20
    else:
        digit["bbox"][0] += 100
    assert not list(exponents.source_candidates(chars))


@pytest.mark.parametrize("digit", [[11, 10, 13, 15], [11, 12, 13, 17], [11, 5, 13, 15], [30, 5, 32, 10]])
def test_baseline_subscript_full_height_and_distant_ink_refuse(digit):
    assert not exponents.raised_ink_agrees([5, 10, 10, 15], digit)


def test_digit_value_is_taken_from_source_not_mapped_from_cubic_unit():
    chars = deepcopy(CAPTURE["source_fragments"][0])
    chars[-1]["text"] = "2"
    candidate, = exponents.source_candidates(chars)
    assert candidate["glyph"]["text"] == "2"
    assert candidate["prefix"] == "2000мм"


def test_real_capture_changes_only_confirmed_glyph_and_preserves_refusal(tmp_path, monkeypatch):
    from docling_core.types.doc import DoclingDocument
    from pipeline.native_text_recovery import apply_native_text_recovery
    from pipeline.treatment_context import chunk_source_context, materialize_treatment_context
    document, pdf, receipt = _prepared(tmp_path, monkeypatch)
    original = document.texts[0].text
    immutable = deepcopy(receipt)
    report = apply_native_text_recovery(document, pdf, receipt)
    assert len(report["repairs"]) == len(report["unresolved"]) == 1
    assert document.texts[0].text == original[:734]+"³"+original[735:]
    assert "0.3-10 мм'" in document.texts[0].text
    assert document.texts[0].orig == original
    note = report["repairs"][0]
    assert note["source_pdf_sha256"] == CAPTURE["source_sha256"]
    assert note["prepared_pdf_sha256"] == hashlib.sha256(pdf.read_bytes()).hexdigest()
    assert note["source_evidence"]["digit_ocr"]["verified"]
    assert report["unresolved"][0]["reason"] == "digit_ocr_disagreement"
    assert receipt == immutable
    artifact = tmp_path / "docling_doc.json"
    document.save_as_json(artifact)
    reloaded = DoclingDocument.load_from_json(artifact)
    metadata = chunk_source_context(reloaded.texts, materialize_treatment_context(reloaded))
    assert {n["status"] for n in metadata["text_integrity"]} == {"repaired", "unresolved"}
    assert not apply_native_text_recovery(reloaded, pdf, receipt)["repairs"]
    assert len(reloaded.texts[0].meta.corpus__native_text_recovery) == 2


@pytest.mark.parametrize("change", ["different_digit", "different_quantity", "distant_prepared",
                                    "ambiguous_owner", "ambiguous_text", "page_geometry", "unconfirmed"])
def test_conflicting_prepared_or_source_evidence_never_changes_text(tmp_path, monkeypatch, change):
    from pipeline import scientific_text
    document, pdf, receipt = _prepared(tmp_path, monkeypatch)
    prepared = deepcopy(CAPTURE["prepared_fragments"])
    glyph = next(c for c in prepared[0] if c["text"] == "?")
    if change == "different_digit":
        glyph["text"] = "2"
    elif change == "different_quantity":
        next(c for c in prepared[0] if c["text"] == "2")["text"] = "4"
    elif change == "distant_prepared":
        glyph["bbox"] = [400, 40, 405, 48]
    elif change == "ambiguous_owner":
        document.texts.append(deepcopy(document.texts[0]))
    elif change == "ambiguous_text":
        document.texts[0].text += " 2000 мм? газа в сравнении"
    elif change == "page_geometry":
        receipt["source_exponents"]["candidates"][0]["page_size"][0] += 10
    else:
        receipt["source_exponents"]["candidates"][0]["status"] = "unresolved"
    monkeypatch.setattr(scientific_text, "_source_lines", lambda page: prepared)
    before = [item.text for item in document.texts]
    report = exponents.apply_source_exponents(document, pdf, receipt)
    assert not report["repairs"]
    assert [item.text for item in document.texts] == before
    assert report["unresolved"]


def test_ink_boxes_use_captured_rasters_and_stop_at_shared_pixel_budget():
    candidate = CAPTURE["receipt"]["source_exponents"]["candidates"][0]
    for role in ("base", "digit"):
        name = "case1-base.png" if role == "base" else "case1-digit-ink.png"
        glyph = candidate["base" if role == "base" else "glyph"]
        data = (FIXTURE / name).read_bytes()
        page = SimpleNamespace(get_pixmap=lambda **kwargs: SimpleNamespace(
            x=math.floor(glyph["bbox"][0]*600/72), y=math.floor(glyph["bbox"][1]*600/72),
            tobytes=lambda kind: data, pil_image=lambda: Image.open(FIXTURE / name)))
        budget = SimpleNamespace(total_pixels=0, MAX_TOTAL_PIXELS=1_000_000)
        result = exponents._ink_box(page, glyph, budget)
        assert result == candidate["source_ink"][role]
        budget.total_pixels = budget.MAX_TOTAL_PIXELS
        assert exponents._ink_box(page, glyph, budget)["reason"] == "source_ink_pixel_budget_exhausted"


def test_preparation_retains_exponent_only_receipt_without_claiming_new_ocr(tmp_path, monkeypatch):
    from pipeline import native_text_recovery
    from pipeline.scan import prepare_pdf
    _, source, receipt = _prepared(tmp_path, monkeypatch)
    receipt.update(candidate_count=0, confirmed_count=0, regions=[], unresolved=[])
    monkeypatch.setattr(native_text_recovery, "inspect_native_text_regions", lambda *args: receipt)
    output = tmp_path / "copied.pdf"
    result = prepare_pdf(source, {"needs_ocr": False}, output)
    assert result["native_text_recovery"] == receipt
    assert source.read_bytes() == output.read_bytes()


def test_original_exponent_policy_and_ocr_identity_invalidate_preparation_and_consumers(monkeypatch):
    from pipeline.build_inputs import config_fingerprints
    before = config_fingerprints({}, panel_mode="ocr")
    monkeypatch.setattr(exponents, "SOURCE_EXPONENT_POLICY", "changed-exponent-policy")
    policy = config_fingerprints({}, panel_mode="ocr")
    changed = {stage for stage in before if before[stage] != policy[stage]}
    assert {"pdf_preparation", "docling_extraction", "text_chunking", "figure_crossref"} <= changed
    import pipeline.source_spaces
    monkeypatch.setattr(pipeline.source_spaces, "source_spacing_producer", lambda: {
        "available": True, "traineddata_sha256": "different-digit-model"})
    model = config_fingerprints({}, panel_mode="ocr")
    assert model["pdf_preparation"] != policy["pdf_preparation"]


@pytest.mark.parametrize("page_limit", [8, 1])
def test_original_inspection_replays_exact_crop_observations_and_bounds_work(tmp_path, monkeypatch, page_limit):
    from pipeline import pdf_cmap_recovery, scientific_text
    _, path, _ = _prepared(tmp_path, monkeypatch)
    monkeypatch.setattr(scientific_text, "_source_lines", lambda page:
                        deepcopy(CAPTURE["source_fragments"]) if page.number == 7 else [])
    monkeypatch.setattr(exponents, "source_exponent_producer", lambda: CAPTURE["receipt"]["source_exponents"]["producer"])
    monkeypatch.setattr(pdf_cmap_recovery._DigitVerifier, "MAX_PAGE", page_limit)
    monkeypatch.setattr(pdf_cmap_recovery.shutil, "which", lambda name: "/recorded/tesseract")

    def render(page, *, clip, dpi, colorspace):
        crop = next(c for c in CAPTURE["crops"] if c["source_bbox"] == list(clip) and c["dpi"] == dpi)
        raw = (FIXTURE / crop["file"]).read_bytes()
        return SimpleNamespace(x=math.floor(clip.x0*dpi/72), y=math.floor(clip.y0*dpi/72),
                               tobytes=lambda kind: raw, pil_image=lambda: Image.open(FIXTURE / crop["file"]))

    monkeypatch.setattr(fitz.Page, "get_pixmap", render)
    calls = []

    def ocr(command, *, input, capture_output, timeout):
        digest = hashlib.sha256(input).hexdigest()
        record = next(c["digit_ocr"] for c in CAPTURE["receipt"]["source_exponents"]["candidates"]
                      if c["digit_ocr"]["crop_sha256"] == digest)
        mode = int(command[-1])
        reading = next(o for o in record["observations"] if o["mode"] == mode)
        assert command[1:] == ["stdin", "stdout", "-l", "eng", "--psm", str(mode)]
        assert 0 < timeout <= 3
        calls.append((digest, mode))
        return SimpleNamespace(returncode=reading["returncode"], stdout=reading["text"].encode())

    monkeypatch.setattr(pdf_cmap_recovery.subprocess, "run", ocr)
    with fitz.open(path) as pdf:
        result = exponents.inspect_source_exponents(pdf)
    assert result["candidate_count"] == 2 and result["confirmed_count"] == 1
    assert result["candidates"][0]["digit_ocr"] == CAPTURE["receipt"]["source_exponents"]["candidates"][0]["digit_ocr"]
    assert result["candidates"][1]["reason"] == ("digit_ocr_disagreement" if page_limit == 8 else "source_exponent_budget_exhausted")
    assert len(calls) == (3 if page_limit == 8 else 2)
