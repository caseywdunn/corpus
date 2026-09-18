"""Source-derived fragments for notation, column order and treatment access.

The library PDFs remain external. Small layout excerpts and three glyph crops
pin the failure mechanisms without adding a built corpuscle fixture.
"""
from copy import deepcopy
import json
import os
from pathlib import Path

import numpy as np
from PIL import Image
import pytest

from pipeline.scientific_text import raster_agrees_with_sign, repair_aligned_text
from pipeline.source_layout import recover_panel_caption_roles, repair_reading_order
from pipeline.treatment_context import materialize_treatment_context

FIXTURES = Path(__file__).parent / "fixtures" / "text_integrity"
CASES = json.loads((FIXTURES / "source_fragments.json").read_text())


def document(case_name):
    from docling_core.types.doc import DoclingDocument, DocItemLabel, ProvenanceItem, Size
    case = CASES[case_name]
    doc = DoclingDocument(name=case_name)
    for number, size in case["pages"].items():
        doc.add_page(page_no=int(number), size=Size(**size))
    for entry in case["items"]:
        kwargs = {"text": entry["text"], "prov": ProvenanceItem.model_validate(entry["prov"][0])}
        if entry["label"] == "section_header":
            item = doc.add_heading(**kwargs)
        else:
            item = doc.add_text(label=DocItemLabel(entry["label"]), **kwargs)
        item.prov = [ProvenanceItem.model_validate(p) for p in entry["prov"]]
    return doc


def evidence(entry):
    result = [{"bbox": [0,0,0,0], "method": None, "verified": False}
              for _ in entry["native_compact"]]
    for item in entry["glyph_evidence"]:
        result[item["offset"]] = item
    return result


def repaired_items(case):
    return [repair_aligned_text(item["text"], item["native_compact"], evidence(item))
            for item in CASES[case]["items"] if "native_compact" in item]


@pytest.mark.parametrize("glyph,agrees", [("sutherland_plusminus", True),
                                        ("hauss_plusminus", True),
                                        ("pakhomov_plusminus", True),
                                        ("kidwai_encoded_plusminus_dash", False)])
def test_rendered_glyph_not_unicode_decides_plusminus(glyph, agrees):
    ink = np.asarray(Image.open(FIXTURES / (glyph+".png"))) < 150
    ys,xs = np.where(ink)
    ink = ink[min(ys):max(ys)+1, min(xs):max(xs)+1]
    assert raster_agrees_with_sign(ink, "±") is agrees


def test_equals_and_plus_are_not_uncertainty_signs():
    equal = np.zeros((20,30), dtype=bool)
    equal[:3,:] = True
    equal[-3:,:] = True
    assert not raster_agrees_with_sign(equal, "±")
    plus = np.zeros((30,30), dtype=bool)
    plus[13:16,:] = True
    plus[:,13:16] = True
    assert not raster_agrees_with_sign(plus, "±")
    assert not raster_agrees_with_sign(plus, "−")


def test_scientific_values_and_raised_units_survive_source_alignment():
    text, edits, unresolved = repaired_items("Sutherland_etal2019b-1-1")[0]
    assert "0.15 ±0.10" in text and "104 ±41 deg. s⁻¹" in text
    assert "1063 ±176 mm s⁻¹" in text
    assert not unresolved
    assert all(e["evidence"] and e["source_bbox"] for e in edits)
    assert "R / L" in text  # Preserve the source's stated ratio.


def test_all_duplicate_equation_copies_preserve_source_minus():
    values = repaired_items("Sutherland_etal2019b-4-4")
    caption = next(text for text, _, _ in values if "0.002" in text)
    assert caption.count("0.002 x −0.0715") == 2
    assert "0.002 x 0.0715" not in caption
    body = values[0][0]
    assert body.count("deg. s⁻¹") == 6
    assert "L / R" in body and "L/R" in body  # Do not correct the source ratio.


def test_volume_exponent_and_scale_bar_micro_unit():
    volume = repaired_items("Hauss_et2016-4-4")[0][0]
    assert "±0.99) m³" in volume
    scale = repaired_items("Haddock_etal2005-1-1")[0][0]
    assert "200 µm" in scale and "1 mm" in scale


def test_native_false_plusminus_date_range_is_unchanged_and_recorded():
    item = CASES["Kidwai_Amjad2000-1-1"]["items"][0]
    text, edits, unresolved = repaired_items("Kidwai_Amjad2000-1-1")[0]
    assert text == item["text"]
    assert "1857-1859" in text and not edits
    assert unresolved[0]["reason"] == "native_glyph_not_confirmed_by_source_raster"


def test_unverified_or_unrelated_native_signs_are_not_applied():
    text = "Abundance 12.5 + 6.1 animals."
    native = "Abundance12.5±6.1animals."
    proof = [{"bbox": [0,0,1,1], "verified": False, "method": "unverified"} for _ in native]
    repaired, edits, unresolved = repair_aligned_text(text, native, proof)
    assert repaired == text and not edits and unresolved
    proof = [{**p, "verified": True} for p in proof]
    repaired, edits, unresolved = repair_aligned_text(text, native, proof)
    assert "12.5 ± 6.1" in repaired and edits and not unresolved
    unrelated = "A different result and specimen."
    assert repair_aligned_text(unrelated, native, proof)[0] == unrelated


def test_notation_repair_is_idempotent():
    for case in CASES.values():
        for item in case["items"]:
            if "native_compact" not in item:
                continue
            once = repair_aligned_text(item["text"], item["native_compact"], evidence(item))[0]
            twice, edits, _ = repair_aligned_text(once, item["native_compact"], evidence(item))
            assert twice == once and not edits


def test_haddock_continuations_stay_in_their_printed_column():
    doc = document("Haddock_etal2005-1-1")
    report = repair_reading_order(doc)
    texts = [ref.resolve(doc).text for ref in doc.body.children]
    first = next(i for i,t in enumerate(texts) if t.endswith("inside two of"))
    assert texts[first+1].startswith("our specimens")
    assert texts[first+3].endswith("yellow to red")
    assert texts[first+4].startswith("(583, 620")
    assert report["reordered_pages"][0]["column_starts"] == pytest.approx([35.8866,213.222856,390.560131])


def test_church_caption_roles_and_minuta_diagnosis_precede_physalis():
    doc = document("Church_etal2025-7-8")
    recovered = recover_panel_caption_roles(doc)
    assert len(recovered) == 5
    report = repair_reading_order(doc)
    assert [p["page"] for p in report["reordered_pages"]] == [8]
    contexts = materialize_treatment_context(doc)
    texts = [ref.resolve(doc) for ref in doc.body.children]
    continuation = next(t for t in texts if t.text.startswith("principal tentacles"))
    next_species = next(t for t in texts if t.text.startswith("Physalia physalis"))
    assert texts.index(continuation) < texts.index(next_species)
    assert contexts[continuation.self_ref]["treatment_context"]["name"] == "Physalia minuta"
    assert contexts[continuation.self_ref]["section_type"] == "diagnosis"
    for item in texts:
        if item.text.startswith("(B) Representative"):
            assert contexts[item.self_ref]["treatment_context"]["status"] == "unknown"


@pytest.mark.parametrize("case,species", [("Siebert_etal2013-4-5", "Apolemia lanosa"),
                                         ("Siebert_etal2013-13-14", "Apolemia rubriversa"),
                                         ("Mapstone2009-213-214", "Muggiaea atlantica")])
def test_inline_and_multipage_diagnoses_have_evidenced_species(case,species):
    doc = document(case)
    contexts = materialize_treatment_context(doc)
    diagnostic = [contexts[item.self_ref] for item in doc.texts if
                  contexts.get(item.self_ref, {}).get("section_type") == "diagnosis"]
    assert diagnostic and all(c["treatment_context"]["name"] == species for c in diagnostic)
    if species == "Muggiaea atlantica":
        assert diagnostic[0]["treatment_context"]["heading"] == "Muggiaea at/antica Cunningham, 1892"
        assert diagnostic[0]["treatment_context"]["name_evidence_ref"]
        assert diagnostic[0]["treatment_context"]["heading_page"] == 213


def test_missing_heading_and_comparative_prose_do_not_imply_species():
    doc = document("Siebert_etal2013-4-5")
    doc.body.children = [r for r in doc.body.children if not r.resolve(doc).text.startswith("Apolemia lanosa")]
    diagnostic = next(t for t in doc.texts if t.text.startswith("Diagnosis."))
    diagnostic.text += " Unlike Apolemia rubriversa, this species lacks diverticula."
    context = materialize_treatment_context(doc)[diagnostic.self_ref]
    assert context["treatment_context"] == {"status": "unknown", "name": None}


def test_fuzzy_heading_alone_is_not_enough_to_classify_diagnosis():
    doc = document("Mapstone2009-213-214")
    heading = next(t for t in doc.texts if t.text == "Diagnosis")
    heading.text = ". ° Diguaosts"
    contexts = materialize_treatment_context(doc)
    body = next(t for t in doc.texts if t.text.startswith("Anterior nectophore"))
    assert contexts[body.self_ref]["section_type"] is None


def test_single_column_and_unconfirmed_columns_keep_input_order():
    doc = document("Siebert_etal2013-4-5")
    before = deepcopy(doc.body.children)
    assert not repair_reading_order(doc)["reordered_pages"]
    assert doc.body.children == before


def test_chunk_output_carries_context_provenance_and_separates_diagnosis(tmp_path, monkeypatch):
    import docling.chunking
    from docling_core.transforms.chunker.tokenizer.base import BaseTokenizer
    from pipeline.chunking import chunk_text
    from pipeline.config import classify_section
    class LocalTokenizer(BaseTokenizer):
        def count_tokens(self,text): return len(text.split())
        def get_max_tokens(self): return 2000
        def get_tokenizer(self): return None
    real = docling.chunking.HybridChunker
    monkeypatch.setattr(docling.chunking, "HybridChunker", lambda **kw: real(tokenizer=LocalTokenizer(), **kw))
    doc = document("Siebert_etal2013-4-5")
    from docling_core.types.doc.common.meta import BaseMeta
    diagnostic_item = next(item for item in doc.texts if item.text.startswith("Diagnosis."))
    diagnostic_item.meta = BaseMeta()
    diagnostic_item.meta.corpus__scientific_text = [{"item_ref": diagnostic_item.self_ref,
                                                  "original": "PROVENANCE_ONLY_SENTINEL", "repairs": []}]
    from pipeline.table_structure import export_source_markdown
    assert "PROVENANCE_ONLY_SENTINEL" not in export_source_markdown(doc)
    source = tmp_path / "docling_doc.json"
    doc.save_as_json(source)
    output = tmp_path / "chunks.json"
    chunk_text(tmp_path / "text.json", chunks_output=output, docling_doc_file=source)
    chunks = json.loads(output.read_text())["chunks"]
    diagnosis = [c for c in chunks if c["section_type"] == "diagnosis"]
    assert len(diagnosis) == 1
    assert diagnosis[0]["treatment_context"]["name"] == "Apolemia lanosa"
    assert diagnosis[0]["section_class"] == "description"
    assert diagnosis[0]["text"].startswith("Diagnosis.")
    assert "Holotype" not in diagnosis[0]["text"]
    assert diagnosis[0]["source_items"][0]["page"] == 5
    assert diagnosis[0]["text_integrity"][0]["original"] == "PROVENANCE_ONLY_SENTINEL"
    assert "PROVENANCE_ONLY_SENTINEL" not in diagnosis[0]["text"]
    assert classify_section(["Re-evaluation of Apolemiidae with description of two new species"]) is None


@pytest.mark.skipif(not os.environ.get("CORPUS_LIBRARY_DIR"), reason="external source PDF replay")
@pytest.mark.parametrize("case", ["Sutherland_etal2019b-1-1", "Sutherland_etal2019b-4-4", "Hauss_et2016-4-4", "Haddock_etal2005-1-1", "Kidwai_Amjad2000-1-1"])
def test_replay_external_source_glyphs(case):
    import hashlib
    from pipeline.scientific_text import prepare_scientific_text
    source = Path(os.environ["CORPUS_LIBRARY_DIR"]) / CASES[case]["pdf"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == CASES[case]["sha256"]
    doc = document(case)
    report = prepare_scientific_text(doc, source)
    joined = "\n".join(t.text for t in doc.texts)
    if case.startswith("Sutherland"):
        assert "s⁻¹" in joined and report["repairs"]
        if "-4-4" in case:
            assert joined.count("0.002 x −0.0715") == 2
    elif case.startswith("Hauss"):
        assert "±0.99) m³" in joined
    elif case.startswith("Haddock"):
        assert "200 µm" in joined
    else:
        assert "1857-1859" in joined and not report["repairs"] and report["unresolved"]


@pytest.mark.skipif(not os.environ.get("CORPUS_LIBRARY_DIR"), reason="external rotated source table")
def test_rotated_table_source_repairs_legacy_plus_but_preserves_correct_cells(monkeypatch):
    from docling_core.types.doc import TableData, TableCell, ProvenanceItem
    from pipeline.scientific_text import prepare_scientific_text
    case = CASES["Pakhomov_etal2000-5-5"]
    doc = document("Pakhomov_etal2000-5-5")
    table = doc.add_table(data=TableData(num_rows=case["table"]["num_rows"],
                                       num_cols=case["table"]["num_cols"],
                                       table_cells=[TableCell.model_validate(c) for c in case["table"]["cells"]]),
                          prov=ProvenanceItem.model_validate(case["table"]["prov"][0]))
    source = Path(os.environ["CORPUS_LIBRARY_DIR"]) / case["pdf"]
    assert not prepare_scientific_text(doc, source)["repairs"]
    for cell in table.data.table_cells:
        cell.text = cell.text.replace("±", "+")  # The historical served corruption.
    report = prepare_scientific_text(doc, source)
    assert [cell.text for cell in table.data.table_cells] == ["12.5 ± 6.1", "35.5 ± 39.7"]
    assert len(report["repairs"]) == 2
    from pipeline import scientific_text
    original_lines = scientific_text._source_lines
    def corrupted_ocr_layer(page):
        lines = original_lines(page)
        for line in lines:
            for char in line:
                if char["text"] == "±":
                    char["text"] = "+"
        return lines
    monkeypatch.setattr(scientific_text, "_source_lines", corrupted_ocr_layer)
    for cell in table.data.table_cells:
        cell.text = cell.text.replace("±", "+")
    report = prepare_scientific_text(doc, source)
    assert len(report["repairs"]) == 2
    assert all(r["evidence"] == ["rendered_plusminus_from_encoded_plus"] for r in report["repairs"])


@pytest.mark.skipif(not os.environ.get("CORPUS_LIBRARY_DIR"), reason="external rendered heading source")
def test_legacy_misread_diagnosis_requires_actual_source_crop():
    from pipeline.treatment_context import recover_section_headings
    case = "Mapstone2009-213-214"
    doc = document(case)
    heading = next(t for t in doc.texts if t.text == "Diagnosis")
    heading.text = ". ° Diguaosts"  # The historical OCR form, not a spelling rule.
    source = Path(os.environ["CORPUS_LIBRARY_DIR"]) / CASES[case]["pdf"]
    result = recover_section_headings(doc, source)
    assert len(result["repairs"]) == 1 and heading.text == "Diagnosis"
    context = materialize_treatment_context(doc)
    body = next(t for t in doc.texts if t.text.startswith("Anterior nectophore"))
    assert context[body.self_ref]["section_type"] == "diagnosis"
    assert context[body.self_ref]["treatment_context"]["name"] == "Muggiaea atlantica"
