"""Source-page regressions: logical cells, key endpoints and word boundaries."""
from __future__ import annotations

import copy
import json
import re
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("docling_core")
from docling_core.types.doc import DoclingDocument, TableCell, TableData
from docling_core.transforms.chunker.hybrid_chunker import HybridChunker
from docling_core.transforms.chunker.tokenizer.base import BaseTokenizer

from pipeline import source_spaces
from pipeline.table_structure import (
    SourceTableSerializerProvider, export_source_markdown, logical_rows,
    prepare_table_structure, restore_source_spaces, serialized_rows, table_chunk_metadata,
)
from pipeline.taxa import extract_lexicon_mentions

FIXTURES = Path(__file__).parent / "fixtures" / "table_structure"


class WordTokenizer(BaseTokenizer):
    limit: int = 80

    def count_tokens(self, text):
        return len(text.split())

    def get_max_tokens(self):
        return self.limit

    def get_tokenizer(self):
        return self.count_tokens


def load(name):
    fixture = json.loads((FIXTURES / f"{name}.json").read_text())
    return DoclingDocument.model_validate(fixture["document"]), fixture


def fake_source(monkeypatch, fixture, *, ocr=False, outputs=None):
    import fitz

    class Page:
        def __init__(self, value):
            self.value = value
            self.rect = fitz.Rect(0, 0, value["width"], value["height"])

        def get_text(self, mode, **kwargs):
            return {"blocks": [{"lines": self.value.get("raw_lines", [])
                                if mode == "rawdict" else self.value["lines"]}]}

        def get_pixmap(self, **kwargs):
            return SimpleNamespace(tobytes=lambda fmt: b"captured-test-crop")

    class Pdf:
        def __enter__(self):
            return self

        def __exit__(self, *_):
            pass

        def __len__(self):
            return max(map(int, fixture["pdf_pages"]))

        def __getitem__(self, index):
            return Page(fixture["pdf_pages"][str(index + 1)])

    monkeypatch.setattr(fitz, "open", lambda *_: Pdf())
    monkeypatch.setattr(source_spaces, "source_spacing_producer", lambda: {
        "available": ocr, "policy": source_spaces.SPACE_POLICY, "dpi": 400,
        "max_candidates_per_page": 8, "timeout_seconds": 15,
        "ocr_modes": [6, 7], "language": "eng", "traineddata_sha256": "fixture-producer"})
    if outputs is not None:
        values = iter(outputs)
        monkeypatch.setattr(source_spaces.shutil, "which", lambda _: "/fixture/tesseract")
        monkeypatch.setattr(source_spaces.subprocess, "run", lambda *a, **k: SimpleNamespace(
            stdout=next(values).encode()))


def chunks(document, limit=80):
    return list(HybridChunker(tokenizer=WordTokenizer(limit=limit), merge_peers=False,
                             serializer_provider=SourceTableSerializerProvider()).chunk(document))


def test_merged_heading_is_one_source_occurrence_and_span_survives():
    doc, _ = load("Daniel1974-20-20")
    old = doc.export_to_markdown()
    assert old.count("Genus ?Epibulia") == 4  # Actual upstream grid expansion.
    markdown = export_source_markdown(doc)
    assert markdown.count("Genus ?Epibulia") == 1
    assert 'colspan="4"' in markdown
    items = chunks(doc)
    text = "\n".join(c.text for c in items)
    assert text.count("Genus ?Epibulia") == 1
    assert ", 1 =" not in text and ", 2 =" not in text
    terms = extract_lexicon_mentions([{"chunk_id": str(i), "text": c.text}
                                     for i, c in enumerate(items)], {"Epibulia": {"synonyms": []}})
    source_count = sum(len(re.findall(r"\bEpibulia\b", c.text))
                       for c in doc.tables[0].data.table_cells)
    assert terms["total_mentions"] == source_count
    metadata = table_chunk_metadata(items[0].meta.doc_items, items[0].text, doc)
    assert any(c["col_span"] == 4 for c in metadata[0]["source_cells"])


def test_equal_text_in_distinct_cells_is_not_deduplicated():
    doc = DoclingDocument(name="repeated observations")
    table = doc.add_table(data=TableData(num_rows=2, num_cols=2, table_cells=[
        TableCell(text="same species", start_row_offset_idx=r, end_row_offset_idx=r+1,
                  start_col_offset_idx=c, end_col_offset_idx=c+1)
        for r in range(2) for c in range(2)]))
    assert export_source_markdown(doc).count("same species") == 4
    assert len([cell for _, cells in logical_rows(table, doc) for cell in cells]) == 4


@pytest.mark.parametrize("name", ["Pugh_Haddock2016-42-42", "Pugh_Haddock2016-43-44"])
def test_key_branches_stay_with_targets_in_bounded_chunks(name):
    doc, _ = load(name)
    rows, is_key = serialized_rows(doc.tables[0], doc)
    assert is_key
    pieces = chunks(doc, limit=80)
    assert len(pieces) > 1
    for _, branch, _ in rows:
        assert sum(branch in piece.text for piece in pieces) == 1
    assert all(len(piece.text.split()) <= 80 for piece in pieces)
    for piece in pieces:
        metadata = table_chunk_metadata(piece.meta.doc_items, piece.text, doc)[0]
        assert metadata["kind"] == "identification_key"
        assert metadata["complete_rows"]
        assert metadata["continuation"]
    all_text = "\n".join(p.text for p in pieces)
    assert ", 1 =" not in all_text
    if "43-44" in name:
        assert any("Largest bracts > 35 mm" in text and text.endswith("2") for _, text, _ in rows)
        assert any("Largest bracts < 35 mm" in text and text.endswith("3") for _, text, _ in rows)
        assert any("complete transverse ridge" in text and "Erenna cornuta Pugh, 2001" in text
                   for _, text, _ in rows)
    else:
        assert any('Long "horn" canals' in text and "Erenna cornuta" in text and "Pugh, 2001" in text
                   for _, text, _ in rows)


def test_oversized_row_is_explicitly_partial():
    doc, _ = load("Pugh_Haddock2016-43-44")
    pieces = chunks(doc, limit=12)
    assert any(not table_chunk_metadata(p.meta.doc_items, p.text, doc)[0]["complete_rows"]
               for p in pieces)


def test_numeric_character_table_is_not_misclassified_as_key():
    doc, _ = load("Hosiaetal2024-5-5")
    source = [cell.text for cell in doc.tables[0].data.table_cells]
    rows, key = serialized_rows(doc.tables[0], doc)
    assert not key
    assert [cell["text"] for _, _, cells in rows for cell in cells] == [
        cell["text"] for _, cells in logical_rows(doc.tables[0], doc) for cell in cells]
    text = "\n".join(p.text for p in chunks(doc, limit=80))
    assert all(value.strip() in text for value in source if value.strip())


def test_geometric_endpoints_keep_four_key_relationships(monkeypatch):
    doc, fixture = load("Daniel1985-271-272")
    fake_source(monkeypatch, fixture)
    before = copy.deepcopy(doc)
    report = prepare_table_structure(doc, "source.pdf")
    associations = report["key_associations"]
    assert len(associations) == 4
    expected = {"sharply conical": "kochi", "wall nearly": "delsmani",
                "beyond  apex": "atlantica", "sausage-shaped": "hargmannae"}
    for phrase, target in expected.items():
        branch = next(x for x in associations if phrase in x["original_lead"])
        assert branch["original_destination"] == target
        item = next(t for t in doc.texts if t.self_ref == branch["lead_ref"])
        assert item.text.endswith(target)
        assert len(item.prov) == 2
        assert getattr(item.meta, "corpus__key_branch")["original_lead"] in [t.text for t in before.texts]
    body = export_source_markdown(doc)
    assert "conical.. kochi" in body and "hydroecium.. . atlantica" in body
    assert "hydroecium.. . hargmannae" in body
    # Source spelling is not inferred from the bibliography or a b/h rule.
    assert all(x["status"] == "geometry_verified_spelling_unverified" for x in associations)
    snapshot = doc.model_dump(mode="json")
    again = prepare_table_structure(doc, "source.pdf")
    assert again["key_associations"] == []
    assert doc.model_dump(mode="json") == snapshot


def test_nearby_endpoint_outside_key_is_not_bound(monkeypatch):
    doc, fixture = load("Daniel1985-271-272")
    doc.texts[0].text = "Descriptive prose"
    fake_source(monkeypatch, fixture)
    assert prepare_table_structure(doc, "source.pdf")["key_associations"] == []


def test_actual_hissmann_source_spaces_restore_independent_anatomy_mentions(monkeypatch):
    doc, fixture = load("Hissmann2005-7-7")
    fake_source(monkeypatch, fixture)
    original = doc.texts[0].text
    lexicon = {"nectophore": {"synonyms": ["nectophores"]}}
    before = extract_lexicon_mentions([{"text": original}], lexicon)["total_mentions"]
    report = prepare_table_structure(doc, "source.pdf")
    repaired = doc.texts[0].text
    assert "The holotype possessed nine nectophores and about nine" in repaired
    assert extract_lexicon_mentions([{"text": repaired}], lexicon)["total_mentions"] == before + 1
    proof = report["space_repairs"][0]
    assert proof["original"] == original
    assert proof["evidence"][0]["route"] == "pdf_line_exact_letters"


def test_mapstone_legacy_space_loss_and_clean_regeneration_agree(monkeypatch):
    clean, fixture = load("Mapstone2009-68-68")
    legacy = copy.deepcopy(clean)
    cell = next(c for c in legacy.tables[0].data.table_cells if "Anterior nectophore alone" in c.text)
    cell.text = cell.text.replace("Anterior nectophore alone developed", "Anteriornectophorealonedeveloped")
    fake_source(monkeypatch, fixture)
    report = prepare_table_structure(legacy, "source.pdf")
    assert report["space_repairs"]
    assert export_source_markdown(legacy) == export_source_markdown(clean)
    assert [p.text for p in chunks(legacy)] == [p.text for p in chunks(clean)]


def test_printed_mapstone_run_is_preserved_as_negative_control(monkeypatch):
    doc, fixture = load("Mapstone2009-200-200")
    fake_source(monkeypatch, fixture)
    original = doc.texts[0].text
    report = prepare_table_structure(doc, "source.pdf")
    assert doc.texts[0].text == original
    assert report["space_repairs"] == []
    assert report["source_gap_observations"] == []
    assert report["unresolved_long_runs"]


def test_duclos_geometry_and_ocr_intersection_ignores_ocr_only_split(monkeypatch):
    doc, fixture = load("DuClos_etal2022-6-6")
    outputs = ["Siphosome length scaled approximately linearly with the nectophore cou nt(L =",
               "Siphosome length scaled approximately linearly with the nectophore count (L ="]
    fake_source(monkeypatch, fixture, ocr=True, outputs=outputs)
    report = prepare_table_structure(doc, "source.pdf")
    assert "Siphosome length scaled approximately linearly with the nectophore count" in doc.texts[0].text
    assert "cou nt" not in doc.texts[0].text
    assert report["source_gap_observations"][0]["ocr_outputs"] == outputs
    assert report["space_repairs"][0]["evidence"][0]["route"] == source_spaces.SPACE_POLICY


@pytest.mark.parametrize("outputs", [None, ["incorrect letters", "incorrect letters"]])
def test_missing_or_disagreeing_ocr_leaves_original_and_reports_reason(monkeypatch, outputs):
    doc, fixture = load("DuClos_etal2022-6-6")
    fake_source(monkeypatch, fixture, ocr=outputs is not None, outputs=outputs)
    original = doc.texts[0].text
    report = prepare_table_structure(doc, "source.pdf")
    assert doc.texts[0].text == original
    assert report["source_gap_observations"][0]["status"] == (
        "ocr_disagreement" if outputs else "ocr_unavailable")


def test_actual_mixed_font_source_retains_disagreeing_ocr_and_taxon(monkeypatch):
    doc, fixture = load("DuClos_etal2022-6-mixed-font")
    fake_source(monkeypatch, fixture, ocr=True, outputs=fixture["ocr_observation"]["outputs"])
    original = doc.texts[0].text
    report = prepare_table_structure(doc, "source.pdf")
    observation = report["source_gap_observations"][0]
    assert observation["status"] == "ocr_disagreement"
    assert observation["accepted"] == []
    assert observation["isolated_font_runs"] is True
    assert observation["choices"][0]["boundaries"] == [4, 13, 15, 21, 27, 39]
    assert observation["choices"][0]["start"] == 23
    assert observation["text"].startswith("offree-swimmingN.bijugawere")
    assert observation["ocr_bbox"][0] > observation["bbox"][0]
    assert observation["crop_bbox"] == fixture["ocr_observation"]["crop_bbox_top_left"]
    assert report["unresolved_long_runs"][0]["runs"] == [
        "werecollectedatFridayHarborLaboratoriesbetween"]
    assert doc.texts[0].text == original
    assert "free-swimming N.bijuga werecollected" in original
    assert report["space_repairs"] == []
    reloaded = DoclingDocument.model_validate_json(doc.model_dump_json())
    assert "N.bijuga werecollectedatFridayHarborLaboratoriesbetween" in " ".join(
        chunk.text for chunk in chunks(reloaded))


@pytest.mark.parametrize("existing_boundary", [True, False])
def test_controlled_subrun_agreement_never_adds_a_font_transition_space(monkeypatch, existing_boundary):
    """Controlled agreement tests application; actual G1 OCR above disagrees."""
    doc, fixture = load("DuClos_etal2022-6-mixed-font")
    if not existing_boundary:
        doc.texts[0].text = doc.texts[0].text.replace("N.bijuga were", "N.bijugawere")
    original = doc.texts[0].text
    phrase = "were collected at Friday Harbor Laboratories between"
    fake_source(monkeypatch, fixture, ocr=True, outputs=[phrase, phrase])
    report = prepare_table_structure(doc, "source.pdf")
    assert doc.texts[0].text == original.replace(phrase.replace(" ", ""), phrase)
    expected_taxon_boundary = "N.bijuga were collected" if existing_boundary else "N.bijugawere collected"
    assert expected_taxon_boundary in doc.texts[0].text
    assert "".join(original.split()) == "".join(doc.texts[0].text.split())
    assert report["source_gap_observations"][0]["status"] == "verified"
    proof = report["space_repairs"][0]["evidence"][0]
    assert proof["route"] == source_spaces.SPACE_POLICY
    assert proof["source_gap_evidence"][0]["original"] == phrase.replace(" ", "")
    assert proof["source_gap_evidence"][0]["source_charspan"] == [23, 69]
    assert doc.texts[0].meta.corpus__source_spaces[0]["original"] == original
    reloaded = DoclingDocument.model_validate_json(doc.model_dump_json())
    assert reloaded.texts[0].text == doc.texts[0].text
    assert reloaded.texts[0].meta.corpus__source_spaces[0]["original"] == original
    if existing_boundary:
        assert prepare_table_structure(reloaded, "source.pdf")["space_repairs"] == []


def _styled_line(parts, *, transition_gap=0):
    """Synthetic native glyphs; spaces describe geometry, never literal text."""
    spans = []
    x = 0
    for index, words in enumerate(parts):
        chars = []
        for char in words:
            if char == " ":
                x += 2
            else:
                chars.append({"c": char, "bbox": [x, 0, x + 5, 10]})
                x += 5
        spans.append({"font": f"style-{index}", "size": 10, "chars": chars})
        x += transition_gap
    return {"dir": [1, 0], "bbox": [0, 0, x, 10], "spans": spans}


def test_subrun_provenance_links_only_overlapping_source_observations(monkeypatch):
    doc, fixture = load("DuClos_etal2022-6-mixed-font")
    other = copy.deepcopy(fixture["pdf_pages"]["6"]["raw_lines"][0])
    for box in [other["bbox"], *[s["bbox"] for s in other["spans"]],
                *[c["bbox"] for s in other["spans"] for c in s["chars"]]]:
        box[1] -= 200
        box[3] -= 200
    fixture["pdf_pages"]["6"]["raw_lines"].append(other)
    phrase = "were collected at Friday Harbor Laboratories between"
    fake_source(monkeypatch, fixture, ocr=True, outputs=[phrase] * 4)
    report = prepare_table_structure(doc, "source.pdf")
    assert len(report["source_gap_observations"]) == 2
    assert len(report["space_repairs"][0]["evidence"][0]["source_gap_evidence"]) == 1


@pytest.mark.parametrize("parts", [
    ["Nesselzellkapsel", "wandverdickung"],
    ["electroencephalograph", "ically"],
    ["ab cd ef gh", "ij kl mn op", "qr st uv wx"],
])
def test_font_transitions_and_short_styled_parts_are_not_word_boundaries(parts):
    assert source_spaces.geometric_space_candidate(_styled_line(parts, transition_gap=2)) is None


def test_multiple_uniform_subruns_share_one_crop_and_candidate_budget(monkeypatch):
    line = _styled_line(["firstlong another thirdword finalword", "secondlong different fourthword lastword"])
    candidate = source_spaces.geometric_space_candidate(line)
    assert len(candidate["choices"]) == 2
    assert candidate["isolated_font_runs"] is True
    # Exact letters and supported internal gaps; the style junction stays joined.
    output = "".join(choice["candidate"] for choice in candidate["choices"])
    calls = []
    monkeypatch.setattr(source_spaces.shutil, "which", lambda _: "tesseract-fixture")
    def ocr(args, **kwargs):
        calls.append(args)
        return SimpleNamespace(stdout=output.encode())
    monkeypatch.setattr(source_spaces.subprocess, "run", ocr)
    import fitz
    crops = []
    def pixmap(**kwargs):
        crops.append(kwargs["clip"])
        return SimpleNamespace(tobytes=lambda _: b"synthetic-region")
    page = SimpleNamespace(rect=fitz.Rect(0, 0, 1000, 20), get_pixmap=pixmap)
    producer = {"available": True, "max_candidates_per_page": 1, "dpi": 400,
                "ocr_modes": [6, 7], "language": "eng", "timeout_seconds": 15}
    lines, observations = source_spaces.corroborate_source_gaps(page, [candidate, candidate], producer)
    assert len(calls) == 2 and len(crops) == 1
    assert lines == [output]
    assert len(observations[0]["accepted"]) == 2
    assert observations[1]["status"] == "candidate_budget_exceeded"
    assert "finalwordsecondlong" in lines[0]


@pytest.mark.parametrize("text", ["Nesselzellkapselwandverdickung", "nectophoral", "nectophores",
                                  "Polymorphismehydrozoaire", "межклеточноевзаимодействие"])
def test_multilingual_compounds_and_morphology_are_not_dictionary_split(text):
    assert restore_source_spaces(text, [text]) == (text, [])


def test_uniform_letter_spacing_is_not_word_boundary_evidence():
    text = "Nesselzellkapselwandverdickung"
    chars = [{"c": c, "bbox": [i * 5.5, 0, i * 5.5 + 5, 10]} for i, c in enumerate(text)]
    line = {"dir": [1, 0], "bbox": [0, 0, len(text)*5.5, 10],
            "spans": [{"font": "fixture", "size": 10, "chars": chars}]}
    assert source_spaces.geometric_space_candidate(line) is None


def test_ocr_failure_is_reviewable_and_does_not_change_native_text(monkeypatch):
    doc, fixture = load("DuClos_etal2022-6-6")
    fake_source(monkeypatch, fixture, ocr=True)
    original = doc.texts[0].text
    def timeout(*args, **kwargs):
        raise source_spaces.subprocess.TimeoutExpired("tesseract", 15)
    monkeypatch.setattr(source_spaces.shutil, "which", lambda _: "/fixture/tesseract")
    monkeypatch.setattr(source_spaces.subprocess, "run", timeout)
    report = prepare_table_structure(doc, "source.pdf")
    assert doc.texts[0].text == original
    assert report["source_gap_observations"][0]["status"] == "ocr_failed"
    assert report["source_gap_observations"][0]["error"] == "TimeoutExpired"


def test_ocr_producer_identity_changes_with_traineddata(monkeypatch, tmp_path):
    executable = tmp_path / "tesseract"
    executable.write_text("fixture executable")
    model = tmp_path / "eng.traineddata"
    model.write_bytes(b"model one")
    monkeypatch.setattr(source_spaces.shutil, "which", lambda _: str(executable))
    monkeypatch.setattr(source_spaces, "_ocr_installation", lambda *args: ("tesseract fixture", str(model)))
    first = source_spaces.source_spacing_producer()
    model.write_bytes(b"a different model")
    second = source_spaces.source_spacing_producer()
    assert first["available"] and second["available"]
    assert first["traineddata_sha256"] != second["traineddata_sha256"]


def test_restore_never_changes_letters_or_ambiguous_boundaries():
    token = "Anteriornectophoredeveloped"
    assert restore_source_spaces(token, ["Anterior nectophore developed"])[0] == "Anterior nectophore developed"
    assert restore_source_spaces(token, ["Anterior nectophore deleted"])[0] == token
    assert restore_source_spaces(token, ["Anterior nectophore developed", "Anteriornectophore developed"])[0] == token


def test_rowspan_metadata_and_html_keep_one_logical_cell():
    doc = DoclingDocument(name="vertical span")
    table = doc.add_table(data=TableData(num_rows=2, num_cols=2, table_cells=[
        TableCell(text="shared taxon", row_span=2, start_row_offset_idx=0, end_row_offset_idx=2,
                  start_col_offset_idx=0, end_col_offset_idx=1),
        TableCell(text="observation A", start_row_offset_idx=0, end_row_offset_idx=1,
                  start_col_offset_idx=1, end_col_offset_idx=2),
        TableCell(text="observation B", start_row_offset_idx=1, end_row_offset_idx=2,
                  start_col_offset_idx=1, end_col_offset_idx=2)]))
    assert export_source_markdown(doc).count("shared taxon") == 1
    assert 'rowspan="2"' in export_source_markdown(doc)
    text = "\n".join(p.text for p in chunks(doc))
    assert text.count("shared taxon") == 1
    assert table_chunk_metadata([table], text, doc)[0]["source_cells"][0]["row_span"] == 2


def test_key_binding_and_raw_observations_survive_save_reload(monkeypatch, tmp_path):
    doc, fixture = load("Daniel1985-271-272")
    fake_source(monkeypatch, fixture)
    prepare_table_structure(doc, "source.pdf")
    expected = export_source_markdown(doc)
    path = tmp_path / "docling.json"
    doc.save_as_json(path)
    restored = DoclingDocument.load_from_json(path)
    assert export_source_markdown(restored) == expected
    assert len([t for t in restored.texts if getattr(t.meta, "corpus__key_branch", None)]) == 4
    assert len([t for t in restored.texts if getattr(t.meta, "corpus__key_destination", None)]) == 4
    assert "original_destination" not in expected  # metadata never becomes source prose
