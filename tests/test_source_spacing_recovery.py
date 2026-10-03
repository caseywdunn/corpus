"""Original spacing lost on the actual OCR route, plus conservative refusals."""
import copy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import fitz
import pytest

from pipeline import source_spacing_recovery as recovery

BASE = Path(__file__).parent / "fixtures/table_structure/original_spacing"
FIXTURE = json.loads((BASE / "mapstone.json").read_text())


def page(lines):
    return SimpleNamespace(get_text=lambda _, **kwargs: {"blocks": [{"lines": lines}]})


@pytest.mark.parametrize("index", range(3))
def test_actual_original_prepared_observations_find_exact_source_spaces(index):
    case = FIXTURE["cases"][index]
    candidates = recovery._candidates(page(case["source_lines"]), page(case["prepared_lines"]), FIXTURE["producer"])
    expected = case["repair"]
    candidate = next(c for c in candidates if c["original"] == expected["original"])
    assert candidate["replacement"] == expected["replacement"]
    assert "reason" not in candidate and len(candidate["anchors"]) >= 2
    assert candidate["replacement"].replace(" ", "") == candidate["original"]
    assert hashlib.sha256((BASE / case["crop_file"]).read_bytes()).hexdigest() == case["recaptured_crop_sha256"]


@pytest.mark.parametrize("change", ["unchanged", "changed_letter", "moved", "rotated", "no_anchors"])
def test_original_spacing_never_authorizes_changed_letters_or_unaligned_pages(change):
    case = copy.deepcopy(FIXTURE["cases"][1])
    original = page(case["source_lines"])
    if change == "unchanged":
        prepared = original
    else:
        lines = case["prepared_lines"]
        if change == "rotated":
            for line in lines:
                line["dir"] = [0, 1]
        for line in lines:
            for span in line["spans"]:
                for char in span["chars"]:
                    if change == "moved":
                        char["bbox"] = [v+100 if i in (0, 2) else v for i, v in enumerate(char["bbox"])]
                    elif change == "changed_letter" and char["c"] == "d":
                        char["c"] = "b"
                    elif change == "no_anchors" and not fitz.Rect(case["repair"]["prepared_bbox"]).contains(fitz.Rect(char["bbox"])):
                        char["c"] = " "
        prepared = page(lines)
    candidates = recovery._candidates(original, prepared, FIXTURE["producer"])
    assert not any(c["original"] == "alonedeveloped" and "reason" not in c for c in candidates)


def document(tmp_path, *, duplicate=False, repeat=False):
    from docling_core.types.doc import (
        BoundingBox, DocItemLabel, DoclingDocument, ProvenanceItem, Size, TableCell, TableData,
    )
    case = FIXTURE["cases"][1]
    observation = case["observation"]
    doc = DoclingDocument(name="actual OCR key cell")
    doc.add_page(page_no=2, size=Size(width=case["page_size"][0], height=case["page_size"][1]))
    text = observation["text"] * (2 if repeat else 1)
    cell = TableCell(text=text, bbox=BoundingBox.model_validate(observation["bbox"]),
                     start_row_offset_idx=0, end_row_offset_idx=1,
                     start_col_offset_idx=0, end_col_offset_idx=1)
    table = doc.add_table(data=TableData(num_rows=1, num_cols=1, table_cells=[cell]),
                          prov=ProvenanceItem(page_no=2, charspan=(0, len(text)), bbox=cell.bbox))
    if duplicate:
        doc.add_text(label=DocItemLabel.TEXT, text=text,
                     prov=ProvenanceItem(page_no=2, charspan=(0, len(text)), bbox=cell.bbox))
    path = tmp_path / "prepared.pdf"
    path.write_bytes(b"identity of actual prepared observation represented in fixture")
    receipt = {"method": recovery.SOURCE_SPACING_RECOVERY_POLICY, "producer": FIXTURE["producer"],
               "source_pdf_sha256": FIXTURE["source_kept_pdf_sha256"],
               "prepared_pdf_sha256": recovery._digest(path),
               "repairs": [copy.deepcopy(case["repair"])], "unresolved": []}
    return doc, table, path, receipt


def test_actual_key_cell_recovers_source_space_and_preserves_chunk_provenance(tmp_path):
    doc, table, path, receipt = document(tmp_path)
    original = table.data.table_cells[0].text
    result = recovery.apply_source_spacing(doc, path, receipt)
    assert table.data.table_cells[0].text == original.replace("alonedeveloped", "alone developed")
    assert len(result["repairs"]) == 1
    assert result["repairs"][0]["original_page"] == 68
    assert result["repairs"][0]["cell_index"] == 0
    assert not recovery.apply_source_spacing(doc, path, receipt)["repairs"]
    from pipeline.treatment_context import chunk_source_context, materialize_treatment_context
    metadata = chunk_source_context([table], materialize_treatment_context(doc))
    assert metadata["text_integrity"][0]["source_pdf_sha256"] == FIXTURE["source_kept_pdf_sha256"]
    assert metadata["text_integrity"][0]["crop_evidence"]["ocr_outputs"] == receipt["repairs"][0]["crop_evidence"]["ocr_outputs"]


@pytest.mark.parametrize("change", ["identity", "duplicate", "repeat", "moved", "different_page"])
def test_application_requires_unique_owner_and_exact_prepared_identity(tmp_path, change):
    doc, table, path, receipt = document(tmp_path, duplicate=change == "duplicate", repeat=change == "repeat")
    before = table.data.table_cells[0].text
    if change == "identity":
        path.write_bytes(b"different prepared PDF")
    elif change == "moved":
        receipt["repairs"][0]["prepared_bbox"] = [0, 0, 20, 20]
    elif change == "different_page":
        receipt["repairs"][0]["page"] = 1
    result = recovery.apply_source_spacing(doc, path, receipt)
    assert not result["repairs"] and result["unresolved"]
    assert table.data.table_cells[0].text == before


def fake_documents(tmp_path, monkeypatch, *, geometry=False, rotated=False):
    case = FIXTURE["cases"][1]
    source, prepared = tmp_path / "source.pdf", tmp_path / "prepared.pdf"
    source.write_bytes(b"original"); prepared.write_bytes(b"prepared")

    class Page:
        rotation = 90 if rotated else 0
        rect = fitz.Rect(0, 0, *case["page_size"])
        cropbox = rect

        def __init__(self, lines, shifted=False):
            self.lines = lines
            if shifted:
                self.cropbox = self.cropbox + (1, 0, 1, 0)

        def get_text(self, _, **kwargs):
            return {"blocks": [{"lines": self.lines}]}

        def get_pixmap(self, **kwargs):
            return SimpleNamespace(tobytes=lambda _: (BASE / case["crop_file"]).read_bytes())

    class Pdf(list):
        def __enter__(self):
            return self

        def __exit__(self, *_):
            pass

    monkeypatch.setattr(fitz, "open", lambda path: Pdf([Page(
        case["source_lines"] if path == source else case["prepared_lines"], geometry and path == prepared)]))
    monkeypatch.setattr(recovery, "source_spacing_recovery_producer", lambda: copy.deepcopy(FIXTURE["producer"]))
    return source, prepared


@pytest.mark.parametrize("mode", ["confirmed", "disagrees", "ocr_only_split", "budget", "geometry", "rotated", "mapping"])
def test_inspection_retains_real_source_boundaries_and_refusal_receipts(tmp_path, monkeypatch, mode):
    from pipeline import source_spaces
    source, prepared = fake_documents(tmp_path, monkeypatch, geometry=mode == "geometry", rotated=mode == "rotated")
    calls = []
    def ocr(args, **kwargs):
        calls.append(args)
        output = "alone developed"
        if mode == "disagrees" and len(calls) % 2 == 0:
            output = "alonedeveloped"
        elif mode == "ocr_only_split":
            output = "alone devel oped"
        return SimpleNamespace(stdout=output.encode())
    monkeypatch.setattr(source_spaces.subprocess, "run", ocr)
    monkeypatch.setattr(source_spaces.shutil, "which", lambda _: "/fixture/tesseract")
    if mode == "budget":
        producer = {**FIXTURE["producer"], "max_candidates_per_document": 0}
        monkeypatch.setattr(recovery, "source_spacing_recovery_producer", lambda: producer)
    result = recovery.inspect_source_spacing(source, prepared, [68, 69] if mode == "mapping" else [68])
    if mode in {"confirmed", "ocr_only_split"}:
        assert result["repairs"][0]["replacement"] == "alone developed"
        assert result["repairs"][0]["original_page"] == 68
        assert [args[args.index("--psm")+1] for args in calls] == ["6", "7"]
        assert all("whitelist" not in " ".join(args) for args in calls)
    else:
        assert not result["repairs"] and result["unresolved"]
    if mode in {"budget", "geometry", "rotated", "mapping"}:
        assert not calls


def test_spacing_identity_invalidates_preparation_and_consumers_only(monkeypatch):
    from pipeline.build_inputs import config_fingerprints
    before = config_fingerprints({}, panel_mode="off")
    monkeypatch.setattr(recovery, "SOURCE_SPACING_RECOVERY_POLICY", "new-policy")
    after = config_fingerprints({}, panel_mode="off")
    for stage in ("pdf_preparation", "docling_extraction", "metadata_extraction", "text_chunking", "figure_materialization"):
        assert before[stage] != after[stage]
    for stage in ("scan_detection", "huge_document_check"):
        assert before[stage] == after[stage]
