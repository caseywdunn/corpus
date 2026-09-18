"""Serve receipts with bounded previews; full source observations stay on disk."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mcpsrv import app
from mcpsrv.chunk_context import (
    ContextProjection, MAX_CONTEXT_ROW_BYTES, MAX_OPTIONAL_CONTEXT_BYTES, context_bytes,
)
from mcpsrv.tools.chunks import get_chunks, get_chunks_by_section
from pipeline.text_encoding import repair_region

HASH = "aaaaaaaaaaaa"
CONTEXT_KEYS = {"treatment_context", "section_type", "source_items", "text_integrity",
                "tables", "key_branches", "context_projection"}


def source_chunk():
    fixture = Path(__file__).parent / "fixtures/text_encoding/source_regions.json"
    region = json.loads(fixture.read_text())["regions"][1]
    glyphs = [{"c": c, "bbox": [0, 0, 0, 0]} for c in region["native_text"]]
    repaired, edits, unresolved = repair_region(region["original"], glyphs)
    assert edits and not unresolved and len(region["original"]) > 1000
    return {
        "chunk_id": "c1", "text": repaired[:80], "section_class": "results",
        "treatment_context": {"status": "unknown", "name": None}, "section_type": None,
        "source_items": [{"item_ref": region["item_ref"], "page": region["physical_page"],
                          "bbox": region["bbox"], "charspan": [0, 80]}],
        "text_integrity": [{"item_ref": region["item_ref"], "page": region["physical_page"],
                            "repairs": edits, "unresolved": unresolved,
                            "method": "source-region-encoding-and-accent-geometry-v1"}],
    }


@pytest.mark.parametrize("tool", [get_chunks, get_chunks_by_section])
@pytest.mark.parametrize("with_text,limit", [(False, 32), (True, 96)])
def test_actual_chinese_paragraph_receipt_is_previewed_without_mutating_source(tmp_path, monkeypatch, tool, with_text, limit):
    chunk = source_chunk()
    artifact = tmp_path / "chunks.json"
    artifact.write_text(json.dumps({"chunks": [chunk], "treatment_context_policy": "source_treatments_v1"}))
    original = artifact.read_bytes()
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(papers={HASH: {"hash_dir": str(tmp_path)}}))
    row = tool(HASH, with_text=with_text)[0]
    receipt = row["text_integrity"][0]
    edit = receipt["repairs"][0]
    original_edit = chunk["text_integrity"][0]["repairs"][0]
    assert edit["charspan"] == original_edit["charspan"]
    assert receipt["item_ref"] == chunk["source_items"][0]["item_ref"]
    assert receipt["page"] == 30
    for field in ("original", "replacement"):
        value = original_edit[field]
        assert edit[field] == {"preview": value[:limit], "char_count": len(value),
                               "preview_charspan": [0, limit], "truncated": True}
        assert value not in json.dumps(row, ensure_ascii=False)
    assert ("text" in row) == with_text
    assert row["context_projection"]["truncated"] is True
    assert receipt["repairs_scope"] == {"available": 1, "returned": 1, "truncated": False}
    assert artifact.read_bytes() == original


def test_shared_optional_budget_preserves_every_requested_row(tmp_path, monkeypatch):
    chunk = source_chunk()
    rows = [dict(chunk, chunk_id=f"c{i}") for i in range(300)]
    (tmp_path / "chunks.json").write_text(json.dumps({"chunks": rows}))
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(papers={HASH: {"hash_dir": str(tmp_path)}}))
    results = get_chunks(HASH, with_text=False)
    assert [r["chunk_id"] for r in results] == [r["chunk_id"] for r in rows]
    optional_bytes = sum(context_bytes(r[field]) for r in results
                         for field in ("source_items", "text_integrity", "tables", "key_branches")
                         if r.get(field))
    assert 0 < optional_bytes <= MAX_OPTIONAL_CONTEXT_BYTES
    assert results[-1]["text_integrity"] == []
    for row in results:
        context = {k: v for k, v in row.items() if k in CONTEXT_KEYS}
        assert context_bytes(context) <= MAX_CONTEXT_ROW_BYTES
        scope = row["context_projection"]["records"]["text_integrity"]
        assert scope == {"available": 1, "returned": len(row["text_integrity"]),
                         "truncated": not bool(row["text_integrity"])}


def test_nested_tables_keys_and_multibyte_prose_are_bounded_with_counts():
    cell = {"cell_index": 4, "text": "海洋" * 10000, "source_bbox": [1, 2, 3, 4]}
    chunk = source_chunk()
    chunk["source_items"] *= 20
    chunk["tables"] = [{"table_ref": "#/tables/0", "column_headers": [cell] * 20,
                        "source_cells": [cell] * 20, "complete_rows": False}] * 20
    chunk["key_branches"] = [{"item_ref": "#/texts/2", "original_destination": "海洋" * 10000,
                              "lead_provenance": [cell] * 20}] * 20
    before = copy.deepcopy(chunk)
    row = ContextProjection(with_text=False).project(chunk)
    assert context_bytes(row) <= MAX_CONTEXT_ROW_BYTES
    assert row["context_projection"]["truncated"]
    for key in ("source_items", "tables", "key_branches"):
        assert row["context_projection"]["records"][key]["available"] == 20
        assert len(row[key]) <= 6
    table = row["tables"][0]
    assert table["column_headers_scope"] == {"available": 20, "returned": 6, "truncated": True}
    assert table["column_headers"][0]["text"]["char_count"] == 20000
    assert chunk == before


def test_untrusted_nested_treatment_metadata_cannot_exceed_row_cap():
    huge = {f"field{i}": ["\U0001f9ec" * 5000] * 20 for i in range(50)}
    chunk = {"treatment_context": {k: huge for k in ("status", "name", "heading", "heading_ref")},
             "section_type": huge, "text_integrity": [huge] * 1000}
    row = ContextProjection(with_text=False).project(chunk)
    assert context_bytes(row) <= MAX_CONTEXT_ROW_BYTES
    assert row["context_projection"]["truncated"]
    assert row["context_projection"]["records"]["text_integrity"]["available"] == 1000
