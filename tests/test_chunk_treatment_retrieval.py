"""Stored treatment retrieval is distinct from literal mentions (#319)."""
import json
from types import SimpleNamespace

import pytest

from mcpsrv import app
from mcpsrv.tools.chunks import get_chunks, get_chunks_by_section

HASH = "aaaaaaaaaaaa"
OTHER = "bbbbbbbbbbbb"


def row(chunk_id, name, section_type, text, *, status=None, page=5):
    return {
        "chunk_id": chunk_id, "text": text, "section_class": "description",
        "headings": ["Diagnosis"], "figure_refs": [],
        "treatment_context": {"status": status or ("resolved" if name else "unknown"),
                              "name": name, "heading_ref": "#/texts/0", "heading_page": 4,
                              "evidence": "standalone_treatment_heading"},
        "section_type": section_type,
        "source_items": [{"item_ref": f"#/texts/{chunk_id}", "page": page,
                          "bbox": {"l": 1, "t": 9, "r": 8, "b": 2, "coord_origin": "BOTTOMLEFT"},
                          "charspan": [0, len(text)]}],
    }


@pytest.fixture
def corpus(tmp_path, monkeypatch):
    rows = [
        row("c0", None, None, "Author affiliations"),
        row("c1", "Apolemia lanosa", "diagnosis", "Tentilla absent; stem densely covered in papillae."),
        row("c2", "Apolemia rubriversa", "diagnosis", "Unlike Apolemia lanosa, this species has tentilla."),
        row("c3", "Apolemia lanosa", "diagnosis", "Caption, uncertain treatment", status="unknown"),
        row("c4", "Apolemia lanosa", "diagnosis", "A name-free continuation on the next page.", page=6),
        row("c5", "Apolemia lanosa", "description", "Additional descriptive text."),
    ]
    rows[1]["text_integrity"] = [{"method": "source_glyph", "status": "repaired", "before": "x2", "after": "x²"}]
    rows[1]["tables"] = [{"table_ref": "#/tables/0", "complete_rows": False, "row_indices": [2]}]
    rows[1]["key_branches"] = [{"lead_ref": "#/texts/1", "original_destination": "kochi",
                                "status": "geometry_verified_spelling_unverified"}]
    artifact = {"treatment_context_policy": "source_treatments_v1", "chunks": rows}
    papers = {}
    for paper_hash in (HASH, OTHER):
        directory = tmp_path / paper_hash
        directory.mkdir()
        (directory / "chunks.json").write_text(json.dumps(artifact))
        papers[paper_hash] = {"hash_dir": str(directory)}
    # No taxonomy, embeddings, or extraction dependency exists in this index.
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(papers=papers))
    return artifact, tmp_path / HASH / "chunks.json"


@pytest.mark.parametrize("tool", [get_chunks, get_chunks_by_section])
def test_stored_resolved_treatment_filter_does_not_search_literal_names(corpus, tool):
    result = tool(HASH, treatment_name="Apolemia lanosa", section_type="diagnosis")
    assert [r["chunk_id"] for r in result] == ["c1", "c4"]
    assert all("Apolemia lanosa" not in r["text"] for r in result)
    # This is exact stored treatment matching, not taxonomy alias resolution.
    assert tool(HASH, treatment_name="A. lanosa") == []
    assert tool(HASH, treatment_name="apolemia lanosa") == []


def test_filtering_precedes_limit_and_retains_class_semantics(corpus):
    result = get_chunks_by_section(HASH, "description", limit=1,
                                   treatment_name="Apolemia lanosa", section_type="diagnosis")
    assert [r["chunk_id"] for r in result] == ["c1"]
    assert result[0]["section_class"] == "description"
    assert get_chunks_by_section(HASH, "results", treatment_name="Apolemia lanosa") == []
    assert len(get_chunks_by_section(HASH, "description")) == 6


@pytest.mark.parametrize("tool", [get_chunks, get_chunks_by_section])
def test_stored_evidence_passes_through_in_metadata_and_text_modes(corpus, tool):
    artifact, _ = corpus
    for with_text in (True, False):
        rows = tool(HASH, with_text=with_text, treatment_name="Apolemia lanosa", section_type="diagnosis")
        result = rows[0]
        for field in ("treatment_context", "section_type", "source_items", "text_integrity", "tables", "key_branches"):
            assert result[field] == artifact["chunks"][1][field]
        assert ("text" in result) == with_text
        assert ("len_chars" in result) != with_text
        if not with_text:
            assert result["len_chars"] == len(artifact["chunks"][1]["text"])
        assert result["tables"][0]["complete_rows"] is False


def test_id_selection_preserves_paper_order_and_intersects_filters(corpus):
    result = get_chunks(HASH, ["c4", "missing", "c1", "c2"], treatment_name="Apolemia lanosa")
    assert [r["chunk_id"] for r in result] == ["c1", "c4"]
    assert get_chunks(HASH, [], treatment_name="Apolemia lanosa") == []


@pytest.mark.parametrize("tool", [get_chunks, get_chunks_by_section])
def test_unknown_context_is_explicit_and_not_a_resolved_match(corpus, tool):
    result = tool(HASH)
    assert result[0]["treatment_context"]["status"] == "unknown"
    assert result[3]["treatment_context"]["status"] == "unknown"
    # A subtype can be known while the enclosing species remains unknown.
    diagnoses = tool(HASH, section_type="diagnosis")
    assert [r["chunk_id"] for r in diagnoses] == ["c1", "c2", "c3", "c4"]


@pytest.mark.parametrize("tool", [get_chunks, get_chunks_by_section])
def test_legacy_stays_readable_but_context_filters_require_rebuild(corpus, tool):
    artifact, path = corpus
    artifact.pop("treatment_context_policy")
    for chunk in artifact["chunks"]:
        for key in ("treatment_context", "section_type", "source_items", "text_integrity", "tables", "key_branches"):
            chunk.pop(key, None)
    path.write_text(json.dumps(artifact))
    rows = tool(HASH)
    assert len(rows) == 6
    assert rows[0]["treatment_context"] == {"status": "unavailable", "name": None}
    assert rows[0]["section_type"] is None and rows[0]["source_items"] == []
    for filters in ({"treatment_name": "Apolemia lanosa"}, {"section_type": "diagnosis"}):
        result = tool(HASH, **filters)
        assert result[0]["code"] == "rebuild_required"
        assert result[0]["capability"] == "treatment_context"
        assert "Rebuild" in result[0]["error"]


@pytest.mark.parametrize("tool", [get_chunks, get_chunks_by_section])
def test_policy_marker_cannot_hide_missing_chunk_context(corpus, tool):
    artifact, path = corpus
    artifact["chunks"][1].pop("treatment_context")
    path.write_text(json.dumps(artifact))
    assert tool(HASH, section_type="diagnosis")[0]["code"] == "rebuild_required"


@pytest.mark.parametrize("tool", [get_chunks, get_chunks_by_section])
@pytest.mark.parametrize("filters", [{"treatment_name": ""}, {"treatment_name": "   "},
                                      {"section_type": 3}, {"section_type": ""}])
def test_invalid_filters_fail_before_paper_lookup(tool, filters, monkeypatch):
    monkeypatch.setattr(app, "_INDEX", None)
    assert tool("missing", **filters)[0]["code"] == "invalid_argument"


def test_retrieval_leaves_materialized_artifact_unchanged(corpus):
    _, path = corpus
    before = path.read_bytes()
    get_chunks_by_section(HASH, treatment_name="Apolemia lanosa", section_type="diagnosis")
    get_chunks(HASH, with_text=False)
    assert path.read_bytes() == before


@pytest.mark.parametrize("case,name,diagnosis_page", [
    ("Siebert_etal2013-4-5", "Apolemia lanosa", 5),
    ("Siebert_etal2013-13-14", "Apolemia rubriversa", 14),
    ("Mapstone2009-213-214", "Muggiaea atlantica", 214),
])
def test_named_source_diagnosis_is_accessible_through_bounded_workflow(tmp_path, monkeypatch, case, name, diagnosis_page):
    # The separate producer tests cover adjacent treatments, captions and
    # missing headings. Here replay its actual source fixture through build,
    # artifact serialization and the documented MCP route.
    pytest.importorskip("docling")
    import docling.chunking
    from docling_core.transforms.chunker.tokenizer.base import BaseTokenizer
    from pipeline.chunking import chunk_text
    from tests.test_source_text_integrity import document

    class LocalTokenizer(BaseTokenizer):
        def count_tokens(self, text):
            return len(text.split())
        def get_max_tokens(self):
            return 2000
        def get_tokenizer(self):
            return None

    real = docling.chunking.HybridChunker
    monkeypatch.setattr(docling.chunking, "HybridChunker", lambda **kw: real(tokenizer=LocalTokenizer(), **kw))
    doc = document(case)
    source = tmp_path / "docling_doc.json"
    doc.save_as_json(source)
    chunk_text(tmp_path / "text.json", chunks_output=tmp_path / "chunks.json", docling_doc_file=source)
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(papers={HASH: {"hash_dir": str(tmp_path)}}))
    scan = get_chunks_by_section(HASH, section_type="diagnosis", treatment_name=name, limit=1, with_text=False)
    assert len(scan) == 1 and "error" not in scan[0]
    assert "text" not in scan[0]
    result = get_chunks(HASH, [scan[0]["chunk_id"]])[0]
    assert result["treatment_context"]["name"] == name
    assert result["section_type"] == "diagnosis"
    assert result["section_class"] == "description"
    assert any(p["page"] == diagnosis_page for p in result["source_items"])
    assert result["text"] and "affiliation" not in result["text"].lower()
    if name == "Muggiaea atlantica":
        assert result["treatment_context"]["heading_page"] == 213
        assert result["treatment_context"]["name_evidence_ref"]
