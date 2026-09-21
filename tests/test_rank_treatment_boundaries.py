"""Source-captured rank transitions and structural table roles (#320)."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from docling_core.types.doc import DocItemLabel, DoclingDocument

from pipeline.treatment_context import materialize_treatment_context

FIXTURE = Path(__file__).parent / "fixtures/text_integrity/totton_treatment_boundaries.json"


@pytest.mark.parametrize("case_index", [0, 1])
@pytest.mark.parametrize("limit", [40, 2000])
def test_captured_rank_boundaries_through_saved_chunks_and_serving(tmp_path, monkeypatch, case_index, limit):
    import docling.chunking
    from mcpsrv import app
    from mcpsrv.tools.chunks import get_chunks, get_chunks_by_section
    from pipeline.chunking import chunk_text
    from pipeline.table_structure import export_source_markdown
    from tests.test_table_structure import WordTokenizer

    case = json.loads(FIXTURE.read_text())["cases"][case_index]
    doc = DoclingDocument.model_validate(case["document"])
    source = tmp_path / "docling_doc.json"
    doc.save_as_json(source)
    before = source.read_bytes()
    (tmp_path / "text.json").write_text(json.dumps({"text": export_source_markdown(doc)}))
    real = docling.chunking.HybridChunker
    monkeypatch.setattr(docling.chunking, "HybridChunker",
                        lambda **kw: real(tokenizer=WordTokenizer(limit=limit), **kw))
    chunk_text(tmp_path / "text.json", chunks_output=tmp_path / "chunks.json", docling_doc_file=source)
    artifact = json.loads((tmp_path / "chunks.json").read_text())
    assert artifact["chunker"] == "hybrid_chunker"
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(papers={"paper": {"hash_dir": str(tmp_path)}}))
    rows = get_chunks("paper", [c["chunk_id"] for c in artifact["chunks"]])
    if case_index == 0:
        table_ref = case["original_refs"]["#/tables/10"]
        assert doc.tables[0].label.value == "document_index"
        tables = [r for r in rows if any(i["item_ref"] == table_ref for i in r["source_items"])]
        assert tables
        assert all(r["treatment_context"] == {"status": "unknown", "name": None} for r in tables)
        assert all(r["section_class"] is None for r in tables)
        assert all({i["item_ref"] for i in r["source_items"]} == {table_ref} for r in tables)
        assert all("PRAYINAE" in " ".join(r["headings"]) for r in tables)
    else:
        diagnosis = get_chunks_by_section("paper", treatment_name="Prayoides", section_type="diagnosis", limit=50)
        assert diagnosis and "Diagnosis: Prayines" in " ".join(r["text"] for r in diagnosis)
        assert all(r["treatment_context"]["heading_page"] == 129 for r in diagnosis)
        assert all(r["treatment_context"]["rank"] == "genus" for r in diagnosis)
        assert get_chunks_by_section("paper", treatment_name="Rosacea cymbiformis", section_type="diagnosis") == []
        species_ref = case["original_refs"]["#/texts/2049"]
        species = [r for r in rows if any(i["item_ref"] == species_ref for i in r["source_items"])]
        assert species and all(r["treatment_context"]["name"] == "Prayoides intermedia" for r in species)
        contexts = materialize_treatment_context(doc)
        context = contexts[case["original_refs"]["#/texts/2047"]]["treatment_context"]
        assert context["rank"] == "genus"
        assert context["evidence"] == "explicit_rank_name_authority_heading"
    assert source.read_bytes() == before


@pytest.mark.parametrize("heading", ["Sub-family ii: PRAYINAE Chun, 1897", "SUBORDER CALYCOPHORAE",
                                     "Key to genera of Prayinae", "Key 10 GENERA OF PRAYINAE"])
def test_higher_rank_and_key_headings_end_species_context(heading):
    doc = DoclingDocument(name="rank reset")
    doc.add_heading(text="Maresearsia praeclara Totton, 1954")
    doc.add_heading(text=heading)
    body = doc.add_text(label=DocItemLabel.TEXT, text="Diagnosis: source wording.")
    assert materialize_treatment_context(doc)[body.self_ref]["treatment_context"]["name"] is None


@pytest.mark.parametrize("text", ["Genus: Prayoides Leloup 1934 contains two species.",
    "Monotypic genus for Prayoides intermedia Leloup.", "Genus concept observed in 1934",
    "Genus: Prayoides sp. nov.", "Genus: Prayoides Leloup",
    "Genus: Prayoides (sp. nov.)", "Genus: Prayoides (Leloup sp. nov.)"])
def test_prose_or_incomplete_genus_authority_cannot_start_treatment(text):
    doc = DoclingDocument(name="rank negative")
    item = doc.add_text(label=DocItemLabel.TEXT, text=text)
    assert materialize_treatment_context(doc)[item.self_ref]["treatment_context"]["name"] is None


@pytest.mark.parametrize("role", [DocItemLabel.CAPTION, DocItemLabel.PAGE_HEADER, DocItemLabel.FOOTNOTE])
def test_genus_heading_in_nonprose_does_not_change_treatment(role):
    doc = DoclingDocument(name="rank role")
    doc.add_heading(text="Rosacea cymbiformis (Chiaje, 1822)")
    doc.add_text(label=role, text="Genus: PRAYOIDES Leloup 1934")
    body = doc.add_text(label=DocItemLabel.TEXT, text="Original prose continues.")
    assert materialize_treatment_context(doc)[body.self_ref]["treatment_context"]["name"] == "Rosacea cymbiformis"
