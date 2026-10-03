"""Captured Daniel key geometry → saved document → real chunks → MCP routes.

Only the tokenizer is substituted with a deterministic local word counter.
Source PDF extraction/OCR is outside this compact acceptance boundary.
"""
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import pytest
from docling_core.transforms.chunker.tokenizer.base import BaseTokenizer
from docling_core.types.doc import DoclingDocument

from mcpsrv import app
from mcpsrv.chunk_context import MAX_CONTEXT_ROW_BYTES, context_bytes
from mcpsrv.indexes import CorpusIndex
from mcpsrv.tools.chunks import get_chunks, get_chunks_by_section
from pipeline.chunking import chunk_text
from pipeline.table_structure import _attach_key_destinations, export_source_markdown


FIXTURE = Path(__file__).parent / "fixtures/table_structure"
HASH = "b4146be29447"
EXPECTED = {"sharply conical": "kochi", "wall nearly": "delsmani",
            "beyond  apex": "atlantica", "sausage-shaped": "hargmannae"}
CONTEXT_FIELDS = {"treatment_context", "section_type", "source_items", "text_integrity",
                  "tables", "key_branches", "context_projection"}


class WordTokenizer(BaseTokenizer):
    budget: int

    def count_tokens(self, text):
        return len(text.split())

    def get_max_tokens(self):
        return self.budget

    def get_tokenizer(self):
        return self.count_tokens


def materialize(tmp_path, monkeypatch, budget, *, associate=True):
    import docling.chunking

    captured = json.loads((FIXTURE / "Daniel1985-271-272.json").read_text())
    assert captured["source"]["sha256"].startswith(HASH)
    crop = captured["rendered_evidence"]
    assert hashlib.sha256((FIXTURE / crop["filename"]).read_bytes()).hexdigest() == crop["sha256"]
    doc = DoclingDocument.model_validate(captured["document"])
    # This is the actual geometry producer called by prepare_table_structure.
    # Spacing/OCR helpers are not needed to associate these captured key items.
    associations = _attach_key_destinations(doc) if associate else []
    directory = tmp_path / "documents" / HASH
    directory.mkdir(parents=True)
    source = directory / "docling_doc.json"
    doc.save_as_json(source)
    restored = DoclingDocument.load_from_json(source)
    assert export_source_markdown(restored) == export_source_markdown(doc)
    (directory / "text.json").write_text(json.dumps({"text": export_source_markdown(restored)}))
    (directory / "metadata.json").write_text(json.dumps({"filename": "Daniel1985.pdf",
                                                         "title": "Daniel key source fragment"}))
    real = docling.chunking.HybridChunker
    monkeypatch.setattr(docling.chunking, "HybridChunker", lambda **kw: real(
        tokenizer=WordTokenizer(budget=budget), **kw))
    chunk_text(directory / "text.json", chunks_output=directory / "chunks.json",
               docling_doc_file=source)
    artifact = json.loads((directory / "chunks.json").read_text())
    assert artifact["chunker"] == "hybrid_chunker"
    from pipeline.key_context import KEY_BRANCH_CONTEXT_POLICY
    assert artifact["key_branch_context_policy"] == KEY_BRANCH_CONTEXT_POLICY
    index = CorpusIndex(tmp_path)
    assert index.load() == 1
    monkeypatch.setattr(app, "_INDEX", index)
    return captured, associations, artifact, directory


@pytest.mark.parametrize("budget", [12, 25, 2000])
def test_source_key_destinations_and_partial_context_survive_bounded_serving(tmp_path, monkeypatch, budget):
    captured, associations, artifact, directory = materialize(tmp_path, monkeypatch, budget)
    assert len(associations) == 4
    for phrase, destination in EXPECTED.items():
        association = next(a for a in associations if phrase in a["original_lead"])
        assert association["original_destination"] == destination
        assert association["status"] == "geometry_verified_spelling_unverified"
    before = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
              for path in directory.iterdir() if path.is_file()}
    scan = get_chunks_by_section(HASH, limit=100, with_text=False)
    assert len(scan) == artifact["total_chunks"]
    assert all("text" not in row for row in scan)
    returned = []
    by_branch = defaultdict(list)
    for summary in scan:
        row = get_chunks(HASH, chunk_ids=[summary["chunk_id"]])[0]
        returned.append(row)
        assert context_bytes({k: v for k, v in row.items() if k in CONTEXT_FIELDS}) <= MAX_CONTEXT_ROW_BYTES
        assert len(row["text"].split()) <= budget
        for branch in row.get("key_branches", []):
            assert branch["status"] == "geometry_verified_spelling_unverified"
            assert branch["original_destination"]["preview"] in EXPECTED.values()
            assert not branch["original_destination"]["truncated"]
            assert branch["lead_provenance"] and branch["destination_provenance"]
            # A split record must explicitly say which context is missing;
            # association confidence is a separate, preserved source status.
            assert "chunk_scope" in branch
            by_branch[branch["original_destination"]["preview"]].append((row, branch))
    assert set(by_branch) == set(EXPECTED.values())
    stored = {row["chunk_id"]: row for row in artifact["chunks"]}
    for destination, fragments in by_branch.items():
        source = next(a for a in associations if a["original_destination"] == destination)
        joined = " ".join(row["text"] for row, _ in fragments)
        assert " ".join(source["original_lead"].split()) in " ".join(joined.split())
        assert destination in joined
        for i, (row, branch) in enumerate(fragments):
            scope = branch["chunk_scope"]
            assert scope["fragment_index"] == i
            assert scope["fragment_count"] == len(fragments)
            assert scope["coverage"] == ("complete" if len(fragments) == 1 else "partial")
            assert scope["continuation"] is (len(fragments) > 1)
            assert scope["complete_branch"] is (len(fragments) == 1)
            assert scope["destination_in_text"] is (destination in row["text"])
            assert scope["previous_chunk_id"] == (fragments[i-1][0]["chunk_id"] if i else None)
            next_id = fragments[i+1][0]["chunk_id"] if i+1 < len(fragments) else None
            assert scope["next_chunk_id"] == next_id
            stored_branch = next(b for b in stored[row["chunk_id"]]["key_branches"]
                                 if b["item_ref"] == branch["item_ref"])
            assert scope == stored_branch["chunk_scope"]
            if next_id:
                next_row = get_chunks(HASH, chunk_ids=[next_id])[0]
                assert any(b["item_ref"] == branch["item_ref"] for b in next_row["key_branches"])
    if budget == 12:
        assert any(len(rows) > 1 for rows in by_branch.values())
        assert any(not b["chunk_scope"]["destination_in_text"] for rows in by_branch.values() for _, b in rows)
    # The source crop reads bargmannae; the captured text reads hargmannae.
    # Preserve that explicitly unverified spelling instead of claiming repair.
    assert "hargmannae" in " ".join(row["text"] for row in returned)
    assert "bargmannae" not in " ".join(row["text"] for row in returned)
    assert before == {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in directory.iterdir() if path.is_file()}
    (tmp_path / "acceptance.json").write_text(json.dumps({
        "source": captured["source"], "rendered_evidence": captured["rendered_evidence"],
        "word_budget": budget, "source_associations": associations,
        "scan": scan, "served_chunks": returned,
        "limits": "Captured source geometry, local word tokenizer; no fresh PDF extraction, spelling repair or embeddings.",
    }, ensure_ascii=False, indent=2) + "\n")


def test_without_geometry_association_names_do_not_create_branch_context(tmp_path, monkeypatch):
    _, associations, _, _ = materialize(tmp_path, monkeypatch, 12, associate=False)
    assert associations == []
    rows = get_chunks(HASH)
    assert set(EXPECTED.values()) <= set(" ".join(row["text"] for row in rows).split())
    assert not any(row.get("key_branches") for row in rows)


def test_rechunking_legacy_output_replaces_missing_scopes_and_is_stable(tmp_path, monkeypatch):
    from pipeline.embedding_state import input_fingerprint, load_document_data, record_payloads

    _, _, fresh, directory = materialize(tmp_path, monkeypatch, 12)
    fresh_doc = load_document_data(directory)
    payloads = record_payloads(HASH, fresh_doc)
    fingerprint = input_fingerprint(HASH, fresh_doc)
    legacy = json.loads(json.dumps(fresh))
    del legacy["key_branch_context_policy"]
    for row in legacy["chunks"]:
        for branch in row.get("key_branches", []):
            del branch["chunk_scope"]
    path = directory / "chunks.json"
    path.write_text(json.dumps(legacy))
    legacy_doc = load_document_data(directory)
    # Scope is served from chunks.json, not stored in vector-row metadata.
    # Rechunking must not waste a model run when its exact inputs are unchanged.
    assert record_payloads(HASH, legacy_doc) == payloads
    assert input_fingerprint(HASH, legacy_doc) == fingerprint
    source = directory / "docling_doc.json"
    before = source.read_bytes()
    for _ in range(2):
        chunk_text(directory / "text.json", chunks_output=path, docling_doc_file=source)
        assert json.loads(path.read_text()) == fresh
        assert source.read_bytes() == before
        assert record_payloads(HASH, load_document_data(directory)) == payloads
        assert input_fingerprint(HASH, load_document_data(directory)) == fingerprint

    changed = load_document_data(directory)
    changed["chunks"]["chunks"][0]["text"] += " changed source text"
    assert input_fingerprint(HASH, changed) != fingerprint


def test_branch_context_policy_invalidates_its_consumers_without_reextracting(monkeypatch):
    from pipeline import key_context
    from pipeline.build_inputs import config_fingerprints

    before = config_fingerprints({}, panel_mode="ocr")
    monkeypatch.setattr(key_context, "KEY_BRANCH_CONTEXT_POLICY", "changed-fragment-policy")
    after = config_fingerprints({}, panel_mode="ocr")
    changed = {stage for stage in before if before[stage] != after[stage]}
    assert changed == {"text_chunking", "taxa_and_lexicon_extraction", "figure_crossref"}


def test_unalignable_or_missing_fragment_text_reports_unknown_coverage():
    from pipeline.key_context import key_branch_context

    captured = json.loads((FIXTURE / "Daniel1985-271-272.json").read_text())
    doc = DoclingDocument.model_validate(captured["document"])
    _attach_key_destinations(doc)
    item = next(t for t in doc.texts if getattr(t.meta, "corpus__key_branch", None))
    before = doc.model_dump(mode="json")
    for text in (None, "An unrelated paragraph with no branch evidence."):
        record = key_branch_context([item], text)[0]
        assert record["chunk_scope"]["coverage"] == "unknown"
        assert record["chunk_scope"]["complete_branch"] is None
        assert record["chunk_scope"]["continuation"] is None
        assert record["status"] == "geometry_verified_spelling_unverified"
    assert doc.model_dump(mode="json") == before
