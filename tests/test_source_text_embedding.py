"""Saved source repairs through real chunking, Stage 2 writes and bounded fetches.

Only token counting and the numerical embedding backend are local doubles. The
production serializers, embedding-input selection, LanceDB writes and tool
functions run unchanged. This tests input preservation, not embedding quality.
"""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mcpsrv import app
from mcpsrv.tools.chunks import get_chunks
from tests.test_embedding_updates import build, cli  # noqa: F401 — pytest fixture

FIXTURE = Path(__file__).parent / "fixtures/text_integrity/materialized"
MANIFEST = json.loads((FIXTURE / "manifest.json").read_text())
EXPECTED = {
    "Sutherland_etal2019b-1-1": [
        "0.15 ±0.10",
        "104 ±41 deg. s⁻¹",
        "1063 ±176 mm s⁻¹",
        "R / L",
    ],
    "Sutherland_etal2019b-4-4": ["0.002 x −0.0715", "deg. s⁻¹", "L / R", "L/R"],
    "Hauss_et2016-4-4": ["±0.99) m³"],
    "Haddock_etal2005-1-1": ["200 µm", "1 mm"],
    "Kidwai_Amjad2000-1-1": ["1857-1859"],
    "Pakhomov_etal2000-5-5": ["12.5 ± 6.1", "35.5 ± 39.7"],
    "Boysen-Ennen1987-12": ["RMT-8-Fänge", "Netzdffnung"],
}


@pytest.fixture
def source_pipeline(build, monkeypatch):  # noqa: F811
    import docling.chunking
    from docling_core.transforms.chunker.tokenizer.base import BaseTokenizer
    from docling_core.types.doc import DoclingDocument
    from pipeline.chunking import chunk_text
    from pipeline.table_structure import export_source_markdown

    class LocalTokenizer(BaseTokenizer):
        def count_tokens(self, text):
            return len(text.split())

        def get_max_tokens(self):
            return 2000

        def get_tokenizer(self):
            return self.count_tokens

    real = docling.chunking.HybridChunker
    monkeypatch.setattr(
        docling.chunking,
        "HybridChunker",
        lambda **kw: real(tokenizer=LocalTokenizer(), **kw),
    )
    calls = []

    def record(texts):
        calls.append(list(texts))
        return [[float(len(text)), 1.0] for text in texts]

    build.backend.embed = record

    def load(name):
        path = FIXTURE / MANIFEST["documents"][name]["file"]
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == MANIFEST["documents"][name]["artifact_sha256"]
        )
        return DoclingDocument.load_from_json(path)

    def save_and_chunk(name, document):
        record = MANIFEST["documents"][name]
        sha = record["source_sha256"][:12]
        folder = build.root / "documents" / sha
        folder.mkdir(parents=True, exist_ok=True)
        document.save_as_json(folder / "docling_doc.json")
        # Actual source-only exporter, so provenance stays structured.
        (folder / "text.json").write_text(
            json.dumps({"text": export_source_markdown(document)})
        )
        (folder / "summary.json").write_text(
            json.dumps({"relative_paths": [record["source"]]})
        )
        chunk_text(folder / "text.json", chunks_output=folder / "chunks.json")
        artifact = json.loads((folder / "chunks.json").read_text())
        assert artifact["chunker"] == "hybrid_chunker"
        monkeypatch.setattr(
            app, "_INDEX", SimpleNamespace(papers={sha: {"hash_dir": str(folder)}})
        )
        return folder, artifact

    return SimpleNamespace(
        build=build, load=load, save_and_chunk=save_and_chunk, calls=calls
    )


def assert_consumed_and_served(replay, folder, artifact):
    # Capture inside embed_document's backend call, not by independently
    # reconstructing what its input ought to be.
    assert replay.build.run(folder) == len(artifact["chunks"])
    texts = [row["text"] for row in artifact["chunks"]]
    assert replay.calls[-1] == texts
    stored = replay.build.table.to_arrow().to_pylist()
    by_id = {row["metadata"]["chunk_id"]: row["text"] for row in stored}
    assert by_id == {row["chunk_id"]: row["text"] for row in artifact["chunks"]}
    before = {
        name: (folder / name).read_bytes()
        for name in ("docling_doc.json", "text.json", "chunks.json")
    }
    for chunk in artifact["chunks"]:
        # One identified chunk per bounded fetch; no corpus-wide text request.
        rows = get_chunks(folder.name, chunk_ids=[chunk["chunk_id"]])
        assert len(rows) == 1
        assert rows[0]["text"] == chunk["text"] == by_id[chunk["chunk_id"]]
        assert len(json.dumps(rows).encode()) < 64 * 1024
    assert all((folder / name).read_bytes() == value for name, value in before.items())
    assert replay.build.problem(folder) is None
    combined = "\n".join(texts)
    assert all(
        token not in combined
        for token in (
            "corpus__",
            "original_native_token",
            "aligned_native_and_rendered_sign",
            "native_glyph_not_confirmed_by_source_raster",
        )
    )
    return " ".join(combined.split())


@pytest.mark.parametrize("name", EXPECTED)
def test_saved_source_text_reaches_actual_embedding_input_and_serving(
    source_pipeline, name
):
    replay = source_pipeline
    document = replay.load(name)
    folder, artifact = replay.save_and_chunk(name, document)
    text = assert_consumed_and_served(replay, folder, artifact)
    for expected in EXPECTED[name]:
        assert expected in text
    if name == "Sutherland_etal2019b-4-4":
        assert text.count("0.002 x −0.0715") == 2
        assert "0.002 x 0.0715" not in text
    if name == "Haddock_etal2005-1-1":
        assert "200 mm" not in text
    if name == "Kidwai_Amjad2000-1-1":
        assert "1857±1859" not in text and "1857 ± 1859" not in text
        assert not MANIFEST["documents"][name]["report"]["repairs"]
        assert any(c.get("text_integrity") for c in artifact["chunks"])
    if name == "Pakhomov_etal2000-5-5":
        assert not MANIFEST["documents"][name]["report"]["repairs"]
    if name == "Boysen-Ennen1987-12":
        assert "RMT-8-Finge" not in text and "RMT-8-FÃ¤ng" not in text
        assert any("RMT-8-Finge" in item.orig for item in document.texts)
        assert any(
            note.get("replacement") == "RMT-8-Fänge"
            for c in artifact["chunks"]
            for note in c.get("text_integrity", [])
        )
        assert any(
            note.get("status") == "unresolved"
            for c in artifact["chunks"]
            for note in c.get("text_integrity", [])
        )


@pytest.mark.parametrize(
    "name, repaired, damaged",
    [
        ("Sutherland_etal2019b-1-1", "1063 ±176 mm s⁻¹", "1063 176 mm s 1"),
        ("Boysen-Ennen1987-12", "RMT-8-Fänge", "RMT-8-Finge"),
    ],
)
def test_source_correction_invalidates_old_vectors_and_unchanged_resume_skips(
    source_pipeline,
    monkeypatch,
    name,
    repaired,
    damaged,
):
    replay = source_pipeline
    corrected = replay.load(name)
    before = deepcopy(corrected)
    # A controlled earlier-generation replay of the exact persisted original
    # extraction, not a new OCR claim or a fabricated source spelling.
    for item in before.texts:
        item.text = item.orig
    folder, old_artifact = replay.save_and_chunk(name, before)
    old_text = assert_consumed_and_served(replay, folder, old_artifact)
    assert damaged in old_text and repaired not in old_text
    folder, corrected_artifact = replay.save_and_chunk(name, corrected)
    assert replay.build.problem(folder) == "embedding inputs changed"
    corrected_text = assert_consumed_and_served(replay, folder, corrected_artifact)
    assert repaired in corrected_text and damaged not in corrected_text
    assert replay.build.table.count_rows() == len(corrected_artifact["chunks"])
    version, count = replay.build.table.version, len(replay.calls)
    assert cli(monkeypatch, replay.build, "--resume") == 0
    assert len(replay.calls) == count
    assert replay.build.table.version == version


def test_scientific_policy_invalidates_extraction_and_its_consumers(monkeypatch):
    from pipeline.build_inputs import config_fingerprints
    from pipeline import scientific_text

    before = config_fingerprints({}, panel_mode="ocr")
    monkeypatch.setattr(
        scientific_text,
        "SCIENTIFIC_TEXT_POLICY",
        "source_glyph_alignment_test_revision",
    )
    after = config_fingerprints({}, panel_mode="ocr")
    assert {stage for stage in before if before[stage] != after[stage]} == {
        "docling_extraction",
        "text_chunking",
        "taxa_and_lexicon_extraction",
        "figure_materialization",
        "figure_crossref",
    }
