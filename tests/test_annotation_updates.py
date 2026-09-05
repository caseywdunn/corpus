"""Annotation changes and deliberate absence converge through both resume gates."""
import json
import sqlite3

import pytest

from pipeline import main, taxon_mentions
from pipeline.annotate import _extract_taxa_and_lexicons, annotation_outputs_problem
from pipeline.stages import _load_pipeline_state, _record_stage_completion
from pipeline.taxa import load_lexicon
from tests import test_metadata_resume

corpus = test_metadata_resume.corpus


def test_removed_category_and_whole_lexicon_match_clean(corpus):
    path = corpus.config.parent / "lexicon.yaml"
    path.write_text("anatomy:\n  text: {}\nmethods:\n  fixed: {}\n")
    corpus.run("--lexicon", str(path), annotations=True)
    hd = corpus.hd()
    assert (hd / "methods.json").exists()
    before = _load_pipeline_state(hd)["stages"]
    path.write_text("anatomy:\n  text: {}\n")
    corpus.run("--lexicon", str(path), annotations=True)
    assert not (hd / "methods.json").exists()
    assert len(list((hd / "annotation_history").glob("*-methods.json"))) == 1
    after = _load_pipeline_state(hd)["stages"]
    assert before["docling_extraction"] == after["docling_extraction"]
    assert before["text_chunking"] == after["text_chunking"]
    corpus.run(annotations=True)
    assert not (hd / "anatomy.json").exists()
    assert json.loads((hd / "annotation_outputs.json").read_text())["outputs"] == {}
    clean = corpus.run(destination=corpus.output.parent / "clean", annotations=True)
    assert (hd / "annotation_outputs.json").read_bytes() == (corpus.hd(destination=clean) / "annotation_outputs.json").read_bytes()
    before = (hd / "pipeline_state.json").read_bytes()
    corpus.run(annotations=True)
    assert (hd / "pipeline_state.json").read_bytes() == before


def test_empty_category_retires_old_output(corpus):
    path = corpus.config.parent / "lexicon.yaml"
    path.write_text("anatomy:\n  text: {}\n")
    corpus.run("--lexicon", str(path), annotations=True)
    path.write_text("anatomy: {}\n")
    corpus.run("--lexicon", str(path), annotations=True)
    assert not (corpus.hd() / "anatomy.json").exists()
    assert annotation_outputs_problem(corpus.hd()) is None


def test_corrupt_or_missing_annotation_is_repaired_on_resume(corpus):
    path = corpus.config.parent / "lexicon.yaml"
    path.write_text("anatomy:\n  text: {}\n")
    corpus.run("--lexicon", str(path), annotations=True)
    target = corpus.hd() / "anatomy.json"
    before = target.read_bytes()
    for corrupt in (True, False):
        if corrupt:
            target.write_text("{")
        else:
            target.unlink()
        assert annotation_outputs_problem(corpus.hd()) is not None
        corpus.run("--lexicon", str(path), annotations=True)
        assert target.read_bytes() == before


def test_missing_configured_lexicon_does_not_retire_evidence(corpus, monkeypatch):
    path = corpus.config.parent / "lexicon.yaml"
    path.write_text("anatomy:\n  text: {}\n")
    corpus.run("--lexicon", str(path), annotations=True)
    before = (corpus.hd() / "anatomy.json").read_bytes()
    path.unlink()
    monkeypatch.setattr("sys.argv", ["extract", str(corpus.source), str(corpus.output),
                                    "--resume", "--no-grobid", "--lexicon", str(path)])
    assert main.main() == 1
    assert (corpus.hd() / "anatomy.json").read_bytes() == before


def test_deliberate_taxonomy_removal_clears_index_but_missing_evidence_fails(tmp_path):
    hd = tmp_path / "documents" / "abc"
    hd.mkdir(parents=True)
    (hd / "chunks.json").write_text('{"chunks": []}')
    taxa = hd / "taxa.json"
    taxa.write_text('{"unique_taxa": 0, "mentions": []}')
    conn = sqlite3.connect(":memory:")
    taxon_mentions.create_schema(conn)
    assert taxon_mentions.build(conn, tmp_path)["papers"] == 1
    taxa.unlink()
    assert taxon_mentions.build(conn, tmp_path)["errors"] == 1
    _extract_taxa_and_lexicons(hd / "chunks.json", hd, None)
    _record_stage_completion(hd, "taxa_and_lexicon_extraction")
    assert taxon_mentions.build(conn, tmp_path)["retired"] == 1
    assert conn.execute("SELECT COUNT(*) FROM papers_processed").fetchone()[0] == 0
    conn.close()


@pytest.mark.parametrize("category", ["../summary", "/tmp/escape", "summary", "taxa", "annotation_outputs"])
def test_category_cannot_overwrite_core_evidence(tmp_path, category):
    path = tmp_path / "lexicon.yaml"
    path.write_text(f'{category!r}: {{}}\n')
    with pytest.raises(ValueError, match="category"):
        load_lexicon(path)
