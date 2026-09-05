"""Per-paper source edits through both resume gates; no models or network."""

import copy
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from bib import BibIndex
from pipeline import main, runner
from pipeline.config import CONFIG
from pipeline.io import calculate_pdf_hash, short_hash
from pipeline.stages import (
    _expected_fingerprints_for_run, _load_pipeline_state,
    _metadata_fingerprint_for_pdf,
)


@pytest.fixture
def corpus(tmp_path, monkeypatch):
    source = tmp_path / "library"
    source.mkdir()
    for name in ("First2001.pdf", "Second2002.pdf"):
        (source / name).write_bytes(b"%PDF stub: " + name.encode())
    bib = tmp_path / "library.bib"
    config = tmp_path / "config.yaml"
    config.write_text("{}")
    output = tmp_path / "output"
    saved_config = copy.deepcopy(CONFIG)
    monkeypatch.setattr(main, "setup_root_logging", lambda: None)
    monkeypatch.setattr(runner, "_pdf_page_count", lambda _: 1)
    monkeypatch.setattr(runner, "detect_scan_type", lambda *a, **kw: {"file_type": "born_digital"})
    monkeypatch.setattr(runner, "prepare_pdf", lambda src, _d, dst: shutil.copyfile(src, dst) and {})

    def extract(pdf, text, figures, figures_dir, **kwargs):
        text.write_text(json.dumps({"text": "fixed extracted text"}))
        figures.write_text(json.dumps({"figures": []}))

    # Real metadata extraction + real fallback chunker; only expensive PDF/
    # figure work is replaced. The two resume gates and disk receipts are real.
    monkeypatch.setattr(runner, "extract_docling_content", extract)
    for name in ("_pass25_annotate_figures", "_pass3a_annotate_rois",
                 "_crossref_chunks_and_figures", "generate_figures_report"):
        monkeypatch.setattr(runner, name, lambda *a, **kw: None)
    monkeypatch.setattr(runner, "_run_quality_gates", lambda _: [])

    class InlineProcess:
        def __init__(self, target):
            self.target = target
            self.exitcode = 0

        def start(self):
            self.target()

        def join(self):
            pass

    monkeypatch.setattr(main.multiprocessing, "get_context",
                        lambda _: SimpleNamespace(Process=InlineProcess))

    def set_bib(first="Original", second="Untouched", first_file="First2001.pdf"):
        bib.write_text(
            f"@article{{first, title={{{first}}}, file={{{first_file}}}, year={{2001}}}}\n"
            f"@article{{second, title={{{second}}}, file={{Second2002.pdf}}, year={{2002}}}}\n")

    def run(*extra, destination=output, grobid=False, annotations=False):
        monkeypatch.setattr("sys.argv", ["extract", str(source), str(destination),
                                       "--resume", *([] if grobid else ["--no-grobid"]),
                                       *([] if annotations else ["--no-taxa"]),
                                       "--config", str(config), "--bib", str(bib), *extra])
        result = main.main()
        assert result in (None, 0)
        return destination

    def hd(name="First2001.pdf", destination=output):
        return destination / "documents" / short_hash(calculate_pdf_hash(source / name))

    set_bib()
    yield SimpleNamespace(source=source, bib=bib, config=config, output=output,
                          run=run, hd=hd, set_bib=set_bib)
    CONFIG.clear()
    CONFIG.update(saved_config)


def test_edit_one_entry_refreshes_only_its_metadata_and_matches_clean(corpus):
    corpus.run()
    a, b = corpus.hd(), corpus.hd("Second2002.pdf")
    before_a = _load_pipeline_state(a)["stages"]
    before_b = (b / "pipeline_state.json").read_bytes()
    unchanged_b_summary = (b / "summary.json").read_bytes()
    corpus.set_bib(first="Corrected")
    corpus.run()
    after_a = _load_pipeline_state(a)["stages"]
    assert json.loads((a / "metadata.json").read_text())["title"] == "Corrected"
    assert before_a["metadata_extraction"] != after_a["metadata_extraction"]
    for stage in ("scan_detection", "pdf_preparation", "docling_extraction", "text_chunking"):
        assert before_a[stage] == after_a[stage]
    assert (b / "pipeline_state.json").read_bytes() == before_b
    assert (b / "summary.json").read_bytes() == unchanged_b_summary
    clean = corpus.run(destination=corpus.output.parent / "clean")
    clean_a = corpus.hd(destination=clean)
    for filename in ("metadata.json", "chunks.json", "references.json"):
        assert json.loads((a / filename).read_text()) == json.loads((clean_a / filename).read_text())
    previous = (a / "pipeline_state.json").read_bytes()
    corpus.run()
    assert (a / "pipeline_state.json").read_bytes() == previous


def test_adding_and_removing_matching_entry_is_not_skipped(corpus):
    corpus.bib.write_text("")
    corpus.run()
    hd = corpus.hd()
    original = json.loads((hd / "metadata.json").read_text())
    corpus.set_bib(first="Curated")
    corpus.run()
    assert json.loads((hd / "metadata.json").read_text())["title"] == "Curated"
    corpus.bib.write_text("")
    corpus.run()
    assert json.loads((hd / "metadata.json").read_text()) == original


@pytest.mark.parametrize("has_entry", [False, True])
def test_rename_refreshes_metadata_with_or_without_a_bib_match(corpus, has_entry):
    if not has_entry:
        corpus.bib.write_text("")
    corpus.run()
    hd = corpus.hd()
    before = _load_pipeline_state(hd)["stages"]
    (corpus.source / "First2001.pdf").rename(corpus.source / "Renamed2010.pdf")
    if has_entry:
        corpus.set_bib(first="Renamed entry", first_file="Renamed2010.pdf")
    corpus.run()
    metadata = json.loads((hd / "metadata.json").read_text())
    assert metadata["filename"] == "Renamed2010.pdf"
    if has_entry:
        assert metadata["title"] == "Renamed entry"
    else:
        assert metadata["year"] == 2010
    assert json.loads((hd / "summary.json").read_text())["relative_paths"] == ["Renamed2010.pdf"]
    assert _load_pipeline_state(hd)["stages"]["docling_extraction"] == before["docling_extraction"]


def test_added_identical_copy_refreshes_paths_without_reextracting(corpus):
    corpus.run()
    hd = corpus.hd()
    before = (hd / "pipeline_state.json").read_bytes()
    (corpus.source / "copies").mkdir()
    shutil.copyfile(corpus.source / "First2001.pdf", corpus.source / "copies" / "First2001.pdf")
    corpus.run()
    summary = json.loads((hd / "summary.json").read_text())
    assert set(summary["relative_paths"]) == {"First2001.pdf", "copies/First2001.pdf"}
    assert summary["total_copies_found"] == 2
    assert (hd / "pipeline_state.json").read_bytes() == before


def test_metadata_edit_preserves_existing_vision_figures(corpus, monkeypatch):
    corpus.run()
    hd = corpus.hd()
    figures = hd / "figures.json"
    figures.write_text(json.dumps({"figures": [{"figure_id": "f1", "pass3_backend": "vision",
                                              "rois": [{"label": "A"}]}]}))
    before = figures.read_bytes()
    for name in ("_pass25_annotate_figures", "_pass3a_annotate_rois",
                 "_pass3b_annotate_rois", "resolve_compound_figures",
                 "_crossref_chunks_and_figures"):
        monkeypatch.setattr(runner, name, lambda *a, **kw: pytest.fail("metadata edit reran figure work"))
    corpus.set_bib(first="Corrected")
    corpus.run()
    assert figures.read_bytes() == before


def test_dry_run_sees_edited_entry_without_mutating_artifacts(corpus, caplog):
    corpus.run()
    hd = corpus.hd()
    before = (hd / "pipeline_state.json").read_bytes()
    corpus.set_bib(first="Changed")
    with caplog.at_level("INFO"):
        corpus.run("--dry-run")
    assert "partial-process 1, skip 1" in caplog.text
    assert (hd / "pipeline_state.json").read_bytes() == before


def test_metadata_fingerprint_is_canonical_and_per_resolved_entry():
    entry = {"file": "Paper.pdf", "title": "Title", "_key": "key"}
    first = _metadata_fingerprint_for_pdf(BibIndex([entry]), "Paper.pdf")
    reordered = dict(reversed(list(entry.items())))
    assert _metadata_fingerprint_for_pdf(BibIndex([reordered, {"file": "Else.pdf", "title": "Other"}]), "Paper.pdf") == first
    assert _metadata_fingerprint_for_pdf(None, "Paper.pdf") != first
    fps = _expected_fingerprints_for_run(config_fingerprints=None, metadata_fingerprint=first,
                                        ocrlang=None, ocrmode=None, keeppages=None)
    assert fps["metadata_extraction"] == first
    assert all(not value for key, value in fps.items() if key != "metadata_extraction")


def test_production_fingerprint_calls_cannot_omit_metadata():
    import ast
    for module in (main, runner):
        tree = ast.parse(Path(module.__file__).read_text())
        calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                 and isinstance(n.func, ast.Name) and n.func.id == "_expected_fingerprints_for_run"]
        assert calls
        for call in calls:
            assert "metadata_fingerprint" in {kw.arg for kw in call.keywords}
