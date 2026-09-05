"""Configuration updates must change materialized data, not just receipts."""

import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from pipeline import main, metadata, runner, stages
from pipeline.build_inputs import config_fingerprints, configuration_drift
from pipeline.figure_materialization import rebuild_figure_base
from pipeline.figures import link_chunks_to_figures
from tests import test_metadata_resume

corpus = test_metadata_resume.corpus  # Shared pipeline fixture, registered here too.


@pytest.mark.parametrize("section,key,value,changed,unchanged", [
    ("ocr", "scan_page_fraction_min", 0.8, "scan_detection", "huge_document_check"),
    ("ocr", "optimize_level", 1, "pdf_preparation", "scan_detection"),
    ("ocr", "jobs", 2, "pdf_preparation", "scan_detection"),
    ("chunking", "max_tokens", 50, "text_chunking", "figure_materialization"),
    ("figures", "images_scale", 4, "docling_extraction", "metadata_extraction"),
    ("grobid", "consolidate_citations", 1, "metadata_extraction", "docling_extraction"),
    ("quality_gates", "empty_text_min_chars", 20, "quality_gates", "scan_detection"),
])
def test_config_consumers_and_nonconsumers(section, key, value, changed, unchanged):
    before = config_fingerprints({}, panel_mode="ocr")
    after = config_fingerprints({section: {key: value}}, panel_mode="ocr")
    assert before[changed] != after[changed]
    assert before[unchanged] == after[unchanged]
    assert config_fingerprints({"server": {"port": 9999}, "output_dir": "elsewhere"}, panel_mode="ocr") == before


def test_defaults_and_explicit_defaults_have_identical_fingerprints():
    from pipeline.config import _DEFAULT_CONFIG
    assert config_fingerprints({}, panel_mode="ocr") == config_fingerprints(
        copy.deepcopy(_DEFAULT_CONFIG), panel_mode="ocr")


def test_chunk_setting_rechunks_without_ocr_metadata_or_figure_detection(corpus, monkeypatch):
    corpus.run()
    hd = corpus.hd()
    before = stages._load_pipeline_state(hd)["stages"]
    corpus.config.write_text("chunking:\n  max_tokens: 4\n")
    with monkeypatch.context() as patch:
        for name in ("detect_scan_type", "prepare_pdf", "extract_docling_content", "extract_metadata", "_pass3a_annotate_rois"):
            patch.setattr(runner, name, lambda *a, **k: pytest.fail("unrelated work reran"))
        corpus.run()
    after = stages._load_pipeline_state(hd)["stages"]
    for stage in ("scan_detection", "pdf_preparation", "metadata_extraction", "docling_extraction", "figure_materialization"):
        assert before[stage] == after[stage]
    chunks = json.loads((hd / "chunks.json").read_text())
    assert len(chunks["chunks"]) > 1
    clean = corpus.run(destination=corpus.output.parent / "clean")
    assert chunks == json.loads((corpus.hd(destination=clean) / "chunks.json").read_text())
    reasons = json.loads((hd / "summary.json").read_text())["processing_summary"]["rerun_reasons"]
    assert "config.chunking.max_tokens" in reasons["text_chunking"]


def _install_figure_stubs(monkeypatch):
    def extract(pdf, text, figures, images, **kwargs):
        images.mkdir(exist_ok=True)
        (images / "raw.png").write_bytes(b"original figure bytes")
        text.write_text(json.dumps({"text": "Figure 1: fixed text"}))
        figures.write_text(json.dumps({"figures": [{"figure_id": "raw", "filename": "raw.png",
                                                   "figure_number": "1", "file_path": str(images / "raw.png")}]}))

    def split(figures, backend=None):
        data = json.loads(figures.read_text())
        image_dir = figures.parent / "figures"
        (image_dir / "raw.png").rename(image_dir / "split.png")
        data["figures"][0].update(filename="split.png", file_path=str(image_dir / "split.png"),
                                  rois=[{"label": getattr(backend, "model", "ocr")}])
        (image_dir / "obsolete.png").write_bytes(b"old derived panel")
        figures.write_text(json.dumps(data))

    monkeypatch.setattr(runner, "extract_docling_content", extract)
    monkeypatch.setattr(runner, "_pass3a_annotate_rois", split)
    monkeypatch.setattr(runner, "_pass3b_annotate_rois", split)
    monkeypatch.setattr(runner, "resolve_compound_figures", lambda _: {})
    monkeypatch.setattr("pipeline.vision.get_vision_backend",
                        lambda name, **kw: SimpleNamespace(name=name, model=kw.get("model", "default")))
    return extract


@pytest.mark.parametrize("start_mode", ["ocr", "vision-local"])
def test_turning_panels_off_resets_renamed_images_and_rois(corpus, monkeypatch, start_mode):
    _install_figure_stubs(monkeypatch)
    corpus.run("--figure-panels", start_mode)
    hd = corpus.hd()
    before = stages._load_pipeline_state(hd)["stages"]
    assert (hd / "figures" / "obsolete.png").exists()
    corpus.config.write_text('figures:\n  panel_detection: "off"\n')
    corpus.run()
    after = stages._load_pipeline_state(hd)["stages"]
    for stage in ("scan_detection", "pdf_preparation", "docling_extraction", "metadata_extraction", "text_chunking"):
        assert before[stage] == after[stage]
    assert sorted(p.name for p in (hd / "figures").iterdir()) == ["raw.png"]
    actual = json.loads((hd / "figures.json").read_text())
    assert "rois" not in actual["figures"][0]
    clean = corpus.run(destination=corpus.output.parent / "clean")
    expected = json.loads((corpus.hd(destination=clean) / "figures.json").read_text())
    for data in (actual, expected):
        data.pop("figures_directory", None)
        for figure in data["figures"]:
            figure["file_path"] = Path(figure["file_path"]).name
    assert actual == expected
    state_bytes = (hd / "pipeline_state.json").read_bytes()
    corpus.run()
    assert (hd / "pipeline_state.json").read_bytes() == state_bytes


def test_vision_model_change_starts_from_unsplit_figures(corpus, monkeypatch):
    _install_figure_stubs(monkeypatch)
    corpus.run("--figure-panels", "vision-local", "--vision-model", "model-a")
    corpus.run("--figure-panels", "vision-local", "--vision-model", "model-b")
    figures = json.loads((corpus.hd() / "figures.json").read_text())["figures"]
    assert figures[0]["rois"] == [{"label": "model-b"}]


def test_failed_figure_reset_is_not_reported_as_success(corpus, monkeypatch):
    _install_figure_stubs(monkeypatch)
    corpus.run()
    hd = corpus.hd()
    before = (hd / "figures.json").read_bytes()
    corpus.config.write_text('figures:\n  panel_detection: "off"\n')

    def fail(*args, **kwargs):
        raise RuntimeError("fresh extraction failed")

    monkeypatch.setattr(runner, "extract_docling_content", fail)
    assert main.main() == 1
    assert (hd / "figures.json").read_bytes() == before
    state = stages._load_pipeline_state(hd)["stages"]
    assert "figure_materialization" not in state
    assert "figure_crossref" not in state


def test_failed_forced_producer_invalidates_its_consumers(tmp_path):
    for stage in ("docling_extraction", "text_chunking", "figure_materialization", "figure_crossref", "metadata_extraction"):
        stages._record_stage_completion(tmp_path, stage)
    with pytest.raises(RuntimeError):
        with stages._stage({}, "docling_extraction", hash_dir=tmp_path):
            raise RuntimeError("interrupted")
    remaining = stages._load_pipeline_state(tmp_path)["stages"]
    assert set(remaining) == {"metadata_extraction"}


def test_keyboard_interrupt_never_records_success(tmp_path):
    stages._record_stage_completion(tmp_path, "docling_extraction")
    with pytest.raises(KeyboardInterrupt):
        with stages._stage({}, "docling_extraction", hash_dir=tmp_path):
            raise KeyboardInterrupt
    assert "docling_extraction" not in stages._load_pipeline_state(tmp_path)["stages"]


def test_raster_change_replaces_old_images_and_missing_docling_sidecar(corpus, monkeypatch):
    _install_figure_stubs(monkeypatch)
    corpus.run("--figure-panels", "off")
    hd = corpus.hd()
    (hd / "figures" / "obsolete.png").write_bytes(b"previous generation")
    (hd / "docling_doc.json").write_text('{"old": true}')
    before = stages._load_pipeline_state(hd)["stages"]
    corpus.config.write_text('figures:\n  panel_detection: "off"\n  images_scale: 4\n')
    corpus.run()
    assert sorted(p.name for p in (hd / "figures").iterdir()) == ["raw.png"]
    assert not (hd / "docling_doc.json").exists()
    after = stages._load_pipeline_state(hd)["stages"]
    assert before["docling_extraction"] != after["docling_extraction"]
    for stage in ("scan_detection", "pdf_preparation", "metadata_extraction"):
        assert before[stage] == after[stage]


@pytest.mark.parametrize("key,value,scan_changes", [
    ("optimize_level", 1, False), ("gibberish_threshold", 0.8, True),
])
def test_ocr_settings_reach_downstream_stages_and_removal_restores_defaults(
        corpus, key, value, scan_changes):
    corpus.run()
    hd = corpus.hd()
    before = stages._load_pipeline_state(hd)["stages"]
    corpus.config.write_text(f"ocr:\n  {key}: {value}\n")
    corpus.run()
    after = stages._load_pipeline_state(hd)["stages"]
    assert (before["scan_detection"] != after["scan_detection"]) == scan_changes
    for stage in ("pdf_preparation", "docling_extraction", "metadata_extraction", "text_chunking",
                  "figure_materialization", "figure_crossref"):
        assert before[stage]["input_fingerprint"] != after[stage]["input_fingerprint"]
    corpus.config.write_text("{}")
    corpus.run()
    restored = stages._load_pipeline_state(hd)["stages"]
    for stage in before:
        assert before[stage]["input_fingerprint"] == restored[stage]["input_fingerprint"]


def test_configuration_drift_names_inputs_without_touching_build(corpus):
    corpus.config.write_text("grobid:\n  disable: true\n")
    corpus.run()
    assert configuration_drift(corpus.output, corpus.config)["documents_with_differences"] == 0
    corpus.config.write_text("grobid:\n  disable: true\nchunking:\n  max_tokens: 4\n")
    state = (corpus.hd() / "pipeline_state.json").read_bytes()
    drift = configuration_drift(corpus.output, corpus.config)
    assert drift["documents_with_differences"] == 2
    changes = drift["differences"][corpus.hd().name]
    assert changes == {"text_chunking": ["config.chunking.max_tokens"],
                       "figure_crossref": ["config.chunking.max_tokens"],
                       "taxa_and_lexicon_extraction": ["config.chunking.max_tokens"]}
    assert (corpus.hd() / "pipeline_state.json").read_bytes() == state
    corpus.run()
    assert configuration_drift(corpus.output, corpus.config)["documents_with_differences"] == 0


def test_status_json_reports_drift_and_missing_config_is_an_error(corpus, monkeypatch, capsys):
    from pipeline import status
    corpus.run()
    capsys.readouterr()
    corpus.config.write_text("chunking:\n  max_tokens: 4\n")
    monkeypatch.setattr("sys.argv", ["status", str(corpus.output), "--config", str(corpus.config), "--json"])
    assert status.main() == 0
    assert json.loads(capsys.readouterr().out)["configuration_drift"]["documents_with_differences"] == 2
    monkeypatch.setattr("sys.argv", ["status", str(corpus.output), "--config", str(corpus.config.parent / "missing.yaml")])
    assert status.main() == 2


def test_status_does_not_apply_unrelated_cwd_config(corpus, monkeypatch):
    from pipeline import cli
    monkeypatch.chdir(corpus.config.parent)
    calls = []
    monkeypatch.setattr(cli, "_passthrough", lambda *a: calls.append(a) or 0)
    args = SimpleNamespace(config=None, output_dir=corpus.output, passthrough=[])
    assert cli._cmd_status(args) == 0
    assert calls[-1] == ("pipeline.status", [str(corpus.output)])
    args.config = corpus.config
    assert cli._cmd_status(args) == 0
    assert calls[-1] == ("pipeline.status", [str(corpus.output), "--config", str(corpus.config)])


def test_cli_panel_override_wins_over_config(corpus, monkeypatch):
    _install_figure_stubs(monkeypatch)
    corpus.config.write_text('figures:\n  panel_detection: "off"\n')
    corpus.run("--figure-panels", "ocr")
    assert (corpus.hd() / "figures" / "split.png").exists()
    assert configuration_drift(corpus.output, corpus.config)["documents_with_differences"] == 2


def test_failed_standalone_vision_cannot_leave_cpu_floor_complete(corpus, monkeypatch):
    extract = _install_figure_stubs(monkeypatch)
    corpus.run()
    monkeypatch.setattr(main, "extract_docling_content", extract)
    monkeypatch.setattr(main, "_pass25_annotate_figures", lambda *a: None)

    def fail(*args):
        raise RuntimeError("vision failed after resetting figures")

    monkeypatch.setattr(main, "_pass3b_annotate_rois", fail)
    monkeypatch.setattr("sys.argv", ["extract", str(corpus.source), str(corpus.output), "--resume",
                                   "--no-grobid", "--no-taxa", "--bib", str(corpus.bib),
                                   "--config", str(corpus.config), "--figure-panels", "vision-local", "--refresh-vision"])
    assert main.main() == 1
    receipts = stages._load_pipeline_state(corpus.hd())["stages"]
    assert "figure_materialization" not in receipts
    assert "vision_refresh" not in receipts
    assert "figure_crossref" not in receipts
    corpus.run()
    assert (corpus.hd() / "figures" / "split.png").exists()


def test_successful_standalone_vision_survives_next_cpu_resume(corpus, monkeypatch):
    extract = _install_figure_stubs(monkeypatch)
    corpus.run()
    baseline = stages._load_pipeline_state(corpus.hd())["stages"]["figure_materialization"]
    monkeypatch.setattr(main, "extract_docling_content", extract)
    monkeypatch.setattr(main, "_pass25_annotate_figures", lambda *a: None)
    monkeypatch.setattr(main, "_pass3b_annotate_rois", runner._pass3b_annotate_rois)
    monkeypatch.setattr(main, "resolve_compound_figures", lambda *a: {})
    monkeypatch.setattr(main, "_crossref_chunks_and_figures", lambda *a: None)
    monkeypatch.setattr(main, "_run_quality_gates", lambda *a: [])
    monkeypatch.setattr("pipeline.figures.generate_figures_report", lambda *a: None)
    corpus.run("--figure-panels", "vision-local", "--vision-model", "model-b", "--refresh-vision")
    receipts = stages._load_pipeline_state(corpus.hd())["stages"]
    assert receipts["figure_materialization"] == baseline
    assert "vision_refresh" in receipts
    vision = (corpus.hd() / "figures.json").read_bytes()
    corpus.run()
    assert (corpus.hd() / "figures.json").read_bytes() == vision
    assert "vision_refresh" in stages._load_pipeline_state(corpus.hd())["stages"]


def test_publication_failure_rolls_back_entire_extraction(tmp_path, monkeypatch):
    (tmp_path / "processed.pdf").write_bytes(b"pdf")
    (tmp_path / "figures").mkdir()
    (tmp_path / "figures" / "old.png").write_bytes(b"old image")
    for name in ("figures.json", "text.json", "docling_doc.json"):
        (tmp_path / name).write_text("old " + name)

    def extract(pdf, text, figures, images, **kwargs):
        text.write_text('{"text": "new"}')
        (images / "new.png").write_bytes(b"new image")
        figures.write_text('{"figures": [{"filename": "new.png"}]}')

    rename = Path.rename

    def fail_text_publication(source, destination):
        if source.name == "text.json" and source.parent.name.startswith(".figure-reset-"):
            raise OSError("simulated publication failure")
        return rename(source, destination)

    monkeypatch.setattr(Path, "rename", fail_text_publication)
    with pytest.raises(OSError, match="publication failure"):
        rebuild_figure_base(tmp_path, extract, figures_only=False)
    assert sorted(p.name for p in (tmp_path / "figures").iterdir()) == ["old.png"]
    for name in ("figures.json", "text.json", "docling_doc.json"):
        assert (tmp_path / name).read_text() == "old " + name


def test_empty_populations_clear_old_crosslinks():
    chunks = [{"text": "Figure 1", "chunk_id": "new", "figure_refs": ["obsolete"]}]
    figures = [{"figure_id": "f1", "figure_number": "1", "referenced_in_chunks": ["old"]}]
    link_chunks_to_figures(chunks, figures)
    assert figures[0]["referenced_in_chunks"] == ["new"]
    link_chunks_to_figures(chunks, [])
    assert chunks[0]["figure_refs"] == []
    link_chunks_to_figures([], figures)
    assert figures[0]["referenced_in_chunks"] == []


def test_missing_fresh_image_leaves_existing_figures_intact(tmp_path):
    (tmp_path / "processed.pdf").write_bytes(b"pdf")
    (tmp_path / "figures").mkdir()
    (tmp_path / "figures" / "old.png").write_bytes(b"old")
    (tmp_path / "figures.json").write_text('{"figures": []}')

    def extract(_pdf, _text, figures, _images, **kwargs):
        figures.write_text('{"figures": [{"filename": "missing.png"}]}')

    with pytest.raises(ValueError, match="image is missing"):
        rebuild_figure_base(tmp_path, extract)
    assert (tmp_path / "figures" / "old.png").read_bytes() == b"old"
    assert (tmp_path / "figures.json").read_text() == '{"figures": []}'


TEI = '<TEI xmlns="http://www.tei-c.org/ns/1.0"><teiHeader><fileDesc><titleStmt><title>Test</title></titleStmt></fileDesc></teiHeader><text><body/><back/></text></TEI>'


@pytest.fixture
def tei_cache(tmp_path, monkeypatch):
    pdf = tmp_path / "processed.pdf"
    pdf.write_bytes(b"original PDF bytes")
    calls = []

    def process(pdf, **kwargs):
        calls.append(kwargs)
        return TEI

    monkeypatch.setitem(metadata.CONFIG, "grobid", {"consolidate_header": 1, "consolidate_citations": 0})
    client = SimpleNamespace(process_fulltext=process)

    def run(client=client, entry=None):
        metadata.extract_metadata(pdf, tmp_path / "metadata.json", grobid_client=client, bib_entry=entry)

    return SimpleNamespace(root=tmp_path, pdf=pdf, calls=calls, run=run)


def test_tei_reuse_tracks_bytes_settings_and_not_bib_header(tei_cache, monkeypatch):
    c = tei_cache
    c.run()
    c.run(entry={"_key": "curated", "title": "Corrected"})
    assert len(c.calls) == 1
    c.pdf.write_bytes(b"changed prepared PDF")
    c.run()
    assert len(c.calls) == 2
    monkeypatch.setitem(metadata.CONFIG["grobid"], "consolidate_citations", 1)
    c.run()
    assert c.calls[-1]["consolidate_citations"] == 1
    assert len(c.calls) == 3
    c.run()
    assert len(c.calls) == 3


def test_raw_citation_policy_invalidates_legacy_tei_receipt(tei_cache):
    c = tei_cache
    c.run()
    assert c.calls[0]["include_raw_citations"] is True
    receipt_path = c.root / "grobid.tei.xml.provenance.json"
    proof = json.loads(receipt_path.read_text())
    del proof["inputs"]["include_raw_citations"]
    receipt_path.write_text(json.dumps(proof))
    c.run()
    assert len(c.calls) == 2
    assert list((c.root / "metadata_cache_history").glob("*-grobid.tei.xml"))
    c.run()
    assert len(c.calls) == 2


@pytest.mark.parametrize("damage", ["missing_receipt", "changed_tei", "changed_pdf"])
def test_invalid_tei_never_remains_active_when_grobid_is_unavailable(tei_cache, damage):
    c = tei_cache
    c.run()
    if damage == "missing_receipt":
        (c.root / "grobid.tei.xml.provenance.json").unlink()
    elif damage == "changed_tei":
        (c.root / "grobid.tei.xml").write_text(TEI + " ")
    else:
        c.pdf.write_bytes(b"new bytes")
    c.run(client=None)
    assert not (c.root / "grobid.tei.xml").exists()
    assert list((c.root / "metadata_cache_history").glob("*-grobid.tei.xml"))


def test_dry_run_does_not_load_requested_vision_backend(corpus, monkeypatch):
    monkeypatch.setattr("pipeline.vision.get_vision_backend",
                        lambda *a, **kw: pytest.fail("dry-run loaded a model"))
    corpus.run("--dry-run", "--figure-panels", "vision-local")
    assert not corpus.output.exists()
