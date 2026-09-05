"""External capability recovery through both real Stage 1 resume gates."""

import json
from types import SimpleNamespace

import pytest
import requests

from pipeline import main, metadata, runner
from pipeline.grobid_client import GrobidClient, GrobidUnavailableError
from pipeline.stages import _load_pipeline_state
from tests import test_metadata_resume

corpus = test_metadata_resume.corpus

TEI = '''<TEI xmlns="http://www.tei-c.org/ns/1.0"><teiHeader><fileDesc>
<titleStmt><title>Extracted title</title></titleStmt></fileDesc></teiHeader>
<text><body/><back><listBibl><biblStruct xml:id="b0"><analytic>
<title level="a">A cited work</title></analytic></biblStruct></listBibl></back></text></TEI>'''


@pytest.fixture
def service(monkeypatch):
    state = SimpleNamespace(alive=True, version="0.8.2", calls=[], probes=[], failures=set(), payload=TEI)

    class Client:
        def __init__(self, **kwargs):
            state.probes.append(kwargs)

        def is_alive(self):
            return state.alive

        def get_version(self):
            return state.version

        def process_fulltext(self, pdf, **kwargs):
            data = pdf.read_bytes()
            state.calls.append((data, kwargs))
            if any(name.encode() in data for name in state.failures):
                raise GrobidUnavailableError("temporary service failure")
            return state.payload

    monkeypatch.setattr(main, "GrobidClient", Client)
    monkeypatch.delenv("GROBID_URL", raising=False)
    monkeypatch.setattr("pipeline.external.STRICT_NETWORK", False)
    return state


def _receipt(corpus, name="First2001.pdf"):
    return _load_pipeline_state(corpus.hd(name))["stages"]["metadata_extraction"]


def _metadata(corpus):
    return json.loads((corpus.hd() / "metadata.json").read_text())


def test_disabled_enabled_and_unchanged_runs(corpus, service, monkeypatch):
    corpus.config.write_text("grobid:\n  disable: true\n")
    corpus.run(grobid=True)  # config, not a CLI opt-out
    assert service.probes == []
    assert _receipt(corpus)["input_fingerprint"]["grobid"]["status"] == "disabled"
    assert _metadata(corpus)["grobid"]["outcome"] == "disabled"
    before = _load_pipeline_state(corpus.hd())["stages"]
    corpus.config.write_text("grobid:\n  disable: false\n")
    with monkeypatch.context() as patch:
        for name in ("detect_scan_type", "prepare_pdf", "extract_docling_content", "chunk_text"):
            patch.setattr(runner, name, lambda *a, **k: pytest.fail("unrelated producer reran"))
        corpus.run(grobid=True)
    assert len(service.calls) == 2
    assert _receipt(corpus)["input_fingerprint"]["grobid"] == {"status": "complete", "service_version": "0.8.2"}
    after = _load_pipeline_state(corpus.hd())["stages"]
    for stage in ("scan_detection", "pdf_preparation", "docling_extraction", "text_chunking", "figure_materialization"):
        assert before[stage] == after[stage]
    old = (corpus.hd() / "pipeline_state.json").read_bytes()
    corpus.run(grobid=True)
    assert (corpus.hd() / "pipeline_state.json").read_bytes() == old
    assert len(service.calls) == 2


@pytest.mark.parametrize("curated", [True, False])
def test_startup_outage_retries_on_recovery_even_with_curated_header(corpus, service, curated):
    if not curated:
        corpus.bib.write_text("")
    service.alive = False
    corpus.run(grobid=True)
    assert _metadata(corpus)["extraction_method"] == ("bib" if curated else "placeholder")
    assert _metadata(corpus)["grobid"]["outcome"] == "unavailable"
    old = (corpus.hd() / "pipeline_state.json").read_bytes()
    corpus.run(grobid=True)
    assert (corpus.hd() / "pipeline_state.json").read_bytes() == old
    service.alive = True
    corpus.run(grobid=True)
    assert len(service.calls) == 2
    assert json.loads((corpus.hd() / "references.json").read_text())["total_references"] == 1
    clean = corpus.run(destination=corpus.output.parent / "clean", grobid=True)
    for name in ("metadata.json", "references.json", "intext_citations.json", "grobid.tei.xml.provenance.json"):
        assert json.loads((corpus.hd() / name).read_text()) == json.loads((corpus.hd(destination=clean) / name).read_text())


def test_request_failure_retries_only_failed_paper(corpus, service):
    service.failures = {"First2001.pdf"}
    corpus.run(grobid=True)
    assert _metadata(corpus)["grobid"]["outcome"] == "request_failed"
    assert _receipt(corpus)["input_fingerprint"]["grobid"]["status"] == "unavailable"
    other = (corpus.hd("Second2002.pdf") / "pipeline_state.json").read_bytes()
    service.failures.clear()
    corpus.run(grobid=True)
    assert len(service.calls) == 3
    assert (corpus.hd("Second2002.pdf") / "pipeline_state.json").read_bytes() == other
    assert _metadata(corpus)["grobid"]["outcome"] == "extracted"


def test_outage_preserves_success_and_cached_bib_refresh(corpus, service):
    corpus.run(grobid=True)
    before = (corpus.hd() / "pipeline_state.json").read_bytes()
    service.alive = False
    corpus.run(grobid=True)
    assert (corpus.hd() / "pipeline_state.json").read_bytes() == before
    corpus.set_bib(first="Revised title")
    corpus.run(grobid=True)
    assert _metadata(corpus)["title"] == "Revised title"
    assert _metadata(corpus)["grobid"]["outcome"] == "cached"
    assert _receipt(corpus)["input_fingerprint"]["grobid"]["status"] == "complete"
    service.alive = True
    corpus.run(grobid=True)
    assert len(service.calls) == 2


def test_disabling_removes_active_tei_but_retains_history(corpus, service):
    corpus.run(grobid=True)
    corpus.run("--no-grobid", grobid=True)
    assert not (corpus.hd() / "grobid.tei.xml").exists()
    assert list((corpus.hd() / "metadata_cache_history").glob("*-grobid.tei.xml"))
    assert json.loads((corpus.hd() / "references.json").read_text())["references"] == []
    assert json.loads((corpus.hd() / "intext_citations.json").read_text())["citations"] == []
    assert _metadata(corpus)["grobid"]["outcome"] == "disabled"
    corpus.run(grobid=True)
    assert len(service.calls) == 4


@pytest.mark.parametrize("change", ["version", "producer_id"])
def test_producer_change_invalidates_tei_and_metadata_only(corpus, service, change):
    corpus.run(grobid=True)
    before = _load_pipeline_state(corpus.hd())["stages"]
    if change == "version":
        service.version = "0.9.0-revision"
    else:
        corpus.config.write_text("grobid:\n  producer_id: image-sha256-models-v2\n")
    corpus.run(grobid=True)
    assert len(service.calls) == 4
    after = _load_pipeline_state(corpus.hd())["stages"]
    assert after["docling_extraction"] == before["docling_extraction"]
    proof = json.loads((corpus.hd() / "grobid.tei.xml.provenance.json").read_text())
    assert proof["inputs"]["service_version"] == service.version
    if change == "producer_id":
        assert proof["inputs"]["producer_id"] == "image-sha256-models-v2"


@pytest.mark.parametrize("payload", ["not xml", "<html>gateway error</html>", '<TEI xmlns="http://www.tei-c.org/ns/1.0"/>'])
def test_invalid_response_is_not_cached_as_success(corpus, service, payload):
    service.payload = payload
    corpus.run(grobid=True)
    assert _metadata(corpus)["grobid"]["outcome"] == "request_failed"
    assert not (corpus.hd() / "grobid.tei.xml").exists()
    service.payload = TEI
    corpus.run(grobid=True)
    assert len(service.calls) == 4
    assert _metadata(corpus)["grobid"]["outcome"] == "extracted"


def test_parser_failure_does_not_reuse_same_bad_cache_forever(corpus, service, monkeypatch):
    with monkeypatch.context() as patch:
        patch.setattr(metadata, "parse_tei_references", lambda *a: (_ for _ in ()).throw(ValueError("bad parse")))
        corpus.run(grobid=True)
    assert _metadata(corpus)["grobid"]["outcome"] == "parse_failed"
    assert not (corpus.hd() / "grobid.tei.xml").exists()
    corpus.run(grobid=True)
    assert len(service.calls) == 4


def test_url_precedence_and_timeout(corpus, service, monkeypatch):
    corpus.config.write_text("grobid:\n  url: http://config:8070\nstage_timeouts:\n  grobid: 123\n")
    corpus.run(grobid=True)
    assert service.probes[-1] == {"base_url": "http://config:8070", "timeout": 123}
    monkeypatch.setenv("GROBID_URL", "http://env:8070")
    corpus.run(grobid=True)
    assert service.probes[-1]["base_url"] == "http://env:8070"
    corpus.run("--grobid-url", "http://cli:8070", grobid=True)
    assert service.probes[-1]["base_url"] == "http://cli:8070"
    assert len(service.calls) == 2  # Placement changes don't imply new models.


def test_dry_run_never_probes_service(corpus, service):
    corpus.run("--dry-run", grobid=True)
    assert service.probes == []
    assert not corpus.output.exists()


def test_unknown_version_is_not_fabricated_and_does_not_loop(corpus, service):
    service.version = None
    corpus.run(grobid=True)
    assert _receipt(corpus)["input_fingerprint"]["grobid"] == {"status": "complete", "service_version": None}
    corpus.run(grobid=True)
    assert len(service.calls) == 2
    service.version = "0.9.0"
    corpus.run(grobid=True)
    assert len(service.calls) == 4


def test_legacy_capability_migration_refreshes_only_metadata(corpus, service):
    corpus.run(grobid=True)
    hd = corpus.hd()
    old = _load_pipeline_state(hd)
    old["stages"]["metadata_extraction"]["input_fingerprint"].pop("grobid")
    (hd / "pipeline_state.json").write_text(json.dumps(old))
    proof_path = hd / "grobid.tei.xml.provenance.json"
    proof = json.loads(proof_path.read_text())
    proof["inputs"].pop("service_version")
    proof["inputs"].pop("producer_id")
    proof_path.write_text(json.dumps(proof))
    corpus.run(grobid=True)
    assert len(service.calls) == 3
    updated = _load_pipeline_state(hd)
    for stage in ("scan_detection", "pdf_preparation", "docling_extraction", "text_chunking", "figure_materialization"):
        assert old["stages"][stage] == updated["stages"][stage]
    assert list((hd / "metadata_cache_history").glob("*-grobid.tei.xml"))


def test_status_reads_persisted_grobid_outcomes_after_unrelated_rerun(corpus, service):
    from pipeline import status
    service.failures = {"First2001.pdf"}
    corpus.run(grobid=True)
    service.alive = False
    corpus.config.write_text("chunking:\n  max_tokens: 4\n")
    corpus.run(grobid=True)
    rollup = status.aggregate(corpus.output / "documents")
    assert rollup["grobid_outcomes"] == {"request_failed": 1, "extracted": 1}
    assert json.loads(status.render_json(rollup))["grobid_outcomes"] == rollup["grobid_outcomes"]
    assert "not live service health" in status.render_text(rollup)


@pytest.mark.parametrize("startup", [True, False])
def test_strict_network_does_not_turn_outages_into_success(corpus, service, monkeypatch, startup):
    service.alive = not startup
    service.failures = {"First2001.pdf"}
    monkeypatch.setattr("sys.argv", ["extract", str(corpus.source), str(corpus.output), "--resume",
                                   "--no-taxa", "--config", str(corpus.config), "--strict-network"])
    assert main.main() == 1
    if not startup:
        assert "metadata_extraction" not in _load_pipeline_state(corpus.hd())["stages"]


@pytest.mark.parametrize("status,text,expected", [(200, "0.8.2\n", "0.8.2"), (404, "missing", None),
                                                (200, "<html>error</html>", None), (200, "", None)])
def test_service_version_probe_is_bounded(monkeypatch, status, text, expected):
    calls = []
    def get(url, **kwargs):
        calls.append((url, kwargs))
        return SimpleNamespace(status_code=status, text=text)
    monkeypatch.setattr(requests, "get", get)
    assert GrobidClient("http://example:8070/").get_version() == expected
    assert calls == [("http://example:8070/api/version", {"timeout": 5})]


def test_unreachable_version_is_explicitly_unknown(monkeypatch):
    monkeypatch.setattr(requests, "get", lambda *a, **k: (_ for _ in ()).throw(requests.Timeout()))
    assert GrobidClient().get_version() is None
