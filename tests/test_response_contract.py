"""Readable shape snapshots of representative success/error responses (#277)."""
import contextlib
import json
from pathlib import Path

from mcpsrv import app
from mcpsrv.tools.bibliography import format_citations
from mcpsrv.tools.chunks import get_chunks
from mcpsrv.tools.figures import get_figure, get_figure_roi_image, get_figure_url
from mcpsrv.tools.papers import get_papers
from mcpsrv.tools.taxonomy import get_taxon_dossier
from tests.test_figure_licensing_states import _make_index, _CLEARED, _UNCLEARED, HASH
from tests.test_format_citation_tool import index as citation_fixture
from tests.test_taxon_dossier_and_get_chunks import corpus as dossier_fixture


def _shape(value):
    if isinstance(value, dict):
        return {key: _shape(val) for key, val in sorted(value.items())}
    if isinstance(value, list):
        shapes = {json.dumps(_shape(item), sort_keys=True) for item in value}
        return [json.loads(item) for item in sorted(shapes)]
    return type(value).__name__


def _sample_results(tmp_path):
    results = {}
    with contextlib.contextmanager(dossier_fixture.__wrapped__)(tmp_path / "dossier"):
        results["chunks_text"] = get_chunks(HASH)
        results["chunks_metadata"] = get_chunks(HASH, with_text=False)
        results["chunks_missing"] = get_chunks("missing")
        results["papers_projected"] = get_papers([HASH], fields=["title", "year", "first_author"])
        results["taxon_dossier"] = get_taxon_dossier("Marrus")
    citation_dir = tmp_path / "citations"
    citation_dir.mkdir()
    with contextlib.contextmanager(citation_fixture.__wrapped__)(citation_dir):
        for tier, work_id in (("bib", "corpus:dunn|2005|marrus"),
                              ("reconciled", "corpus:siebert|2011|hapless"),
                              ("unresolved", "corpus:lensia|1965|unknown"),
                              ("missing", "absent")):
            results["citations_" + tier] = format_citations(work_ids=[work_id])
    original = app._INDEX
    try:
        idx = _make_index(tmp_path / "figures", work=dict(_CLEARED))
        idx.figure_url_base = "https://corpus.example.test"
        for profile in ("report", "manuscript"):
            results["figure_" + profile] = get_figure(HASH, "docling_1", profile=profile)
            results["figure_url_" + profile] = get_figure_url(HASH, "docling_1", "A", profile=profile)
            results["figure_roi_" + profile] = get_figure_roi_image(HASH, "docling_1", "A", profile=profile)
        idx.biblio_db._work = dict(_UNCLEARED)
        results["figure_url_refused"] = get_figure_url(HASH, "docling_1", profile="manuscript")
    finally:
        app.set_index(original)
    return _shape(results)


def test_representative_response_shapes_are_frozen(tmp_path):
    expected = json.loads((Path(__file__).parent / "fixtures/mcp_response_contract.json").read_text())
    assert _sample_results(tmp_path) == expected
