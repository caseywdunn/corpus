"""Explicit image exclusions survive build stages and all delivery gates (#302)."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from pipeline.figure_rights import materialize_figure_rights
from pipeline.figure_passes import _pass25_annotate_figures, _crossref_chunks_and_figures
from mcpsrv.tools.figures import get_figure, get_figure_image, get_figure_roi_image, get_figure_url
from tests.test_figure_licensing_states import _make_index, _CLEARED, HASH
from tests.test_signed_figure_urls import request
from mcpsrv.figure_http import make_figure_app

# Independently transcribed Hosia et al. 2024 Figure 1, PDF page 3, already
# present in the reference library's gold set (paper hash 1648cd91e973).
EXCLUSION = (
    "These images are not covered by the terms of the Creative Commons license "
    "of this publication."
)
CAPTION = (
    "Figure 1. Range of Nanomia nectophore shapes. "
    "Reproduced with permission of The Trustees of the Natural History Museum, "
    "London; Canadian Science Publishing; and the Zoological Society of Japan, "
    "respectively. " + EXCLUSION + " For permission to reuse, please contact "
    "the relevant rights holder."
)


def test_materializes_explicit_source_exclusion_with_provenance():
    fig = {"figure_id": "docling_3", "caption_text": CAPTION, "page": 3,
           "caption_source": "docling", "caption_status": "bound"}
    materialize_figure_rights([fig])
    rights = fig["figure_rights"]
    assert rights["status"] == "excluded_from_publication_license"
    assert rights["evidence"] == [{"source": "caption_explicit_exclusion",
                                  "figure_id": "docling_3", "page": 3,
                                  "caption_source": "docling", "caption_status": "bound",
                                  "text": EXCLUSION}]


def test_source_caption_reconstruction_reaches_strict_delivery(tmp_path):
    from pipeline.figures import extract_caption_info
    from tests.test_caption_fragment_recovery import source_document

    document = source_document("hosia")
    extracted = extract_caption_info(document.pictures[2], document)
    assert EXCLUSION in extracted["caption_text"]
    _idx, hd = _build_figure(tmp_path, caption=extracted["caption_text"])
    built = json.loads((hd / "figures.json").read_text())["figures"][0]
    assert built["figure_rights"]["status"] == "excluded_from_publication_license"
    refusal = get_figure_image(HASH, "docling_1", profile="manuscript")
    assert refusal.is_error
    assert refusal.structured_content["license_source"] == "figure_caption_exclusion"
    assert get_figure_url(HASH, "docling_1", profile="manuscript")["code"] == "forbidden"


@pytest.mark.parametrize("caption", [
    "Figure 1. A colony.",
    "Reproduced with permission of the Natural History Museum.",
    "These images are covered by the Creative Commons license of this publication.",
    "A CC-BY-4.0 licensed illustration.",
])
def test_credit_or_normal_license_does_not_invent_an_exclusion(caption):
    fig = {"figure_id": "a", "caption_text": caption}
    materialize_figure_rights([fig])
    assert fig["figure_rights"]["status"] == "inherit"
    assert fig["figure_rights"]["evidence"] == []


def test_shared_image_children_inherit_and_refresh_is_order_independent():
    figures = [
        {"figure_id": "host", "filename": "mixed.png", "caption_text": CAPTION},
        {"figure_id": "child", "image_shared_with": "host"},
        {"figure_id": "sibling", "filename": "mixed.png"},
        {"figure_id": "separate", "filename": "other.png"},
    ]
    materialize_figure_rights(figures)
    baseline = deepcopy(figures)
    materialize_figure_rights(figures)
    assert figures == baseline
    figures.reverse()
    materialize_figure_rights(figures)
    assert figures == list(reversed(baseline))
    for fig in figures:
        assert (fig["figure_rights"]["status"] == "inherit") == (fig["figure_id"] == "separate")
    next(f for f in figures if f["figure_id"] == "host")["caption_text"] = "Corrected caption."
    materialize_figure_rights(figures)
    assert all(f["figure_rights"]["status"] == "inherit" for f in figures)


def test_later_caption_splitting_does_not_erase_source_exclusion():
    figures = [{"figure_id": "host", "caption_text": CAPTION}]
    materialize_figure_rights(figures)
    figures[0]["caption_text"] = "Only the first caption fragment."
    figures.append({"figure_id": "child", "image_shared_with": "host"})
    materialize_figure_rights(figures, preserve_existing=True)
    assert all(f["figure_rights"]["status"] == "excluded_from_publication_license" for f in figures)
    before = deepcopy(figures)
    materialize_figure_rights(figures, preserve_existing=True)
    assert figures == before


def _build_figure(tmp_path, *, caption=CAPTION):
    idx = _make_index(tmp_path, work=dict(_CLEARED))
    idx.figure_url_base = "https://corpus.example.test"
    hd = Path(idx.papers[HASH]["hash_dir"])
    figures_file = hd / "figures.json"
    payload = json.loads(figures_file.read_text())
    payload["figures"][0]["caption_text"] = caption
    figures_file.write_text(json.dumps(payload))
    (hd / "text.json").write_text(json.dumps({"text": "Figure 1. A colony."}))
    (hd / "chunks.json").write_text(json.dumps({"chunks": []}))
    _pass25_annotate_figures(hd / "text.json", figures_file)
    _crossref_chunks_and_figures(figures_file, hd / "chunks.json")
    return idx, hd


def test_build_and_bundle_scrubbing_preserve_rights(tmp_path):
    from mcpsrv.bundle import _scrub_figures

    _, hd = _build_figure(tmp_path)
    path = hd / "figures.json"
    before = json.loads(path.read_text())["figures"][0]["figure_rights"]
    _scrub_figures(path, tmp_path)
    assert json.loads(path.read_text())["figures"][0]["figure_rights"] == before
    assert before["status"] == "excluded_from_publication_license"


@pytest.mark.parametrize("label", [None, "A", "missing-panel"])
def test_all_strict_delivery_paths_use_figure_exclusion(tmp_path, label):
    idx, _ = _build_figure(tmp_path)
    meta = get_figure(HASH, "docling_1", include_licensing=True)
    assert meta["publication_clearance"] == "undetermined"
    assert meta["license"] is None
    assert meta["license_source"] == "figure_caption_exclusion"
    assert meta["figure_rights"]["evidence"][0]["text"] == EXCLUSION
    assert "figure_rights" not in get_figure(HASH, "docling_1")
    refused = get_figure_image(HASH, "docling_1", label, profile="manuscript")
    assert refused.is_error
    assert refused.structured_content["code"] == "forbidden"
    assert "explicitly excluded" in refused.structured_content["error"]
    assert get_figure_url(HASH, "docling_1", label, profile="manuscript")["code"] == "forbidden"
    if label:
        assert get_figure_roi_image(HASH, "docling_1", label, profile="manuscript")["code"] == "forbidden"
    app = make_figure_app(idx)
    url = f"/figures/{HASH}/docling_1?profile=manuscript"
    if label:
        url += "&label=" + label
    reply = request(app, url)
    assert reply.status_code == 403
    assert reply.json()["code"] == "forbidden"
    assert reply.json()["figure_rights"]["status"] == "excluded_from_publication_license"
    assert request(app, url.replace("manuscript", "report")).status_code == 200
    assert get_figure_image(HASH, "docling_1", label, profile="report")


def test_permitted_control_stays_usable(tmp_path):
    idx, _ = _build_figure(tmp_path, caption="Figure 1. A colony.")
    assert get_figure(HASH, "docling_1", include_licensing=True)["publication_clearance"] == "licensed_open"
    assert get_figure_image(HASH, "docling_1", profile="manuscript")
    assert "url" in get_figure_url(HASH, "docling_1", profile="manuscript")
    assert request(make_figure_app(idx), f"/figures/{HASH}/docling_1?profile=manuscript").status_code == 200


def test_server_does_not_interpret_legacy_caption_text(tmp_path):
    idx = _make_index(tmp_path, work=dict(_CLEARED))
    path = Path(idx.papers[HASH]["hash_dir"]) / "figures.json"
    data = json.loads(path.read_text())
    data["figures"][0]["caption_text"] = CAPTION
    path.write_text(json.dumps(data))
    # Repair is a rebuild requirement, never a new query-time inference.
    assert get_figure(HASH, "docling_1", include_licensing=True)["publication_clearance"] == "licensed_open"
