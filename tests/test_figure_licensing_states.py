"""Figure licensing: states, advisory-metadata suppression, gate coverage (#154).

Three defects, fixed together because §1 and §3 pull in opposite
directions:

* **§1** Advisory license metadata was injected regardless of profile, so
  a model *just permitted* to display a figure read ``publishable: false``
  beside it and self-censored. On the served corpus 86% of works were
  ``publishable=0, license_source=unknown`` and not one was asserted
  ``all-rights-reserved`` — the flag was ~100% false-positive-shaped for
  its apparent meaning.
* **§2** ``publishable`` collapsed "explicitly restricted", "unrecognized
  license string", and "nothing on file" into a single ``0``.
* **§3** ``get_figure_roi_image`` never called the gate, so a client
  refused by ``get_figure_image`` under ``manuscript`` could obtain the
  same pixels — including the whole uncropped figure via the no-pixel-ROI
  fallback.
"""
from __future__ import annotations

import json
import types
from pathlib import Path

import pytest
from PIL import Image as PILImage

from bib.authority import clearance_state, derive_publishable
from mcpsrv import app as mcp_app
from mcpsrv.tools.figures import (
    get_figure,
    get_figure_image,
    get_figure_roi_image,
)

HASH = "aaaaaaaaaaaa"


# --- §2 the five states -------------------------------------------------------


def test_no_record_is_not_restricted():
    """The whole point: an absence of evidence must not read as refusal."""
    assert clearance_state({
        "publishable": 0, "license": None, "license_source": "unknown",
    }) == "no_record"


def test_missing_work_is_no_record():
    assert clearance_state(None) == "no_record"
    assert clearance_state({}) == "no_record"


def test_asserted_restriction_is_restricted():
    assert clearance_state({
        "publishable": 0, "license": "all-rights-reserved",
        "license_source": "bibtex",
    }) == "restricted"


def test_age_based_pd_is_public_domain():
    assert clearance_state({
        "publishable": 1, "license": None, "license_source": "age_based_pd",
    }) == "public_domain"


def test_open_license_is_licensed_open():
    assert clearance_state({
        "publishable": 1, "license": "CC-BY-4.0", "license_source": "bibtex",
    }) == "licensed_open"


def test_unrecognized_license_is_undetermined():
    """The typo case: `license = {CC-BY 4.0}` (space, not hyphen)."""
    pub, source = derive_publishable("CC-BY 4.0", 2010)
    assert pub is None and source == "bibtex"
    assert clearance_state({
        "publishable": None, "license": "CC-BY 4.0", "license_source": source,
    }) == "undetermined"


def test_unrecognized_license_warns_at_build_time(caplog):
    """It used to be NULLed silently, so a typo blocked figures under
    strict profiles with nothing in the log to explain it."""
    import logging
    with caplog.at_level(logging.WARNING, logger="bib.authority"):
        derive_publishable("CC-BY 4.0 (see publisher)", 2010)
    joined = "\n".join(r.getMessage() for r in caplog.records)
    assert "unrecognized license string" in joined, joined


def test_recognized_licenses_do_not_warn(caplog):
    import logging
    with caplog.at_level(logging.WARNING, logger="bib.authority"):
        derive_publishable("CC-BY-4.0", 2010)
        derive_publishable("all-rights-reserved", 2010)
    assert [r for r in caplog.records if r.levelno >= logging.WARNING] == []


# --- fixtures for the tool-level tests ----------------------------------------


class _StubBiblio:
    def __init__(self, work):
        self._work = work

    def get_work_by_corpus_hash(self, _h):
        return self._work

    def get_authors(self, _w):
        return [{"surname": "Dunn"}]


def _make_index(tmp_path, *, work, default_profile="report"):
    hash_dir = tmp_path / "documents" / HASH
    (hash_dir / "figures").mkdir(parents=True)
    PILImage.new("RGB", (40, 40), "white").save(hash_dir / "figures" / "fig_1.png")
    (hash_dir / "figures.json").write_text(json.dumps({"figures": [{
        "figure_id": "docling_1",
        "filename": "fig_1.png",
        "figure_type": "figure",
        "caption_text": "Figure 1. Panels A and B.",
        "panels_from_caption": [{"label": "A", "description": "left"}],
        "rois": [{"label": "A", "roi_px": [0, 0, 20, 20]}],
    }]}))
    idx = types.SimpleNamespace(
        papers={HASH: {"hash_dir": str(hash_dir), "title": "A paper"}},
        biblio_db=_StubBiblio(work),
        default_profile=default_profile,
        figure_url_base=None,
        figure_auth_token=None,
    )
    mcp_app.set_index(idx)
    return idx


_UNCLEARED = {
    "publishable": 0, "license": None, "license_source": "unknown",
    "title": "A paper", "year": 2010, "work_id": "w1",
}
_CLEARED = {
    "publishable": 1, "license": "CC-BY-4.0", "license_source": "bibtex",
    "title": "A paper", "year": 2010, "work_id": "w1",
}


# --- §1 advisory metadata is not leaked into permissive use -------------------


def test_get_figure_omits_the_determination_under_report(tmp_path):
    _make_index(tmp_path, work=_UNCLEARED)
    out = get_figure(HASH, "docling_1")            # default profile = report
    assert "publishable" not in out, (
        "a permission-shaped flag under the permissive profile is what made "
        "clients withhold cleared figures (#154 §1)"
    )
    assert "publication_clearance" not in out
    # ...but attribution fields a caption needs are still there.
    assert "license" in out and "attribution" in out


def test_get_figure_includes_the_determination_under_strict(tmp_path):
    _make_index(tmp_path, work=_UNCLEARED)
    out = get_figure(HASH, "docling_1", profile="manuscript")
    assert out["publication_clearance"] == "no_record"
    assert out["license_source"] == "unknown"


def test_get_figure_includes_the_determination_on_request(tmp_path):
    _make_index(tmp_path, work=_UNCLEARED)
    out = get_figure(HASH, "docling_1", include_licensing=True)
    assert out["publication_clearance"] == "no_record"


def test_get_figure_never_reports_the_old_boolean(tmp_path):
    """`publishable` was the field name that read as 'may I display this'."""
    _make_index(tmp_path, work=_UNCLEARED)
    for kwargs in ({}, {"profile": "manuscript"}, {"include_licensing": True}):
        assert "publishable" not in get_figure(HASH, "docling_1", **kwargs)


def test_get_figure_rejects_an_unknown_profile(tmp_path):
    _make_index(tmp_path, work=_CLEARED)
    out = get_figure(HASH, "docling_1", profile="thesis")
    assert out["code"] == "invalid_argument"


# --- §3 get_figure_roi_image honors the gate ----------------------------------


def test_roi_image_refused_under_strict_profile(tmp_path):
    """The bypass: this returned pixels a strict client had been refused."""
    _make_index(tmp_path, work=_UNCLEARED)
    out = get_figure_roi_image(HASH, "docling_1", "A", profile="manuscript")
    assert out.get("code") == "forbidden", out
    assert "image_path" not in out


def test_roi_image_allowed_under_report(tmp_path):
    _make_index(tmp_path, work=_UNCLEARED)
    out = get_figure_roi_image(HASH, "docling_1", "A", profile="report")
    assert out.get("image_path"), out


def test_roi_image_allowed_when_cleared(tmp_path):
    _make_index(tmp_path, work=_CLEARED)
    out = get_figure_roi_image(HASH, "docling_1", "A", profile="manuscript")
    assert out.get("image_path"), out


def test_roi_image_whole_figure_fallback_is_also_gated(tmp_path):
    """The widest part of the bypass: with no pixel ROI the tool returns the
    *whole uncropped figure*, so the gate has to come first."""
    idx = _make_index(tmp_path, work=_UNCLEARED)
    hash_dir = Path(idx.papers[HASH]["hash_dir"])
    figs = json.loads((hash_dir / "figures.json").read_text())
    figs["figures"][0]["rois"] = []          # force the fallback path
    (hash_dir / "figures.json").write_text(json.dumps(figs))
    out = get_figure_roi_image(HASH, "docling_1", "A", profile="manuscript")
    assert out.get("code") == "forbidden", out
    assert "image_path" not in out


def test_roi_image_gate_matches_get_figure_image(tmp_path):
    """The three pixel paths must agree — disagreement *is* the bug."""
    _make_index(tmp_path, work=_UNCLEARED)
    with pytest.raises(ValueError):
        get_figure_image(HASH, "docling_1", profile="manuscript")
    assert get_figure_roi_image(
        HASH, "docling_1", "A", profile="manuscript",
    ).get("code") == "forbidden"


def test_roi_image_rejects_an_unknown_profile(tmp_path):
    _make_index(tmp_path, work=_CLEARED)
    out = get_figure_roi_image(HASH, "docling_1", "A", profile="thesis")
    assert out["code"] == "invalid_argument"


# --- refusals explain themselves ---------------------------------------------


def test_refusal_names_the_state_and_says_it_is_not_a_prohibition(tmp_path):
    _make_index(tmp_path, work=_UNCLEARED)
    with pytest.raises(ValueError) as exc:
        get_figure_image(HASH, "docling_1", profile="manuscript")
    msg = str(exc.value)
    assert "no_record" in msg, msg
    assert "ABSENCE of evidence" in msg, msg


def test_refusal_for_a_real_restriction_says_so(tmp_path):
    _make_index(tmp_path, work={
        "publishable": 0, "license": "all-rights-reserved",
        "license_source": "bibtex", "title": "A paper", "year": 2010,
        "work_id": "w1",
    })
    with pytest.raises(ValueError) as exc:
        get_figure_image(HASH, "docling_1", profile="manuscript")
    msg = str(exc.value)
    assert "restricted" in msg, msg
    assert "forbids republication" in msg, msg
