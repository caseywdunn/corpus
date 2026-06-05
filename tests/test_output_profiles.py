"""Output-profile figure-licensing gate (#101).

The gate is keyed to the per-call ``profile=`` (client/session property)
with a server ``--default-profile`` fallback. A ``strict`` profile
(manuscript / presentation) refuses an unpublishable figure; the default
permissive ``report`` allows it. The same index serving two calls with
different ``profile=`` must yield independent outcomes — the core
per-session regression for the shared-SSE deploy.
"""
from __future__ import annotations

import json
import types
from pathlib import Path

import pytest
from PIL import Image as PILImage

from mcpsrv import app as mcp_app
from mcpsrv.profiles import BUILTIN_PROFILES, get_profile, resolve_profile
from mcpsrv.tools.figures import get_figure_image, get_figure_url
from mcpsrv.tools.profiles import get_active_profile, list_output_profiles


class _StubBiblio:
    def __init__(self, publishable: bool):
        self._pub = publishable

    def get_work_by_corpus_hash(self, _h):
        return {
            "publishable": 1 if self._pub else 0,
            "license": "CC-BY-4.0" if self._pub else "unknown",
            "license_source": "bibtex" if self._pub else "unknown",
            "title": "A paper", "year": 2010,
        }


def _make_index(tmp_path: Path, *, publishable: bool, default_profile="report"):
    h = "aaaaaaaaaaaa"
    hash_dir = tmp_path / "documents" / h
    (hash_dir / "figures").mkdir(parents=True)
    PILImage.new("RGB", (4, 4), "white").save(hash_dir / "figures" / "fig_1.png")
    (hash_dir / "figures.json").write_text(json.dumps({"figures": [
        {"figure_id": "docling_1", "filename": "fig_1.png",
         "figure_type": "figure", "caption_text": "Figure 1. A thing."},
    ]}))
    return types.SimpleNamespace(
        papers={h: {"hash_dir": str(hash_dir), "title": "A paper"}},
        biblio_db=_StubBiblio(publishable),
        default_profile=default_profile,
        figure_url_base="http://127.0.0.1:9999",
        figure_auth_token="tok",
    ), h


@pytest.fixture(autouse=True)
def _restore():
    original = mcp_app._INDEX
    yield
    mcp_app.set_index(original)


# --- profile vocabulary ----------------------------------------------------


def test_builtin_vocabulary_and_resolution():
    assert set(BUILTIN_PROFILES) == {"report", "manuscript", "presentation"}
    assert BUILTIN_PROFILES["report"].figure_licensing == "permissive"
    assert BUILTIN_PROFILES["manuscript"].figure_licensing == "strict"
    # per-call wins; unknown falls through to server default then report.
    assert resolve_profile("manuscript", "report").name == "manuscript"
    assert resolve_profile(None, "manuscript").name == "manuscript"
    assert resolve_profile(None, None).name == "report"
    assert get_profile("nope") is None


# --- get_figure_image gate -------------------------------------------------


def test_strict_profile_refuses_unpublishable(tmp_path):
    idx, h = _make_index(tmp_path, publishable=False)
    mcp_app.set_index(idx)
    with pytest.raises(ValueError, match="not publishable under profile 'manuscript'"):
        get_figure_image(h, "docling_1", profile="manuscript")


def test_permissive_profile_allows_unpublishable(tmp_path):
    idx, h = _make_index(tmp_path, publishable=False)
    mcp_app.set_index(idx)
    img = get_figure_image(h, "docling_1", profile="report")
    assert img is not None  # returned, not raised


def test_default_fallback_is_permissive(tmp_path):
    idx, h = _make_index(tmp_path, publishable=False, default_profile="report")
    mcp_app.set_index(idx)
    assert get_figure_image(h, "docling_1") is not None  # no per-call profile


def test_server_default_manuscript_gates_calls_without_profile(tmp_path):
    idx, h = _make_index(tmp_path, publishable=False, default_profile="manuscript")
    mcp_app.set_index(idx)
    with pytest.raises(ValueError, match="not publishable"):
        get_figure_image(h, "docling_1")  # falls back to strict default


def test_publishable_passes_strict(tmp_path):
    idx, h = _make_index(tmp_path, publishable=True)
    mcp_app.set_index(idx)
    assert get_figure_image(h, "docling_1", profile="manuscript") is not None


def test_unknown_profile_raises(tmp_path):
    idx, h = _make_index(tmp_path, publishable=True)
    mcp_app.set_index(idx)
    with pytest.raises(ValueError, match="unknown profile"):
        get_figure_image(h, "docling_1", profile="thesis")


def test_per_call_profiles_are_independent_on_one_index(tmp_path):
    # The shared-SSE invariant: same server/index, two sessions passing
    # different profile= get independent gate outcomes.
    idx, h = _make_index(tmp_path, publishable=False)
    mcp_app.set_index(idx)
    assert get_figure_image(h, "docling_1", profile="report") is not None
    with pytest.raises(ValueError):
        get_figure_image(h, "docling_1", profile="manuscript")


# --- get_figure_url --------------------------------------------------------


def test_url_carries_resolved_profile_and_gates(tmp_path):
    idx, h = _make_index(tmp_path, publishable=True)
    mcp_app.set_index(idx)
    out = get_figure_url(h, "docling_1", profile="manuscript")
    assert out["profile"] == "manuscript"
    assert "profile=manuscript" in out["url"]
    assert out["attribution"] is not None or out["attribution"] is None  # field present


def test_url_strict_refuses_unpublishable(tmp_path):
    idx, h = _make_index(tmp_path, publishable=False)
    mcp_app.set_index(idx)
    out = get_figure_url(h, "docling_1", profile="manuscript")
    assert "not publishable" in out["error"]
    assert "url" not in out


def test_url_unknown_profile_error(tmp_path):
    idx, h = _make_index(tmp_path, publishable=True)
    mcp_app.set_index(idx)
    out = get_figure_url(h, "docling_1", profile="thesis")
    assert out["error"].startswith("unknown profile")


# --- discovery tools -------------------------------------------------------


def test_list_output_profiles(tmp_path):
    idx, _ = _make_index(tmp_path, publishable=True, default_profile="report")
    mcp_app.set_index(idx)
    out = list_output_profiles()
    assert out["server_default"] == "report"
    names = {p["name"] for p in out["profiles"]}
    assert names == {"report", "manuscript", "presentation"}


def test_get_active_profile_reports_server_default(tmp_path):
    idx, _ = _make_index(tmp_path, publishable=True, default_profile="manuscript")
    mcp_app.set_index(idx)
    out = get_active_profile()
    assert out["name"] == "manuscript"
    assert out["is_server_default"] is True
