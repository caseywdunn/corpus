"""Figure-DPI config + backfill (#121).

`figures.images_scale` controls docling's render scale (saved DPI =
72 * scale); the default moved from docling's effective 1.0 (72 dpi) to
2.0 (144 dpi). `backfill_figure_dpi.py` re-renders an *existing* bundle's
figure PNGs at a target scale from the stored bbox + source PDF, without
re-running docling.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

import backfill_figure_dpi as bf
from pipeline.config_schema import CorpuscleConfig, validate_config

FIXTURE_PDF = Path(__file__).parent / "fixtures" / "round2_paper" / "Siebert_etal2011.pdf"


# --- config -------------------------------------------------------------------


def test_images_scale_default_and_bounds():
    assert CorpuscleConfig().figures.images_scale == 2.0
    # in-range value accepted
    assert validate_config({"figures": {"images_scale": 4.0}}).figures.images_scale == 4.0
    # out-of-range rejected (ge=1.0, le=8.0)
    with pytest.raises(Exception):
        validate_config({"figures": {"images_scale": 0.5}})
    with pytest.raises(Exception):
        validate_config({"figures": {"images_scale": 99.0}})


def test_pipeline_config_default_carries_images_scale():
    from pipeline.config import _DEFAULT_CONFIG
    assert _DEFAULT_CONFIG["figures"]["images_scale"] == 2.0


# --- backfill -----------------------------------------------------------------


def _make_paper(tmp_path: Path, bbox, coord_system, page=1, record_dims=False):
    """Build a minimal documents/<hash>/ with processed.pdf + figures.json."""
    hd = tmp_path / "documents" / "abc123def456"
    (hd / "figures").mkdir(parents=True)
    shutil.copy(FIXTURE_PDF, hd / "processed.pdf")
    fig = {
        "figure_id": "fig_1", "filename": "fig_1.png", "figure_type": "figure",
        "page": page, "bbox": bbox, "bbox_coord_system": coord_system,
    }
    if record_dims:
        fig["width"], fig["height"] = 100, 80
    (hd / "figures.json").write_text(json.dumps({"figures": [fig]}))
    return hd


@pytest.mark.skipif(not FIXTURE_PDF.is_file(), reason="fixture PDF missing")
def test_backfill_renders_at_target_scale(tmp_path):
    # A 200×150 pt top-left region; at scale 2.0 → ~400×300 px.
    hd = _make_paper(tmp_path, [50, 50, 250, 200], "pdf_pts_top_left",
                     record_dims=True)
    rc = bf.main([str(tmp_path), "--scale", "2.0"])
    assert rc == 0

    png = hd / "figures" / "fig_1.png"
    assert png.is_file() and png.stat().st_size > 0
    from PIL import Image
    w, h = Image.open(png).size
    assert abs(w - 400) <= 2 and abs(h - 300) <= 2  # 200pt*2, 150pt*2

    # width/height refreshed in figures.json + provenance stamped
    rec = json.loads((hd / "figures.json").read_text())["figures"][0]
    assert rec["width"] == w and rec["height"] == h
    assert rec["images_scale"] == 2.0


@pytest.mark.skipif(not FIXTURE_PDF.is_file(), reason="fixture PDF missing")
def test_backfill_bottom_left_coords_render(tmp_path):
    """A docling-style bottom-left bbox should y-flip and still render a
    non-empty crop."""
    hd = _make_paper(tmp_path, [50, 100, 250, 300], "pdf_pts_bottom_left")
    assert bf.main([str(tmp_path), "--scale", "3.0"]) == 0
    png = hd / "figures" / "fig_1.png"
    from PIL import Image
    w, h = Image.open(png).size
    assert w > 0 and h > 0
    assert abs(w - 600) <= 2  # 200pt * scale 3.0


@pytest.mark.skipif(not FIXTURE_PDF.is_file(), reason="fixture PDF missing")
def test_backfill_dry_run_writes_nothing(tmp_path):
    hd = _make_paper(tmp_path, [50, 50, 250, 200], "pdf_pts_top_left")
    before = json.loads((hd / "figures.json").read_text())
    assert bf.main([str(tmp_path), "--dry-run"]) == 0
    assert not (hd / "figures" / "fig_1.png").exists()
    assert json.loads((hd / "figures.json").read_text()) == before


def test_backfill_skips_figure_without_bbox(tmp_path):
    hd = tmp_path / "documents" / "deadbeef0001"
    (hd / "figures").mkdir(parents=True)
    shutil.copy(FIXTURE_PDF, hd / "processed.pdf")
    (hd / "figures.json").write_text(json.dumps({"figures": [
        {"figure_id": "f1", "filename": "f1.png", "figure_type": "figure"},  # no bbox
    ]}))
    assert bf.main([str(tmp_path)]) == 0
    assert not (hd / "figures" / "f1.png").exists()
