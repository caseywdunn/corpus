"""#281 — a vision refresh must not re-extract documents with a clean base.

`--only vision` passes `--refresh-vision`, which used to reset the figure base
unconditionally. Resetting runs a full docling conversion per document, so on
the 1775-paper reference library the GPU vision phase went from ~1h27m (which
logged zero docling conversions) to a projected 35 hours that could not finish
inside its 24-hour allocation. Only documents a prior Pass 3c actually split
need the reset; there were 10 of them in 1775.
"""
import json

from pipeline.figure_materialization import has_split_figure_state


def write(tmp_path, figures):
    p = tmp_path / "figures.json"
    p.write_text(json.dumps({"figures": figures}))
    return p


def test_a_clean_base_needs_no_reset(tmp_path):
    p = write(tmp_path, [
        {"figure_id": "f1", "pass3_status": "completed"},
        {"figure_id": "f2", "pass3_status": "no_labels_found"},
        {"figure_id": "f3", "pass3_status": "partial_ocr", "rois": [{"label": "A"}]},
    ])
    assert has_split_figure_state(p) is False


def test_compound_status_forces_a_reset(tmp_path):
    p = write(tmp_path, [{"figure_id": "f1", "pass3_status": "completed_compound"}])
    assert has_split_figure_state(p) is True
    p = write(tmp_path, [{"figure_id": "f1", "pass3_status": "partial_vision_compound"}])
    assert has_split_figure_state(p) is True


def test_a_shared_image_forces_a_reset(tmp_path):
    """Pass 3c emits sub-figure records sharing the host figure's renamed PNG."""
    p = write(tmp_path, [
        {"figure_id": "fig_3", "previous_filenames": ["fig_3.png"]},
        {"figure_id": "fig_4", "image_shared_with": "fig_3"},
    ])
    assert has_split_figure_state(p) is True


def test_each_split_marker_alone_is_enough(tmp_path):
    assert has_split_figure_state(write(tmp_path, [{"image_shared_with": "fig_1"}])) is True
    assert has_split_figure_state(write(tmp_path, [{"previous_filenames": ["a.png"]}])) is True


def test_an_unreadable_or_absent_base_resets_conservatively(tmp_path):
    """A base we cannot vouch for is one we rebuild: guessing wrong here
    silently annotates split images."""
    assert has_split_figure_state(tmp_path / "absent.json") is True
    bad = tmp_path / "figures.json"
    bad.write_text("{not json")
    assert has_split_figure_state(bad) is True


def test_no_figures_key_is_treated_as_clean(tmp_path):
    p = tmp_path / "figures.json"
    p.write_text(json.dumps({"total_figures": 0}))
    assert has_split_figure_state(p) is False
