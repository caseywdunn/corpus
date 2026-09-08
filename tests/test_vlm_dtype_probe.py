"""The dtype probe's comparison arithmetic (#258).

`tools/qc/vlm_dtype_probe.py` is what decides whether the local VLM's
default dtype moves, so its verdict has to be trustworthy on a machine
nobody can check afterwards. The probe itself needs Apple Silicon; this
does not.
"""
from __future__ import annotations

import json

import pytest

from tools.qc.vlm_dtype_probe import _iou, compare, select_figures


# --- IoU ----------------------------------------------------------------


def test_identical_boxes_agree_completely():
    assert _iou([0, 0, 10, 10], [0, 0, 10, 10]) == 1.0


def test_a_small_shift_is_detected_but_not_overstated():
    """A 3 px shift on a ~1000 px panel is the scale of difference this is
    looking for: real, and not a reason to call the panel lost."""
    assert 0.94 < _iou([0, 0, 1000, 1000], [3, 3, 1003, 1003]) < 1.0


def test_disjoint_boxes_score_zero():
    assert _iou([0, 0, 10, 10], [100, 100, 110, 110]) == 0.0


@pytest.mark.parametrize("bad", [None, [], [1, 2], "nonsense"])
def test_a_missing_or_malformed_box_scores_zero_rather_than_raising(bad):
    """A dtype that returns a malformed ROI must show up as disagreement,
    not as a traceback that loses the whole run."""
    assert _iou([0, 0, 10, 10], bad) == 0.0
    assert _iou(bad, [0, 0, 10, 10]) == 0.0


# --- per-figure comparison ---------------------------------------------


def _run(dtype, figures):
    return {"dtype": dtype, "loaded": True,
            "figures": {k: {"rois": v, "error": None} for k, v in figures.items()}}


def test_an_exact_match_reports_as_such():
    rois = [{"label": "A", "roi_px": [0, 0, 10, 10]},
            {"label": "B", "roi_px": [10, 0, 20, 10]}]
    cmp = compare(_run("float32", {"d/f1": rois}), _run("bfloat16", {"d/f1": rois}))
    assert cmp["same_roi_count"] == 1
    assert cmp["tight_boxes"] == 1
    assert cmp["rows"][0]["min_iou"] == 1.0


def test_a_lost_panel_is_visible_as_a_count_difference():
    """The failure that matters most: half precision quietly finding fewer
    panels. The corpuscle still builds and no gate fires."""
    base = [{"label": "A", "roi_px": [0, 0, 10, 10]},
            {"label": "B", "roi_px": [10, 0, 20, 10]}]
    cmp = compare(_run("float32", {"d/f1": base}),
                  _run("bfloat16", {"d/f1": base[:1]}))
    row = cmp["rows"][0]
    assert (row["baseline_rois"], row["other_rois"]) == (2, 1)
    assert cmp["same_roi_count"] == 0


def test_boxes_are_matched_by_label_not_by_position():
    """"Panel B moved" is a different statement from "there is one fewer
    box", and ordering must not turn one into the other."""
    base = [{"label": "A", "roi_px": [0, 0, 10, 10]},
            {"label": "B", "roi_px": [50, 50, 60, 60]}]
    reordered = [base[1], base[0]]
    cmp = compare(_run("float32", {"d/f1": base}), _run("bfloat16", {"d/f1": reordered}))
    assert cmp["rows"][0]["min_iou"] == 1.0


def test_an_unlabelled_panel_is_not_silently_matched():
    base = [{"label": "A", "roi_px": [0, 0, 10, 10]}]
    other = [{"label": None, "roi_px": [0, 0, 10, 10]}]
    row = compare(_run("float32", {"d/f1": base}), _run("bfloat16", {"d/f1": other}))["rows"][0]
    assert row["matched_labels"] == 0
    assert row["min_iou"] is None


def test_a_failed_figure_carries_its_error_into_the_row():
    base = {"d/f1": [{"label": "A", "roi_px": [0, 0, 10, 10]}]}
    other = _run("bfloat16", {})
    other["figures"]["d/f1"] = {"rois": None, "error": "RuntimeError: NaN"}
    row = compare(_run("float32", base), other)["rows"][0]
    assert "NaN" in row["error"]
    assert row["other_rois"] == 0


# --- fixture selection --------------------------------------------------


def _fixture(tmp_path, figures):
    d = tmp_path / "documents" / "aaaaaaaaaaaa"
    (d / "figures").mkdir(parents=True)
    for f in figures:
        (d / "figures" / f["filename"]).write_bytes(b"\x89PNG\r\n\x1a\n")
    (d / "figures.json").write_text(json.dumps({"figures": figures}))
    return tmp_path


def test_figures_whose_caption_names_panels_are_preferred(tmp_path):
    """Those have an expected label set, so a lost panel reads as a miss
    rather than a judgement call."""
    root = _fixture(tmp_path, [
        {"filename": "plain.png", "figure_id": "a", "caption_text": "Fig 1.",
         "panels_from_caption": []},
        {"filename": "panelled.png", "figure_id": "b", "caption_text": "Fig 2. A, B.",
         "panels_from_caption": ["A", "B"]},
    ])
    assert select_figures(root, 1)[0]["figure_id"] == "b"


def test_it_falls_back_rather_than_returning_too_few(tmp_path):
    root = _fixture(tmp_path, [
        {"filename": "plain.png", "figure_id": "a", "caption_text": "Fig 1.",
         "panels_from_caption": []},
        {"filename": "panelled.png", "figure_id": "b", "caption_text": "Fig 2.",
         "panels_from_caption": ["A"]},
    ])
    assert len(select_figures(root, 2)) == 2


def test_selection_is_deterministic_so_every_dtype_sees_the_same_figures(tmp_path):
    root = _fixture(tmp_path, [
        {"filename": f"f{i}.png", "figure_id": f"id{i}", "caption_text": "c",
         "panels_from_caption": ["A"]} for i in range(6)
    ])
    first = [e["figure_id"] for e in select_figures(root, 3)]
    assert first == [e["figure_id"] for e in select_figures(root, 3)]


def test_recorded_rois_are_carried_as_a_reference(tmp_path):
    """The fixture ships production's H200 ROIs, which is the comparison
    that does not depend on float32 loading at all."""
    root = _fixture(tmp_path, [
        {"filename": "f.png", "figure_id": "a", "caption_text": "c",
         "panels_from_caption": ["A"],
         "rois": [{"label": "A", "roi_px": [0, 0, 5, 5],
                   "source": "vision:qwen2.5-vl-7b-instruct"}]},
    ])
    entry = select_figures(root, 1)[0]
    assert entry["reference_rois"] == [{"label": "A", "roi_px": [0, 0, 5, 5]}]
    assert entry["reference_source"] == "vision:qwen2.5-vl-7b-instruct"


def test_a_directory_with_no_corpuscle_says_what_to_do(tmp_path):
    with pytest.raises(SystemExit) as exc:
        select_figures(tmp_path, 4)
    assert "corpus run --only extract" in str(exc.value)
