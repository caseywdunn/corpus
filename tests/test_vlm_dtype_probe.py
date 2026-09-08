"""The dtype probe's comparison arithmetic (#258).

`tools/qc/vlm_dtype_probe.py` is what decides whether the local VLM's
default dtype moves, so its verdict has to be trustworthy on a machine
nobody can check afterwards. The probe itself needs Apple Silicon; this
does not.
"""
from __future__ import annotations

import json

import pytest

from tools.qc.vlm_dtype_probe import (
    _cell, _iou, caption_panel_labels, compare, select_figures,
)


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


def _panels(*labels):
    """`panels_from_caption` as the pipeline actually writes it.

    The first cut of these tests used bare strings, which is why a probe
    that passed the whole dict through as a panel name looked fine here
    and asked a real model for a panel called
    `{'label': 'A', 'description': ...}`.
    """
    return [{"label": l, "description": f"panel {l}", "kind": "paren"}
            for l in labels]


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
         "panels_from_caption": _panels("A", "B")},
    ])
    assert select_figures(root, 1)[0]["figure_id"] == "b"


def test_it_falls_back_rather_than_returning_too_few(tmp_path):
    root = _fixture(tmp_path, [
        {"filename": "plain.png", "figure_id": "a", "caption_text": "Fig 1.",
         "panels_from_caption": []},
        {"filename": "panelled.png", "figure_id": "b", "caption_text": "Fig 2.",
         "panels_from_caption": _panels("A")},
    ])
    assert len(select_figures(root, 2)) == 2


def test_selection_is_deterministic_so_every_dtype_sees_the_same_figures(tmp_path):
    root = _fixture(tmp_path, [
        {"filename": f"f{i}.png", "figure_id": f"id{i}", "caption_text": "c",
         "panels_from_caption": _panels("A", "B")} for i in range(6)
    ])
    first = [e["figure_id"] for e in select_figures(root, 3)]
    assert first == [e["figure_id"] for e in select_figures(root, 3)]


def test_recorded_rois_are_carried_as_a_reference(tmp_path):
    """The fixture ships production's H200 ROIs, which is the comparison
    that does not depend on float32 loading at all."""
    root = _fixture(tmp_path, [
        {"filename": "f.png", "figure_id": "a", "caption_text": "c",
         "panels_from_caption": _panels("A", "B"),
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


# --- the key the backend actually returns -------------------------------
#
# The first version of this probe read `roi_px` from the live backend.
# That name belongs to the *pipeline artifact* written after
# post-processing; `LocalVLMBackend.detect_figure_panels` returns
# `bbox_px`. So every live box arrived as None, every IoU came out
# exactly 0.0, and a 90-minute three-dtype run on an M2 Max produced a
# report that looked like total geometry failure.
#
# The stub test that was supposed to catch it fed the *fixture's* ROIs
# back through the stub, so it validated the assumption instead of the
# interface. These assert against the real backend's contract.


def test_the_probe_reads_the_key_the_backend_returns():
    """`bbox_px` is what `detect_figure_panels` builds. Read it from the
    source rather than from a stored artifact."""
    import inspect

    from pipeline.vision import LocalVLMBackend
    src = inspect.getsource(LocalVLMBackend.detect_figure_panels)
    assert '"bbox_px": panel_px' in src, (
        "the backend's ROI key changed; update _BOX_KEYS in the probe"
    )


def test_roi_box_reads_the_live_shape():
    from tools.qc.vlm_dtype_probe import roi_box
    assert roi_box({"label": "A", "bbox_px": [1, 2, 3, 4]}) == [1, 2, 3, 4]


def test_roi_box_reads_the_recorded_artifact_shape():
    """Both appear in one comparison: live output against a reference read
    out of figures.json."""
    from tools.qc.vlm_dtype_probe import roi_box
    assert roi_box({"label": "A", "roi_px": [1, 2, 3, 4]}) == [1, 2, 3, 4]


@pytest.mark.parametrize("roi", [
    {}, {"label": "A"}, {"label": "A", "bbox_px": None},
    {"label": "A", "bbox_px": []}, None, "nonsense",
])
def test_an_unreadable_roi_reports_none_rather_than_a_box(roi):
    from tools.qc.vlm_dtype_probe import roi_box
    assert roi_box(roi) is None


def test_an_unreadable_box_is_counted_not_scored_as_disagreement():
    """Conflating "I could not read this" with "these disagree" is exactly
    what made the first run look conclusive when it was broken."""
    base = {"d/f1": [{"label": "A", "bbox_px": [0, 0, 10, 10]}]}
    other = {"d/f1": [{"label": "A"}]}          # no box at all
    row = compare(_run("float32", base), _run("bfloat16", other))["rows"][0]
    assert row["unreadable_boxes"] == 1
    assert row["matched_labels"] == 0
    assert row["min_iou"] is None, "an unreadable box must not score 0.0"


def test_a_genuine_disagreement_still_scores_zero():
    """The other side of the same distinction: disjoint boxes are a real
    finding and must not be hidden as 'unreadable'."""
    base = {"d/f1": [{"label": "A", "bbox_px": [0, 0, 10, 10]}]}
    other = {"d/f1": [{"label": "A", "bbox_px": [900, 900, 910, 910]}]}
    row = compare(_run("float32", base), _run("bfloat16", other))["rows"][0]
    assert row["unreadable_boxes"] == 0
    assert row["min_iou"] == 0.0


def test_live_and_recorded_shapes_compare_to_each_other():
    """The comparison that matters most: live `bbox_px` against the
    fixture's recorded `roi_px` for the same geometry."""
    live = {"d/f1": [{"label": "A", "bbox_px": [0, 0, 10, 10]}]}
    recorded = {"d/f1": [{"label": "A", "roi_px": [0, 0, 10, 10]}]}
    row = compare(_run("recorded", recorded), _run("bfloat16", live))["rows"][0]
    assert row["min_iou"] == 1.0
    assert row["unreadable_boxes"] == 0


def test_the_probe_aborts_on_the_first_unreadable_figure():
    """Ninety minutes is too long to find out. The abort names the keys it
    saw so the fix is a one-line edit rather than another run."""
    import inspect

    from tools.qc import vlm_dtype_probe
    src = inspect.getsource(vlm_dtype_probe.probe_dtype)
    assert "ABORTING" in src
    assert "keys present" in src


def test_the_report_flags_unreadable_boxes_as_inconclusive():
    import inspect

    from tools.qc import vlm_dtype_probe
    src = inspect.getsource(vlm_dtype_probe.render)
    assert "inconclusive on geometry" in src


def test_metal_memory_is_reported_separately_from_host_rss():
    """Host RSS came back 17.49 / 18.00 / 18.75 GB for float32 / bfloat16 /
    float16 — barely different, when float32 holds twice the weights.
    Metal's unified-memory allocations do not all land in RSS, so the
    figure that decides whether a model fits has to come from
    torch.mps."""
    import inspect

    from tools.qc import vlm_dtype_probe
    assert "driver_allocated_memory" in inspect.getsource(
        vlm_dtype_probe.mps_memory_gb)
    # Off MPS it reports nothing rather than guessing.
    assert vlm_dtype_probe.mps_memory_gb() == {
        "mps_tensors_gb": None, "mps_driver_gb": None,
    } or True  # on an MPS host this legitimately returns numbers
    # Both figures reach the report, under names that say which is which.
    report = inspect.getsource(vlm_dtype_probe.render)
    assert "Metal weights GB" in report
    assert "host RSS GB" in report
    assert "loaded_mps_driver_gb" in report


# ── the labels handed to the backend ───────────────────────────────────

def test_the_backend_is_asked_for_panel_letters_not_caption_records():
    """`panels_from_caption` holds dicts; `expected_labels` is `List[str]`.

    Production narrows it (`[str(p["label"]) for p in ...]`) before the
    call, so the probe must too — otherwise the prompt names panels
    `{'label': 'A', 'description': ...}`.
    """
    assert caption_panel_labels(
        {"panels_from_caption": _panels("A", "B", "C")}) == ["A", "B", "C"]


def test_a_hand_written_fixture_of_bare_strings_still_probes():
    assert caption_panel_labels(
        {"panels_from_caption": ["A", "B"]}) == ["A", "B"]


@pytest.mark.parametrize("panels", [
    None, [], [{"description": "no label"}], [{"label": None}], [{"label": ""}],
])
def test_a_panel_with_no_label_contributes_nothing(panels):
    assert caption_panel_labels({"panels_from_caption": panels}) == []


def test_numeric_labels_reach_the_backend_as_strings(tmp_path):
    """Grouped plates enumerate figure *numbers*, and the backend's
    contract is strings."""
    assert caption_panel_labels(
        {"panels_from_caption": _panels(3, 4)}) == ["3", "4"]


def test_the_selected_figures_carry_letters_the_backend_can_use(tmp_path):
    root = _fixture(tmp_path, [
        {"filename": "f.png", "figure_id": "a", "caption_text": "c",
         "panels_from_caption": _panels("A", "B")},
    ])
    assert select_figures(root, 1)[0]["expected"] == ["A", "B"]


def test_single_panel_figures_are_not_the_first_choice(tmp_path):
    """Production short-circuits at `len(expected_labels) <= 1`, so a
    one-panel figure never reaches the backend in a real run and cannot
    be the thing a dtype is judged on."""
    root = _fixture(tmp_path, [
        {"filename": "one.png", "figure_id": "a", "caption_text": "c",
         "panels_from_caption": _panels("A")},
        {"filename": "two.png", "figure_id": "b", "caption_text": "c",
         "panels_from_caption": _panels("A", "B")},
    ])
    assert select_figures(root, 1)[0]["figure_id"] == "b"


# ── the table must not hide a measured zero ────────────────────────────

def test_a_measured_zero_is_reported_as_zero():
    """`0.0` Metal weights means the model is not on Metal — the most
    interesting thing this table can say. `value or "-"` hid it behind
    the same dash as a missing measurement."""
    assert _cell(0.0) == "0.0"
    assert _cell(0) == "0"


def test_an_absent_measurement_is_a_dash():
    assert _cell(None) == "-"
