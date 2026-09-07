"""Model bboxes arrive as pixels, not the 0–1 floats the prompt demands (#253).

Both backends multiplied by the image dimensions on the prompt's assurance
that "each coordinate is a float in 0.0 .. 1.0". Qwen2.5-VL frequently
ignores it: 130 of 142 observable responses on one cluster carried absolute
pixels, and 100% of them since 2026-05-30.

The old conversion turned that into loss two different ways, neither logged:

    [17, 808, 150, 1365]  x0 = 17*w = 8432, x1 = min(1.0,150)*w = w
                          x1 <= x0 -> dropped, silently.
    [0, 0, 487, 606]      x0 = 0, x1 = w, y1 = h -> the guard passes and
                          you get a "panel" spanning the whole figure.
                          Wrong data rather than missing data, and it lands
                          in the `completed` counter where nothing looks.

Both are covered here with the values observed in the wild. A coordinate
above 1.0 cannot be normalized, so the units are recoverable without
guessing, and repairing strictly dominates dropping.
"""
from __future__ import annotations

import pytest

from pipeline.vision import (
    _BBOX_DEGENERATE,
    _BBOX_MALFORMED,
    _BBOX_NORMALIZED,
    _BBOX_OUT_OF_RANGE,
    _BBOX_PIXELS,
    _bbox_to_px,
    VisionBackendError,
)

W, H = 496, 1400


# --- the values actually observed --------------------------------------------


@pytest.mark.parametrize("bbox", [
    [17, 808, 150, 1365],
    [50, 30, 250, 500],
    [6, 10, 496, 400],
    [294, 70, 328, 104],
    [20, 811, 45, 836],
])
def test_pixel_bboxes_are_repaired_not_dropped(bbox):
    px, units = _bbox_to_px(bbox, W, H)
    assert units == _BBOX_PIXELS
    assert px is not None, "these were all silently dropped before"
    assert px == [int(v) for v in bbox]


def test_the_whole_figure_impostor_is_gone():
    """[0, 0, 487, 606] used to pass the guard and become a panel covering
    the entire figure — the worst case, because it was counted as success."""
    px, units = _bbox_to_px([0, 0, 487, 606], W, H)
    assert units == _BBOX_PIXELS
    assert px == [0, 0, 487, 606]
    assert px != [0, 0, W, H], "this is the impostor the guard let through"


@pytest.mark.parametrize("bbox", [
    [0.0, 0.0, 0.357, 0.205],
    [0.06, 0.06, 0.22, 0.4],
    [0.0, 0.0, 0.51, 0.5],
])
def test_normalized_bboxes_still_scale_by_the_image(bbox):
    """The May 2026 responses were genuinely normalized; that path must not
    regress."""
    px, units = _bbox_to_px(bbox, W, H)
    assert units == _BBOX_NORMALIZED
    assert px == [int(bbox[0] * W), int(bbox[1] * H),
                  int(bbox[2] * W), int(bbox[3] * H)]


# --- and the things that genuinely cannot be used ----------------------------


def test_a_degenerate_normalized_box_is_still_rejected():
    px, units = _bbox_to_px([0.5, 0.5, 0.5, 0.9], W, H)
    assert px is None and units == _BBOX_DEGENERATE


def test_coordinates_beyond_the_image_are_rejected():
    px, units = _bbox_to_px([0, 0, 99999, 99999], W, H)
    assert px is None and units == _BBOX_OUT_OF_RANGE


@pytest.mark.parametrize("bbox", [None, [], [1, 2, 3], "nope", [1, 2, 3, "x"]])
def test_malformed_input_is_rejected_not_raised(bbox):
    px, units = _bbox_to_px(bbox, W, H)
    assert px is None and units == _BBOX_MALFORMED


def test_a_repaired_box_is_clamped_to_the_image():
    px, _ = _bbox_to_px([0, 0, W, H], W, H)
    assert px == [0, 0, W, H]


# --- one converter, not two ---------------------------------------------------


def test_both_backends_use_the_shared_converter():
    """It was duplicated, so this defect had to be found twice. It should not
    be possible to fix one and miss the other again."""
    import inspect
    from pipeline import vision
    src = inspect.getsource(vision)
    assert src.count("def _bbox_to_px") == 1
    assert src.count("_bbox_to_px(bbox_norm, w, h)") == 2, (
        "each backend should delegate to the shared converter"
    )


def test_the_dispositions_are_counted_for_the_log():
    """The drop paths had no logging at all, which is why the units question
    could only be answered from responses that failed to parse."""
    import inspect
    from pipeline import vision
    src = inspect.getsource(vision)
    assert src.count("_log_bbox_dispositions(bbox_counts") == 2


# --- end to end: the real response shape through a real backend ---------------


PIXEL_RESPONSE = """{"panels": [
  {"label": "A", "parent_figure_index": 0,
   "panel_bbox_norm": [17, 808, 150, 1365],
   "label_bbox_norm": [20, 811, 45, 836],
   "confidence": 1.0, "description": "left panel"},
  {"label": "B", "parent_figure_index": 0,
   "panel_bbox_norm": [294, 70, 328, 104],
   "label_bbox_norm": [294, 70, 310, 86],
   "confidence": 1.0, "description": "right panel"}
], "embedded_figures": []}"""


def _claude_backend_returning(text, tmp_path, monkeypatch, stop_reason="end_turn"):
    """A real ClaudeVisionBackend with only the network call stubbed, so the
    conversion under test is the shipped one."""
    from types import SimpleNamespace
    from PIL import Image
    from pipeline.vision import ClaudeVisionBackend

    img = tmp_path / "fig_219.png"
    Image.new("RGB", (496, 1400), "white").save(img)

    be = object.__new__(ClaudeVisionBackend)
    be.model = "stub"
    be.max_tokens = 4096
    stop_reasons = (list(stop_reason) if isinstance(stop_reason, (list, tuple))
                    else [stop_reason, stop_reason])
    budgets = []

    def create(**kw):
        reason = stop_reasons[min(len(budgets), len(stop_reasons) - 1)]
        budgets.append(kw["max_tokens"])
        return SimpleNamespace(
            content=[SimpleNamespace(text=text, type="text")],
            stop_reason=reason)

    be.client = SimpleNamespace(messages=SimpleNamespace(create=create))
    be._test_budgets = budgets
    return be, img


def test_a_pixel_coordinate_response_now_yields_rois(tmp_path, monkeypatch):
    """Before #253 this response produced zero ROIs and the figure was
    recorded no_labels_found — with no warning on any path."""
    be, img = _claude_backend_returning(PIXEL_RESPONSE, tmp_path, monkeypatch)
    rois = be.detect_figure_panels(img, "(A) left (B) right", ["A", "B"])
    panels = [r for r in rois if r.get("type") == "panel"]
    assert len(panels) == 2, "both pixel-coordinate panels should survive now"
    assert {p["label"] for p in panels} == {"A", "B"}
    for p in panels:
        x0, y0, x1, y1 = p["bbox_px"]
        assert x1 > x0 and y1 > y0
        assert [x0, y0, x1, y1] != [0, 0, 496, 1400], "not the whole figure"


def test_claude_max_tokens_stop_is_failure_even_when_json_parses(
    tmp_path, monkeypatch,
):
    be, img = _claude_backend_returning(
        PIXEL_RESPONSE, tmp_path, monkeypatch, stop_reason="max_tokens",
    )
    with pytest.raises(VisionBackendError, match="token-truncated"):
        be.detect_figure_panels(img, "(A) left (B) right", ["A", "B"])
    assert be._test_budgets == [4096, 8192]


def test_claude_retries_one_truncated_response_at_twice_the_budget(
    tmp_path, monkeypatch,
):
    be, img = _claude_backend_returning(
        PIXEL_RESPONSE, tmp_path, monkeypatch,
        stop_reason=["max_tokens", "end_turn"],
    )
    rois = be.detect_figure_panels(img, "(A) left (B) right", ["A", "B"])
    assert len([r for r in rois if r["type"] == "panel"]) == 2
    assert be._test_budgets == [4096, 8192]
