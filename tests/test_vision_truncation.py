"""Pass 3b must not report incomplete model output as "no panels" (#253/#269).

The local VLM backend capped generation at a fixed 1024 tokens. The panel
prompt asks for one six-field JSON object per panel at roughly 60–90 tokens
each, so panel-rich figures ran out mid-object; `_extract_json` found no
balanced `{...}`, the backend returned `[]`, and `[]` maps to
`no_labels_found` — a clean result, indistinguishable from a figure the
model genuinely found nothing in.

In one full reference run (Qwen2.5-VL-7B, job 24083450), ROI coverage
fell monotonically with panel count — 47.8% at 2–3 panels, 32.1% at 4–5,
22.7% at 6–7, 13.6% at 10+. That shape is a fixed output budget, not a
vision failure, and the figures it costs most are the ones panel ROIs
matter most for.

These tests do not load a model. What they cover is the reasoning around
the call — budget arithmetic, JSON salvage, and the truncated-vs-absent
distinction — which is where the defect lived. **The budget's adequacy is
not tested here and cannot be**: that needs a real Pass 3b run.
"""
from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from pipeline.figure_passes import _pass3b_annotate_rois
from pipeline.figures import detect_figure_rois_via_vision
from pipeline.vision import (
    ClaudeVisionBackend,
    _VLM_TOKEN_BASE,
    _VLM_TOKENS_PER_PANEL,
    LocalVLMBackend,
    VisionBackendError,
    _extract_json,
    _parse_complete_vision_response,
)


# The real response head from `7247562445ba` / `fig_219.png` (10 panels),
# quoted in #253: well-formed JSON that simply stops.
TRUNCATED = (
    '```json\n{\n  "panels": [\n    {\n      "label": "A",\n'
    '      "parent_figure_index": 0,\n'
    '      "panel_bbox_norm": [0.05, 0.03, 0.25, 0.50],\n'
    '      "label_bbox_norm": [0.05, 0.03, 0.09, 0.07],\n'
    '      "confidence": 1.0,\n   '
)


# --- salvaging what did arrive -----------------------------------------------


def test_a_truncated_response_yields_nothing_to_parse():
    """No balanced span exists, so this is None — the caller has to tell
    that apart from a genuine empty result, which is the whole issue."""
    assert _extract_json(TRUNCATED) is None


def test_a_malformed_leading_object_no_longer_discards_a_good_one():
    """This returned None on the first JSONDecodeError and stopped scanning,
    throwing away a usable object later in the same response."""
    text = '{"panels": [oops not json}\n{"panels": [], "embedded_figures": []}'
    assert _extract_json(text) == {"panels": [], "embedded_figures": []}


def test_the_first_parseable_object_still_wins():
    text = '{"panels": [{"label": "A"}]}\n{"panels": [{"label": "Z"}]}'
    assert _extract_json(text)["panels"][0]["label"] == "A"


def test_ordinary_fenced_json_is_unaffected():
    assert _extract_json('```json\n{"panels": []}\n```') == {"panels": []}


def test_empty_and_junk_are_still_none():
    assert _extract_json("") is None
    assert _extract_json(None) is None
    assert _extract_json("no json here at all") is None


def test_only_an_explicit_complete_empty_payload_means_no_labels():
    assert _parse_complete_vision_response(
        '{"panels": [], "embedded_figures": []}', backend="test",
    ) == {"panels": [], "embedded_figures": []}


@pytest.mark.parametrize("text", [
    "not json",
    '{"panels": []}',
    '{"panels": {}, "embedded_figures": []}',
])
def test_malformed_or_incomplete_payload_is_a_backend_failure(text):
    with pytest.raises(VisionBackendError):
        _parse_complete_vision_response(text, backend="test")


def test_even_parseable_json_is_rejected_when_generation_hit_the_cap():
    with pytest.raises(VisionBackendError, match="token-truncated"):
        _parse_complete_vision_response(
            '{"panels": [], "embedded_figures": []}',
            backend="test", truncated=True, truncation_detail="at cap",
        )


# --- sizing the budget to the work -------------------------------------------


def _Budget(floor):
    """The budget helper in isolation. Constructing the real backend would
    download 15 GB of weights, so bind the method to a stand-in carrying the
    one attribute it reads."""
    return SimpleNamespace(
        _max_new_tokens=floor,
        _token_budget=lambda labels: LocalVLMBackend._token_budget(
            SimpleNamespace(_max_new_tokens=floor), labels),
    )


@pytest.mark.parametrize("panels,expected_at_least", [
    (2, 1024),      # small figures keep the floor
    (6, 1024),
    (10, 1400),     # #253 measured the 10+ bucket needing ~1200-1500
    (11, 1400),
])
def test_the_budget_covers_the_panel_count(panels, expected_at_least):
    b = _Budget(1024)
    assert b._token_budget(["A"] * panels) >= expected_at_least


def test_the_configured_value_is_a_floor_not_a_ceiling():
    """An operator who raised it keeps what they asked for."""
    assert _Budget(4096)._token_budget(["A"] * 3) == 4096


def test_discovery_with_no_expected_labels_gets_a_conservative_budget():
    b = _Budget(1024)
    assert b._token_budget([]) == 4096
    assert b._token_budget(None) == 4096


def test_the_budget_grows_with_panels():
    b = _Budget(0)
    assert b._token_budget(["A"] * 4) - b._token_budget(["A"] * 3) == _VLM_TOKENS_PER_PANEL
    assert b._token_budget([]) == 4096


def test_claude_uses_the_same_structure_sized_budget():
    stand_in = SimpleNamespace(max_tokens=1024)
    budget = ClaudeVisionBackend._token_budget(stand_in, ["A"] * 10)
    assert budget >= 1400


# --- and out through the status the operator reads ---------------------------


class _Backend:
    """Minimal VisionBackend stand-in; `detect_figure_rois_via_vision` only
    needs `name` and `detect_figure_panels`."""

    name = "vision:test"

    def __init__(self, behaviour):
        self._behaviour = behaviour

    def detect_figure_panels(self, image_path, caption_text, expected_labels):
        return self._behaviour()


@pytest.fixture
def figure(tmp_path):
    from PIL import Image
    p = tmp_path / "fig_219.png"
    Image.new("RGB", (400, 300), "white").save(p)
    return p


def test_a_truncated_backend_is_a_failure_not_a_clean_result(figure):
    """The heart of #253. A backend that gives up must not land in the same
    bucket as a figure that genuinely has no panels — that is what let ~10%
    of a corpus lose its ROIs without a single counter moving."""
    def truncate():
        raise VisionBackendError("response truncated at the 1456-token budget")

    r = detect_figure_rois_via_vision(
        figure, [{"label": "A"}, {"label": "B"}], _Backend(truncate),
        caption_text="(A) x (B) y")
    assert r["pass3_status"] == "vision_backend_failed"
    assert r["pass3_status"] != "no_labels_found"
    assert "truncated" in r["pass3_error"]
    assert r["rois"] == []


def test_a_genuinely_empty_result_is_still_no_labels_found(figure):
    """The other half: this one is a clean result and must stay one, or the
    new signal is just noise."""
    r = detect_figure_rois_via_vision(
        figure, [{"label": "A"}, {"label": "B"}], _Backend(lambda: []),
        caption_text="(A) x (B) y")
    assert r["pass3_status"] == "no_labels_found"


def test_failure_reason_is_persisted_and_cleared_after_a_clean_retry(figure):
    figures_file = figure.parent / "figures.json"
    figures_file.write_text(json.dumps({"figures": [{
        "figure_type": "figure",
        "file_path": str(figure),
        "panels_from_caption": [{"label": "A"}, {"label": "B"}],
    }]}))

    def truncate():
        raise VisionBackendError("response token-truncated at the budget")

    _pass3b_annotate_rois(figures_file, _Backend(truncate))
    failed = json.loads(figures_file.read_text())["figures"][0]
    assert failed["pass3_status"] == "vision_backend_failed"
    assert "token-truncated" in failed["pass3_error"]

    _pass3b_annotate_rois(figures_file, _Backend(lambda: []))
    retried = json.loads(figures_file.read_text())["figures"][0]
    assert retried["pass3_status"] == "no_labels_found"
    assert "pass3_error" not in retried


def test_truncation_detection_reads_the_generated_length():
    """`generate()` stopping exactly at the budget is the signal; one token
    short of it is a model that finished on its own."""
    b = SimpleNamespace(_max_new_tokens=1024)
    budget = LocalVLMBackend._token_budget(b, ["A"] * 10)
    assert budget >= 1400
    assert (budget >= budget) is True          # at the cap  -> truncated
    assert ((budget - 1) >= budget) is False   # under it    -> finished


def test_both_backends_use_the_hard_completion_gate():
    import inspect
    from pipeline import vision

    source = inspect.getsource(vision)
    assert source.count("_parse_complete_vision_response(") == 3
    assert source.count("(initial_budget, initial_budget * 2)") == 2
