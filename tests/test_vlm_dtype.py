"""The local VLM's weight dtype is selectable (#258).

`LocalVLMBackend` picked its dtype from whether the device was CUDA, not
from what the device supports:

    dtype = torch.bfloat16 if self._device == "cuda" else torch.float32

so MPS got float32 — ~30.4 GB of weights for Qwen2.5-VL-7B against
~15.2 GB at bfloat16, which is the difference between "does not fit on a
32 GB Mac" and "fits". Apple Silicon supports bf16 and fp16, and the
branch keying on `== "cuda"` is the shape of "CUDA is the one I tested"
rather than a statement about MPS.

**The default is deliberately unchanged.** A mocked dtype-selection test
does not establish that Qwen2.5-VL is numerically sound in half precision
on Metal, and shipping a silently worse panel detector is the failure
class this cycle spent itself on. So the dtype is *selectable* and the
default moves when someone with the hardware reports numbers — which also
settles the open sub-question of whether MPS prefers float16 to bfloat16
for this model, since all three become one run each rather than three
source edits.
"""
from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from pipeline.config_schema import FiguresConfig          # noqa: E402
from pipeline.vision import (                              # noqa: E402
    VisionBackendError,
    estimated_vlm_weight_gb,
    resolve_vlm_dtype,
)


# --- today's behaviour is preserved -------------------------------------


def test_auto_is_unchanged_on_every_device():
    """The whole point of not moving the default yet."""
    assert resolve_vlm_dtype("cuda", "auto") is torch.bfloat16
    assert resolve_vlm_dtype("mps", "auto") is torch.float32
    assert resolve_vlm_dtype("cpu", "auto") is torch.float32


def test_auto_is_the_default():
    assert FiguresConfig().vision_dtype == "auto"
    assert resolve_vlm_dtype("mps") is torch.float32


# --- what an Apple Silicon owner can now try ----------------------------


@pytest.mark.parametrize("setting,expected", [
    ("float32", torch.float32),
    ("float16", torch.float16),
    ("bfloat16", torch.bfloat16),
])
def test_an_explicit_dtype_is_honoured_on_mps(setting, expected):
    assert resolve_vlm_dtype("mps", setting) is expected


def test_an_explicit_dtype_is_honoured_on_cuda_too():
    """So a CUDA host can reproduce a Mac's dtype when comparing output."""
    assert resolve_vlm_dtype("cuda", "float32") is torch.float32


def test_case_and_whitespace_are_forgiven():
    assert resolve_vlm_dtype("mps", " BFloat16 ") is torch.bfloat16


@pytest.mark.parametrize("bad", ["fp16", "half", "int8", "float64"])
def test_a_dtype_it_cannot_load_is_refused_by_name(bad):
    """Named in the error, so a typo does not silently fall back to
    float32 and quietly need 30 GB."""
    with pytest.raises(VisionBackendError) as exc:
        resolve_vlm_dtype("mps", bad)
    assert "vision_dtype" in str(exc.value)
    assert bad in str(exc.value)


def test_the_config_rejects_an_unknown_dtype():
    with pytest.raises(Exception):
        FiguresConfig(vision_dtype="fp16")


# --- the footprint that decides whether it fits -------------------------


def test_the_reported_footprint_matches_the_measured_figures():
    """#258 measured ~30.4 GB at float32 and ~15.2 GB at bfloat16 for
    Qwen2.5-VL-7B. A round 7B would report 28 and 14, so the parameter
    count has to be the real one for the number to be worth printing."""
    assert estimated_vlm_weight_gb(torch.float32) == 30.4
    assert estimated_vlm_weight_gb(torch.bfloat16) == 15.2
    assert estimated_vlm_weight_gb(torch.float16) == 15.2


def test_half_precision_halves_it():
    assert (estimated_vlm_weight_gb(torch.bfloat16)
            == pytest.approx(estimated_vlm_weight_gb(torch.float32) / 2, rel=0.01))


# --- precedence and wiring ---------------------------------------------


def test_precedence_is_flag_then_env_then_config():
    """The env var exists so all three dtypes can be tried in one sitting
    without editing a corpuscle's config.yaml."""
    import inspect

    from pipeline.vision import LocalVLMBackend
    src = inspect.getsource(LocalVLMBackend.__init__)
    order = [src.index(x) for x in
             ("dtype\n", 'os.environ.get("CORPUS_VLM_DTYPE")',
              '"vision_dtype"')]
    assert order == sorted(order), "flag must beat env, env must beat config"


def test_the_loader_uses_the_resolver_not_a_device_comparison():
    """Guard the thing that was wrong: a `== "cuda"` dtype branch."""
    import inspect

    from pipeline import vision
    src = inspect.getsource(vision.LocalVLMBackend)
    assert "resolve_vlm_dtype(self._device, self._dtype_setting)" in src
    assert 'torch.bfloat16 if self._device == "cuda"' not in src


def test_the_footprint_is_logged_at_load():
    """A machine that cannot hold the model should say so with a number
    rather than dying in the allocator."""
    import inspect

    from pipeline import vision
    src = inspect.getsource(vision.LocalVLMBackend)
    assert "estimated_vlm_weight_gb(dtype)" in src
    assert "GB of weights" in src


def test_device_map_still_only_applies_to_cuda():
    """Not part of the dtype question: on MPS the loader correctly passes
    device_map=None and then .to(device), so that branch stays."""
    import inspect

    from pipeline import vision
    src = inspect.getsource(vision.LocalVLMBackend)
    assert 'device_map=(self._device if self._device == "cuda" else None)' in src
