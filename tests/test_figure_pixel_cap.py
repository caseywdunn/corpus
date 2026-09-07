"""A pixel ceiling on saved figures (#184).

Figures are ~88% of a served bundle, so their size is most of what a
colleague downloads. Measured on the 1,775-document reference tree:
21,521 figures holding 11.88 GiB, median longest side 1,351 px against a
p99 of 5,697 px and a maximum of 16,237 px — the mass is a small tail.

**`max_dpi` cannot reach it, and that is why this is a separate knob.**
The byte-heavy figures are full plate pages at an ordinary 400 dpi —
11-17 MB each — so capping density to 300 would leave a 16,000 px figure
at 12,000 px. There is a second population at 2,000+ dpi, but those files
are 0.4-0.5 MB and hold no bytes.

**The flat cap wins on measurement, which is what #184 asked to settle.**
A 3000 px cap recovers 2.74 GiB (23%) and costs the median
panel-detected figure 0%, because 966 of 1,015 such figures (95%) are
already under it. Capping *only* figures with no detected panels recovers
2.71 GiB — 0.03 GiB more — while depending on `rois == 0`, which on this
tree means "ROI detection never ran" for 95% of figures rather than "has
no panels". So: flat, and dependent on nothing.
"""
from __future__ import annotations

import pytest

from pipeline.config_schema import FiguresConfig
from pipeline.extract import _cap_image_pixels
from pipeline.figures import cap_scale_to_pixels


# --- the scale arithmetic -----------------------------------------------


def test_an_oversized_figure_is_scaled_down():
    # A 2925 pt plate rendered at 400 dpi is ~16,240 px.
    scale, capped = cap_scale_to_pixels(2925, 1000, 5.55, 3000)
    assert capped is True
    assert 2925 * scale <= 3000


def test_a_figure_already_inside_the_cap_is_untouched():
    scale, capped = cap_scale_to_pixels(500, 400, 2.0, 3000)
    assert (scale, capped) == (2.0, False)


@pytest.mark.parametrize("cap", [None, 0])
def test_no_cap_means_no_change(cap):
    assert cap_scale_to_pixels(2925, 1000, 5.55, cap) == (5.55, False)


def test_the_shorter_side_does_not_trigger_the_cap():
    """The bound is on the longest side, so a tall narrow figure is
    judged on its height."""
    scale, capped = cap_scale_to_pixels(100, 2925, 5.55, 3000)
    assert capped is True
    assert 2925 * scale <= 3000


def test_the_cap_scales_the_whole_figure_not_just_one_axis():
    """Aspect ratio has to survive, or a plate's panels distort."""
    scale, _ = cap_scale_to_pixels(2000, 1000, 4.0, 3000)
    assert pytest.approx(2000 * scale / (1000 * scale), rel=1e-9) == 2.0


# --- the PIL path (fixed mode / the docling save) -----------------------


class _FakeImage:
    def __init__(self, size):
        self.size = size
        self.resized_to = None

    def resize(self, size, _filter):
        out = _FakeImage(size)
        self.resized_to = size
        return out


def test_the_docling_save_path_caps_too():
    """`native` mode re-renders from bboxes and caps there, but `fixed`
    mode keeps docling's render as-is — so without this the cap would
    silently not apply in that mode."""
    img = _FakeImage((8000, 4000))
    out = _cap_image_pixels(img, 3000, "fig_1.png")
    assert out.size == (3000, 1500)


def test_the_pil_path_leaves_a_small_figure_alone():
    img = _FakeImage((800, 400))
    assert _cap_image_pixels(img, 3000, "fig_1.png") is img
    assert img.resized_to is None


@pytest.mark.parametrize("cap", [None, 0])
def test_the_pil_path_respects_no_cap(cap):
    img = _FakeImage((8000, 4000))
    assert _cap_image_pixels(img, cap, "fig_1.png") is img


def test_a_resize_failure_returns_the_original_rather_than_losing_it():
    """A figure that cannot be downscaled should still be saved. Losing
    it to save bytes is the wrong trade."""
    class _Broken(_FakeImage):
        def resize(self, *_a):
            raise RuntimeError("no")

    img = _Broken((8000, 4000))
    assert _cap_image_pixels(img, 3000, "fig_1.png") is img


# --- configuration ------------------------------------------------------


def test_the_cap_is_on_by_default_at_3000():
    """Chosen because it is free where it matters: 95% of
    panel-detected figures are already under it."""
    assert FiguresConfig().max_pixels_long_side == 3000


def test_it_can_be_disabled():
    assert FiguresConfig(max_pixels_long_side=None).max_pixels_long_side is None


def test_a_nonsense_cap_is_refused():
    with pytest.raises(Exception):
        FiguresConfig(max_pixels_long_side=0)
    with pytest.raises(Exception):
        FiguresConfig(max_pixels_long_side=-1)


def test_it_is_distinct_from_max_dpi():
    """Conflating the two is the trap: the byte mass is at 400 dpi, so a
    density cap does not reach it."""
    cfg = FiguresConfig(max_dpi=300, max_pixels_long_side=3000)
    assert cfg.max_dpi == 300 and cfg.max_pixels_long_side == 3000


# --- wiring -------------------------------------------------------------


def test_both_write_paths_and_the_backfill_are_wired():
    import inspect

    from pipeline import extract, figures
    from tools import backfill_figure_dpi

    assert 'pixel_cap=fig_cfg.get("max_pixels_long_side")' in inspect.getsource(extract)
    assert "_cap_image_pixels(image, pixel_cap, filename)" in inspect.getsource(extract)
    assert "cap_scale_to_pixels(" in inspect.getsource(figures)
    # The backfill exists so an existing bundle can be shrunk without
    # re-running docling.
    src = inspect.getsource(backfill_figure_dpi)
    assert "--max-pixels-long-side" in src
    assert "pixel_cap=args.max_pixels_long_side" in src


def test_a_capped_render_is_recorded_as_such():
    """`resolution_mode` is the audit trail for why a figure is the size
    it is; a silently shrunk figure would be indistinguishable from a
    small source."""
    import inspect

    from pipeline import figures
    assert 'pixel_capped' in inspect.getsource(figures)
