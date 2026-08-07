"""Lossless PNG size reduction for figures (#184).

The pass exists because figures are ~97% of the served bundle and, on a
2,527-figure sample of the production corpuscle, 47.6% of them are
greyscale stored as RGB. Dropping the two redundant channels measured
-33.2% across the whole figure set.

What these tests pin down is the part that makes it *safe* rather than the
part that makes it small:

* the channel drop is bit-exact, and a figure carrying real colour is
  never touched;
* the original is always a candidate, so a figure that PIL would re-encode
  *larger* keeps its original bytes — 52.1% of sampled figures did;
* nothing throws. The pass is an optimization, so every failure mode has
  to degrade to "leave the file alone".
"""
from __future__ import annotations

import numpy as np
import pytest
from PIL import Image

from pipeline.figures import (
    _is_greyscale_rgb,
    shrink_figure_png,
    shrink_figures_dir,
)


def _rgb_from(arr: np.ndarray) -> Image.Image:
    return Image.fromarray(arr, mode="RGB")


def _grey_as_rgb(w: int, h: int, seed: int = 0) -> Image.Image:
    """A greyscale image stored in RGB — the case worth optimizing."""
    rng = np.random.default_rng(seed)
    lum = rng.integers(0, 256, size=(h, w), dtype=np.uint8)
    return _rgb_from(np.stack([lum, lum, lum], axis=-1))


def _colour(w: int, h: int, seed: int = 0) -> Image.Image:
    rng = np.random.default_rng(seed)
    return _rgb_from(rng.integers(0, 256, size=(h, w, 3), dtype=np.uint8))


# ── channel-equality detection ──────────────────────────────────────────

def test_detects_greyscale_stored_as_rgb():
    assert _is_greyscale_rgb(_grey_as_rgb(40, 30)) is True


def test_rejects_genuine_colour():
    assert _is_greyscale_rgb(_colour(40, 30)) is False


def test_a_single_coloured_pixel_defeats_the_drop():
    """Exact equality, not a tolerance. One coloured pixel — a scale bar, a
    stamp, JPEG ringing — must keep the figure in RGB, because dropping
    channels would then lose real information."""
    a = np.zeros((20, 20, 3), dtype=np.uint8)
    a[10, 10] = (255, 0, 0)
    assert _is_greyscale_rgb(_rgb_from(a)) is False


# ── the shrink itself ───────────────────────────────────────────────────

def test_greyscale_figure_shrinks_and_is_bit_exact(tmp_path):
    p = tmp_path / "fig_1.png"
    src = _grey_as_rgb(300, 200, seed=1)
    src.save(p)
    before = p.stat().st_size

    saved = shrink_figure_png(p)

    assert saved > 0, "a greyscale-in-RGB figure should get smaller"
    assert p.stat().st_size == before - saved
    with Image.open(p) as out:
        assert out.mode == "L"
        # The invariant a consumer observes: converting back to RGB
        # reproduces the original pixels exactly.
        assert np.array_equal(np.asarray(out.convert("RGB")), np.asarray(src))


def test_colour_figure_is_left_alone_or_only_shrunk_losslessly(tmp_path):
    p = tmp_path / "fig_2.png"
    src = _colour(200, 150, seed=2)
    src.save(p)

    shrink_figure_png(p)

    with Image.open(p) as out:
        assert out.mode == "RGB"
        assert np.array_equal(np.asarray(out), np.asarray(src))


def test_never_grows_a_file(tmp_path):
    """Random colour noise is the adversarial case: PIL's optimize pass
    reliably produces a *larger* PNG than the source encoder. The original
    must stay on disk. This is the rule that separates -33.2% from -18.1%
    over the real corpus."""
    p = tmp_path / "fig_3.png"
    src = _colour(400, 400, seed=3)
    src.save(p, optimize=False, compress_level=9)
    before = p.stat().st_size

    saved = shrink_figure_png(p)

    assert saved >= 0
    assert p.stat().st_size <= before


def test_already_greyscale_mode_is_handled(tmp_path):
    p = tmp_path / "fig_4.png"
    rng = np.random.default_rng(4)
    src = Image.fromarray(rng.integers(0, 256, (80, 80), dtype=np.uint8), mode="L")
    src.save(p)

    shrink_figure_png(p)

    with Image.open(p) as out:
        assert out.mode == "L"
        assert np.array_equal(np.asarray(out), np.asarray(src))


def test_rgba_is_skipped_rather_than_flattened(tmp_path):
    """Dropping channels on RGBA would silently discard alpha. No figure in
    the corpus is RGBA; the point is that the pass declines rather than
    guessing."""
    p = tmp_path / "fig_5.png"
    src = Image.new("RGBA", (30, 30), (10, 10, 10, 128))
    src.save(p)
    before = p.stat().st_size

    assert shrink_figure_png(p) == 0
    assert p.stat().st_size == before
    with Image.open(p) as out:
        assert out.mode == "RGBA"


# ── failure modes must not propagate ────────────────────────────────────

def test_corrupt_png_is_skipped_silently(tmp_path):
    p = tmp_path / "broken.png"
    p.write_bytes(b"\x89PNG\r\n\x1a\n" + b"garbage")
    assert shrink_figure_png(p) == 0
    assert p.exists()


def test_missing_file_is_skipped_silently(tmp_path):
    assert shrink_figure_png(tmp_path / "absent.png") == 0


def test_no_temp_files_left_behind(tmp_path):
    p = tmp_path / "fig_6.png"
    _grey_as_rgb(120, 90, seed=6).save(p)
    shrink_figure_png(p)
    assert list(tmp_path.glob("*.tmp")) == []


# ── directory sweep ─────────────────────────────────────────────────────

def test_sweep_reports_only_files_it_changed(tmp_path):
    d = tmp_path / "figures"
    d.mkdir()
    _grey_as_rgb(200, 200, seed=7).save(d / "a.png")
    _grey_as_rgb(180, 160, seed=8).save(d / "b.png")
    # Not a PNG, and not ours to touch.
    (d / "notes.txt").write_text("ignore me")

    changed, saved = shrink_figures_dir(d)

    assert changed == 2
    assert saved > 0
    assert (d / "notes.txt").read_text() == "ignore me"


def test_sweep_on_missing_dir_is_a_no_op(tmp_path):
    assert shrink_figures_dir(tmp_path / "nope") == (0, 0)


def test_sweep_is_idempotent(tmp_path):
    d = tmp_path / "figures"
    d.mkdir()
    _grey_as_rgb(200, 200, seed=9).save(d / "a.png")

    shrink_figures_dir(d)
    size_after_first = (d / "a.png").stat().st_size
    changed, saved = shrink_figures_dir(d)

    assert (changed, saved) == (0, 0)
    assert (d / "a.png").stat().st_size == size_after_first
