"""`dedupe_figures` — merging redundant crops without merging real figures (#205).

The function groups same-numbered figures and then decides, from bounding
boxes, whether they are redundant crops, coequal panels, or a whole figure
plus its subpanels. It had no direct tests, which is how the defect below
survived: **a bbox carries no page**, so the grouping key had to carry one.

Grouped on figure number alone, two figures at similar coordinates on
different pages read as redundant crops of each other and one was dropped.
That is routine in this material — a document that is its own translation
prints "Fig. 4" once in the original and again in the translation. On
`Carre1969_Nanomia_tr` it cost nine of twenty-two figures; keying on
(page, number) took that document from 13 extracted figures to 22, which is
exactly what docling detects and exactly what the gold transcription records.
"""
from __future__ import annotations

from pipeline.figures import (
    FIGURE_TYPE_FIGURE,
    FIGURE_TYPE_GRAPHICAL,
    FIGURE_TYPE_SUBPANEL,
    dedupe_figures,
)


def fig(idx, page, bbox, num="4", ftype=FIGURE_TYPE_FIGURE):
    return {"docling_idx": idx, "page": page, "bbox": list(bbox),
            "figure_number": num, "figure_type": ftype}


# --- the defect this file exists for -----------------------------------------


def test_same_number_on_different_pages_is_two_figures():
    """The translation case. Identical coordinates, different leaves."""
    items = [fig(0, 3, (0, 0, 100, 100)), fig(1, 24, (0, 0, 100, 100))]
    kept = dedupe_figures(items)
    assert len(kept) == 2
    assert {k["page"] for k in kept} == {3, 24}


def test_same_number_on_different_pages_keeps_both_even_when_boxes_nest():
    """An encompassing bbox on another page is not a whole-figure overview —
    it is a different figure that happens to be larger."""
    items = [fig(0, 3, (0, 0, 200, 200)), fig(1, 24, (10, 10, 100, 100))]
    kept = dedupe_figures(items)
    assert len(kept) == 2
    assert all(k["figure_type"] == FIGURE_TYPE_FIGURE for k in kept)
    assert not any("panel_letter" in k for k in kept)


# --- what dedupe was built for, which must keep working ----------------------


def test_overlapping_crops_on_one_page_still_merge():
    """The case the function exists for: docling emitting a whole-figure crop
    and a near-identical duplicate of it."""
    items = [fig(0, 5, (0, 0, 100, 100)), fig(1, 5, (2, 2, 98, 98))]
    kept = dedupe_figures(items)
    assert len(kept) == 1
    assert kept[0]["docling_idx"] == 0        # the larger crop survives


def test_coequal_panels_on_one_page_get_letters():
    """PLOS-style: Figure 3 is really four side-by-side panel images. They are
    all real figures and are lettered in reading order."""
    items = [fig(0, 5, (0, 0, 40, 40)), fig(1, 5, (50, 0, 90, 40)),
             fig(2, 5, (0, 50, 40, 90)), fig(3, 5, (50, 50, 90, 90))]
    kept = dedupe_figures(items)
    assert len(kept) == 4
    assert all(k["figure_type"] == FIGURE_TYPE_FIGURE for k in kept)
    assert sorted(k["panel_letter"] for k in kept) == ["a", "b", "c", "d"]


def test_whole_figure_plus_subpanels_is_now_reachable():
    """The branch that could never fire (#207).

    Both stages used to share `_bbox_overlap_fraction`, which divides by the
    *smaller* box and is therefore symmetric: a fully contained panel scored
    1.0, tripping stage 1's 0.5 redundancy threshold and being discarded
    before stage 2's 0.8 containment test could classify it. Anything stage 2
    would have accepted, stage 1 had already thrown away — so
    `FIGURE_TYPE_SUBPANEL` was never assigned, and no figure among the 420 in
    the reference corpuscle carried it.

    The two stages ask different questions and now use different measures.
    "Are these the same box?" is intersection over *union*, which punishes a
    size difference. "Does this box contain that one?" keeps the original
    formula, which is right for containment.
    """
    items = [fig(0, 5, (0, 0, 100, 100)),
             fig(1, 5, (5, 5, 45, 45)), fig(2, 5, (55, 55, 95, 95))]
    kept = dedupe_figures(items)
    primary = [k for k in kept if k["figure_type"] == FIGURE_TYPE_FIGURE]
    subs = [k for k in kept if k["figure_type"] == FIGURE_TYPE_SUBPANEL]
    assert len(primary) == 1 and primary[0]["docling_idx"] == 0
    assert len(subs) == 2
    assert sorted(s["panel_letter"] for s in subs) == ["a", "b"]
    assert all(s["primary_figure_docling_idx"] == 0 for s in subs)


def test_the_two_measures_answer_different_questions():
    """The distinction the fix rests on, stated directly."""
    from pipeline.figures import _bbox_iou, _bbox_overlap_fraction
    big, nested = [0, 0, 100, 100], [5, 5, 45, 45]
    # A nested panel is *contained* but is not the same box.
    assert _bbox_overlap_fraction(big, nested) == 1.0
    assert _bbox_iou(big, nested) < 0.5
    # A near-identical crop is both.
    dup = [2, 2, 98, 98]
    assert _bbox_iou(big, dup) > 0.5


# --- what must never be grouped ----------------------------------------------


def test_graphical_elements_never_participate():
    """A journal logo that latched onto a figure number via the caption
    proximity heuristic must not merge with, or relabel, a real figure."""
    items = [fig(0, 5, (0, 0, 100, 100)),
             fig(1, 5, (0, 0, 100, 100), ftype=FIGURE_TYPE_GRAPHICAL)]
    kept = dedupe_figures(items)
    assert len(kept) == 2
    assert {k["figure_type"] for k in kept} == {FIGURE_TYPE_FIGURE,
                                               FIGURE_TYPE_GRAPHICAL}


def test_unnumbered_figures_pass_through_untouched():
    items = [dict(fig(0, 5, (0, 0, 100, 100)), figure_number=None),
             dict(fig(1, 5, (0, 0, 100, 100)), figure_number=None)]
    kept = dedupe_figures(items)
    assert len(kept) == 2


def test_different_numbers_on_one_page_are_not_merged():
    items = [fig(0, 5, (0, 0, 100, 100), num="4"),
             fig(1, 5, (2, 2, 98, 98), num="5")]
    assert len(dedupe_figures(items)) == 2


def test_a_missing_page_does_not_collapse_everything():
    """Items whose page could not be resolved share a `None` key. They must
    still be compared on their bboxes rather than silently vanishing."""
    a = dict(fig(0, None, (0, 0, 100, 100)))
    b = dict(fig(1, None, (500, 500, 600, 600)))
    kept = dedupe_figures([a, b])
    assert len(kept) == 2
