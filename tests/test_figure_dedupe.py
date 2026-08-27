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


def test_contained_panels_are_dropped_not_labelled_as_subpanels():
    """Pins *actual* behaviour, which is not what the docstring describes.

    `dedupe_figures` documents a "whole figure plus subpanels" mode that
    relabels contained items `subpanel`. That branch is unreachable:
    `_bbox_overlap_fraction` divides by the *smaller* box, so a fully
    contained panel scores 1.0 against its parent — over step 1's 0.5
    overlap-duplicate threshold, which drops it before step 2's 0.8
    "encompassed" test is ever consulted. Since the measure is symmetric,
    anything step 2 would accept step 1 has already removed.

    Consistent with the reference corpuscle, where none of 420 figures
    carries `figure_type: subpanel`. Tracked separately; this test exists so
    that whoever fixes it sees a failure here rather than a silent change.
    """
    items = [fig(0, 5, (0, 0, 100, 100)),
             fig(1, 5, (5, 5, 45, 45)), fig(2, 5, (55, 55, 95, 95))]
    kept = dedupe_figures(items)
    assert len(kept) == 1
    assert kept[0]["docling_idx"] == 0
    assert kept[0]["figure_type"] == FIGURE_TYPE_FIGURE
    assert not any(k["figure_type"] == FIGURE_TYPE_SUBPANEL for k in kept)


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
