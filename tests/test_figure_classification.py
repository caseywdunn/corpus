"""`classify_figure` and the furniture-recurrence guard (#204).

`figure_type` is not cosmetic: `_REAL_FIGURE_TYPES` in
`mcpsrv/tools/figures.py` excludes `graphical_element` from **every** tool
that returns figures. Misclassifying a real figure does not mislabel it, it
makes it unreachable.

The size floor was condemning small real figures on their dimensions alone —
Vanhoeffen 1906 Fig. 11 is an engraved nectophore 49 pts wide against a 50 pt
floor, captioned and numbered. A caption carrying a parseable figure number is
evidence the page provides, and it now overrides the floor.

That evidence is not unconditionally safe, because caption proximity sometimes
attaches a real caption to a logo. The guard is recurrence: running furniture
sits at the same place on page after page, and no real figure in the
35-document reference corpus repeats a position even twice, while one paper's
masthead repeats on 24 of 25 pages.
"""
from __future__ import annotations

from pipeline.figures import (
    FIGURE_TYPE_FIGURE,
    FIGURE_TYPE_GRAPHICAL,
    FIGURE_TYPE_UNCLASSIFIED,
    _FURNITURE_RECURRENCE_PAGES,
    classify_figure,
    furniture_positions,
)


def item(page, bbox, caption="", num=None):
    return {"page": page, "bbox": list(bbox), "caption_text": caption,
            "figure_number": num}


# --- the defect ---------------------------------------------------------------


def test_a_small_captioned_numbered_figure_is_not_furniture():
    """Vanhoeffen 1906 Fig. 11: 49 x 96 pts, one point under the floor."""
    e = item(7, (243.3, 226.1, 292.0, 321.8), "Fig. 11.", "11")
    assert classify_figure(e) == FIGURE_TYPE_FIGURE


def test_a_small_uncaptioned_item_is_still_furniture():
    """The floor still does its job where there is no evidence against it."""
    assert classify_figure(item(1, (0, 0, 40, 40))) == FIGURE_TYPE_GRAPHICAL


def test_a_caption_without_a_number_does_not_rescue_a_small_item():
    """A caption alone is weak — proximity attaches them to all sorts of
    things. It takes a parseable figure number to overrule the floor."""
    assert classify_figure(item(1, (0, 0, 40, 40), "some nearby text")) \
        == FIGURE_TYPE_GRAPHICAL


# --- the recurrence guard -----------------------------------------------------


def test_a_masthead_repeating_across_pages_is_furniture():
    """Even with a caption and number wrongly attached to it.

    This is the real case: one paper's PLOS One logo sits at the same bbox on
    24 of its 25 pages and caption proximity gave it "Fig 3. Syntenic
    transitions ...". Without this guard the size relaxation above promotes it
    ten times over.
    """
    logo = (34.4, 734.5, 118.1, 775.1)
    items = [item(p, logo, "Fig 3.  Syntenic transitions", "3") for p in range(1, 26)]
    rec = furniture_positions(items)
    assert classify_figure(items[0], recurring=rec) == FIGURE_TYPE_GRAPHICAL


def test_a_figure_appearing_once_is_not_caught_by_the_guard():
    items = [item(7, (243.3, 226.1, 292.0, 321.8), "Fig. 11.", "11")]
    rec = furniture_positions(items)
    assert classify_figure(items[0], recurring=rec) == FIGURE_TYPE_FIGURE


def test_the_guard_needs_several_pages_not_two():
    """Two figures happening to share a position is coincidence, not a running
    header. The threshold sits in the wide gap between 1 and 24."""
    assert _FURNITURE_RECURRENCE_PAGES > 2
    logo = (0, 700, 80, 740)
    few = [item(p, logo) for p in range(1, 3)]
    assert furniture_positions(few) == frozenset()


def test_positions_are_matched_with_tolerance_not_exactly():
    """Sub-point jitter between pages must not defeat the match."""
    base = (34.0, 734.0, 118.0, 775.0)
    jittered = [item(p, (base[0] + p * 0.3, base[1], base[2], base[3]))
                for p in range(1, 10)]
    assert len(furniture_positions(jittered)) == 1


def test_items_without_a_bbox_are_ignored_by_the_guard():
    assert furniture_positions([{"page": 1, "bbox": None}] * 10) == frozenset()


# --- unchanged behaviour ------------------------------------------------------


def test_a_large_uncaptioned_item_is_still_unclassified():
    assert classify_figure(item(1, (0, 0, 400, 400))) == FIGURE_TYPE_UNCLASSIFIED


def test_a_large_captioned_numbered_item_is_still_a_figure():
    assert classify_figure(item(1, (0, 0, 400, 400), "Fig. 2. A caption.", "2")) \
        == FIGURE_TYPE_FIGURE
