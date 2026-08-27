"""Plate legends: one record per figure the legend names (#203).

A historical plate carries several separately-numbered engravings under one
legend, printed as a run of lines:

    Fig. 31.   Agalmopsis elegans Sars.
    Fig. 32.   Schwimmglocke.
    Fig. 33.   Deckstück.
    ...

docling emits those as separate text items and extracts the plate itself as a
*single* picture. Caption association then binds whichever labelled line is
vertically nearest — on the page this was written for, a bare "Fig. 36." at the
foot of the page — and the other five figures existed nowhere in the output.
Asking for Figure 33 returned nothing.

The records share the plate's image, because locating an individual engraving
needs OCR of the lettering printed on the plate and that is a separate problem.
What this fixes is retrieval and counting, and it is worth 0.849 -> 0.917
figure recall across the reference corpus, with no document regressing.
"""
from __future__ import annotations

from pipeline.figures import (
    FIGURE_TYPE_FIGURE,
    _MIN_PLATE_LEGEND_ENTRIES,
    expand_plate_figures,
    plate_legend_entries,
)

LEGEND = [
    {"text": "Fig. 31. Agalmopsis elegans Sars.", "bbox": [246, 518, 400, 528]},
    {"text": "Fig. 32. Schwimmglocke.", "bbox": [246, 507, 400, 517]},
    {"text": "Fig. 33. Deckstück.", "bbox": [246, 495, 400, 505]},
]


def plate(idx=1, page=17, bbox=(0, 0, 500, 700), num=None, cap=""):
    return {"docling_idx": idx, "page": page, "bbox": list(bbox),
            "figure_number": num, "caption_text": cap, "image": object()}


# --- recognising a legend -----------------------------------------------------


def test_a_run_of_numbered_lines_is_a_legend():
    entries = plate_legend_entries(LEGEND)
    assert [e["figure_number"] for e in entries] == ["31", "32", "33"]
    assert entries[0]["caption_text"].startswith("Fig. 31.")


def test_one_caption_is_not_a_legend():
    """The ordinary page: one figure and its caption. Expanding that would
    invent figures."""
    assert plate_legend_entries([LEGEND[0]]) == []
    assert _MIN_PLATE_LEGEND_ENTRIES >= 2


def test_the_same_number_twice_counts_once():
    """A figure whose caption is repeated — a running head, a continuation —
    must not inflate the count."""
    dup = [LEGEND[0], dict(LEGEND[0]), LEGEND[1]]
    assert [e["figure_number"] for e in plate_legend_entries(dup)] == ["31", "32"]


def test_unnumbered_prose_is_ignored():
    texts = [{"text": "Siphonophoren.", "bbox": None},
             {"text": "n. Fewkes.", "bbox": None}] + LEGEND
    assert len(plate_legend_entries(texts)) == 3


# --- expanding ----------------------------------------------------------------


def test_a_plate_gains_one_record_per_named_figure():
    items = [plate(num="36", cap="Fig.  36.")]
    out = expand_plate_figures(items, {17: plate_legend_entries(LEGEND)})
    assert len(out) == 4                       # the original plus 31, 32, 33
    nums = sorted(str(o["figure_number"]) for o in out)
    assert nums == ["31", "32", "33", "36"]


def test_the_new_records_share_the_plate_image_and_say_so():
    items = [plate(idx=7, num="36")]
    out = expand_plate_figures(items, {17: plate_legend_entries(LEGEND)})
    siblings = [o for o in out if o.get("shares_image_with") is not None]
    assert len(siblings) == 3
    assert all(s["shares_image_with"] == 7 for s in siblings)
    assert all(s["bbox"] == items[0]["bbox"] for s in siblings)
    assert all(s["figure_type"] == FIGURE_TYPE_FIGURE for s in siblings)
    assert all(s["caption_source"] == "plate_legend" for s in siblings)


def test_a_figure_that_already_has_its_own_record_is_not_duplicated():
    items = [plate(idx=1, num="31", cap="Fig. 31. Agalmopsis elegans Sars.")]
    out = expand_plate_figures(items, {17: plate_legend_entries(LEGEND)})
    assert sorted(str(o["figure_number"]) for o in out) == ["31", "32", "33"]


def test_a_page_docling_already_separated_is_left_alone():
    """Three legend lines and three pictures: nothing to add. Expanding here
    would double the page."""
    items = [plate(idx=i, num=str(30 + i)) for i in (1, 2, 3)]
    out = expand_plate_figures(items, {17: plate_legend_entries(LEGEND)})
    assert len(out) == 3


def test_a_page_with_no_pictures_gains_nothing():
    """A legend with no plate to attach to is not evidence of a figure the
    extractor can serve."""
    assert expand_plate_figures([], {17: plate_legend_entries(LEGEND)}) == []


def test_the_largest_picture_on_the_page_is_taken_as_the_plate():
    small = plate(idx=2, bbox=(0, 0, 40, 40))
    big = plate(idx=3, bbox=(0, 0, 500, 700))
    out = expand_plate_figures([small, big], {17: plate_legend_entries(LEGEND)})
    siblings = [o for o in out if o.get("shares_image_with") is not None]
    assert siblings and all(s["shares_image_with"] == 3 for s in siblings)


def test_other_pages_are_untouched():
    items = [plate(idx=1, page=17), plate(idx=2, page=18, num="99")]
    out = expand_plate_figures(items, {17: plate_legend_entries(LEGEND)})
    assert [o for o in out if o["page"] == 18] == [items[1]]


# --- the wiring ---------------------------------------------------------------


def test_siblings_reuse_the_plate_file_instead_of_copying_it():
    """Six records pointing at one 987 KB plate must not be six copies of it.
    Before this, one document's figures directory was 30 MB instead of 13 MB.
    """
    import inspect
    from pipeline import extract
    src = inspect.getsource(extract.extract_docling_content)
    assert "saved_filenames[shares]" in src
    assert "shares in saved_filenames" in src
