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
What this fixes is retrieval and counting: figure recall across the reference
corpus goes 0.849 -> 0.910, and on the served surface 0.833 -> 0.894.

The first version of this claimed the recall with no cost, which the gold set
did not bear out — it also read the cross-references in running prose as legend
entries, putting 34 wrong figures into one monograph. What the line has to
open with is the point, and the tests at the foot of this file are about that.
"""
from __future__ import annotations

import pytest

from pipeline.figures import (
    FIGURE_TYPE_FIGURE,
    _LEGEND_OPENER,
    _MIN_PLATE_LEGEND_ENTRIES,
    expand_plate_figures,
    caption_figure_entries,
    plate_legend_entries,
    reconcile_plate_legend_numbers,
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


def test_one_enumerating_caption_yields_separate_figure_entries():
    entries = caption_figure_entries(
        "Fig. 10. Colony. Fig. 11. Eudoxid. Fig. 12. Gonophore.",
    )
    assert [e["figure_number"] for e in entries] == ["10", "11", "12"]
    assert entries[1]["caption_text"] == "Fig. 11. Eudoxid."


def test_number_list_shares_caption_without_becoming_panels():
    entries = caption_figure_entries(
        "Figur 8 und 9: Muggiaea, A-C lateral sections.",
    )
    assert [e["figure_number"] for e in entries] == ["8", "9"]
    assert entries[0]["caption_text"] == entries[1]["caption_text"]


def test_comma_list_with_conjunction_yields_every_number():
    entries = caption_figure_entries("Figs. 1, 2, and 3. Developmental stages.")
    assert [entry["figure_number"] for entry in entries] == ["1", "2", "3"]


def test_reference_inside_caption_prose_is_not_an_enumerated_entry():
    entries = caption_figure_entries(
        "Fig. 1. Young stage (see Fig. 9). Fig. 2. Mature stage.",
    )
    assert [entry["figure_number"] for entry in entries] == ["1", "2"]


def test_numeric_range_expands_but_panel_range_does_not():
    entries = caption_figure_entries("Fig. 58-63. Velella, panels A-C.")
    assert [e["figure_number"] for e in entries] == [
        "58", "59", "60", "61", "62", "63",
    ]


def test_single_grouped_text_item_can_form_a_plate_legend():
    entries = plate_legend_entries([{
        "text": "Fig. 10. Colony. Fig. 11. Eudoxid. Fig. 12. Gonophore.",
        "bbox": [0, 0, 200, 20],
    }])
    assert [e["figure_number"] for e in entries] == ["10", "11", "12"]


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
    assert all("_full_caption_text" not in s["caption_candidates"][0]
               for s in siblings)


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


def test_equal_counts_reconcile_a_lower_confidence_duplicate():
    entries = plate_legend_entries([{
        "text": "Fig. 10. Colony. Fig. 11. Eudoxid. Fig. 12. Gonophore.",
        "bbox": [100, 120, 300, 160],
    }])
    items = [
        {**plate(idx=1, bbox=(100, 300, 500, 700), num="10", cap="Fig. 10."),
         "caption_source": "docling_caption_link", "caption_confidence": "high"},
        {**plate(idx=2, bbox=(130, 180, 200, 290), num="11", cap="Fig. 11."),
         "caption_source": "docling_caption_link", "caption_confidence": "high"},
        {**plate(idx=3, bbox=(240, 180, 300, 280), num="11", cap="Fig. 11."),
         "caption_source": "heuristic_proximity", "caption_confidence": "medium"},
    ]

    out = expand_plate_figures(items, {17: entries})

    assert len(out) == 3
    assert sorted(str(item["figure_number"]) for item in out) == ["10", "11", "12"]
    repaired = next(item for item in out if item["figure_number"] == "12")
    assert repaired["caption_text"] == "Fig. 12. Gonophore."
    assert repaired["caption_source"] == "plate_legend_reconciled"
    assert repaired["caption_confidence"] == "medium"
    assert repaired["figure_number_raw"] == "11"
    assert repaired["figure_number_source"] == "plate_legend_reconciled"
    assert repaired["caption_text_raw"] == "Fig. 11."


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


def test_collected_legend_enriches_previous_page_instead_of_wrong_plate():
    prior = {
        **plate(idx=1, page=14, num="26", cap="Figur 26."),
        "caption_candidates": [],
    }
    current = plate(idx=2, page=15, num="29", cap="Figur 29.")
    entries = plate_legend_entries([
        {"text": "Figur 26. Hippopodius after Huxley.", "bbox": [0, 80, 200, 90]},
        {"text": "Figur 28. Young stage.", "bbox": [0, 60, 200, 70]},
        {"text": "Figur 29 und 30. Older stages.", "bbox": [0, 40, 200, 50]},
    ])

    out = expand_plate_figures([prior, current], {15: entries})

    assert prior["caption_text"] == "Figur 26. Hippopodius after Huxley."
    assert prior["caption_source"] == "facing_page_plate_legend"
    assert prior["caption_page_distance"] == 1
    assert prior["caption_status"] == "bound"
    page15 = [item for item in out if item["page"] == 15]
    assert sorted(item["figure_number"] for item in page15) == ["28", "29", "30"]
    assert not any(item["figure_number"] == "26" for item in page15)


def test_caption_only_page_can_enrich_previous_page_figure():
    prior = plate(idx=1, page=14, num="26", cap="Figur 26.")
    entries = plate_legend_entries([
        {"text": "Figur 26. Hippopodius after Huxley.", "bbox": [0, 80, 200, 90]},
        {"text": "Figur 27. Longitudinal section.", "bbox": [0, 60, 200, 70]},
    ])

    out = expand_plate_figures([prior], {15: entries})

    assert out == [prior]
    assert prior["caption_text"] == "Figur 26. Hippopodius after Huxley."
    assert prior["caption_page"] == 15


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


# ---------------------------------------------------------------------------
# A legend line opens with its figure's label; a cross-reference doesn't.
#
# The first version of #203 scanned every text item on a page for a figure
# number anywhere in it. On a plate page that is right. On a page of running
# prose it reads the monograph's own cross-references as legend entries, and
# the expansion then serves *some other figure's image* under the referenced
# number — 34 of them in Totton 1965, which is worse than not finding the
# figure at all. Measured against the gold set, anchoring the label recovers
# the precision that cost (0.892 -> 0.962 on the served surface) while
# keeping essentially all of the recall it bought (0.901 -> 0.894).
# ---------------------------------------------------------------------------

# Lines lifted verbatim from Totton 1965 pages the unanchored rule expanded.
TOTTON_CROSS_REFERENCES = [
    "Plate XX, figures 1, 2",                    # a species heading
    "Plate XXXIX",
    "Plate XIII, figures 1-3",
    "figured by Bigelow (1911b, PL. 21, as Anthophysa rosea) and by Leloup "
    "(1941, Pl. 11). Bigelow shows the giant cells in the septa",
]

# Lines from plates the rule is meant to expand — Vanhoeffen 1906 sets three
# separately-numbered engravings under one legend block.
REAL_LEGEND_LINES = [
    "Figur 55. Physalia arethusa Browne nach Agassiz.",
    "Figur 56. Physalialarve nach Huxley.",
    "Figur 57. Ältere Physalialarve nach Chun.",
]


@pytest.mark.parametrize("line", TOTTON_CROSS_REFERENCES)
def test_a_cross_reference_is_not_a_legend_entry(line):
    assert plate_legend_entries([{"text": line, "bbox": [0, 0, 100, 10]},
                                 {"text": line, "bbox": [0, 20, 100, 30]}]) == []


def test_a_real_legend_block_still_expands():
    page = [{"text": t, "bbox": [0, i * 20, 200, i * 20 + 10]}
            for i, t in enumerate(REAL_LEGEND_LINES)]
    entries = plate_legend_entries(page)
    assert [e["figure_number"] for e in entries] == ["55", "56", "57"]


def test_prose_mentioning_figures_does_not_become_a_legend():
    """The failing page: ordinary text, one text-figure, two cross-references."""
    page = [
        {"text": "Forskalia edwardsi Kolliker, 1853", "bbox": [0, 0, 200, 10]},
        {"text": "Plate XX, figures 1, 2", "bbox": [0, 20, 200, 30]},
        {"text": "Nectophores (text-figs. 52, 53): The nectosome consists of "
                 "a conical or more cylindrical structure", "bbox": [0, 40, 200, 60]},
        {"text": "FIG. 53. Forskalia edwardsi Kolliker", "bbox": [0, 80, 200, 90]},
    ]
    entries = plate_legend_entries(page)
    # Only the actual caption qualifies, which is one entry — below the
    # threshold, so the page is left alone entirely.
    assert entries == []


def test_the_label_must_be_followed_by_the_number():
    """`figured by Bigelow` opens with "fig" and is not a figure label."""
    assert not _LEGEND_OPENER.match("figured by Bigelow (1911b, PL. 21)")
    assert _LEGEND_OPENER.match("Fig. 33.")
    assert _LEGEND_OPENER.match("FiG. 128. Eudoxoides spiralis")
    assert _LEGEND_OPENER.match("Abb. 7. Physalia")
    assert _LEGEND_OPENER.match("Text-figure 24")
    assert _LEGEND_OPENER.match("Figuren 55 und 56.")


def test_plate_number_ocr_loss_requires_sequence_and_missing_evidence():
    host = {
        "figure_id": "docling_18", "figure_number": "36",
        "caption_text": "Fig. 36. Taster.",
    }
    siblings = [
        {
            "figure_id": f"docling_18_fig{number}",
            "figure_number": number,
            "caption_text": f"Fig. {number}. Description.",
            "caption_source": "plate_legend",
            "shares_image_with": 18,
            "page": 17,
            "caption_page": 17,
        }
        for number in ("31", "32", "33", "3", "35")
    ]

    repaired = reconcile_plate_legend_numbers(
        [host, *siblings], [{"figure_number": "34"}],
    )

    assert repaired == 1
    corrected = next(f for f in siblings if f["figure_number"] == "34")
    assert corrected["figure_number_raw"] == "3"
    assert corrected["caption_text"].startswith("Fig. 34.")
    assert corrected["caption_text_raw"].startswith("Fig. 3.")
    assert corrected["caption_source"] == "plate_legend_sequence_reconciled"


def test_plate_number_sequence_is_not_rewritten_without_independent_missing():
    host = {"figure_id": "docling_18", "figure_number": "36"}
    siblings = [
        {
            "figure_id": f"docling_18_fig{number}",
            "figure_number": number,
            "caption_text": f"Fig. {number}.",
            "caption_source": "plate_legend",
            "shares_image_with": 18,
        }
        for number in ("31", "32", "33", "3", "35")
    ]

    assert reconcile_plate_legend_numbers([host, *siblings], []) == 0
    assert any(f["figure_number"] == "3" for f in siblings)
