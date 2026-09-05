"""Caption association scoring (#195).

Detection asks whether the figure objects are right. This asks whether the
caption is attached to the figure it belongs to — what `dev_docs/OVERVIEW.md`
calls "the highest-value annotation per figure and the hardest in historical
layouts".

Two design decisions are load-bearing and both were arrived at by getting them
wrong first:

* **Bind on the figure number, not the caption text.** A caption-similarity
  matcher reports 44% and is mostly artifact: one paper prints every caption
  twice (Chinese and English) and scored 0 of 10 while every figure was in
  fact bound correctly; plate pages carry `FIG. 1` and nothing else, which no
  threshold can match; and a document that is its own translation prints
  `Fig. 1` twice, legitimately.

* **Count only numbers printed in a figure block (plus its adjacent plate
  heading).** Gold pages are full of figure numbers that are references —
  "see Fig. 18", "figured by Bigelow (op. cit., fig. 34)". Counting those
  measures objects that are not on that page. Restricted to the structural
  block and its plate identity, the measure means what it says.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from pipeline.figures import parse_panels_from_caption

_PATH = Path(__file__).parent.parent / "tools" / "qc" / "caption_binding.py"
_spec = importlib.util.spec_from_file_location("qc_caption_binding", _PATH)
cb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cb)

FIXTURE = Path(__file__).parent / "fixtures" / "gold_fidelity"
PANEL_FIXTURE = Path(__file__).parent / "fixtures" / "caption_panels"


# --- what a gold block is offering -------------------------------------------


def test_a_prose_caption_is_recognised():
    kind, nums, first = cb.classify_block(
        "Fig. 99. A, Lensia conoidea (Kefferstein & Ehlers) lateral view of "
        "the anterior nectophore from Thor St. 36.")
    assert kind == "prose_caption"
    assert nums == {"99"}


def test_a_bare_label_is_not_a_prose_caption():
    """Scoring a bare label for caption text reports a failure where the page
    prints no prose to bind."""
    kind, nums, _ = cb.classify_block("FIG. 1")
    assert kind == "bare_label"
    assert nums == {"1"}


def test_plate_lettering_with_no_label_is_its_own_class():
    """An engraved plate whose only text is numerals beside the drawings."""
    kind, nums, _ = cb.classify_block("1 7 Stéphanomie triangulaire. 8 12")
    assert kind == "lettering_only"
    assert nums == set()


def test_standalone_numbers_are_labels_only_inside_a_plate_block():
    body = "PLATE I\nF. 1.\n3\nF.2\n3a"

    kind, nums, _ = cb.classify_block(body, structural_kind="PLATE")

    assert kind == "bare_label"
    assert nums == {"1", "2", "3", "3a"}
    # The same shape inside an ordinary figure can be chart ticks or scale
    # values and must not manufacture figure identities.
    assert cb.classify_block("1\n3\n2", structural_kind="FIGURE")[1] == set()


def test_a_block_with_no_printed_text_is_counted_not_scored():
    kind, nums, _ = cb.classify_block("[NOTE: engraved ornament, no lettering]")
    assert kind == "nothing_printed"
    assert nums == set()


def test_a_legend_naming_several_figures_yields_all_of_them():
    """The plate case from #203 — six engravings under one legend."""
    _kind, nums, _ = cb.classify_block(
        "Fig. 31. Agalmopsis elegans Sars.\n"
        "Fig. 32. Schwimmglocke.\n"
        "Fig. 33. Deckstück.")
    assert nums == {"31", "32", "33"}


def test_one_shared_caption_label_yields_every_number():
    _kind, nums, _ = cb.classify_block(
        "Figur 8 und 9: Muggiaea atlantica Cunningham nach Cunningham.",
    )
    assert nums == {"8", "9"}


def test_several_captions_on_one_line_are_all_counted():
    _kind, nums, _ = cb.classify_block(
        "Fig. 10. Colony. Fig. 11. Eudoxid. Fig. 12. Gonophore.",
    )
    assert nums == {"10", "11", "12"}


def test_gold_parser_handles_lists_ranges_and_roman_labels():
    assert cb.gold_caption_entries("Figs. 1, 2, and 3. Stages.") == ["1", "2", "3"]
    assert cb.gold_caption_entries("Figg. 2-5. Sezioni.") == ["2", "3", "4", "5"]
    assert cb.gold_caption_entries("Fig. 58-60. Stages.") == ["58", "59", "60"]
    assert cb.gold_caption_entries("Plate IV. Colony.") == ["4"]


def test_gold_parser_does_not_read_mouth_plate_key_as_plate_1000():
    assert cb.gold_caption_entries("Pl.M.") == []
    assert cb.gold_caption_entries("Pl.M.\tSom.") == []
    assert cb.gold_caption_entries(
        "Fig. 4: R.r. R.l. Pl.m. Ost. To.d. Vel."
    ) == ["4"]


def test_gold_panel_parser_starts_after_engraved_mouth_plate_key():
    body = "Pl.M.\nA.\nB.\nFIG. 1. Whole colony without separate panels."
    assert cb.gold_panel_labels(body) == []


def test_gold_parser_rejects_caption_internal_source_citations():
    assert cb.gold_caption_entries(
        "Fig. 77. Nectopyramis, from Totton, 1954, fig. 36)",
    ) == ["77"]


def test_gold_parser_does_not_accept_production_only_ocr_damage():
    """A hand-transcription oracle must not inherit extractor heuristics."""
    assert cb.gold_caption_entries("F16. 1. OCR-damaged extraction") == []


def test_transcriber_commentary_does_not_contribute_numbers():
    """A note mentioning Fig. 4 is not a figure on the page."""
    _kind, nums, _ = cb.classify_block(
        "[NOTE: this plate is reproduced from Fig. 4 of Bigelow 1911]\n"
        "Fig. 12. Nectophore.")
    assert nums == {"12"}


# --- caption-enumerated panel splits -----------------------------------------


@pytest.mark.parametrize(("caption", "labels"), [
    ("Fig. 1. A. colony. B. nectophore. C. bract.", ["A", "B", "C"]),
    ("Fig. 2. (A) upper view; (B) lower view.", ["A", "B"]),
    ("Fig. 3. A-C. successive developmental stages.", ["A", "B", "C"]),
    ("Fig. 4. A, anterior; B-D, lateral views.", ["A", "B", "C", "D"]),
])
def test_gold_panel_parser_covers_declared_styles(caption, labels):
    """The gold parser is intentionally independent of production's parser."""
    assert cb.gold_panel_labels(caption) == labels


def test_gold_panel_parser_ignores_image_lettering_before_caption():
    body = "A\nB\nC\nFig. 5. Whole colony without caption-enumerated panels."
    assert cb.gold_panel_labels(body) == []


def test_gold_panel_parser_rejects_initials_but_preserves_printed_gaps():
    assert cb.gold_panel_labels(
        "Fig. 6. Whole colony. Photo credit to C. Munro."
    ) == []
    assert cb.gold_panel_labels("Fig. 7. A. upper view; C. lower view.") == ["A", "C"]


def test_gold_panel_parser_preserves_a_deliberately_gapped_set():
    caption = (
        "Fig. 74. Bracts. A, first; B, second; C, third; D, fourth; "
        "E, fifth; F, sixth; G, seventh; H, eighth; K, ninth; L, tenth."
    )
    assert cb.gold_panel_labels(caption) == [
        "A", "B", "C", "D", "E", "F", "G", "H", "K", "L",
    ]


def test_gold_panel_parser_ignores_credit_and_taxonomic_initials():
    assert cb.gold_panel_labels(
        "Fig. 83. A, male; B, female. From sketches by the late Miss M. J. Delap."
    ) == ["A", "B"]
    assert cb.gold_panel_labels(
        "Fig. 132. A, first. Chuniphyes multidentata L. & v. R. B, second."
    ) == ["A", "B"]


def test_panel_case_requires_same_page_and_figure_number():
    figures = {
        3: [{
            "figure_number": "8",
            "panels_from_caption": [{"label": "A"}, {"label": "B"}],
        }],
        4: [{
            "figure_number": "9",
            "panels_from_caption": [{"label": "A"}, {"label": "B"}],
        }],
    }
    exact = cb.score_panel_case(3, "8", ["A", "B"], figures)
    missing = cb.score_panel_case(3, "9", ["A", "B"], figures)
    assert exact["exact"] is True
    assert exact["number_bound"] is True
    assert missing["reported"] is False
    assert missing["number_bound"] is False


def test_panel_case_reports_missing_and_surplus_labels():
    case = cb.score_panel_case(
        5,
        "10",
        ["A", "B", "C"],
        {5: [{
            "figure_number": "10",
            "panels_from_caption": [
                {"label": "A"}, {"label": "B"}, {"label": "D"},
            ],
        }]},
    )
    assert case["matched_labels"] == ["A", "B"]
    assert case["missing_labels"] == ["C"]
    assert case["surplus_labels"] == ["D"]
    assert case["exact"] is False


def test_panel_report_scores_positive_and_negative_fixture_cases():
    report = cb.build_report(
        PANEL_FIXTURE / "gold", PANEL_FIXTURE / "corpuscle",
    )
    totals = report["surfaces"]["default MCP types"]["totals"]
    assert totals["reported_pair_capacity"] == 2
    assert totals["reported_pair_capacity_rate"] == 1.0
    assert totals["gold_panel_figures"] == 1
    assert totals["reported_panel_figures"] == 1
    assert totals["exact_panel_figures"] == 1
    assert totals["panel_label_recall"] == 1.0
    assert totals["panel_label_precision"] == 1.0
    assert totals["number_matched_unpanelled_figures"] == 1
    assert totals["unexpected_panel_figures"] == 1
    assert totals["panel_figure_precision"] == 0.5
    assert report["segments"]["era"]["2000-"]["panel_exact_rate"] == 1.0


def test_production_panel_parser_handles_comma_style():
    panels = parse_panels_from_caption(
        "FIG. 99. A, Lensia conoidea, lateral view; "
        "B, Lensia multicristata, dorsal view."
    )
    assert [panel["label"] for panel in panels] == ["A", "B"]
    assert [panel["kind"] for panel in panels] == ["comma", "comma"]


def test_production_panel_parser_handles_comma_lists_and_ranges():
    panels = parse_panels_from_caption(
        "FIG. 104. A, anterior nectophore; B-D, other views."
    )
    assert [panel["label"] for panel in panels] == ["A", "B", "C", "D"]


def test_production_panel_parser_rejects_noncontiguous_comma_initials():
    assert parse_panels_from_caption(
        "FIG. 1. Whole animal, after A, Smith and C, Jones."
    ) == []


def test_production_panel_parser_ignores_letters_beyond_supported_panel_range():
    panels = parse_panels_from_caption(
        "Fig. 4. A. first. B. second. C. third. D. fourth. E. fifth. "
        "Map made in R."
    )
    assert [panel["label"] for panel in panels] == ["A", "B", "C", "D", "E"]


def test_production_panel_parser_keeps_a_strong_gapped_set_through_l():
    panels = parse_panels_from_caption(
        "Fig. 74. A, one; B, two; C, three; D, four; E, five; F, six; "
        "G, seven; H, eight; K, nine; L, ten."
    )
    assert [panel["label"] for panel in panels] == [
        "A", "B", "C", "D", "E", "F", "G", "H", "K", "L",
    ]


def test_production_panel_parser_stops_at_abbreviation_glossary():
    panels = parse_panels_from_caption(
        "Fig. 28. A, upper side; B, lateral view. "
        "C.ped = pedicular canal; C. rad.lat = lateral radial canal."
    )
    assert [panel["label"] for panel in panels] == ["A", "B"]


def test_one_scientific_equality_does_not_hide_later_panels():
    panels = parse_panels_from_caption(
        "Fig. 1. A. Sample size n = 16. B. Replicate measurement."
    )
    assert [panel["label"] for panel in panels] == ["A", "B"]


def test_repeated_scientific_equalities_do_not_look_like_a_glossary():
    panels = parse_panels_from_caption(
        "Fig. 1. A. First sample, n = 16 and p = 0.04. "
        "B. Replicate measurement."
    )
    assert [panel["label"] for panel in panels] == ["A", "B"]


def test_production_panel_parser_ignores_taxonomic_author_initial():
    panels = parse_panels_from_caption(
        "Fig. 65. Rosacea plicata Q. & G. Polygastric phase. "
        "A, lateral; B, dorsal; C, ventral."
    )
    assert [panel["label"] for panel in panels] == ["A", "B", "C"]


def test_production_panel_parser_ignores_parenthesized_author_initial():
    panels = parse_panels_from_caption(
        "Fig. 90. Sulculeolaria chuni (L. & v. R.). "
        "A, anterior; B, posterior; C, ventral."
    )
    assert [panel["label"] for panel in panels] == ["A", "B", "C"]


# --- the denominator ----------------------------------------------------------


def test_only_numbers_inside_a_figure_block_are_counted():
    """A body-text reference must not become a figure the pipeline missed.

    This is the difference between a recall of 0.296 and one that means
    something.
    """
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = Path(d)
        (p / "page_001.txt").write_text(
            "[PAGE 1]\n"
            "As figured by Bigelow (op. cit., fig. 34) and see Fig. 18 below.\n"
            "[FIGURE]\nFig. 12. Nectophore of Nanomia.\n[/FIGURE]\n",
            encoding="utf-8")
        blocks = list(cb.gold_blocks(p))
        assert len(blocks) == 1
        _page, _kind, body = blocks[0]
        _k, nums, _f = cb.classify_block(body)
        assert nums == {"12"}, "references outside the block must not count"


def test_explicit_gold_plate_inventory_counts_every_engraved_number():
    body = (
        "Figure numbers on the plate, in reading order: "
        "1, 2, 3, 4, 5, 6, 7a, 7b, 8, 9, 10, 11, 12, 13.\n\n"
        "Lettering on the individual figures, in reading order:\n"
        "Fig. 2: Sb H\nFig. 7a: W\nFig. 13: H\n"
    )

    kind, numbers, first_line = cb.classify_block(body)

    assert kind == "bare_label"
    assert numbers == {
        "1", "2", "3", "4", "5", "6", "7a", "7b", "8", "9",
        "10", "11", "12", "13",
    }
    assert first_line == "Fig. 2: Sb H"


def test_plate_inventory_lettering_notes_do_not_become_caption_prose():
    body = (
        "Figure numbers on the plate, in reading order: 1, 2.\n"
        "Lettering, figure by figure:\n"
        "Fig. 1 [illegible] (single italic letters scattered over the drawing)\n"
    )

    kind, numbers, _first_line = cb.classify_block(body, structural_kind="PLATE")

    assert kind == "bare_label"
    assert numbers == {"1", "2"}


def test_plate_inventory_rejects_non_number_prose():
    assert cb.gold_plate_number_list(
        "Figure numbers on the plate: 1, two, 3."
    ) == []


def test_a_nested_figure_marker_is_not_a_second_block():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = Path(d)
        (p / "page_001.txt").write_text(
            "[PAGE 1]\n[NOTE: lettering inside [FIGURE] is transcribed]\n"
            "[FIGURE]\nFig. 3. A real one.\n[/FIGURE]\n", encoding="utf-8")
        assert len(list(cb.gold_blocks(p))) == 1


def test_plate_heading_immediately_outside_block_is_identity_evidence():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = Path(d)
        (p / "page_017.txt").write_text(
            "[PAGE 17]\nPLATE XVII\n[PLATE]\n1\n2\n3\n[/PLATE]\n",
            encoding="utf-8",
        )

        page, structural_kind, body = next(cb.gold_blocks(p))
        kind, numbers, _ = cb.classify_block(body, structural_kind)

        assert page == 17
        assert kind == "bare_label"
        assert numbers == {"1", "2", "3", "17"}


def test_nonadjacent_plate_reference_is_not_host_identity(tmp_path):
    (tmp_path / "page_017.txt").write_text(
        "[PAGE 17]\nPLATE IX mentioned in an index\nIntervening prose.\n"
        "[PLATE]\n1\n2\n[/PLATE]\n",
        encoding="utf-8",
    )

    _page, structural_kind, body = next(cb.gold_blocks(tmp_path))
    _kind, numbers, _first = cb.classify_block(body, structural_kind)

    assert numbers == {"1", "2"}


def test_plate_and_same_numbered_child_are_distinct_scoring_identities(tmp_path):
    gold = tmp_path / "gold"
    corpus = tmp_path / "corpus"
    gold.mkdir()
    corpus.mkdir()
    (gold / "page_010.txt").write_text(
        "[PAGE 10]\nPLATE X\n[PLATE]\n10.\n[/PLATE]\n",
        encoding="utf-8",
    )
    (corpus / "scan_detection.json").write_text("{}", encoding="utf-8")
    figures = [{
        "figure_id": "plate_10",
        "figure_type": "plate",
        "figure_number": "10",
        "page": 10,
    }]
    (corpus / "figures.json").write_text(
        json.dumps({"figures": figures}), encoding="utf-8",
    )

    plate_only = cb.score_document(gold, corpus)

    assert plate_only["gold_numbers"] == 2
    assert plate_only["found_numbers"] == 1
    assert plate_only["matched_numbers"] == 1
    assert plate_only["pages"]["10"]["missing_identities"] == ["figure:10"]

    figures.append({
        "figure_id": "figure_10",
        "figure_type": "figure",
        "figure_number": "10",
        "page": 10,
    })
    (corpus / "figures.json").write_text(
        json.dumps({"figures": figures}), encoding="utf-8",
    )

    complete = cb.score_document(gold, corpus)
    assert complete["gold_numbers"] == complete["found_numbers"] == 2
    assert complete["matched_numbers"] == 2


# --- reporting ----------------------------------------------------------------


@pytest.fixture(scope="module")
def report():
    return cb.build_report(FIXTURE / "gold", FIXTURE / "corpuscle")


def test_documents_bind_and_unmatched_are_reported(report):
    assert report["schema_version"] == 7
    assert report["documents_bound"] == 2
    assert report["documents_unmatched"] == ["DocTwo"]


def test_block_kinds_are_reported_so_no_rate_hides_them(report):
    assert report["gold_block_kinds"] == {"bare_label": 1}


def test_precision_is_none_rather_than_zero_when_nothing_was_reported(report):
    """The fixture's entries carry no figure_number. Precision is undefined,
    not zero — printing 0.000 would read as "every number was wrong"."""
    assert report["totals"]["found_numbers"] == 0
    assert report["totals"]["number_precision"] is None
    assert report["totals"]["number_recall"] == 0.0
    assert report["totals"]["reported_pair_capacity"] == 0
    assert report["totals"]["reported_pair_capacity_rate"] == 0.0


def test_f1_is_zero_when_rates_are_defined_but_nothing_matches():
    packed = cb._pack_counts({
        "gold_numbers": 2,
        "found_numbers": 3,
        "matched_numbers": 0,
        "captions_compared": 0,
    })
    assert packed["f1"] == 0.0


def test_reported_pair_capacity_rate_is_derived_from_aggregate_counts():
    packed = cb._pack_counts({
        "gold_numbers": 4,
        "found_numbers": 3,
        "matched_numbers": 2,
        "reported_pair_capacity": 3,
        "captions_compared": 0,
    })
    assert packed["reported_pair_capacity_rate"] == 0.75


def test_the_report_states_why_it_does_not_use_caption_similarity(report):
    assert "artifact" in report["caveat"]
    assert "not a ceiling" in report["caveat"]


def test_report_states_independent_panel_measurement(report):
    assert "independently" in report["panel_method"]
    totals = report["surfaces"]["default MCP types"]["totals"]
    assert totals["gold_panel_figures"] == 0
    assert totals["panel_exact_rate"] is None
    assert "panel_figure_precision" in totals


def test_segments_include_rates_by_era_file_type_and_layout(report):
    assert set(report["segments"]) == {"era", "file_type", "layout"}
    plate = report["segments"]["layout"]["plate_blocks"]
    assert plate["gold_numbers"] > 0
    assert "number_recall" in plate and "number_precision" in plate


def test_raw_and_default_mcp_surfaces_are_reported_separately(report):
    assert set(report["surfaces"]) == {"all entries", "default MCP types"}
    raw = report["surfaces"]["all entries"]["totals"]
    served = report["surfaces"]["default MCP types"]["totals"]
    assert served["found_numbers"] <= raw["found_numbers"]
    assert set(report["surfaces"]["default MCP types"]["segments"]) == {
        "era", "file_type", "layout",
    }


def test_page_diagnostics_name_missing_and_reported_evidence(report):
    page = report["documents"]["DocOne"]["pages"]["3"]
    assert page["gold_numbers"]
    assert page["missing_numbers"] == page["gold_numbers"]
    assert page["reported_figures"]
    assert page["reported_figures"][0]["caption_status"] == "bound"


def test_summary_renders(report, capsys):
    cb.print_summary(report)
    out = capsys.readouterr().out
    assert "binding recall" in out and "binding precision" in out
    assert "reported-pair capacity" in out
    assert "by era" in out and "by layout" in out
