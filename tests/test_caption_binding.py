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

* **Count only numbers printed *inside* a figure block.** Gold pages are full
  of figure numbers that are references — "see Fig. 18", "figured by Bigelow
  (op. cit., fig. 34)". Counting those measures objects that are not on that
  page. Restricted to numbers inside a `[FIGURE]`/`[PLATE]` block the measure
  means what it says.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_PATH = Path(__file__).parent.parent / "tools" / "qc" / "caption_binding.py"
_spec = importlib.util.spec_from_file_location("qc_caption_binding", _PATH)
cb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cb)

FIXTURE = Path(__file__).parent / "fixtures" / "gold_fidelity"


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


def test_transcriber_commentary_does_not_contribute_numbers():
    """A note mentioning Fig. 4 is not a figure on the page."""
    _kind, nums, _ = cb.classify_block(
        "[NOTE: this plate is reproduced from Fig. 4 of Bigelow 1911]\n"
        "Fig. 12. Nectophore.")
    assert nums == {"12"}


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


def test_a_nested_figure_marker_is_not_a_second_block():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = Path(d)
        (p / "page_001.txt").write_text(
            "[PAGE 1]\n[NOTE: lettering inside [FIGURE] is transcribed]\n"
            "[FIGURE]\nFig. 3. A real one.\n[/FIGURE]\n", encoding="utf-8")
        assert len(list(cb.gold_blocks(p))) == 1


# --- reporting ----------------------------------------------------------------


@pytest.fixture(scope="module")
def report():
    return cb.build_report(FIXTURE / "gold", FIXTURE / "corpuscle")


def test_documents_bind_and_unmatched_are_reported(report):
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


def test_f1_is_zero_when_rates_are_defined_but_nothing_matches():
    packed = cb._pack_counts({
        "gold_numbers": 2,
        "found_numbers": 3,
        "matched_numbers": 0,
        "captions_compared": 0,
    })
    assert packed["f1"] == 0.0


def test_the_report_states_why_it_does_not_use_caption_similarity(report):
    assert "artifact" in report["caveat"]


def test_segments_include_rates_by_era_file_type_and_layout(report):
    assert set(report["segments"]) == {"era", "file_type", "layout"}
    plate = report["segments"]["layout"]["plate_blocks"]
    assert plate["gold_numbers"] > 0
    assert "number_recall" in plate and "number_precision" in plate


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
    assert "by era" in out and "by layout" in out
