"""Figure detection scoring against the gold set (#194).

Two questions, and they are not the same one: is every figure found (recall),
and is publisher furniture being called a figure (precision)? Caption
association is #195, and OCR of text *inside* a figure is explicitly not a
target.

`tools/qc/figure_detection.py` answers both by counting per page. The gold set
records no bounding boxes, so an individual gold block cannot be matched to an
individual `figures.json` entry — but per-page counts can be, and the
corpus-wide totals are actively misleading: 424 gold blocks against 420
entries looks like agreement and is a coincidence, hiding one paper with 6
gold blocks against 31 entries and another with 67 against 34.

These tests use the committed fixture and synthetic counts, with no corpus
dependency.
"""
from __future__ import annotations

import importlib.util
from collections import Counter
from pathlib import Path

import pytest

_PATH = Path(__file__).parent.parent / "tools" / "qc" / "figure_detection.py"
_spec = importlib.util.spec_from_file_location("qc_figure_detection", _PATH)
fd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fd)

FIXTURE = Path(__file__).parent / "fixtures" / "gold_fidelity"


@pytest.fixture(scope="module")
def report():
    return fd.build_report(FIXTURE / "gold", FIXTURE / "corpuscle")


# --- what counts as a gold figure ---------------------------------------------


def test_tables_are_not_figures(report):
    """`[TABLE]` is body content — `fidelity.py` scores it as prose for the
    same reason, and docling reports tables separately from pictures, so
    counting them here would compare two different populations.

    The fixture's DocOne page 2 carries a `[TABLE]` and nothing else.
    """
    assert report["gold_block_kinds"] == {"PLATE": 1}


def test_a_figure_quoted_inside_a_note_is_not_a_second_figure():
    """Notes discuss figures — `[NOTE: ... inside [FIGURE] ...]` occurs in the
    real set — and the shared span scanner already declines to yield a nested
    marker. Counting it would inflate the gold."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = Path(d)
        (p / "page_001.txt").write_text(
            "[PAGE 1]\n[NOTE: lettering inside [FIGURE] is transcribed]\n"
            "[FIGURE]\nFig. 1. A real one.\n[/FIGURE]\n", encoding="utf-8")
        assert fd.gold_blocks_by_page(p)[1] == Counter({"FIGURE": 1})


# --- the counting ------------------------------------------------------------


def test_score_counts_missed_and_surplus_per_page():
    """Recall and precision failures are per page and can coexist in one
    document — a page missing a figure does not offset a page inventing one."""
    gold = {1: Counter({"FIGURE": 2}), 2: Counter({"PLATE": 1}), 3: Counter()}
    entries = [
        {"page": 1, "figure_type": "figure"},        # page 1: 2 gold, 1 found
        {"page": 3, "figure_type": "figure"},        # page 3: 0 gold, 2 found
        {"page": 3, "figure_type": "figure"},
    ]
    t = fd.score(gold, entries, lambda f: True)
    assert t["gold"] == 3 and t["found"] == 3
    assert t["matched"] == 1          # only page 1 matches, and only one of two
    assert t["missed"] == 2           # one on page 1, one on page 2
    assert t["surplus"] == 2          # both on page 3
    # The totals agree at 3 and 3 while nothing is right — which is exactly why
    # this is counted per page rather than corpus-wide.
    assert t["pages_exact"] == 0


def test_a_page_with_the_right_count_is_exact():
    gold = {1: Counter({"FIGURE": 1})}
    t = fd.score(gold, [{"page": 1, "figure_type": "figure"}], lambda f: True)
    assert (t["matched"], t["missed"], t["surplus"], t["pages_exact"]) == (1, 0, 0, 1)


def test_entries_without_a_page_are_dropped_not_guessed_at():
    gold = {1: Counter({"FIGURE": 1})}
    t = fd.score(gold, [{"page": None, "figure_type": "figure"}], lambda f: True)
    assert t["found"] == 0 and t["missed"] == 1


# --- the furniture filters ----------------------------------------------------


def test_filters_are_scored_side_by_side(report):
    """The furniture question is a filter question, and the point of the
    scorer is to answer it with evidence rather than by which field name
    sounds most like furniture."""
    assert set(report["filters"]) == set(fd.FILTERS)


def test_dropping_graphical_element_raises_precision(report):
    """On the fixture, three of the four surplus entries are
    `graphical_element` — journal-furniture marks on pages the gold says carry
    no figure."""
    base = report["filters"]["all entries"]
    filtered = report["filters"]["drop graphical_element"]
    assert base["precision"] == 0.2
    assert filtered["precision"] == 0.5
    assert filtered["recall"] == base["recall"], "no real figure may be lost"


def test_dropping_unclassified_too_is_over_reach(report):
    """`unclassified` is a mixed population — it holds real figures as well as
    furniture — so a filter that takes it costs recall. The fixture keeps one
    `unclassified` surplus entry so this stays visible."""
    assert "unclassified" in str(report["surplus_page_profile"])


def test_surplus_profile_names_what_sits_on_over_counted_pages(report):
    assert report["surplus_page_profile"]["graphical_element/uncaptioned"] == 3
    assert report["surplus_page_profile"]["unclassified/uncaptioned"] == 1


# --- reporting ----------------------------------------------------------------


def test_a_document_with_no_gold_figures_has_no_recall_not_zero_recall(report):
    """DocThree's gold marks no figures. Recall is undefined there, and
    printing 0.000 would read as a failure to find figures that do not exist —
    while precision is genuinely 0, because something was found."""
    d = report["documents"]["DocThree"]
    assert d["gold_figures"] == 0
    assert d["recall"] is None
    assert d["precision"] == 0.0


def test_segments_split_the_axes_that_behave_differently(report):
    """Born-digital papers carry publisher furniture that historical scans do
    not, so one corpus-wide precision number mixes two unlike populations."""
    assert set(report["segments"]) == {"era", "file_type"}
    assert set(report["segments"]["file_type"]) == {"scanned", "born_digital"}


def test_unmatched_gold_documents_are_reported(report):
    assert report["documents_unmatched"] == ["DocTwo"]


def test_report_states_what_page_level_counting_cannot_do(report):
    """The caveat is part of the output, not folded knowledge — a page-level
    count cannot tell a page that found the right two figures from one that
    found two wrong ones."""
    assert "cannot" in report["caveat"]


def test_summary_renders_without_a_rate(report, capsys):
    """Documents with no gold figures make rates legitimately None; the
    formatter must print a dash rather than crash or show 0.000."""
    fd.print_summary(report)
    out = capsys.readouterr().out
    assert "recall" in out and "precis" in out


# --- panel siblings are one figure, not several (#211) -------------------------
#
# A figure with panels produces one image per panel — `fig_99_a.png`,
# `fig_99_b.png` — while the gold records one `[FIGURE]` block carrying both,
# with the panels enumerated inside its caption. Counting entries against
# blocks penalises the pipeline for splitting panels correctly, and it flips
# the sign of the error: `Totton1965a` reads as over-counting by 3 entries
# when it is in fact under-counting by 4 figures.


def _entry(page, num, fid):
    return {"page": page, "figure_number": num, "figure_id": fid,
            "figure_type": "figure", "caption_text": "x"}


def test_panel_siblings_collapse_to_one_figure():
    figures, groups = fd.collapse_panels([
        _entry(168, "99", "a"), _entry(168, "99", "b"), _entry(169, "100", "c"),
    ])
    assert len(figures) == 2
    assert groups == {(168, "99"): 2}


def test_same_number_on_different_pages_stays_two_figures():
    """The translation case again — collapsing must key on page as well."""
    figures, groups = fd.collapse_panels([_entry(3, "4", "a"), _entry(24, "4", "b")])
    assert len(figures) == 2
    assert groups == {}


def test_unnumbered_entries_are_never_merged():
    """Without a number there is nothing to group on, and merging two real
    figures is worse than counting one twice."""
    figures, groups = fd.collapse_panels([
        {"page": 5, "figure_number": None, "figure_type": "figure"},
        {"page": 5, "figure_number": "", "figure_type": "figure"},
    ])
    assert len(figures) == 2 and groups == {}


def test_a_correctly_split_panel_is_not_scored_as_surplus():
    """The defect in one assertion: one gold block, two panel images, and the
    page must score as matched rather than one match plus one false positive."""
    gold = {168: Counter({"FIGURE": 1})}
    entries = [_entry(168, "99", "a"), _entry(168, "99", "b")]
    t = fd.score(gold, entries, lambda f: True)
    assert t["found"] == 1
    assert t["matched"] == 1
    assert t["surplus"] == 0
    assert t["missed"] == 0


def test_panel_splitting_is_reported_beside_detection(report):
    """Not discarded — a split panel is a correct extra image, and #195 scores
    whether the split matches the panels the gold caption enumerates."""
    for d in report["documents"].values():
        assert "entries" in d and "panel_groups" in d and "panel_images" in d
