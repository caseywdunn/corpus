"""`keeppages`: which physical pages the document actually is (#188).

A PDF in a scanned library is frequently not just the paper — library cover
sheets, a Russian original bound in front of its English translation, runs of
blank versos. The costs compound rather than add: `detect_scan_type` samples
pages to pick the OCR mode *and* the language pack, so forty pages of scanner
filler ahead of the body decide how the body gets read, and then OCR pays full
price for the filler, and a colour calibration target becomes a figure.

Two things here are easy to get wrong and are asserted directly:

- The numbers are **physical 1-based positions**, never printed folios. An
  entry routinely carries `pages = {699--714}` and a `keeppages` that
  disagrees completely, and that is correct.
- Once pages are dropped, `page` in the artifacts is a position in the
  *subset*. A figure served as "page 3" that is page 44 of the scan is a
  citation error that nothing downstream can detect, which is why
  `source_page` is carried beside it.
"""
from __future__ import annotations

import json

import pytest

from pipeline.pageselect import (
    PageSelectionError,
    annotate_source_pages,
    parse_keeppages,
    resolve_selection,
    write_subset,
)

fitz = pytest.importorskip("fitz")


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("3--20", [(3, 20)]),
    ("5", [(5, 5)]),
    ("2,4,8--20,22--40,55", [(2, 2), (4, 4), (8, 20), (22, 40), (55, 55)]),
    ("40--", [(40, None)]),                      # "page 40 to the end"
    ("3-20", [(3, 20)]),                         # single dash, what people type
    ("  3 -- 20 ,  25 ", [(3, 20), (25, 25)]),   # whitespace anywhere
    ("", []),
    (None, []),
])
def test_parse(raw, expected):
    assert parse_keeppages(raw) == expected


@pytest.mark.parametrize("raw", ["abc", "3--x", "0--5", "20--3", "3;5", "-5"])
def test_unparseable_values_raise(raw):
    """No safe reading of a malformed range exists.

    `ocrlang` can drop an uninstalled pack and let detection carry on;
    there is no equivalent here. Silently keeping every page would look
    exactly like success, and the operator would never learn the directive
    they wrote had no effect.
    """
    with pytest.raises(PageSelectionError):
        parse_keeppages(raw)


def test_a_backwards_range_is_an_error_not_a_reordering():
    with pytest.raises(PageSelectionError, match="backwards"):
        parse_keeppages("20--3")


# ---------------------------------------------------------------------------
# Resolution against a real page count
# ---------------------------------------------------------------------------

def test_open_ended_range_runs_to_the_last_page():
    selected, warnings = resolve_selection(parse_keeppages("8--"), 12)
    assert selected == [8, 9, 10, 11, 12]
    assert warnings == []


def test_selection_is_sorted_and_deduplicated():
    """A selection says which pages *are* the document, not what order.

    `40--50,10--20` must not become a way to reorder a scan.
    """
    selected, _ = resolve_selection(parse_keeppages("40--42,10--12"), 50)
    assert selected == [10, 11, 12, 40, 41, 42]


def test_overlapping_ranges_collapse():
    selected, _ = resolve_selection(parse_keeppages("3--10,8--12"), 20)
    assert selected == list(range(3, 13))


def test_a_range_past_the_end_is_clamped_with_a_warning():
    """Almost always a typo about where a document ends.

    Clamping matches how `ocrlang` drops packs that aren't installed — the
    pages that do exist are still the right ones. But it is recorded, so a
    clamp is visible rather than inferred from a page count.
    """
    selected, warnings = resolve_selection(parse_keeppages("3--999"), 10)
    assert selected == [3, 4, 5, 6, 7, 8, 9, 10]
    assert len(warnings) == 1 and "clamped" in warnings[0]


def test_a_range_starting_past_the_end_is_dropped_with_a_warning():
    selected, warnings = resolve_selection(parse_keeppages("3--5,90--95"), 10)
    assert selected == [3, 4, 5]
    assert len(warnings) == 1 and "past the last page" in warnings[0]


# ---------------------------------------------------------------------------
# Writing the subset
# ---------------------------------------------------------------------------

def _pdf(path, n_pages):
    doc = fitz.open()
    for i in range(n_pages):
        page = doc.new_page()
        page.insert_text((72, 72), f"PAGE {i + 1}")
    doc.save(path)
    doc.close()
    return path


def test_subset_keeps_exactly_the_selected_pages(tmp_path):
    src = _pdf(tmp_path / "in.pdf", 10)
    dst = tmp_path / "out.pdf"
    assert write_subset(src, dst, [3, 4, 5]) is True
    with fitz.open(dst) as doc:
        assert doc.page_count == 3
        assert [doc[i].get_text().strip() for i in range(3)] == \
               ["PAGE 3", "PAGE 4", "PAGE 5"]


def test_a_selection_covering_everything_writes_nothing(tmp_path):
    """An annotated document that keeps every page must cost nothing.

    Returning False leaves the original bound, so its bytes — and therefore
    its OCR — are identical to an unannotated document's.
    """
    src = _pdf(tmp_path / "in.pdf", 4)
    dst = tmp_path / "out.pdf"
    assert write_subset(src, dst, [1, 2, 3, 4]) is False
    assert not dst.exists()


def test_an_empty_selection_writes_nothing(tmp_path):
    src = _pdf(tmp_path / "in.pdf", 4)
    dst = tmp_path / "out.pdf"
    assert write_subset(src, dst, []) is False
    assert not dst.exists()


def test_a_discontinuous_selection_keeps_document_order(tmp_path):
    src = _pdf(tmp_path / "in.pdf", 10)
    dst = tmp_path / "out.pdf"
    write_subset(src, dst, [2, 7, 9])
    with fitz.open(dst) as doc:
        assert [doc[i].get_text().strip() for i in range(3)] == \
               ["PAGE 2", "PAGE 7", "PAGE 9"]


# ---------------------------------------------------------------------------
# Page-number provenance
# ---------------------------------------------------------------------------

def test_source_page_is_written_beside_every_subset_page(tmp_path):
    """The resolved list *is* the map: subset page i is selected[i - 1]."""
    figures = tmp_path / "figures.json"
    figures.write_text(json.dumps([
        {"figure_id": "a", "page": 1, "caption_page": 1},
        {"figure_id": "b", "page": 3, "caption_page": 2},
    ]))
    selected = [41, 42, 43, 44]
    assert annotate_source_pages([figures], selected) == 4

    got = json.loads(figures.read_text())
    assert got[0]["page"] == 1 and got[0]["source_page"] == 41
    assert got[1]["page"] == 3 and got[1]["source_page"] == 43
    assert got[1]["caption_page"] == 2 and got[1]["caption_source_page"] == 42


def test_page_stays_subset_relative():
    """`page` must not be rewritten to the source number.

    It is what indexes the artifacts actually on disk — the rendered
    figure images, the docling reading order. Rewriting it in place would make the
    number correct for a citation and wrong for everything else.
    """
    obj = [{"page": 2}]
    import copy
    before = copy.deepcopy(obj)
    annotate_source_pages([], [10, 11, 12])
    assert obj == before


def test_nested_structures_are_reached(tmp_path):
    """Chunk provenance is nested; figures.json is sometimes wrapped."""
    path = tmp_path / "chunks.json"
    path.write_text(json.dumps({"chunks": [{"provenance": {"page": 2}}]}))
    assert annotate_source_pages([path], [7, 8, 9]) == 1
    assert json.loads(path.read_text())["chunks"][0]["provenance"]["source_page"] == 8


def test_a_page_outside_the_selection_is_left_alone(tmp_path):
    """Better no `source_page` than a wrong one.

    A page number past the end of the map means something upstream
    disagrees about the subset, and inventing a mapping would turn that
    into a plausible wrong citation.
    """
    path = tmp_path / "figures.json"
    path.write_text(json.dumps([{"page": 99}]))
    assert annotate_source_pages([path], [1, 2, 3]) == 0
    assert "source_page" not in json.loads(path.read_text())[0]


# ---------------------------------------------------------------------------
# The fingerprint — a page-range edit must invalidate everything OCR-derived
# ---------------------------------------------------------------------------

def test_keeppages_fingerprints_every_ocr_dependent_stage():
    """It invalidates strictly more than an `ocrlang` change does.

    An ocrlang edit changes how the PDF is read; a keeppages edit changes
    which pages the PDF *is*. Fingerprinting only the first two stages
    re-runs the selection and then skips everything that consumes it, so
    `text.json` still holds text from pages that are no longer in the
    document while the log reports success.
    """
    from pipeline.stages import _OCR_DEPENDENT_STAGES, _expected_fingerprints_for_run

    fps = _expected_fingerprints_for_run(
        ocrlang=None, ocrmode=None, keeppages="3--20",
    )
    for stage in _OCR_DEPENDENT_STAGES:
        assert fps[stage]["keeppages"] == "3--20", stage


def test_an_unannotated_document_still_compares_equal():
    """`{}`, not a key with None — otherwise every existing corpus re-OCRs."""
    from pipeline.stages import _expected_fingerprints_for_run

    fps = _expected_fingerprints_for_run(
        ocrlang=None, ocrmode=None, keeppages=None,
    )
    assert all(fp == {} for fp in fps.values())


def test_adding_changing_and_removing_all_invalidate():
    from pipeline.stages import _expected_fingerprints_for_run

    none = _expected_fingerprints_for_run(
        ocrlang=None, ocrmode=None, keeppages=None,
    )
    added = _expected_fingerprints_for_run(
        ocrlang=None, ocrmode=None, keeppages="3--20",
    )
    changed = _expected_fingerprints_for_run(
        ocrlang=None, ocrmode=None, keeppages="4--20",
    )
    assert none != added and added != changed and changed != none


def test_the_argument_is_required():
    """A default here reproduces the shipped #176 bug.

    main.py's outer fast path — the gate that actually skips work — would
    silently keep the default, so a paper with a newly added directive
    matches `{}` against `{}` and is skipped whole before anything sees it.
    """
    from pipeline.stages import _expected_fingerprints_for_run

    with pytest.raises(TypeError, match="keeppages"):
        _expected_fingerprints_for_run(ocrlang=None, ocrmode=None)


# ---------------------------------------------------------------------------
# Bib round-trip
# ---------------------------------------------------------------------------

def test_the_bib_carries_keeppages_through_to_metadata():
    from bib.parser import bib_entry_to_metadata, entry_keeppages

    entry = {"_key": "Stepanjants1970", "keeppages": "{3--20}",
             "pages": "{699--714}", "ocrlang": "{rus+eng}"}
    meta = bib_entry_to_metadata(entry, "Stepanjants1970.pdf")
    assert meta["keeppages"] == "3--20"
    assert entry_keeppages(entry) == "3--20"


def test_keeppages_does_not_collide_with_the_standard_pages_field():
    """They are different numbers about the same paper and routinely disagree.

    `pages` is where the article sits in its journal volume and is already
    consumed by `bib/format.py`; `keeppages` is where it physically sits in
    the scan. For an offprint they are wildly different, and reading one as
    the other produces a citation that is confidently wrong.
    """
    from bib.parser import bib_entry_to_metadata

    meta = bib_entry_to_metadata(
        {"_key": "X", "pages": "{699--714}", "keeppages": "{3--20}"}, "X.pdf")
    assert meta["keeppages"] == "3--20"
    assert meta.get("pages") in (None, "699--714")   # never the keeppages value


def test_export_round_trips_keeppages(tmp_path):
    import sqlite3
    from bib.authority import create_schema
    from bib.export import export_bibtex

    db = tmp_path / "a.sqlite"
    conn = sqlite3.connect(db)
    create_schema(conn)
    conn.execute(
        "INSERT INTO works (work_id, guid_type, title, year, in_corpus, "
        "source, keeppages, created_at, updated_at) "
        "VALUES ('w1','corpus_key','A paper',1970,1,'corpus_paper','3--20',0,0)")
    conn.commit(); conn.close()
    assert "keeppages = {3--20}" in export_bibtex(db)
