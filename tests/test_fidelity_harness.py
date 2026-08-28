"""The gold-transcription fidelity harness, on a committed fixture (#193).

`tools/qc/fidelity.py` scores a built corpuscle's extraction against an
independently transcribed gold set. The full run is manual and release-time
against a real corpuscle; this file covers the scorer's own arithmetic with no
corpus dependency, so a change that silently inverts a measure fails in T0
rather than at the next release.

The fixture under `tests/fixtures/gold_fidelity/` is three pages of one
document plus one page of a second, each page chosen to pin one behaviour:

  DocOne p1   verbatim agreement            -> coverage 1.0
  DocOne p2   columns emitted out of order  -> coverage 1.0, similarity ~0.41
              plus a docling table, whose text lives only in
              `data.table_cells[]` and is invisible to a `text`-only walk
  DocOne p3   nothing extracted at all      -> scored 0.0, NOT skipped
  DocThree p1 Cyrillic page, Latin output   -> `script_missing`
  DocTwo      transcribed, not in corpuscle -> reported unmatched

Run locally with::

    python -m pytest tests/test_fidelity_harness.py
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

# tools/ is not an importable package (and is excluded from the wheel), so the
# script is loaded by path — the same approach the other tools/ tests take.
_FIDELITY_PATH = Path(__file__).parent.parent / "tools" / "qc" / "fidelity.py"
_spec = importlib.util.spec_from_file_location("qc_fidelity", _FIDELITY_PATH)
fid = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fid)

FIXTURE = Path(__file__).parent / "fixtures" / "gold_fidelity"
GOLD = FIXTURE / "gold"
CORPUSCLE = FIXTURE / "corpuscle"


@pytest.fixture(scope="module")
def report():
    return fid.build_report(GOLD, CORPUSCLE)


def _page(report, document, n):
    for d in report["detail"]:
        if d["document"] == document:
            for p in d["pages"]:
                if p["page"] == n:
                    return p
    raise AssertionError(f"no page {n} of {document} in the report")


# --- binding ------------------------------------------------------------------


def test_documents_bind_on_checksum_not_stem(report):
    """Directory names are sha256 prefixes bearing no relation to the stem.

    Filenames are not unique across a library — two editions of the same 1594
    travel narrative were both once named `Lery1594.pdf`, setting the same
    passage on different folios, with only one transcribed. Binding on the stem
    would score one edition's gold against the other's extraction and report
    the disagreement as an extractor defect.
    """
    assert report["documents_bound"] == 2
    assert set(report["documents"]) == {"DocOne", "DocThree"}
    # The fixture's directories are `a1b2c3d4e5f6` and `c3d4e5f6a1b2`; nothing
    # in either name could be derived from "DocOne" or "DocThree".
    assert not (CORPUSCLE / "documents" / "DocOne").exists()


def test_gold_document_missing_from_corpuscle_is_reported_not_skipped(report):
    """A transcribed document that was never built has to surface.

    Silently dropping it would let a corpuscle score well on the subset it
    happened to contain.
    """
    assert report["documents_unmatched"] == ["DocTwo"]


# --- the measures -------------------------------------------------------------


def test_verbatim_page_scores_one(report):
    p = _page(report, "DocOne", 1)
    assert p["status"] == "ok"
    assert p["coverage"] == 1.0
    assert p["recall"] == 1.0
    assert p["similarity"] == 1.0
    assert p["excess_novel"] == 0.0


def test_reading_order_difference_spares_coverage_but_not_similarity(report):
    """The single most important property of the pair of measures.

    Page 2's two prose columns are extracted in the reverse order. Every word
    is present, so `coverage` — the headline measure — is 1.0, while the
    order-sensitive `similarity` collapses. A harness reporting only
    `similarity` would send someone to fix an extraction that is correct.
    """
    p = _page(report, "DocOne", 2)
    assert p["status"] == "ok"
    assert p["coverage"] == 1.0
    assert p["similarity"] < 0.6
    assert p["excess_novel"] == 0.0


def test_table_cell_text_is_recovered_from_the_docling_artifact(report):
    """A docling table has no `text` field: its content is in
    `data.table_cells[].text` only.

    Page 2's gold holds a three-row table, and the extraction carries it as a
    table item. A walk reading only `text` drops every cell, and the page then
    scores as lost prose. Measured on the real gold set before this was fixed,
    it cost up to 26% of a document's tokens.
    """
    p = _page(report, "DocOne", 2)
    # 39 gold tokens, 8 of which are the table's cells; coverage is only 1.0 if
    # all of them were found.
    assert p["gold_words"] == p["extracted_words"] == 39
    assert p["coverage"] == 1.0

    doc = fid._read_json(
        CORPUSCLE / "documents" / "a1b2c3d4e5f6" / "docling_doc.json")
    text = fid.page_texts_from_docling(doc)[2]
    assert "Nanomia bijuga" in text and "Physalia physalis" in text
    # Cells are stored unsorted in the fixture; the walk must re-impose
    # row-major order rather than trusting the JSON's sequence.
    assert text.index("Species") < text.index("Nanomia")


def test_empty_extraction_is_scored_as_failure_not_excluded(report):
    """Which side is on trial, in one assertion.

    The library-side cross-check treats an empty comparison layer as `no_signal`
    and drops the page, because there it says nothing about the gold. Here the
    corpuscle is what is on trial, so an empty extraction of a page the gold
    says carries text is the finding — and scoring it None would delete the
    pipeline's worst pages from every median. Measured on the real 35-document
    set, the exclusion policy would hide 57 of 761 pages and lift the median
    coverage from 0.891 to 0.908.
    """
    p = _page(report, "DocOne", 3)
    assert p["status"] == "extraction_empty"
    assert p["coverage"] == 0.0
    assert p["similarity"] == 0.0
    # It counts toward `scored`, so the medians include it.
    assert report["documents"]["DocOne"]["scored"] == 3
    assert report["documents"]["DocOne"]["pages_below_0.5_coverage"] == 1


def test_script_missing_is_flagged_and_still_scored(report):
    """A page whose writing system never reached the extraction.

    This is an OCR language-pack question (#176) rather than a text-quality
    one, so it gets its own status — but it is a real failure and stays in the
    scores, unlike on the library side where the same condition means "no
    signal".
    """
    p = _page(report, "DocThree", 1)
    assert p["script"] == "cyrillic"
    assert p["status"] == "script_missing"
    assert p["coverage"] == 0.0
    # Latin salad where Cyrillic belongs: word-like, but none of it near a gold
    # word.
    assert p["excess_novel"] == 1.0


# --- markup handling ----------------------------------------------------------


def test_printed_folio_and_pdf_index_are_kept_apart(report):
    """`page_001.txt` is PDF page 1; `[PAGE 3]` is the folio printed on it.

    Page-number provenance on served figures is an open question (#188), and
    this set is the only place both numbers are known for the same leaf, so
    conflating them here would waste the evidence.
    """
    p1 = _page(report, "DocOne", 1)
    assert p1["page"] == 1 and p1["printed_page"] == 3
    p3 = _page(report, "DocOne", 3)
    assert p3["printed_page"] is None and p3["printed_page_unnumbered"] is True


def test_transcriber_commentary_is_stripped_but_printed_text_is_kept():
    """`[NOTE:]` is commentary; `[RUNNING HEAD]` labels text that was printed.

    Keeping the first would score the pipeline against words that are not on
    the page; dropping the second would penalise it for finding words that are.
    """
    raw = (GOLD / "DocOne" / "page_001.txt").read_text(encoding="utf-8")
    stripped = fid.strip_markup(raw)
    assert "printer's ornament" not in stripped.text     # [NOTE:] dropped whole
    assert "Histoire Naturelle" in stripped.text         # [RUNNING HEAD] label only
    assert "PAGE 3" not in stripped.text


def test_uncertain_reading_keeps_the_transcribers_best_guess():
    """`[?reading:gastrozooids]` is what the transcriber saw, so it is compared."""
    raw = (GOLD / "DocOne" / "page_003.txt").read_text(encoding="utf-8")
    stripped = fid.strip_markup(raw)
    assert "gastrozooids" in stripped.text
    assert "?reading" not in stripped.text
    assert _page(fid.build_report(GOLD, CORPUSCLE), "DocOne", 3)["uncertain"] == 1


def test_structural_share_separates_plate_lettering_from_prose(report):
    """A page of engraved plate lettering fails differently from a page of prose.

    Page 3 is entirely inside a `[PLATE]` block and page 1 entirely outside one.
    Reporting the share lets a low coverage be read as "missed the plate" rather
    than "missed the text", which are different pieces of work.
    """
    assert _page(report, "DocOne", 3)["gold_structural_share"] == 1.0
    assert _page(report, "DocOne", 1)["gold_structural_share"] == 0.0
    assert 0.0 < _page(report, "DocOne", 2)["gold_structural_share"] < 1.0


def test_nested_markers_do_not_leak_commentary():
    """Notes quote other markers inside themselves, so the span scanner nests.

    A non-nesting regex ends a `[NOTE: ... [FIGURE] ... ]` at the wrong bracket
    and leaks the remainder of the commentary into the compared text.
    """
    text = "before [NOTE: see the [FIGURE] on the facing page] after"
    assert fid.strip_markup(text).text.split() == ["before", "after"]


# --- brackets that are not markup ---------------------------------------------
#
# The gold convention uses `[` for markers, but the pages themselves print
# brackets too, and notes talk *about* brackets. Each case below was found by
# parsing all 761 gold pages and asking which notes failed to close; together
# they moved one document from 0.767 to 0.927 median coverage, and brought the
# structural tag counts to exactly the 348 / 76 / 65 the gold set documents.


def test_a_note_quoting_a_bracket_character_still_closes():
    """`[NOTE: ... "[" is my marker.]` — a mention of the character itself.

    A scanner counting every `[` as a level never closes this note, so the
    whole of it leaks into the compared page as though it had been printed
    there. That is what made one document post 0.767 coverage against 0.998
    recall: gold apparently holding text no extractor could find, because the
    page does not contain it. 13 of that document's 17 pages were affected.
    """
    text = 'printed [NOTE: everywhere else on this page "[" is my marker.] words'
    assert fid.strip_markup(text).text.split() == ["printed", "words"]


def test_a_bracket_pair_inside_a_note_does_not_close_it_early():
    """`[sic]` and `[21]` are quoted from the page; their `]` is not the note's.

    Treating it as the note's own ends the note early, and every marker after
    it is then read at the wrong nesting level — which is how a `[FIGURE]` two
    paragraphs later came to be counted as an unclosed block.
    """
    text = 'printed [NOTE: the bracketed "[sic]" is the translator\'s own] words'
    assert fid.strip_markup(text).text.split() == ["printed", "words"]
    text2 = "printed [NOTE: continuing from the citation [21]. on page 9] words"
    assert fid.strip_markup(text2).text.split() == ["printed", "words"]


def test_an_unterminated_bracket_does_not_swallow_the_marker_after_it():
    """A transcription typo — `[continued opposite` with no `]`.

    Consuming it takes the `[/FIGURE]` that follows with it, so the figure
    block never ends and the rest of the page is scored as plate lettering.
    """
    text = "[FIGURE]\nFig. 1. A caption.\n[continued opposite\n[/FIGURE]\nbody text"
    stripped = fid.strip_markup(text)
    # The figure block closed, so the body text is outside it.
    assert "body text" in stripped.text
    assert "body" not in stripped.structural
    # The tag itself never reaches the compared text.
    assert "/FIGURE" not in stripped.text


def test_unterminated_brackets_are_counted_not_absorbed():
    """A gold-integrity signal: only inspection separates a typo from a note
    mentioning the bracket character, so the count is surfaced per page."""
    assert fid.unterminated_brackets("[FIGURE]\nclean\n[/FIGURE]") == 0
    assert fid.unterminated_brackets("text [continued opposite\n[/FIGURE]") == 1
    rec = fid.score_page("[PAGE 1]\n[continued opposite\n" + "word " * 40, "x " * 40)
    assert rec["unterminated_brackets"] == 1


def test_only_known_keywords_open_a_marker():
    """The vocabulary is what separates markup from a printed bracket."""
    assert fid._marker_at("[FIGURE]", 0)
    assert fid._marker_at("[NOTE: something]", 0)
    assert fid._marker_at("[?reading:gastrozooids]", 0)
    assert fid._marker_at("[/PLATE]", 0)
    assert not fid._marker_at("[21]", 0)          # a printed citation
    assert not fid._marker_at("[ind. m^-3]", 0)   # a printed unit


# --- tokenisation and segmentation --------------------------------------------


def test_tokens_are_unicode_aware():
    """An ASCII-only filter deletes the Cyrillic, Greek and CJK documents'
    actual content and leaves stray Latin species names to compare — which
    yields a confident-looking score computed from a few percent of the page.
    """
    assert fid.tokens("Сифонофоры представляют") == ["сифонофоры", "представляют"]
    # Unspaced scripts are compared character by character; whitespace
    # tokenisation would make one enormous token per line.
    assert fid.tokens("管水母") == ["管", "水", "母"]
    # Case, punctuation and diacritics are normalised away: they are typography,
    # not content, and OCR mangles all three constantly.
    assert fid.tokens("Physalia, arethusa!") == fid.tokens("physalia arethusa")
    assert fid.tokens("Mańko") == fid.tokens("Manko")


def test_dominant_script_is_measured_per_page_not_per_document():
    """At least one document in the gold set is an English translation followed
    by the vertical Japanese original, so a document-level label mis-segments
    both halves.
    """
    assert fid.dominant_script("the siphonophore floats") == "latin"
    assert fid.dominant_script("Сифонофоры представляют собой") == "cyrillic"
    assert fid.dominant_script("管水母は海に浮かぶ") == "cjk"
    # A Latin page quoting one Greek word is still a Latin page.
    assert fid.dominant_script(
        "the nematocyst, from Greek νῆμα, is a stinging organelle") == "latin"


def test_segments_report_the_axes_the_gold_set_was_built_to_span(report):
    """A corpus-wide mean over 13 languages and five centuries cannot tell a
    pipeline that handles born-digital PDFs well and Fraktur not at all from one
    that is mediocre everywhere — and those call for different work.
    """
    assert set(report["segments"]) == {"script", "era", "file_type"}
    assert set(report["segments"]["script"]) == {"latin", "cyrillic"}
    assert set(report["segments"]["era"]) == {"1800-1899", "2000-"}
    assert set(report["segments"]["file_type"]) == {"scanned", "born_digital"}
    latin = report["segments"]["script"]["latin"]
    assert latin["pages"] == 3 and latin["median_coverage"] == 1.0


def test_era_buckets_split_on_typography_not_round_numbers():
    assert fid.era_bucket(1594) == "pre-1800"
    assert fid.era_bucket(1799) == "pre-1800"      # long-s and Fraktur side
    assert fid.era_bucket(1800) == "1800-1899"
    assert fid.era_bucket(1965) == "1950-1999"
    assert fid.era_bucket(2026) == "2000-"
    assert fid.era_bucket(None) == "unknown"


def test_blank_gold_page_is_counted_but_not_scored():
    """There is nothing to recover from a blank leaf, so scoring it would
    reward the extractor for the page being empty."""
    rec = fid.score_page("[PAGE 7]\n\n[BLANK PAGE]\n", "")
    assert rec["status"] == "gold_blank"
    assert rec["coverage"] is None


def test_ocr_damage_is_told_apart_from_genuinely_novel_text():
    """`excess_novel` is what makes a low `recall` actionable.

    Long-s typography guarantees a mismatch on eighteenth-century pages — an
    extractor reads "ſ" as "f", so nearly every word differs while the content
    is right. Those near-misses must not read as the pipeline inventing text.
    """
    gold = "Vesica integra brachia basi ramosa aequalia membrana tenuissima"
    damaged = "Vdica integra bnchia bafi ramofa aequllia membrana tenuiffima"
    invented = ("Vesica integra brachia basi ramosa aequalia membrana tenuissima "
                "helicopter submarine telephone")
    assert fid.score_page(gold, damaged)["excess_novel"] < 0.5
    assert fid.score_page(gold, invented)["excess_novel"] > 0.9


def test_character_salad_is_told_apart_from_real_words():
    """`excess_wordlike` separates the two causes of low recall.

    Text hallucinated out of image texture — a marbled endpaper reads as
    `¿rit 1H äL «•f ^H` — is not the pipeline having found something.
    """
    gold = "the pneumatophore is inflated with gas and floats upon the surface"
    salad = gold + " rit 1h al f h xzq bnk"
    words = gold + " nematocyst gastrozooid siphosome nectophore tentacle"
    assert fid.score_page(gold, salad)["excess_wordlike"] < 0.5
    assert fid.score_page(gold, words)["excess_wordlike"] > 0.9


def test_median_ignores_unscored_pages():
    assert fid._median([0.2, 0.4, 0.6]) == 0.4
    assert fid._median([0.2, None, 0.6]) == 0.4      # even count after filtering
    assert fid._median([None, None]) is None


# ---------------------------------------------------------------------------
# keeppages and the gold set use different page coordinates (#188)
# ---------------------------------------------------------------------------
#
# The gold transcriptions were made over the *whole* PDF — `Beklemishev1969`'s
# gold page 1 is an ownership endpaper carrying a collector's stamp,
# `Kawamura1911a`'s is twelve pages of English typescript bound in front of
# the Japanese original. A `keeppages` document's extracted page numbers are
# positions in the subset. 20 of the 35 gold documents carry a selection, so
# binding by raw page number shifts every one of them and scores near zero —
# a numbering artifact indistinguishable from a catastrophic regression.

def test_no_selection_leaves_page_numbers_alone():
    """The ordinary case must be untouched."""
    from tools.qc.fidelity import keeppages_map, rebase_gold_pages

    assert keeppages_map({}) == {}
    assert keeppages_map({"file_type": "scanned"}) == {}
    gold = {1: "a", 2: "b", 3: "c"}
    assert rebase_gold_pages(gold, {}) is gold


def test_the_recorded_selection_inverts_to_a_gold_to_subset_map():
    from tools.qc.fidelity import keeppages_map

    scan = {"keeppages_selected": [2, 3, 4, 5]}
    assert keeppages_map(scan) == {2: 1, 3: 2, 4: 3, 5: 4}


def test_gold_pages_are_restated_in_subset_coordinates():
    from tools.qc.fidelity import rebase_gold_pages

    # Beklemishev1969's shape: keeppages = {2--45}, so gold page 2 is the
    # extraction's page 1.
    gold = {1: "endpaper", 2: "title", 3: "body"}
    assert rebase_gold_pages(gold, {"keeppages_selected": [2, 3]}) == {1: "title", 2: "body"}


def test_excluded_pages_are_dropped_not_scored_as_misses():
    """Counting them would penalise the selection for doing its job.

    Corpus was *told* the ownership endpaper is not the paper. Scoring it as
    a page whose text corpus failed to recover would make `keeppages` look
    like a regression on exactly the documents it helps most.
    """
    from tools.qc.fidelity import rebase_gold_pages

    gold = {1: "endpaper", 2: "body", 3: "body", 9: "appended translation"}
    out = rebase_gold_pages(gold, {"keeppages_selected": [2, 3]})
    assert set(out.values()) == {"body"}
    assert len(out) == 2


def test_a_discontinuous_selection_maps_correctly():
    """Eschscholtz1825's shape: keeppages = {2--9,11}."""
    from tools.qc.fidelity import keeppages_map

    assert keeppages_map({"keeppages_selected": [2, 3, 4, 5, 6, 7, 8, 9, 11]})[11] == 9
