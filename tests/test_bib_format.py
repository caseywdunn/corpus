"""Tests for citation formatters in bib.format (#79).

Pins the contract that downstream code (format_citation MCP tool)
relies on: each formatter takes a fields dict and returns
``{formatted, inline}``. Author-list shape rules — APA-style
surname + initials, Oxford comma + ``&`` before final — are tested
explicitly because the MCP tool's whole value proposition is that
its output is what the LLM emits verbatim. Drift here is a
hallucination vector.
"""
from __future__ import annotations

import pytest

from bib.format import (
    SUPPORTED_STYLES,
    _format_author_list,
    _initials,
    _inline_author_phrase,
    format_author_year,
    format_citation,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("forename,expected", [
    ("Anne Marie", "A. M."),
    ("A. K.", "A. K."),
    ("Jane", "J."),
    ("",     ""),
    (None,   ""),
    ("jean-luc", "J."),  # only the first letter of the first token
])
def test_initials(forename, expected):
    assert _initials(forename) == expected


def test_format_author_list_one():
    assert _format_author_list([
        {"surname": "Totton", "forename": "A. K."},
    ]) == "Totton, A. K."


def test_format_author_list_two_uses_oxford_ampersand():
    """APA convention: two authors get '<a>, & <b>' (Oxford comma
    before the ampersand). Drift here changes every two-author
    citation in the corpus."""
    assert _format_author_list([
        {"surname": "Smith", "forename": "Jane"},
        {"surname": "Jones", "forename": "Karim"},
    ]) == "Smith, J., & Jones, K."


def test_format_author_list_three():
    assert _format_author_list([
        {"surname": "Smith", "forename": "Jane"},
        {"surname": "Jones", "forename": "Karim"},
        {"surname": "Lee",   "forename": "Quentin"},
    ]) == "Smith, J., Jones, K., & Lee, Q."


def test_format_author_list_surname_only():
    """An author with no forename surfaces as just the surname,
    with no trailing comma."""
    assert _format_author_list([
        {"surname": "Totton", "forename": ""},
    ]) == "Totton"


def test_format_author_list_skips_blank_surnames():
    assert _format_author_list([
        {"surname": "", "forename": "Anonymous"},
        {"surname": "Smith", "forename": "Jane"},
    ]) == "Smith, J."


def test_format_author_list_empty():
    assert _format_author_list([]) == ""
    assert _format_author_list(None) == ""


@pytest.mark.parametrize("authors,expected", [
    ([{"surname": "Smith"}], "Smith"),
    ([{"surname": "Smith"}, {"surname": "Jones"}], "Smith & Jones"),
    ([{"surname": "Smith"}, {"surname": "Jones"}, {"surname": "Lee"}],
     "Smith et al."),
    ([{"surname": "Smith"}, {"surname": "Jones"}, {"surname": "Lee"}, {"surname": "Wu"}],
     "Smith et al."),
    ([], ""),
])
def test_inline_author_phrase(authors, expected):
    assert _inline_author_phrase(authors) == expected


# ---------------------------------------------------------------------------
# format_author_year — end-to-end shapes
# ---------------------------------------------------------------------------


def test_full_record_renders_apa_style():
    out = format_author_year({
        "authors": [{"surname": "Totton", "forename": "A. K."}],
        "year": 1965,
        "title": "A synopsis of the Siphonophora",
        "journal": "Bulletin of the British Museum (Natural History)",
        "volume": "7",
        "pages": "1-200",
        "doi": "10.5962/bhl.part.20269",
    })
    assert out["formatted"] == (
        "Totton, A. K. (1965). A synopsis of the Siphonophora. "
        "*Bulletin of the British Museum (Natural History)*, *7*, 1-200. "
        "https://doi.org/10.5962/bhl.part.20269"
    )
    assert out["inline"] == "(Totton, 1965)"


def test_two_author_inline_uses_ampersand():
    out = format_author_year({
        "authors": [
            {"surname": "Smith", "forename": "Jane"},
            {"surname": "Jones", "forename": "Karim"},
        ],
        "year": 2010,
        "title": "Some paper",
    })
    assert out["inline"] == "(Smith & Jones, 2010)"


def test_three_author_inline_uses_et_al():
    out = format_author_year({
        "authors": [
            {"surname": "Smith"}, {"surname": "Jones"}, {"surname": "Lee"},
        ],
        "year": 2010,
        "title": "Some paper",
    })
    assert out["inline"] == "(Smith et al., 2010)"


def test_missing_year_uses_nd():
    out = format_author_year({
        "authors": [{"surname": "Smith", "forename": "Jane"}],
        "title": "Some paper",
    })
    assert out["formatted"].startswith("Smith, J. (n.d.).")
    assert out["inline"] == "(Smith, n.d.)"


def test_missing_doi_journal_volume_pages_omits_gracefully():
    """Sparse records (e.g. a Grobid-extracted reference with only
    author + year + title) shouldn't surface empty fields or
    placeholder text."""
    out = format_author_year({
        "authors": [{"surname": "Totton", "forename": "A. K."}],
        "year": 1965,
        "title": "A synopsis of the Siphonophora",
    })
    assert out["formatted"] == (
        "Totton, A. K. (1965). A synopsis of the Siphonophora."
    )
    assert "*" not in out["formatted"]            # no journal italics
    assert "doi.org" not in out["formatted"]


def test_no_authors_year_only_inline():
    out = format_author_year({
        "year": 1900,
        "title": "Anonymous treatise",
    })
    assert out["inline"] == "(1900)"
    assert out["formatted"].startswith("(1900). Anonymous treatise.")


def test_title_terminal_period_not_doubled():
    """A title that already ends with a period shouldn't get an
    extra one appended."""
    out = format_author_year({
        "authors": [{"surname": "Smith"}],
        "year": 2010,
        "title": "Some paper.",
    })
    assert "Some paper.." not in out["formatted"]
    assert "Some paper. " in out["formatted"] or out["formatted"].endswith("Some paper.")


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


def test_format_citation_dispatches_to_author_year():
    fields = {"authors": [{"surname": "Smith"}], "year": 2010}
    assert format_citation(fields, style="author-year") == format_author_year(fields)


def test_format_citation_rejects_unknown_style():
    with pytest.raises(ValueError, match="unknown citation style"):
        format_citation({}, style="vancouver")


def test_supported_styles_set_drives_dispatcher():
    """SUPPORTED_STYLES is the public list; the dispatcher must
    handle every entry. Today that's just author-year; this test
    guards against a future entry being added to the set without a
    corresponding dispatch branch."""
    for style in SUPPORTED_STYLES:
        result = format_citation({"year": 2010}, style=style)
        assert "formatted" in result and "inline" in result
