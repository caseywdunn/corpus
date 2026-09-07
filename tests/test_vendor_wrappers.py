"""Vendor wrapper detection (#216).

`_VENDOR_BOILERPLATE` drives one decision: a text layer that clears the volume
gate but contains nothing except a vendor banner belongs to a document whose
real content is raster underneath, so it is re-routed to OCR.

The list held three markers and missed the two most common wrappers in a
scanned library. Measured in one full reference-library audit, pages
1-2:

    20  JSTOR cover        'Your use of the JSTOR archive', 'links.jstor.org'
     6  ResearchGate
     6  BHL
     5  blank notice
     1  ProQuest
     1  Google Books       incl. the German 'Über dieses Buch'
    34  documents with any wrapper

WRAPPERS ARE NOT IMPRINTS, and the separation is the load-bearing part.
An imprint is publisher branding on the paper's *own* pages — a ScienceDirect
header, a Springer footer. In the same library **373** documents carry one,
against 34 carrying a wrapper. Here a false match only costs a wasteful
re-OCR; for #188, where a marker is evidence to *drop* a page, matching a
Springer footer would delete the article's first page.
"""
from __future__ import annotations

import pytest

from pipeline.scan import _VENDOR_BOILERPLATE


@pytest.mark.parametrize("marker", [
    "links.jstor.org",
    "Your use of the JSTOR archive",
    "books.google.com",
    "digitized by Google",
    "Über dieses Buch",
    "researchgate.net",
])
def test_the_common_wrappers_are_present(marker):
    assert marker in _VENDOR_BOILERPLATE


@pytest.mark.parametrize("marker", [
    "ProQuest ebrary",
    "biodiversitylibrary.org",
    "This page intentionally left blank",
])
def test_the_original_markers_are_kept(marker):
    assert marker in _VENDOR_BOILERPLATE


@pytest.mark.parametrize("imprint", [
    "ScienceDirect",
    "Springer",
    "Wiley",
    "This content downloaded",
    "Downloaded by",
    "Downloaded from",
])
def test_publisher_imprints_are_not_in_the_wrapper_list(imprint):
    """The distinction this issue exists to encode.

    373 documents in the reference library carry an imprint against 34 a
    wrapper. Adding one here would re-OCR a tenth of the library for nothing —
    and would be actively destructive if #188 ever shares this list, since it
    would mark the article's own first page as droppable.
    """
    assert not any(imprint.lower() in m.lower() for m in _VENDOR_BOILERPLATE)


def test_a_jstor_cover_matches_the_way_detect_scan_type_checks():
    """The real page-1 text of a JSTOR cover sheet, abridged."""
    cover = (
        "Analysis of Locomotion in a Siphonophore Colony / G. 0. Mackie\n"
        "Proceedings of the Royal Society of London. Series B, Volume 159, "
        "Issue 975 (Jan. 14, 1964), 366-391.\n"
        "Stable URL: http://links.jstor.org/sici?sici=0080-4649\n"
        "Your use of the JSTOR archive indicates your acceptance of JSTOR's "
        "Terms and Conditions of Use.\n"
    )
    assert any(m in cover for m in _VENDOR_BOILERPLATE)
    assert len(cover) < 5000, "must also clear the char gate to re-route"


def test_a_google_books_german_notice_matches():
    """`Über dieses Buch` is the localised variant on German scans, and is the
    one that actually fires — the English notice sits on a 201-character page
    while this one is 3.7 KB, still inside the gate."""
    assert "Über dieses Buch" in _VENDOR_BOILERPLATE


def test_an_ordinary_paper_does_not_match():
    body = ("The siphonophore floats upon the surface of the sea, its "
            "pneumatophore inflated with gas. Downloaded from a Springer "
            "journal via ScienceDirect.")
    assert not any(m in body for m in _VENDOR_BOILERPLATE), (
        "imprint words in ordinary text must not trip the wrapper check")
