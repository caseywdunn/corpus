"""Authority-string parsing + original-description resolution.

Covers the linker fixes that let held original descriptions link to their
Darwin Core taxa despite awkward `scientificNameAuthorship` strings:

- comma-separated author lists ("Gershwin, Zeidler & Davie, 2010")
- the zoological "X in Y" idiom ("Brandt in Mertens, 1833"), where the
  containing work Y physically holds the description
- multi-letter forename initials ("O.F. Müller") and stray connectors
  ("of Lesson")
- preference for the in-corpus paper over a same-key cited-reference ghost
- disambiguation of two same-author/same-year works by co-author overlap
- a conservative fuzzy fallback for authorship spelling variants
  (Eschscholz→Eschscholtz)
"""
from __future__ import annotations

import time

import sqlite3
import pytest

from bib.authority import (
    create_schema,
    insert_authors,
    insert_work,
    parse_authority,
    resolve_authority_work,
)


# ── parse_authority (pure) ───────────────────────────────────────────

@pytest.mark.parametrize("authority, first, year", [
    ("Eschscholtz, 1829", "Eschscholtz", 1829),
    ("(Huxley, 1859)", "Huxley", 1859),
    ("Quoy & Gaimard, 1833", "Quoy", 1833),
    ("L. Agassiz, 1862", "Agassiz", 1862),
    ("Gershwin, Zeidler & Davie, 2010", "Gershwin", 2010),
    ("Harbison, Matsumoto & Robison, 2001", "Harbison", 2001),
    ("O.F. Müller, 1776", "Müller", 1776),
    ("of Lesson, 1843", "Lesson", 1843),
])
def test_first_surname(authority, first, year):
    surnames, y = parse_authority(authority)
    assert surnames[0] == first
    assert y == year


def test_comma_list_splits_all_authors():
    surnames, _ = parse_authority("Gershwin, Zeidler & Davie, 2010")
    assert surnames == ["Gershwin", "Zeidler", "Davie"]


def test_in_idiom_puts_containing_work_first():
    # The description of Idya mertensii was published by Brandt *inside*
    # Mertens 1833 — Mertens holds the PDF, so it must be tried first.
    surnames, year = parse_authority("Brandt in Mertens, 1833")
    assert surnames[0] == "Mertens"
    assert "Brandt" in surnames
    assert year == 1833


def test_in_idiom_with_initials():
    surnames, year = parse_authority("(A. Agassiz in L. Agassiz, 1860)")
    assert surnames == ["Agassiz"]
    assert year == 1860


def test_particle_surname_preserved():
    surnames, _ = parse_authority("Lens & van Riemsdijk, 1908")
    assert surnames == ["Lens", "van Riemsdijk"]


def test_unparseable_returns_none():
    assert parse_authority("no year here") is None
    assert parse_authority("") is None


# ── resolve_authority_work (DB-backed) ───────────────────────────────

def _work(conn, work_id, year, surnames, *, in_corpus):
    insert_work(conn, work_id, "corpus_key", title=f"t {work_id}", year=year,
                journal="J", doi="", corpus_hash=(work_id if in_corpus else None),
                in_corpus=in_corpus, source="corpus_paper" if in_corpus
                else "cited_reference")
    insert_authors(conn, work_id, [(s, "") for s in surnames])


@pytest.fixture()
def conn():
    c = sqlite3.connect(":memory:")
    create_schema(c)
    return c


def test_prefers_in_corpus_over_ghost(conn):
    # Same (surname, year) as an empty cited-reference ghost — must return the
    # held paper, not the stub.
    _work(conn, "ghost", 1829, ["Eschscholtz"], in_corpus=False)
    _work(conn, "held", 1829, ["Eschscholtz"], in_corpus=True)
    wid, conf = resolve_authority_work(conn, ["Eschscholtz"], 1829)
    assert wid == "held"


def test_coauthor_disambiguates_same_author_year(conn):
    # Two Harbison 2001 papers; the authority's later surnames pick the right.
    _work(conn, "invasion", 2001, ["Harbison"], in_corpus=True)
    _work(conn, "lampocteis", 2001,
          ["Harbison", "Matsumoto", "Robison"], in_corpus=True)
    wid, conf = resolve_authority_work(
        conn, ["Harbison", "Matsumoto", "Robison"], 2001)
    assert wid == "lampocteis"
    assert conf >= 0.85


def test_in_idiom_resolves_to_containing_work(conn):
    _work(conn, "mertens1833", 1833, ["Mertens"], in_corpus=True)
    # surnames as parse_authority would return them: containing first.
    wid, _ = resolve_authority_work(conn, ["Mertens", "Brandt"], 1833)
    assert wid == "mertens1833"


def test_fuzzy_spelling_variant(conn):
    # DwC typo 'Eschscholz' (missing 't') should still find the held work.
    _work(conn, "held", 1829, ["Eschscholtz"], in_corpus=True)
    wid, conf = resolve_authority_work(conn, ["Eschscholz"], 1829)
    assert wid == "held"
    assert conf <= 0.6


def test_no_match_returns_none(conn):
    _work(conn, "other", 1900, ["Nobody"], in_corpus=True)
    wid, conf = resolve_authority_work(conn, ["Ghostwriter"], 1850)
    assert wid is None
    assert conf == 0.0
