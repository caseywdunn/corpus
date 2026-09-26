"""Author-list syntax must agree across retrieval and citation formatting (#310)."""
import sqlite3
from types import SimpleNamespace

import pytest

from bib.authority import create_schema, insert_authors, insert_work
from mcpsrv import app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import format_citations, resolve_reference


@pytest.fixture
def query_index(tmp_path, monkeypatch):
    path = tmp_path / "biblio.sqlite"
    conn = sqlite3.connect(path)
    create_schema(conn)

    def add(wid, title, year, names):
        insert_work(conn, wid, "corpus_key", title, year, None,
                    None, None, False, "cited_reference")
        insert_authors(conn, wid, [(name, "") for name in names])

    add("manko", "Footprints of Atlantification", 2020, ["Mańko", "Pugh"])
    add("totton", "A Synopsis of the Siphonophora", 1965, ["Totton", "Bargmann"])
    # Twelve independent paired-query controls, with ambiguous short forms.
    # These exercise the reported syntax; they are not the unsaved audit calls.
    names = ["Pugh", "Haddock", "Dunn", "Siebert", "Mapstone", "Hosia",
             "Sutherland", "van Soest", "De Haan", "Lo Bianco", "Müller", "D’Urville"]
    for n, name in enumerate(names):
        for variant in ("alpha", "beta"):
            add(f"pair:{n}:{variant}", f"Query control {variant} {n}", 2010,
                [name, "Colleague"])
    # Full compound surnames must win over a distinct surname's suffix.
    add("soest", "A different work", 2010, ["Soest"])
    conn.commit()
    authority = BiblioAuthority(path)
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(biblio_db=authority))
    yield names
    authority.conn.close()
    conn.close()


def _formatted(query):
    return format_citations(queries=[query])["citations"][0]


@pytest.mark.parametrize("query,expected", [
    ("Manko 2020", "manko"),
    ("Manko et al 2020 Footprints of Atlantification", "manko"),
    ("Mańko et al. (2020) Footprints of Atlantification", "manko"),
    ("Man\u0301ko et al. 2020 Footprints of Atlantification", "manko"),
    ("Totton and Bargmann 1965 A Synopsis of the Siphonophora", "totton"),
    ("Totton & Bargmann 1965", "totton"),
    ("Totton, A. K. & Bargmann, H. E., 1965 A Synopsis", "totton"),
])
def test_named_author_forms_agree_across_tools(query_index, query, expected):
    assert resolve_reference(query)["work_id"] == expected
    assert _formatted(query)["work_id"] == expected


def test_twelve_paired_queries_retain_candidates_with_more_information(query_index):
    for n, surname in enumerate(query_index):
        short = f"{surname} 2010"
        expected = {f"pair:{n}:alpha", f"pair:{n}:beta"}
        assert {row["work_id"] for row in resolve_reference(short)["matches"]} == expected
        assert _formatted(short)["code"] == "ambiguous"
        for separator in ("et al.", "and Colleague", "& Colleague"):
            long = f"{surname} {separator} 2010 Query control beta {n}"
            assert resolve_reference(long)["work_id"] == f"pair:{n}:beta"
            assert _formatted(long)["work_id"] == f"pair:{n}:beta"


def test_particles_and_initials_survive_explicit_author_list(query_index):
    query = "R. W. M. van Soest & Colleague 2010 Query control alpha 7"
    assert resolve_reference(query)["work_id"] == "pair:7:alpha"
    assert _formatted(query)["work_id"] == "pair:7:alpha"


def test_explicit_author_and_year_still_accept_title_fragment(query_index):
    assert resolve_reference("Query control alpha 7", author="van Soest", year=2010)["work_id"] == "pair:7:alpha"


def test_true_query_miss_does_not_claim_a_corpus_census(query_index):
    query = "AbsentAuthor et al. 1890 A nonexistent control"
    missing = resolve_reference(query)
    assert missing["not_found"]
    formatted = _formatted(query)
    assert formatted["code"] == "not_found"
    assert "does not establish publication absence" in formatted["error"]
    assert formatted["authors_tried"] == missing["authors_tried"]


@pytest.mark.parametrize("query", ["Footprints of Atlantification", "123 !!! 2020 Title"])
def test_unsupported_parse_is_distinct_from_a_lookup_miss(query_index, query):
    assert resolve_reference(query)["code"] == "invalid_argument"
    assert _formatted(query)["code"] == "invalid_argument"
