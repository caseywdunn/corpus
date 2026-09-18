"""Curated source authority and publication locator regressions (#296, #301)."""
from __future__ import annotations

import json
import sqlite3
from types import SimpleNamespace

import pytest

from bib.authority import create_schema, insert_authors, insert_work, phase1_corpus_papers
from bib.documents import document_metadata, find_work
from bib.export import export_bibtex
from bib.importer import apply_entry, import_bibtex
from bib.parser import bib_entry_to_metadata, parse_bibtex
from bib.reconcile import merge_phase1_into_ghost
from mcpsrv import app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import format_citations


CURATED = {"_key": "MankoPugh2018", "title": "A curated species description",
           "author": "Mańko, Maciej K. and Pugh, Philip R.", "year": "2018",
           "doi": "10.11646/zootaxa.4441.2.7", "journal": "Zootaxa"}


def build(tmp_path, entries):
    for sha, entry in entries.items():
        folder = tmp_path / "documents" / sha
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "metadata.json").write_text(json.dumps(bib_entry_to_metadata(entry, sha + ".pdf")))
    db = tmp_path / "biblio.sqlite"
    conn = sqlite3.connect(db)
    create_schema(conn)
    phase1_corpus_papers(conn, tmp_path)
    return conn, db


def authors(conn, wid):
    return list(conn.execute("SELECT surname,forename FROM work_authors WHERE work_id=? ORDER BY position", (wid,)))


def ghost(conn, wid="ghost"):
    insert_work(conn, wid, "corpus_key", CURATED["title"], 2018, "Zootaxa", "", None, False, "cited_reference")
    insert_authors(conn, wid, [("Pugh", "Incorrect"), ("Manko", "Wrong"), ("Extra", "Parsed")])
    return wid


def test_build_with_bib_preserves_origin_order_unicode_refresh_and_formatter(tmp_path, monkeypatch):
    conn, db = build(tmp_path, {"paper": CURATED})
    wid = find_work(conn, "paper")
    assert document_metadata(conn, "paper")["bib_key"] == "MankoPugh2018"
    assert document_metadata(conn, "paper")["extraction_method"] == "bib"
    expected = [("Mańko", "Maciej K."), ("Pugh", "Philip R.")]
    assert authors(conn, wid) == expected
    before = list(conn.iterdump())
    assert phase1_corpus_papers(conn, tmp_path) == 0
    assert list(conn.iterdump()) == before
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(biblio_db=BiblioAuthority(db)))
    result = format_citations(paper_hashes=["paper"])["citations"][0]
    assert result["provenance"] == "bib"
    assert result["bib_key"] == "MankoPugh2018"
    assert [(a["surname"], a["forename"]) for a in result["fields"]["authors"]] == expected
    # A changed metadata artifact refreshes the current source, not just its stamp.
    revised = dict(CURATED, author="Pugh, Philip R. and Mańko, Maciej K.")
    p = tmp_path / "documents/paper/metadata.json"
    p.write_text(json.dumps(bib_entry_to_metadata(revised, "paper.pdf")))
    phase1_corpus_papers(conn, tmp_path)
    assert authors(conn, wid) == expected[::-1]


@pytest.mark.parametrize("import_first", [True, False])
def test_import_merge_order_preserves_entire_curated_author_list(tmp_path, import_first):
    entry = dict(CURATED)
    conn, _ = build(tmp_path, {"paper": entry})
    source = find_work(conn, "paper")
    target = ghost(conn)
    curated = dict(entry, author="Pugh, Philip R. and Mańko, Maciej K.")
    if import_first:
        apply_entry(conn, source, curated)
    merge_phase1_into_ghost(conn, source, target, "paper")
    if not import_first:
        apply_entry(conn, target, curated)
    assert authors(conn, target) == [("Pugh", "Philip R."), ("Mańko", "Maciej K.")]
    assert apply_entry(conn, target, curated) == 0
    conn.commit()
    # Unchanged phase 1 must retain a valid reconciliation and curated fields.
    phase1_corpus_papers(conn, tmp_path)
    assert find_work(conn, "paper") == target
    assert authors(conn, target) == [("Pugh", "Philip R."), ("Mańko", "Maciej K.")]


@pytest.mark.parametrize("reverse", [False, True])
def test_conflicting_bib_entries_have_deterministic_explicit_policy(tmp_path, reverse):
    conn, db = build(tmp_path, {"paper": CURATED})
    wid = find_work(conn, "paper")
    entries = [dict(CURATED, _key="A", author="Mańko, Maciej K."),
               dict(CURATED, _key="Z", author="Mańko, Different Name and Extra, Author")]
    for entry in reversed(entries) if reverse else entries:
        apply_entry(conn, wid, entry)
    conn.commit()
    work = BiblioAuthority(db).get_work(wid)
    assert authors(conn, wid) == [("Mańko", "Maciej K.")]
    conflict = next(item for item in work["bibliographic_conflicts"] if item["field"] == "authors")
    assert conflict["selected_source"] == "import:A"
    assert conflict["policy"] == "explicit_import_then_source_id"
    assert any(item["value"][0]["forename"] == "Different Name" for item in conflict["sources"])


@pytest.mark.parametrize("fields", [
    {"volume": "324", "number": "5", "pages": "435--449"},  # Church et al. 2015, #301
    {"eid": "e1000505", "volume": "6", "number": "7"},  # synthetic article-number branch
    {"booktitle": "Collected studies", "volume": "2", "chapter": "4", "pages": "39--57", "publisher": "Example Press"},
])
def test_locator_build_export_import_and_unchanged_rebuild(tmp_path, monkeypatch, fields):
    entry = dict(CURATED, **fields, keeppages="2--9")
    conn, db = build(tmp_path, {"paper": entry})
    wid = find_work(conn, "paper")
    ba = BiblioAuthority(db)
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(biblio_db=ba))
    def assert_citation():
        result = format_citations(paper_hashes=["paper"])["citations"][0]
        for key, value in fields.items():
            assert result["fields"][key] == value
        assert "keeppages" not in result["fields"]
        if fields.get("pages"):
            assert fields["pages"].replace("--", "–") in result["formatted"]
        if fields.get("number"):
            assert f"({fields['number']})" in result["formatted"]
        return result
    first = assert_citation()
    exported = tmp_path / "roundtrip.bib"
    exported.write_text(export_bibtex(db))
    parsed = parse_bibtex(exported.read_text())[0]
    assert all(parsed[key] == value for key, value in fields.items())
    assert parsed["keeppages"] == "2--9"
    import_bibtex(db, exported)
    import_bibtex(db, exported)
    phase1_corpus_papers(conn, tmp_path)
    assert assert_citation()["fields"] == first["fields"]
    assert all(ba.get_work(wid)[key] == value for key, value in fields.items())
