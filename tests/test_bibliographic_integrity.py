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


@pytest.mark.parametrize("reverse", [False, True])
def test_seven_book_parts_and_duplicate_scans_retain_identity(tmp_path, monkeypatch, reverse):
    from bib.authority import _resolve_reference
    from mcpsrv.tools.bibliography import get_citation_graph
    entries = {f"part{part}": dict(CURATED, _key="Book" + part, title="A book. Volume " + part,
                                    doi="10.5962/bhl.title.10031", volume=part)
               for part in ("1", "2", "3", "4", "5", "8", "6-7 atlas")}
    entries["duplicate"] = entries["part1"]
    if reverse:
        entries = dict(reversed(list(entries.items())))
    conn, db = build(tmp_path, entries)
    ids = {sha: find_work(conn, sha) for sha in entries}
    assert len(set(ids.values())) == 7
    assert ids["duplicate"] == ids["part1"]
    ba = BiblioAuthority(db)
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(biblio_db=ba))
    for sha, entry in entries.items():
        item = format_citations(paper_hashes=[sha])["citations"][0]
        assert item["fields"]["title"] == entry["title"]
        assert item["shared_identifier"] == entry["doi"]
        assert item["corpus_hash"] == sha
        graph = get_citation_graph(paper_hash=sha)
        assert graph["root"]["work_id"] == ids[sha]
        assert graph["root"]["title"] == entry["title"]
    # A specific cited part can be resolved; a DOI-only citation cannot pick one.
    ref = {"doi": entries["part5"]["doi"], "title": entries["part5"]["title"], "year": 2018,
           "authors": ["Mańko M K", "Pugh P R"]}
    assert _resolve_reference(conn, ref)[0] == ids["part5"]
    vague = _resolve_reference(conn, {"doi": ref["doi"]})
    assert vague[0] == ref["doi"]
    assert vague[1] == "shared_doi_unresolved_part"
    assert _resolve_reference(conn, {"doi": ref["doi"]})[1] == "shared_doi_unresolved_part"
    conn.commit()
    assert phase1_corpus_papers(conn, tmp_path) == 0
    assert {sha: find_work(conn, sha) for sha in entries} == ids


def test_same_short_key_distinguishes_moore_parts_and_full_title_suffixes(tmp_path):
    from bib.authority import _resolve_reference
    # Exact reported titles/journals; extra author slots are synthetic.
    entries = {
        "part2": {"_key": "Moore1953", "author": "Moore, Hilary B.", "year": "1953",
                  "title": "Plankton of the Florida Current",
                  "journal": "II. Siphonophora. Bulletin of Marine Science"},
        "part3": {"_key": "Moore_etal1953", "author": "Moore, Hilary B. and Second, A. and Third, R. and Fourth, W.",
                  "year": "1953", "title": "Plankton of the Florida Current",
                  "journal": "III. The control of the vertical distribution of zooplankton in the daytime by light and temperature. Bulletin of Marine Science"},
        "long1": dict(CURATED, doi="", title="A very long identical prefix of a title that exceeds forty characters: first distinct study"),
        "long2": dict(CURATED, doi="", title="A very long identical prefix of a title that exceeds forty characters: second distinct study"),
    }
    conn, _ = build(tmp_path, entries)
    assert len({find_work(conn, sha) for sha in entries}) == 4
    assert len(authors(conn, find_work(conn, "part2"))) == 1
    assert len(authors(conn, find_work(conn, "part3"))) == 4
    for sha in ("part2", "part3"):
        meta = bib_entry_to_metadata(entries[sha], "")
        ref = dict(meta, authors=[a["surname"] for a in meta["authors"]])
        assert _resolve_reference(conn, ref)[0] == find_work(conn, sha)


@pytest.mark.parametrize("reverse", [False, True])
def test_new_shared_part_rekeys_existing_member_like_clean_build(tmp_path, reverse):
    from bib.authority import phase2_references
    from bib.documents import work_map
    from tools.qc.index_reference import index_snapshot
    first = dict(CURATED, title="A book. Volume 1", volume="1")
    second = dict(CURATED, title="A book. Volume 5", volume="5")
    if reverse:
        first, second = second, first
    incremental = tmp_path / "incremental"
    fresh = tmp_path / "fresh"
    conn, _ = build(incremental, {"one": first})
    old = find_work(conn, "one")
    folder = incremental / "documents/five"
    folder.mkdir()
    (folder / "metadata.json").write_text(json.dumps(bib_entry_to_metadata(second, "five.pdf")))
    phase1_corpus_papers(conn, incremental)
    fresh_conn, _ = build(fresh, {"five": second, "one": first})
    assert work_map(conn) == work_map(fresh_conn)
    assert find_work(conn, "one") != old
    assert conn.execute("SELECT in_corpus,bib_imported_at FROM works WHERE work_id=?", (old,)).fetchone() == (0, None)
    assert conn.execute("SELECT COUNT(*) FROM work_identity_decisions WHERE reason='rematerialized_document_identity'").fetchone()[0] == 1
    phase2_references(conn, incremental)
    phase2_references(fresh_conn, fresh)
    # Compare every current bibliography table, including curated sources.
    # Build history/audit decisions intentionally are not current content.
    conn.close()
    fresh_conn.close()
    (incremental / "biblio.sqlite").rename(incremental / "biblio_authority.sqlite")
    (fresh / "biblio.sqlite").rename(fresh / "biblio_authority.sqlite")
    assert index_snapshot(incremental)["bibliography"] == index_snapshot(fresh)["bibliography"]


def test_curated_chun_identity_rejects_popular_incompatible_title(tmp_path):
    from bib.reconcile import reconcile
    title = "Über den Excretionsporus an der Pneumatophore von Physophora"
    entry = {"_key": "Chun1898b", "title": title, "author": "Chun, Carl", "year": "1898", "journal": "Zoologischer Anzeiger"}
    conn, _ = build(tmp_path, {"chun": entry})
    folder = tmp_path / "documents/chun"
    # A long opening text may mention another paper; title-page similarity is
    # not permission to substitute that paper for the curated identity.
    wrong = "Die Ctenophoren der Plankton Expedition"
    (folder / "text.json").write_text(json.dumps({"text": wrong * 20}))
    p = folder / "metadata.json"
    meta = json.loads(p.read_text())
    meta["filename"] = "Chun1898b.pdf"
    p.write_text(json.dumps(meta))
    phase1_corpus_papers(conn, tmp_path)
    actual = find_work(conn, "chun")
    insert_work(conn, "wrong-chun", "corpus_key", wrong, 1898, "", "", None, False, "cited_reference")
    insert_authors(conn, "wrong-chun", [("Chun", "C.")])
    result = reconcile(conn, tmp_path)
    assert result["identity_conflict"] == 1
    assert find_work(conn, "chun") == actual
    assert conn.execute("SELECT title FROM works WHERE work_id=?", (actual,)).fetchone()[0] == title
    decision = conn.execute("SELECT evidence_json FROM work_identity_decisions WHERE reason='rejected_curated_identity_conflict'").fetchone()
    assert "authoritative_title_disagrees" in json.loads(decision[0])["reasons"]


def test_legacy_wrong_membership_is_rederived_without_destroying_observations(tmp_path):
    from bib.authority import phase2_references
    conn, _ = build(tmp_path, {"paper": CURATED})
    original = find_work(conn, "paper")
    folder = tmp_path / "documents/paper"
    (folder / "references.json").write_text(json.dumps({"references": [{"title": "An independent observation", "authors": ["Other A"], "year": 1910, "raw": "Raw source bytes", "xml_id": "b0"}]}))
    phase2_references(conn, tmp_path)
    evidence = list(conn.execute("SELECT * FROM reference_observations"))
    target = ghost(conn)
    merge_phase1_into_ghost(conn, original, target, "paper")
    # Simulate a pre-v1.5 membership receipt, which cannot assert the new
    # identity policy even when metadata bytes have not changed.
    conn.execute("UPDATE work_documents SET source_sha256='old-producer'")
    phase1_corpus_papers(conn, tmp_path)
    phase2_references(conn, tmp_path)
    assert find_work(conn, "paper") == original
    assert list(conn.execute("SELECT * FROM reference_observations")) == evidence
    assert conn.execute("SELECT COUNT(*) FROM work_reconciliation_decisions").fetchone()[0] == 1


@pytest.mark.parametrize("curated", [False, True])
def test_reconcile_still_accepts_uncurated_repairs_and_supported_curated_duplicates(tmp_path, curated):
    from bib.reconcile import reconcile
    target_title = "A synopsis of the marine colonial animals"
    entry = dict(CURATED, title=target_title if curated else "Table (continued)")
    conn, _ = build(tmp_path, {"paper": entry})
    folder = tmp_path / "documents/paper"
    metadata = bib_entry_to_metadata(entry, "Manko2018.pdf")
    if not curated:
        metadata.pop("extraction_method")
        metadata.pop("bib_key")
    (folder / "metadata.json").write_text(json.dumps(metadata))
    (folder / "text.json").write_text(json.dumps({"text": target_title}))
    phase1_corpus_papers(conn, tmp_path)
    insert_work(conn, "supported-ghost", "corpus_key", target_title, 2018, "Zootaxa", "", None, False, "cited_reference")
    insert_authors(conn, "supported-ghost", [("Mańko", "M.K.")])
    result = reconcile(conn, tmp_path)
    assert result["matched"] == 1
    assert find_work(conn, "paper") == "supported-ghost"


def test_curated_volume_number_in_title_is_identity_evidence(tmp_path):
    from bib.reconcile import _curated_candidate_conflict
    conn, _ = build(tmp_path, {"paper": dict(CURATED, title="A long book about colonial animals. Volume 1")})
    wid = find_work(conn, "paper")
    insert_work(conn, "wrong-volume", "corpus_key", "A long book about colonial animals. Volume 5", 2018, "Zootaxa", "", None, False, "cited_reference")
    conflict = _curated_candidate_conflict(conn, wid, "wrong-volume")
    assert "authoritative_title_part_disagrees" in conflict["reasons"]
