"""One work can have many PDFs without losing evidence or sharing permissions."""
import json
import shutil
import sqlite3

from bib import authority
from bib.documents import document_metadata, work_map
from bib.export import export_bibtex
from bib.importer import import_bibtex
from bib.parser import parse_bibtex
from mcpsrv.bundle import _load_skipped_hashes
from mcpsrv.indexes import BiblioAuthority


def paper(root, sha, *, doi="10.1234/shared", **extra):
    hd = root / "documents" / sha
    hd.mkdir(parents=True, exist_ok=True)
    meta = dict(title="Collected studies", year=2000, journal="Journal", doi=doi,
                authors=[{"surname": "Author", "forename": "A"}], **extra)
    (hd / "metadata.json").write_text(json.dumps(meta))
    (hd / "references.json").write_text(json.dumps({"references": [{
        "xml_id": "b0", "title": "An independent source", "year": 1990,
        "authors": ["Other"], "raw": "Other 1990. An independent source."}]}))
    return hd


def database(root):
    path = root / "authority.sqlite"
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys=ON")
    authority.create_schema(conn)
    return conn


def run(conn, root):
    authority.phase1_corpus_papers(conn, root)
    authority.phase2_references(conn, root)
    authority.apply_publishable_derivation(conn)


def logical(conn):
    return {"members": work_map(conn),
            "works": conn.execute("SELECT work_id,title,year,doi,corpus_hash,in_corpus,license,serve FROM works ORDER BY work_id").fetchall(),
            "authors": conn.execute("SELECT work_id,position,surname,forename FROM work_authors ORDER BY work_id,position").fetchall(),
            "mappings": conn.execute("SELECT observation_id,work_id,match_method FROM observation_work ORDER BY observation_id").fetchall(),
            "citations": conn.execute("SELECT citing_work_id,cited_work_id,citing_corpus_hash FROM citations ORDER BY citing_corpus_hash,cited_work_id").fetchall()}


def test_shared_work_keeps_both_citing_documents_and_noop(tmp_path):
    paper(tmp_path, "aaa")
    paper(tmp_path, "bbb")
    conn = database(tmp_path)
    run(conn, tmp_path)
    assert work_map(conn) == {"aaa": "10.1234/shared", "bbb": "10.1234/shared"}
    assert conn.execute("SELECT COUNT(*) FROM observation_work").fetchone()[0] == 2
    before = conn.total_changes
    authority.phase1_corpus_papers(conn, tmp_path)
    assert authority.phase2_references(conn, tmp_path) == (0, 0)
    assert conn.total_changes == before
    conn.close()


def test_incremental_addition_in_reverse_order_matches_clean(tmp_path):
    incremental, clean = tmp_path / "incremental", tmp_path / "clean"
    incremental.mkdir()
    clean.mkdir()
    paper(incremental, "bbb", license="CC-BY-4.0", serve=0)
    a = database(incremental)
    run(a, incremental)
    paper(incremental, "aaa", license="unknown", serve=1)
    run(a, incremental)
    shutil.copytree(incremental / "documents", clean / "documents")
    b = database(clean)
    run(b, clean)
    assert logical(a) == logical(b)
    a.close()
    b.close()


def test_identity_edit_moves_only_one_member_and_matches_clean(tmp_path):
    paper(tmp_path, "aaa")
    paper(tmp_path, "bbb")
    conn = database(tmp_path)
    run(conn, tmp_path)
    paper(tmp_path, "bbb", doi="10.1234/revised")
    run(conn, tmp_path)
    assert work_map(conn) == {"aaa": "10.1234/shared", "bbb": "10.1234/revised"}
    clean = tmp_path / "clean"
    shutil.copytree(tmp_path / "documents", clean / "documents")
    other = database(clean)
    run(other, clean)
    assert logical(conn) == logical(other)
    conn.close()
    other.close()


def test_removed_member_keeps_historical_observations_but_no_live_edges(tmp_path):
    paper(tmp_path, "aaa")
    hd = paper(tmp_path, "bbb")
    conn = database(tmp_path)
    run(conn, tmp_path)
    shutil.rmtree(hd)
    run(conn, tmp_path)
    assert work_map(conn) == {"aaa": "10.1234/shared"}
    assert conn.execute("SELECT COUNT(*) FROM reference_observations").fetchone()[0] == 2
    assert conn.execute("SELECT citing_corpus_hash FROM citations").fetchall() == [("aaa",)]
    conn.close()


def test_permissions_and_directives_are_document_local_on_export_and_import(tmp_path):
    paper(tmp_path, "aaa", license="CC-BY-4.0", serve=0, keeppages="2--8", ocrlang="eng")
    paper(tmp_path, "bbb", license="unknown", serve=1, keeppages="9--20", ocrlang="deu")
    conn = database(tmp_path)
    run(conn, tmp_path)
    path = tmp_path / "authority.sqlite"
    assert _load_skipped_hashes(path) == {"aaa"}
    reader = BiblioAuthority(path)
    assert reader.get_work_by_corpus_hash("bbb")["publishable"] == 0
    assert reader.get_work_by_corpus_hash("aaa")["publishable"] == 1
    assert reader.get_work_by_corpus_hash("bbb")["corpus_hash"] == "bbb"
    reader.conn.close()
    entries = parse_bibtex(export_bibtex(path, documents_dir=tmp_path / "documents"))
    assert len(entries) == 2
    assert len({e["_key"] for e in entries}) == 2
    by_hash = {e["corpus_hash"]: e for e in entries}
    assert by_hash["aaa"]["keeppages"] == "2--8"
    assert by_hash["bbb"]["keeppages"] == "9--20"
    update = tmp_path / "changes.bib"
    update.write_text('@article{edit, corpus_hash={bbb}, ocrlang={fra}, serve={false}}')
    assert import_bibtex(path, update)["matched_corpus_hash"] == 1
    assert document_metadata(conn, "aaa")["ocrlang"] == "eng"
    assert document_metadata(conn, "bbb")["ocrlang"] == "fra"
    assert _load_skipped_hashes(path) == {"aaa", "bbb"}
    conn.close()


def test_removing_a_directive_clears_old_value_even_with_unchanged_mtime(tmp_path):
    hd = paper(tmp_path, "aaa", ocrlang="deu", keeppages="2--5", serve=0)
    conn = database(tmp_path)
    run(conn, tmp_path)
    import os
    stamp = (hd / "metadata.json").stat().st_mtime_ns
    paper(tmp_path, "aaa")
    os.utime(hd / "metadata.json", ns=(stamp, stamp))
    run(conn, tmp_path)
    assert conn.execute("SELECT ocrlang,keeppages,serve FROM works WHERE in_corpus=1").fetchone() == (None, None, 1)
    conn.close()


def test_reconciliation_moves_all_members_without_overwriting_permissions(tmp_path):
    from bib.reconcile import merge_phase1_into_ghost
    paper(tmp_path, "aaa", license="CC-BY-4.0", serve=0)
    paper(tmp_path, "bbb", license="unknown", serve=1)
    conn = database(tmp_path)
    run(conn, tmp_path)
    authority.insert_work(conn, "10.1234/ghost", "doi", "A better identity", 2000, "", "10.1234/ghost", None, False, "cited_reference")
    merge_phase1_into_ghost(conn, "10.1234/shared", "10.1234/ghost", "aaa")
    conn.commit()
    assert work_map(conn) == {"aaa": "10.1234/ghost", "bbb": "10.1234/ghost"}
    assert conn.execute("SELECT corpus_hash FROM works WHERE work_id='10.1234/ghost'").fetchone()[0] == "bbb"
    assert _load_skipped_hashes(tmp_path / "authority.sqlite") == {"aaa"}
    conn.close()


def test_legacy_scalar_membership_migrates_without_reverting_reconciliation(tmp_path):
    paper(tmp_path, "aaa")
    conn = database(tmp_path)
    run(conn, tmp_path)
    conn.execute("DELETE FROM work_documents")
    conn.execute("UPDATE works SET title='Curated title' WHERE in_corpus=1")
    authority.phase1_corpus_papers(conn, tmp_path)
    assert conn.execute("SELECT title FROM works WHERE in_corpus=1").fetchone()[0] == "Curated title"
    assert work_map(conn) == {"aaa": "10.1234/shared"}
    conn.close()


def test_work_only_curation_cannot_change_multiple_documents_implicitly(tmp_path):
    import pytest
    paper(tmp_path, "aaa")
    paper(tmp_path, "bbb")
    conn = database(tmp_path)
    run(conn, tmp_path)
    bib = tmp_path / "ambiguous.bib"
    bib.write_text('@article{x,doi={10.1234/shared},serve={false}}')
    with pytest.raises(ValueError, match="require corpus_hash"):
        import_bibtex(tmp_path / "authority.sqlite", bib)
    assert _load_skipped_hashes(tmp_path / "authority.sqlite") == set()
    conn.close()


def test_work_only_curation_updates_the_only_document(tmp_path):
    paper(tmp_path, "aaa")
    conn = database(tmp_path)
    run(conn, tmp_path)
    bib = tmp_path / "update.bib"
    bib.write_text('@article{x,doi={10.1234/shared},serve={false},license={CC-BY-4.0}}')
    import_bibtex(tmp_path / "authority.sqlite", bib)
    assert _load_skipped_hashes(tmp_path / "authority.sqlite") == {"aaa"}
    assert document_metadata(conn, "aaa")["publishable"] == 1
    conn.close()


def test_unknown_hash_cannot_apply_policy_via_doi_fallback(tmp_path):
    import pytest
    paper(tmp_path, "aaa")
    conn = database(tmp_path)
    run(conn, tmp_path)
    bib = tmp_path / "update.bib"
    bib.write_text('@article{x,corpus_hash={missing},doi={10.1234/shared},serve={false}}')
    with pytest.raises(ValueError, match="not a member"):
        import_bibtex(tmp_path / "authority.sqlite", bib)
    assert _load_skipped_hashes(tmp_path / "authority.sqlite") == set()
    conn.close()


def test_member_removal_preserves_curated_work_header(tmp_path):
    hd = paper(tmp_path, "aaa")
    paper(tmp_path, "bbb")
    conn = database(tmp_path)
    run(conn, tmp_path)
    bib = tmp_path / "update.bib"
    bib.write_text('@article{x,doi={10.1234/shared},title={Curated title},author={Curator, C}}')
    import_bibtex(tmp_path / "authority.sqlite", bib)
    shutil.rmtree(hd)
    run(conn, tmp_path)
    assert conn.execute("SELECT title FROM works WHERE in_corpus=1").fetchone()[0] == "Curated title"
    assert conn.execute("SELECT surname FROM work_authors WHERE work_id='10.1234/shared'").fetchone()[0] == "Curator"
    conn.close()
