"""Taxonomy and bibliography edits replace derived original-description links."""
import sqlite3

from bib import authority as a


def taxonomy(path, rows):
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE IF NOT EXISTS taxa(taxon_id TEXT, scientific_name_authorship TEXT)")
    conn.execute("DELETE FROM taxa")
    conn.executemany("INSERT INTO taxa VALUES (?,?)", rows)
    conn.commit()
    conn.close()


def database():
    conn = sqlite3.connect(":memory:")
    conn.execute("PRAGMA foreign_keys=ON")
    a.create_schema(conn)
    return conn


def work(conn, work_id, surname, year, second=None):
    a.insert_work(conn, work_id, "corpus_key", "A real work", year, "", None, work_id, True, "corpus_paper")
    a.insert_authors(conn, work_id, [(surname, ""), *([(second, "")] if second else [])])


def logical(conn):
    return (conn.execute("SELECT * FROM taxon_work_links ORDER BY 1,2,3").fetchall(),
            conn.execute("SELECT work_id,title,year,source FROM works ORDER BY 1").fetchall())


def test_authority_edit_removal_and_noop_match_clean(tmp_path):
    tx = tmp_path / "taxonomy.sqlite"
    taxonomy(tx, [("1", "Author, 1900"), ("2", "Other, 1910")])
    conn = database()
    a.phase3_authority_links(conn, tx)
    changes = conn.total_changes
    assert a.phase3_authority_links(conn, tx) == 0
    assert conn.total_changes == changes
    taxonomy(tx, [("1", "Correct, 1920")])
    a.phase3_authority_links(conn, tx)
    clean = database()
    a.phase3_authority_links(clean, tx)
    assert logical(conn) == logical(clean)
    conn.close()
    clean.close()


def test_adding_real_work_replaces_taxon_stub_and_matches_clean(tmp_path):
    tx = tmp_path / "taxonomy.sqlite"
    taxonomy(tx, [("1", "Author, 1900")])
    conn = database()
    a.phase3_authority_links(conn, tx)
    work(conn, "paper", "Author", 1900)
    a.phase3_authority_links(conn, tx)
    clean = database()
    work(clean, "paper", "Author", 1900)
    a.phase3_authority_links(clean, tx)
    assert logical(conn) == logical(clean)
    assert conn.execute("SELECT work_id FROM taxon_work_links").fetchone()[0] == "paper"
    conn.close()
    clean.close()


def test_missing_snapshot_retires_derived_but_keeps_curated_links(tmp_path):
    tx = tmp_path / "taxonomy.sqlite"
    taxonomy(tx, [("1", "Author, 1900")])
    conn = database()
    work(conn, "paper", "Author", 1900)
    a.phase3_authority_links(conn, tx)
    conn.execute("INSERT INTO taxon_work_links VALUES ('2','paper','curator',1.0)")
    tx.unlink()
    a.phase3_authority_links(conn, tx)
    assert conn.execute("SELECT taxon_id,link_type FROM taxon_work_links").fetchall() == [("2", "curator")]
    conn.close()


def test_ambiguous_second_author_does_not_choose_by_insertion_order(tmp_path):
    tx = tmp_path / "taxonomy.sqlite"
    taxonomy(tx, [("1", "Author & Other, 1900")])
    conn = database()
    work(conn, "first", "Author", 1900, "Other")
    work(conn, "second", "Author", 1900, "Other")
    a.phase3_authority_links(conn, tx)
    work_id, confidence = conn.execute("SELECT work_id,confidence FROM taxon_work_links").fetchone()
    assert work_id not in {"first", "second"}
    assert confidence == 0.5
    conn.close()


def test_retiring_taxon_does_not_remove_a_work_still_cited(tmp_path):
    tx = tmp_path / "taxonomy.sqlite"
    taxonomy(tx, [("1", "Author, 1900")])
    conn = database()
    a.phase3_authority_links(conn, tx)
    stub = conn.execute("SELECT work_id FROM taxon_work_links").fetchone()[0]
    work(conn, "paper", "Other", 2000)
    a.insert_citation(conn, "paper", stub, "paper", "b0", "Author, 1900", "exact", 1.0)
    taxonomy(tx, [])
    a.phase3_authority_links(conn, tx)
    assert conn.execute("SELECT cited_work_id FROM citations").fetchone()[0] == stub
    assert conn.execute("SELECT 1 FROM works WHERE work_id=?", (stub,)).fetchone()
    conn.close()
