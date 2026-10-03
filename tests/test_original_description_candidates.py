"""Source-supported candidates must remain distinct from curator verdicts (#311)."""
import json
import sqlite3
from pathlib import Path
from types import SimpleNamespace

import pytest

from bib.authority import (create_schema, insert_authors, insert_work, parse_authority,
                           phase1_corpus_papers, phase3_authority_links)
from bib.parser import bib_entry_to_metadata, parse_bibtex
from mcpsrv import app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import get_original_description


@pytest.mark.parametrize("authority,expected", [
    ("Siebert, Pugh, Haddock & Dunn, 2013", ["Siebert", "Pugh", "Haddock", "Dunn"]),
    ("Siebert, S., Pugh, P.R., Haddock, S.H.D. and Dunn, C.W., 2013", ["Siebert", "Pugh", "Haddock", "Dunn"]),
    ("S. Siebert, P.R.Pugh, S.H.D. Haddock & C.W. Dunn, 2013", ["Siebert", "Pugh", "Haddock", "Dunn"]),
    ("Lens, van Riemsdijk and de Blainville, 1908", ["Lens", "van Riemsdijk", "de Blainville"]),
    ("(Linnaeus, 1758)", ["Linnaeus"]),
    ("É. Perrier, 1882", ["Perrier"]),
])
def test_complete_zoological_author_lists(authority, expected):
    assert parse_authority(authority)[0] == expected


def test_botanical_conjunction_is_not_one_combined_surname():
    assert parse_authority("Rehder and Hesse") == (["Rehder", "Hesse"], None)


def taxonomy(path, rows):
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE IF NOT EXISTS taxa(taxon_id TEXT,scientific_name TEXT,scientific_name_authorship TEXT)")
    conn.execute("DELETE FROM taxa")
    conn.executemany("INSERT INTO taxa VALUES (?,?,?)", rows)
    conn.commit()
    conn.close()


def install(monkeypatch, db, name, taxon_id):
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(
        biblio_db=BiblioAuthority(db),
        taxonomy_db=SimpleNamespace(lookup=lambda query: {"accepted_taxon_id": taxon_id, "accepted_name": name})))


def test_real_apolemia_entry_and_description_excerpt_are_supported_candidate(tmp_path, monkeypatch):
    fixture = Path(__file__).parent / "fixtures/bibliographic_integrity"
    excerpt = json.loads((fixture / "description_excerpt.json").read_text())
    entry = next(e for e in parse_bibtex((fixture / "source.bib").read_text()) if e['_key']=='Siebertetal2013')
    doc = tmp_path / "documents" / excerpt["corpus_hash"]
    doc.mkdir(parents=True)
    (doc / "metadata.json").write_text(json.dumps(bib_entry_to_metadata(entry, entry['file'])))
    (doc / "text.json").write_text(json.dumps({"text": excerpt['text']}))
    db = tmp_path / "biblio_authority.sqlite"
    conn = sqlite3.connect(db)
    create_schema(conn)
    phase1_corpus_papers(conn, tmp_path)
    tx = tmp_path / "taxonomy.sqlite"
    taxonomy(tx, [(excerpt["taxon_id"], excerpt["taxon_name"], excerpt["authorship"])])
    # A pre-fix malformed-author stub is derived history, not an authority.
    insert_work(conn, 'old-stub', 'corpus_key', '', 2013, '', '', None, False, 'taxon_authority')
    insert_authors(conn, 'old-stub', [('Siebert, Pugh, Haddock', ''), ('Dunn', '')])
    conn.execute("INSERT INTO taxon_work_links VALUES (?,'old-stub','authority_match',0.5)", (excerpt['taxon_id'],))
    phase3_authority_links(conn, tx)
    assert conn.execute("SELECT 1 FROM works WHERE work_id='old-stub'").fetchone() is None
    before = conn.total_changes
    assert phase3_authority_links(conn, tx) == 0
    assert conn.total_changes == before
    install(monkeypatch, db, excerpt["taxon_name"], excerpt["taxon_id"])
    result = get_original_description(excerpt["taxon_name"])
    assert result["original_description"] is None
    candidate, = result["candidate_works"]
    assert candidate["work_id"] == '10.11646/zootaxa.3702.3.1'
    assert candidate["basis"]["complete_author_list_match"] is True
    evidence = candidate['basis']['source_evidence']
    assert evidence["kind"] == "opening_text_new_species_marker"
    assert "A. rubriversa sp. nov" in evidence["excerpt"]
    assert candidate["confidence"] == 0.95
    assert result["authority_stubs"] is None
    # A curator verdict is a distinct input and survives re-derivation.
    conn.execute("INSERT INTO taxon_work_links VALUES (?,?,'curator',1)", (excerpt['taxon_id'], candidate['work_id']))
    conn.commit()
    (doc / 'text.json').write_text(json.dumps({'text': 'A taxon mention alone: Apolemia rubriversa.'}))
    phase3_authority_links(conn, tx)
    reviewed = get_original_description(excerpt['taxon_name'])
    assert reviewed['original_description']['work_id'] == candidate['work_id']
    assert reviewed['original_description']['basis']['kind'] == 'curator_reviewed'


@pytest.mark.parametrize("name,authority", [("Physalia physalis", "(Linnaeus, 1758)"), ("Nectopyramis thetis", "Bigelow, 1911")])
def test_ambiguous_historical_candidates_are_not_promoted_without_source_review(tmp_path, monkeypatch, name, authority):
    db = tmp_path / 'biblio_authority.sqlite'
    conn = sqlite3.connect(db)
    create_schema(conn)
    surnames,year = parse_authority(authority)
    for wid,title in [('edition1','First historical edition'),('edition2','Second historical edition')]:
        insert_work(conn,wid,'corpus_key',title,year,'','',wid,True,'corpus_paper')
        insert_authors(conn,wid,[(surnames[0],'')])
        doc=tmp_path/'documents'/wid
        doc.mkdir(parents=True)
        (doc/'text.json').write_text(json.dumps({'text': name + ' is mentioned in a later discussion.'}))
    tx=tmp_path/'taxonomy.sqlite'
    taxonomy(tx,[('t',name,authority)])
    phase3_authority_links(conn,tx)
    install(monkeypatch,db,name,'t')
    result=get_original_description(name)
    assert result['original_description'] is None
    assert {w['work_id'] for w in result['candidate_works']} == {'edition1','edition2'}
    assert result['authority_stubs'][0]['title'] == ''
    assert all('source_evidence' not in w['basis'] for w in result['candidate_works'])
    assert all(w['basis']['requires_source_review'] for w in result['candidate_works'])


def test_empty_legacy_stub_is_not_a_located_description(tmp_path, monkeypatch):
    db=tmp_path/'biblio_authority.sqlite'
    conn=sqlite3.connect(db)
    create_schema(conn)
    tx=tmp_path/'taxonomy.sqlite'
    taxonomy(tx,[('t','Example species','Author, 1900')])
    phase3_authority_links(conn,tx)
    install(monkeypatch,db,'Example species','t')
    result=get_original_description('Example species')
    assert result['original_description'] is None
    assert result['candidate_works'] is None
    assert result['authority_stubs'][0]['title'] == ''


def test_reference_rematerialization_retains_curator_reviewed_ghost(tmp_path):
    from bib.authority import _clear_derived_reference_materialization
    conn = sqlite3.connect(':memory:')
    conn.execute('PRAGMA foreign_keys=ON')
    create_schema(conn)
    insert_work(conn,'reviewed','corpus_key','Historical description',1900,'','',None,False,'cited_reference')
    insert_authors(conn,'reviewed',[('Author','')])
    conn.execute("INSERT INTO taxon_work_links VALUES ('t','reviewed','curator',1)")
    _clear_derived_reference_materialization(conn)
    assert conn.execute("SELECT work_id FROM taxon_work_links WHERE link_type='curator'").fetchone()[0] == 'reviewed'
    assert conn.execute("SELECT title FROM works WHERE work_id='reviewed'").fetchone()[0] == 'Historical description'
