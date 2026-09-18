"""Reference fragments stay inspectable without becoming publications (#313)."""
import json
import sqlite3
from pathlib import Path
from types import SimpleNamespace

from bib.authority import (create_schema, insert_work, phase1_corpus_papers,
                           phase2_references)
from bib.reference_quality import classify
from mcpsrv import app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import (format_citations, get_bibliography,
    get_missing_references, get_works_by_author, resolve_reference)


def build(tmp_path, references):
    sha='citing'
    doc=tmp_path/'documents'/sha
    doc.mkdir(parents=True)
    (doc/'metadata.json').write_text(json.dumps({'title':'An actual publication','authors':[{'surname':'Real','forename':'A'}],'year':2000}))
    (doc/'references.json').write_text(json.dumps({'references':references}))
    db=tmp_path/'biblio_authority.sqlite'
    conn=sqlite3.connect(db)
    create_schema(conn)
    phase1_corpus_papers(conn,tmp_path)
    phase2_references(conn,tmp_path)
    return conn,db,doc


def test_reported_fragments_across_all_authority_routes(tmp_path,monkeypatch):
    fixture=Path(__file__).parent/'fixtures/bibliographic_integrity/reference_fragments.json'
    refs=json.loads(fixture.read_text())['references']
    # The newer issue includes this independent figure label; the retained
    # older parsed-only title does not justify quarantine by itself.
    refs[0]['raw']='R.gp-lat. Nempect Fic. 37. '+refs[0]['title']
    # Also exercise the newer report's raw-evidence shape, independently of
    # the retained older bundle's parsed values and absent raw strings.
    refs.extend([
        {'xml_id':'b3','raw':'Plate XIV, figures 1-6','title':'','authors':[],'year':None,'journal':'Plate XIV'},
        {'xml_id':'b4','raw':'C.rad.lat = lateral radial canal; Nem.p.ect = ectodermal patch; R.ap.lat = ridge;',
         'title':'','authors':['C Rad'],'year':None},
    ])
    conn,db,doc=build(tmp_path,refs)
    before=list(conn.execute('SELECT * FROM reference_observations ORDER BY observation_id'))
    assert len(before)==5
    assert conn.execute('SELECT COUNT(*) FROM observation_work').fetchone()[0]==0
    assert conn.execute('SELECT COUNT(*) FROM citations').fetchone()[0]==0
    assert conn.execute("SELECT COUNT(*) FROM reference_observation_quality WHERE disposition='quarantined_fragment'").fetchone()[0]==5
    ba=BiblioAuthority(db)
    monkeypatch.setattr(app,'_INDEX',SimpleNamespace(biblio_db=ba,papers={'citing':{'hash_dir':str(doc)}}))
    for resolved in (True,False):
        results=get_bibliography('citing',resolved=resolved)
        assert len(results)==5
        for expected,actual in zip(refs,results):
            assert actual['raw']==expected.get('raw','')
            assert actual['quality']['disposition']=='quarantined_fragment'
            assert actual['quality']['grobid_xml_id']==expected['xml_id']
            if resolved:
                assert actual['work_id'] is None
    assert get_works_by_author('Gp-Lat')==[]
    assert get_works_by_author('Rad')==[]
    assert get_missing_references(min_citations=1)==[]
    assert resolve_reference('Gp-Lat 1900').get('not_found')
    assert format_citations(queries=['Gp-Lat 1900'])['citations'][0]['code']=='not_found'
    changes=conn.total_changes
    assert phase2_references(conn,tmp_path)==(0,0)
    assert conn.total_changes==changes
    assert list(conn.execute('SELECT * FROM reference_observations ORDER BY observation_id'))==before


def test_sparse_historical_and_long_titles_remain_usable(tmp_path):
    refs=[{'xml_id':'b0','authors':['A Author'],'title':'','year':None,'raw':'A. Author. An old treatise.'},
          {'xml_id':'b1','authors':['B Writer'],'title':'A study of '+('historical observations '*30),'year':None,'raw':''},
          {'xml_id':'b2','authors':[],'title':'Plate XIV','year':1870,'raw':'Plate XIV. 1870.'}]
    conn,_,_=build(tmp_path,refs)
    assert conn.execute('SELECT COUNT(*) FROM observation_work').fetchone()[0]==3
    assert conn.execute("SELECT COUNT(*) FROM reference_observation_quality WHERE disposition='quarantined_fragment'").fetchone()[0]==0


def test_uncurated_panel_wording_is_reviewable_and_visible_across_routes(tmp_path,monkeypatch):
    refs=[{'xml_id':'b0','authors':['A Author'],
           'title':'A, B, C views of historical observations','year':None,'raw':''}]
    conn,db,doc=build(tmp_path,refs)
    work_id=conn.execute('SELECT work_id FROM observation_work').fetchone()[0]
    assert conn.execute('SELECT disposition FROM reference_observation_quality').fetchone()[0]=='review_needed'
    assert conn.execute('SELECT COUNT(*) FROM citations').fetchone()[0]==1
    ba=BiblioAuthority(db)
    monkeypatch.setattr(app,'_INDEX',SimpleNamespace(biblio_db=ba,papers={'citing':{'hash_dir':str(doc)}}))
    assert get_bibliography('citing',resolved=True)[0]['work_id']==work_id
    assert get_bibliography('citing',resolved=True)[0]['quality']['disposition']=='review_needed'
    for work in [ba.get_work(work_id),ba.search_works('Author')[0],
                 get_works_by_author('Author')[0],get_missing_references(min_citations=1)[0],
                 format_citations(work_ids=[work_id])['citations'][0]]:
        assert work['reference_quality_warnings'][0]['reasons'][0]['code']=='possible_panel_description'
    changes=conn.total_changes
    assert phase2_references(conn,tmp_path)==(0,0)
    assert conn.total_changes==changes


def test_old_caption_title_without_independent_anchor_is_only_reviewable():
    conn=sqlite3.connect(':memory:')
    create_schema(conn)
    fixture=Path(__file__).parent/'fixtures/bibliographic_integrity/reference_fragments.json'
    ref=json.loads(fixture.read_text())['references'][0]
    assert classify(conn,ref)[0]=='review_needed'


def test_curated_title_is_counterevidence_for_unusual_publication():
    conn=sqlite3.connect(':memory:')
    create_schema(conn)
    title='Plate XIV'
    insert_work(conn,'curated','corpus_key',title,None,'','',None,False,'cited_reference')
    conn.execute("UPDATE works SET bib_imported_at=1 WHERE work_id='curated'")
    disposition,reasons=classify(conn,{'title':title,'authors':[]})
    assert disposition=='usable'
    assert reasons[0]['code']=='curated_title_counterevidence'


def test_quarantine_is_rederived_when_evidence_changes_and_raw_history_survives(tmp_path):
    conn,_,doc=build(tmp_path,[{'xml_id':'b0','title':'Plate XIV','authors':[],'raw':'Plate XIV'}])
    old_observation=conn.execute('SELECT observation_id FROM reference_observations').fetchone()[0]
    (doc/'references.json').write_text(json.dumps({'references':[{'xml_id':'b0','title':'An actual historical work','authors':['A Author'],'year':1900,'raw':'Author 1900. An actual historical work.'}]}))
    phase2_references(conn,tmp_path)
    assert conn.execute('SELECT COUNT(*) FROM reference_observations').fetchone()[0]==2
    assert conn.execute('SELECT 1 FROM reference_observations WHERE observation_id=?',(old_observation,)).fetchone()
    assert conn.execute('SELECT COUNT(*) FROM observation_work').fetchone()[0]==1
    assert conn.execute('SELECT disposition FROM reference_observation_quality').fetchone()[0]=='usable'


def test_truncated_surnames_are_reviewed_and_never_guessed_or_merged(tmp_path,monkeypatch):
    refs=[{'xml_id':f'b{i}','authors':[f'A {surname}'],'title':f'A distinct historical observation {i}',
           'year':1853,'raw':f'A {surname} 1853. Observation {i}.'}
          for i,surname in enumerate(['¨lliker','˚mstedt','¨stman','˜o'])]
    conn,db,doc=build(tmp_path,refs)
    assert conn.execute("SELECT COUNT(*) FROM reference_observation_quality WHERE disposition='review_needed'").fetchone()[0]==4
    assert conn.execute("SELECT COUNT(*) FROM observation_work WHERE match_method='unresolved_author'").fetchone()[0]==4
    for author, in conn.execute("SELECT authors_json FROM reference_observations"):
        assert json.loads(author)[0] in [ref['authors'][0] for ref in refs]
    assert not conn.execute("SELECT 1 FROM work_aliases WHERE work_id LIKE 'corpus:unresolved-author|%'").fetchone()
    monkeypatch.setattr(app,'_INDEX',SimpleNamespace(biblio_db=BiblioAuthority(db),papers={'citing':{'hash_dir':str(doc)}}))
    assert all(row['quality']['reasons'][0]['code']=='suspect_truncated_surname' for row in get_bibliography('citing',resolved=True))


def test_orphan_accents_and_complete_multilingual_names_have_distinct_quality():
    from bib.reference_quality import author_quality_reasons
    for author in ['D. ¨ Ursprung','¨rsprung, D.','D. \u0308 Ursprung']:
        assert author_quality_reasons({'authors':[author]})
    for author in ['D. Ursprung','R. Kölliker','A. Alvariño','J. Åmstedt','S. Östman','A. Niño','F. Pacifici','B. Mu\u0308ller']:
        assert author_quality_reasons({'authors':[author]})==[]


def test_orphan_author_identity_is_separate_in_both_ingestion_orders():
    from bib.authority import _resolve_reference
    from bib.reconcile import find_candidates
    for reverse in (False,True):
        conn=sqlite3.connect(':memory:');create_schema(conn)
        refs=[{'authors':['D ¨Ursprung'],'title':'A substantial study of plankton and animals','year':1965},
              {'authors':['D Ursprung'],'title':'A substantial study of plankton and animals','year':1965}]
        if reverse:refs.reverse()
        results={ref['authors'][0]:_resolve_reference(conn,ref,fallback_key='bad-observation') for ref in refs}
        assert results['D ¨Ursprung'][0]!=results['D Ursprung'][0]
        assert results['D ¨Ursprung'][1]=='unresolved_author'
        assert all(not row[0].startswith('corpus:unresolved-author|') for row in find_candidates(conn,'Ursprung',1965))


def test_independent_doi_still_resolves_a_suspect_author_without_an_alias():
    from bib.authority import _resolve_reference
    conn=sqlite3.connect(':memory:');create_schema(conn)
    work,method,_=_resolve_reference(conn,{'authors':['A ¨lliker'],'title':'Some observations','year':1853,'doi':'10.1234/source'})
    assert method=='doi_exact' and work=='10.1234/source'
    assert not conn.execute('SELECT 1 FROM work_aliases WHERE work_id=?',(work,)).fetchone()
