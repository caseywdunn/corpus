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
