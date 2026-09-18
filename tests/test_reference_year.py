"""Per-observation publication-year evidence; title dates never imply a merge."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from bib.authority import phase1_corpus_papers, phase2_references
from bib.documents import find_work
from bib.parser import parse_bibtex, bib_entry_to_metadata
from bib.reference_year import adjudicate, candidate_index
from mcpsrv import app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import get_bibliography, get_missing_references
from tests.test_bibliographic_integrity import build

FIXTURE = Path(__file__).parent / 'fixtures/bibliographic_integrity'
ENTRY = next(e for e in parse_bibtex((FIXTURE / 'source.bib').read_text()) if e['_key']=='Pugh1974')
REF = json.loads((FIXTURE / 'publication_year_conflict.json').read_text())['reference']


def references(tmp_path, sha, refs):
    path=tmp_path/'documents'/sha/'references.json'
    path.write_text(json.dumps({'references':refs}))


def test_reported_raw_year_repairs_old_mapping_and_deduplicates_citing_edges(tmp_path,monkeypatch):
    from bib import reference_year
    conn,db=build(tmp_path,{'b5a7af6140ca':ENTRY,'db7338ea1867':{'title':'A citing study','author':'Writer, A','year':'2000'},
                            'second':{'title':'Another citing study','author':'Reader, B','year':'2001'}})
    target=find_work(conn,'b5a7af6140ca')
    correct=dict(REF,xml_id='b41',year=1974,title=ENTRY['title'])
    missing={'xml_id':'b42','authors':['Absent A'],'title':'A genuinely unavailable historical publication','year':1850,'raw':'Absent 1850. A genuinely unavailable historical publication.'}
    references(tmp_path,'db7338ea1867',[REF,correct,missing])
    references(tmp_path,'second',[dict(REF,xml_id='b2')])
    with monkeypatch.context() as old:
        old.setattr(reference_year,'adjudicate',lambda ref,index:(None,[]))
        phase2_references(conn,tmp_path)
    assert conn.execute('SELECT COUNT(*) FROM works WHERE year=1965').fetchone()[0]==1
    raw_before=list(conn.execute('SELECT * FROM reference_observations ORDER BY observation_id'))
    conn.execute("UPDATE observation_work SET producer_version='legacy'")
    conn.commit()
    phase2_references(conn,tmp_path)
    assert list(conn.execute('SELECT * FROM reference_observations ORDER BY observation_id'))==raw_before
    repaired=list(conn.execute("SELECT ow.work_id,ow.match_method,q.reasons_json FROM observation_work ow JOIN reference_observations o USING(observation_id) JOIN reference_observation_quality q USING(observation_id) WHERE o.year=1965"))
    assert len(repaired)==2
    assert all(wid==target and method=='raw_publication_year_title_authors' for wid,method,_ in repaired)
    assert all(json.loads(reason)[0]['publication_year']==1974 for _,_,reason in repaired)
    assert conn.execute('SELECT COUNT(*) FROM citations WHERE cited_work_id=?',(target,)).fetchone()[0]==2
    assert conn.execute('SELECT COUNT(*) FROM works WHERE year=1965').fetchone()[0]==0
    ba=BiblioAuthority(db)
    monkeypatch.setattr(app,'_INDEX',SimpleNamespace(biblio_db=ba,papers={'db7338ea1867':{'hash_dir':str(tmp_path/'documents/db7338ea1867')}}))
    assert get_bibliography('db7338ea1867',resolved=True)[0]['work_id']==target
    suggestions=get_missing_references(min_citations=1)
    assert len(suggestions)==1 and suggestions[0]['year']==1850
    changes=conn.total_changes
    assert phase2_references(conn,tmp_path)==(0,0)
    assert conn.total_changes==changes


@pytest.mark.parametrize('raw',[
    '',
    'Pugh. '+ENTRY['title'],
    'Pugh 1974, 1975. '+ENTRY['title'],
    'Different 1974. '+ENTRY['title'],
    'Pugh 1965. '+ENTRY['title'],
])
def test_missing_or_ambiguous_raw_publication_evidence_stays_reviewable(tmp_path,monkeypatch,raw):
    conn,db=build(tmp_path,{'article':ENTRY,'citing':{'title':'A citing study','author':'Writer, A','year':'2000'}})
    references(tmp_path,'citing',[dict(REF,raw=raw)])
    phase2_references(conn,tmp_path)
    wid=conn.execute('SELECT work_id FROM observation_work').fetchone()[0]
    assert wid!=find_work(conn,'article')
    assert conn.execute('SELECT disposition FROM reference_observation_quality').fetchone()[0]=='review_needed'
    monkeypatch.setattr(app,'_INDEX',SimpleNamespace(biblio_db=BiblioAuthority(db)))
    missing=get_missing_references(min_citations=1)
    assert missing[0]['reference_quality_warnings'][0]['reasons'][0]['code']=='possible_publication_year_conflict'


def test_date_in_title_correct_year_does_not_trigger_repair(tmp_path):
    conn,_=build(tmp_path,{'article':ENTRY})
    assert adjudicate(dict(REF,title=ENTRY['title'],year=1974),candidate_index(conn))==(None,[])
    # Same surname and date with a different substantive title is not a lead.
    assert adjudicate(dict(REF,title='An unrelated study of currents during 1965'),candidate_index(conn))==(None,[])
    # Incomplete or contradictory authors/DOI must not repair this observation.
    assert adjudicate(dict(REF,authors=['P R Pugh','B Another']),candidate_index(conn))==(None,[])
    assert adjudicate(dict(REF,doi='10.9999/different'),candidate_index(conn))==(None,[])


@pytest.mark.parametrize('shared_doi',[False,True])
def test_indistinguishable_editions_or_parts_are_not_selected(tmp_path,shared_doi):
    second=dict(ENTRY,_key='OtherEdition',edition='2',doi=ENTRY['doi'] if shared_doi else '10.9999/edition2')
    conn,_=build(tmp_path,{'first':dict(ENTRY,edition='1'),'second':second})
    target,reasons=adjudicate(REF,candidate_index(conn))
    assert target is None
    assert reasons[0]['candidate_count']==2


def test_curated_metadata_change_rederives_year_decision(tmp_path):
    conn,_=build(tmp_path,{'article':ENTRY,'citing':{'title':'A citing study','author':'Writer, A','year':'2000'}})
    references(tmp_path,'citing',[REF])
    phase2_references(conn,tmp_path)
    assert conn.execute('SELECT match_method FROM observation_work').fetchone()[0]=='raw_publication_year_title_authors'
    path=tmp_path/'documents/article/metadata.json'
    path.write_text(json.dumps(bib_entry_to_metadata(dict(ENTRY,year='1975'),'article.pdf')))
    phase1_corpus_papers(conn,tmp_path)
    phase2_references(conn,tmp_path)
    assert conn.execute('SELECT disposition FROM reference_observation_quality').fetchone()[0]=='review_needed'
    assert conn.execute('SELECT year FROM reference_observations').fetchone()[0]==1965


def test_year_evidence_does_not_override_a_damaged_author(tmp_path):
    conn,_=build(tmp_path,{'article':ENTRY})
    assert adjudicate(dict(REF,authors=['P R ¨Pugh']),candidate_index(conn))==(None,[])
