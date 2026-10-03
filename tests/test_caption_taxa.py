"""Caption abbreviation evidence is built once and shared by figure routes."""
import json
import sqlite3
from pathlib import Path
from types import SimpleNamespace

import pytest

from pipeline.caption_taxa import caption_evidence, materialize
from pipeline.taxa import TaxonomyDB
from pipeline.taxonomy_ingest import create_schema as taxonomy_schema, insert_records, make_record
from pipeline.taxon_mentions import build, create_schema
from mcpsrv import app
from mcpsrv.indexes import TaxonMentionDB
from mcpsrv.tools.figures import get_figure, get_figure_dossier_for_taxon, get_figures_for_taxon


@pytest.fixture
def taxonomy(tmp_path):
    conn=sqlite3.connect(tmp_path/'taxonomy.sqlite')
    taxonomy_schema(conn)
    insert_records(conn,[
        make_record(taxon_id='1',scientific_name='Nanomia',taxon_rank='Genus'),
        make_record(taxon_id='2',scientific_name='Nanomia septata',taxon_rank='Species',extra_names=['Nanomia oldname']),
        make_record(taxon_id='3',scientific_name='Nanomia cara',taxon_rank='Species'),
        make_record(taxon_id='4',scientific_name='Nanomia bijuga',taxon_rank='Species'),
        make_record(taxon_id='5',scientific_name='Notheria',taxon_rank='Genus'),
        make_record(taxon_id='6',scientific_name='Notheria septata',taxon_rank='Species'),
    ])
    conn.commit();conn.close()
    db=TaxonomyDB(tmp_path/'taxonomy.sqlite')
    yield db
    db.close()


def evidence(taxonomy,text,context=None):
    from collections import Counter
    return caption_evidence(text,taxonomy,taxonomy.name_set(),Counter(context or {}))


def test_real_hosia_caption_uses_explicit_local_genus(taxonomy):
    caption=json.loads((Path(__file__).parent/'fixtures/hosia_caption_taxa.json').read_text())['caption']
    result=evidence(taxonomy,caption,{'Notheria':1})
    hits=[m for m in result['matches'] if m['accepted_taxon_id']=='2']
    assert len(hits)==1
    hit=hits[0]
    assert hit['mention_text']=='N. septata'
    assert hit['method']=='contextual_abbreviation'
    assert hit['context_scope']=='caption'
    assert hit['context_genera']==['Nanomia']
    assert caption[slice(*hit['text_span'])]=='N. septata'
    assert not result['unresolved']


@pytest.mark.parametrize('text',[
    'N. septata.',
    'Nanomia and Notheria. N. septata.',
])
def test_ambiguous_initial_is_explicit_and_never_guessed(taxonomy,text):
    result=evidence(taxonomy,text)
    assert not any(m['method']=='contextual_abbreviation' for m in result['matches'])
    assert result['unresolved'][0]['reason']=='ambiguous_abbreviation'
    assert result['unresolved'][0]['candidate_names']==['Nanomia septata','Notheria septata']


@pytest.mark.parametrize('text',[
    '(NANOMIA SEPTATA); nANOmia oldname.',
    'Nanomia: n. septata; N, septata.',
])
def test_case_punctuation_and_full_aliases(taxonomy,text):
    result=evidence(taxonomy,text)
    assert len([m for m in result['matches'] if m['accepted_taxon_id']=='2'])==2


def test_epithet_overlap_and_word_substrings_never_link(taxonomy):
    result=evidence(taxonomy,'PseudoNanomia septata; Nanomia septatax; septata alone; Notheria septata.')
    assert not any(m['accepted_taxon_id']=='2' for m in result['matches'])
    assert any(m['accepted_taxon_id']=='6' for m in result['matches'])


def prepare(tmp_path,caption,chunk_text=''):
    doc=tmp_path/'documents'/'paper'
    doc.mkdir(parents=True,exist_ok=True)
    (doc/'figures.json').write_text(json.dumps({'figures':[{
        'figure_id':'docling_3','figure_number':1,'figure_type':'figure','page':3,
        'caption_text':caption,'caption_status':'bound'}]}))
    (doc/'chunks.json').write_text(json.dumps({'chunks':[{'chunk_id':'chunk_0','text':chunk_text}]}))
    db=tmp_path/'taxon_mentions.sqlite'
    conn=sqlite3.connect(db)
    create_schema(conn)
    return conn,db,doc


def test_post_build_routes_direct_provenance_and_caption_only_discovery(tmp_path,taxonomy,monkeypatch):
    caption=json.loads((Path(__file__).parent/'fixtures/hosia_caption_taxa.json').read_text())['caption']
    conn,db,doc=prepare(tmp_path,caption)
    assert build(conn,tmp_path)['captions']['papers']==1  # no taxa.json needed
    stored=TaxonMentionDB(db)
    idx=SimpleNamespace(taxonomy_db=taxonomy,taxon_mention_db=stored,biblio_db=None,
        papers={'paper':{'hash_dir':str(doc),'title':'Hosia et al.'}},taxon_to_papers={},taxon_mention_counts={})
    monkeypatch.setattr(app,'_INDEX',idx)
    hits=get_figures_for_taxon('Nanomia septata',caption_only=True)
    assert len(hits)==1 and hits[0]['caption_has_taxon'] is True
    dossier=get_figure_dossier_for_taxon('Nanomia septata')['figures']
    direct=get_figure('paper','docling_3')
    assert hits[0]['caption_taxa']==dossier[0]['caption_taxa']==direct['caption_taxa']
    # The explicit genus remains queryable just as in older full-name matching.
    assert get_figures_for_taxon('Nanomia',caption_only=True)
    assert direct['caption_taxa']['input_fingerprint']['taxonomy_sha256']
    changes=conn.total_changes
    assert build(conn,tmp_path)['captions']['skipped']==1
    assert conn.total_changes==changes


def test_caption_document_taxonomy_changes_each_rederive_links(tmp_path,taxonomy):
    conn,db,doc=prepare(tmp_path,'N. septata.','Nanomia')
    materialize(conn,tmp_path)
    reader=TaxonMentionDB(db)
    assert reader.caption_papers('2')==['paper']
    hit=next(m for m in reader.caption_evidence('paper','docling_3')['matches'] if m['accepted_taxon_id']=='2')
    assert hit['context_scope']=='document'
    assert hit['context_mentions']==[{'chunk_id':'chunk_0','mention_text':'Nanomia','text_span':[0,7]}]
    (doc/'chunks.json').write_text(json.dumps({'chunks':[{'text':'Nanomia and Notheria'}]}))
    materialize(conn,tmp_path)
    assert reader.caption_papers('2')==[]
    assert reader.caption_evidence('paper','docling_3')['unresolved']
    (doc/'figures.json').write_text(json.dumps({'figures':[{'figure_id':'docling_3','caption_text':'Nanomia septata.'}]}))
    materialize(conn,tmp_path)
    assert reader.caption_papers('2')==['paper']
    tx=sqlite3.connect(tmp_path/'taxonomy.sqlite')
    tx.execute("UPDATE taxa SET accepted_name_usage_id='6',accepted_name='Notheria septata' WHERE taxon_id='2'")
    tx.commit();tx.close()
    assert materialize(conn,tmp_path)['papers']==1
    assert reader.caption_papers('2')==[]
    assert reader.caption_papers('6')==['paper']
    # An intentionally emptied figure artifact retires old links.
    (doc/'figures.json').write_text('{"figures":[]}')
    materialize(conn,tmp_path)
    assert reader.caption_papers('6')==[]


def test_legacy_bundle_keeps_full_name_lookup_without_abbreviation_inference(tmp_path,taxonomy,monkeypatch):
    _,db,doc=prepare(tmp_path,'Nanomia: N. septata.')
    # Explicitly simulate a legacy schema rather than an unbuilt current one.
    conn=sqlite3.connect(db)
    for table in ('caption_taxon_evidence','caption_taxon_links','caption_taxon_receipts'):
        conn.execute(f'DROP TABLE {table}')
    conn.commit();conn.close()
    reader=TaxonMentionDB(db)
    idx=SimpleNamespace(taxonomy_db=taxonomy,taxon_mention_db=reader,biblio_db=None,
        papers={'paper':{'hash_dir':str(doc)}},taxon_to_papers={'2':['paper']},taxon_mention_counts={})
    monkeypatch.setattr(app,'_INDEX',idx)
    assert get_figures_for_taxon('Nanomia septata',caption_only=True)==[]
    assert get_figure('paper','docling_3')['caption_taxa']['availability']=='legacy_unavailable'
    (doc/'figures.json').write_text(json.dumps({'figures':[{'figure_id':'docling_3','figure_type':'figure','caption_text':'Nanomia septata.'}]}))
    assert get_figures_for_taxon('Nanomia septata',caption_only=True)


def test_corrupt_update_preserves_prior_evidence_and_blocks_success(tmp_path,taxonomy):
    conn,db,doc=prepare(tmp_path,'Nanomia septata.')
    materialize(conn,tmp_path)
    (doc/'figures.json').write_text('{broken')
    assert build(conn,tmp_path)['errors']==1
    assert TaxonMentionDB(db).caption_papers('2')==['paper']


def test_clean_and_incremental_caption_evidence_agree(tmp_path,taxonomy):
    import shutil
    conn,_,doc=prepare(tmp_path,'N. septata.','Nanomia')
    materialize(conn,tmp_path)
    (doc/'figures.json').write_text(json.dumps({'figures':[{
        'figure_id':'replacement','caption_text':'Nanomia oldname.'}]}))
    materialize(conn,tmp_path)
    fresh=tmp_path/'fresh'
    shutil.copytree(tmp_path/'documents',fresh/'documents')
    shutil.copyfile(tmp_path/'taxonomy.sqlite',fresh/'taxonomy.sqlite')
    clean=sqlite3.connect(fresh/'taxon_mentions.sqlite')
    create_schema(clean)
    materialize(clean,fresh)
    for table in ('caption_taxon_links','caption_taxon_evidence','caption_taxon_receipts'):
        assert list(conn.execute(f'SELECT * FROM {table} ORDER BY 1,2'))==list(clean.execute(f'SELECT * FROM {table} ORDER BY 1,2'))
    shutil.rmtree(doc)
    materialize(conn,tmp_path)
    assert conn.execute('SELECT COUNT(*) FROM caption_taxon_links').fetchone()[0]==0
