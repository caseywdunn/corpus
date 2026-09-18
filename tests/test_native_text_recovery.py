"""Original damaged-layer evidence must survive full-page OCR (#312)."""
import json
import os
from pathlib import Path
from types import SimpleNamespace

import fitz
import pytest

from pipeline import native_text_recovery as recovery

SOURCE = json.loads((Path(__file__).parent / 'fixtures/text_integrity/german_native_layer.json').read_text())


@pytest.mark.parametrize('raw', ['Ã', 'Ãngstrom', 'literal Ã here', 'Ã¤', 'valid Fänge', 'FÃ¤ng 12²'])
def test_literal_or_nonreversible_text_is_not_a_candidate(raw):
    assert recovery.decoded_hint(raw) is None


def test_reversible_decoding_only_selects_a_candidate():
    assert recovery.decoded_hint('RMT-8-FÃ¤ng') == 'RMT-8-Fäng'
    assert recovery.supported_reading('RMT-8-Fäng', ['RMT-8-Fänge', 'RMT-8-Finge']) is None
    assert recovery.supported_reading('RMT-8-Fäng', ['RMT-9-Fänge']*2) is None
    assert recovery.supported_reading('RMT-8-Fäng', ['RMT-8-Fänge']) is None
    assert recovery.supported_reading('RMT-8-Fäng', ['RMT-8-Fänge RMT-8-Fängen']*2) is None
    assert recovery.supported_reading('RMT-8-Fäng', ['RMT-8-Fänge.']*2) == 'RMT-8-Fänge'


def make_document(tmp_path, text='RMT-8-Finge', duplicate=False):
    from docling_core.types.doc import BoundingBox, CoordOrigin, DocItemLabel, DoclingDocument, ProvenanceItem, Size
    path = tmp_path / 'prepared.pdf'
    with fitz.open() as pdf:
        page = pdf.new_page(width=300, height=100)
        page.insert_text((20, 35), text)
        pdf.save(path)
        box = page.get_text('words')[0][:4]
    doc = DoclingDocument(name='prepared')
    doc.add_page(page_no=1, size=Size(width=300, height=100))
    for _ in range(2 if duplicate else 1):
        doc.add_text(label=DocItemLabel.TEXT, text=text, prov=ProvenanceItem(
            page_no=1, charspan=(0, len(text)), bbox=BoundingBox(
                l=box[0],t=box[1],r=box[2],b=box[3],coord_origin=CoordOrigin.TOPLEFT)))
    receipt = {'method':recovery.NATIVE_TEXT_RECOVERY_POLICY,'languages':['deu','eng'],
               'source_pdf_sha256':'original-source-digest', 'unresolved':[], 'regions':[
                   {'page':1,'page_size':[300,100],'bbox':list(box),'candidates':[
                       {'original':'RMT-8-FÃ¤ng','decoded_hint':'RMT-8-Fäng',
                        'confirmed_text':'RMT-8-Fänge','bbox':list(box)}]}]}
    return doc, path, receipt


def test_geometry_aligned_repair_preserves_original_and_is_idempotent(tmp_path):
    doc, pdf, receipt = make_document(tmp_path)
    result = recovery.apply_native_text_recovery(doc, pdf, receipt)
    assert doc.texts[0].text == 'RMT-8-Fänge'
    assert doc.texts[0].orig == 'RMT-8-Finge'
    note = result['repairs'][0]
    assert note['original_native_token'] == 'RMT-8-FÃ¤ng'
    assert note['source_pdf_sha256'] == 'original-source-digest'
    assert note['dpi'] == [300,600] and note['languages'] == ['deu','eng']
    assert not recovery.apply_native_text_recovery(doc,pdf,receipt)['repairs']
    assert doc.texts[0].meta.corpus__native_text_recovery == [note]
    # Provenance is retained as structured metadata and never source prose.
    from pipeline.table_structure import export_source_markdown
    assert export_source_markdown(doc).strip() == 'RMT-8-Fänge'


@pytest.mark.parametrize('change', ['different_page_size','different_word','distant_box','ambiguous_owner','unconfirmed'])
def test_missing_or_conflicting_evidence_does_not_change_text(tmp_path,change):
    doc,pdf,receipt=make_document(tmp_path, 'RMT-9-Finge' if change=='different_word' else 'RMT-8-Finge',
                                 duplicate=change=='ambiguous_owner')
    before=[t.text for t in doc.texts]
    region=receipt['regions'][0]
    if change=='different_page_size':
        region['page_size'][0]=301
    if change=='distant_box':
        region['candidates'][0]['bbox']=[250,80,290,95]
    if change=='unconfirmed':
        region['candidates'][0]['confirmed_text']=None
    result=recovery.apply_native_text_recovery(doc,pdf,receipt)
    assert [t.text for t in doc.texts] == before and not result['repairs']


def test_unresolved_source_observation_survives_to_chunk_metadata(tmp_path):
    doc,pdf,receipt=make_document(tmp_path)
    candidate=receipt['regions'][0]['candidates'][0]
    candidate['confirmed_text']=None
    receipt['unresolved']=[{'page':1,'bbox':candidate['bbox'],'original':candidate['original'],
                            'reason':'regional_ocr_did_not_agree_with_decoded_accent'}]
    recovery.apply_native_text_recovery(doc,pdf,receipt)
    from pipeline.treatment_context import chunk_source_context,materialize_treatment_context
    metadata=chunk_source_context(doc.texts,materialize_treatment_context(doc))
    assert metadata['text_integrity'][0]['status']=='unresolved'
    assert metadata['text_integrity'][0]['original']=='RMT-8-FÃ¤ng'
    assert doc.texts[0].text=='RMT-8-Finge'
    recovery.apply_native_text_recovery(doc,pdf,receipt)
    assert len(doc.texts[0].meta.corpus__native_text_recovery)==1


def test_native_inspection_requires_two_raster_observations_and_retains_prep_receipt(tmp_path,monkeypatch):
    _,source,_=make_document(tmp_path,'RMT-8-FÃ¤ng')
    monkeypatch.setattr(recovery.shutil,'which',lambda name:'/fake/tesseract')
    monkeypatch.setattr(recovery,'native_text_recovery_producer',lambda languages:{'models':dict.fromkeys(languages,'digest')})
    calls=[]
    def ocr(args,**kwargs):
        calls.append((args,kwargs['input']))
        return SimpleNamespace(returncode=0,stdout=b'RMT-8-F\xc3\xa4nge')
    monkeypatch.setattr(recovery.subprocess,'run',ocr)
    from pipeline.scan import prepare_pdf
    output=tmp_path/'copied.pdf'
    outcome=prepare_pdf(source,{'needs_ocr':False,'tesseract_packs':['deu','eng']},output)
    receipt=outcome['native_text_recovery']
    assert receipt['candidate_count']==receipt['confirmed_count']==1
    assert receipt['source_pdf_sha256'] and output.read_bytes()==source.read_bytes()
    assert len(calls)==2 and calls[0][1]!=calls[1][1]
    assert all(args[-4:]==['-l','deu+eng','--psm','6'] for args,_ in calls)


def test_literal_damaged_byte_example_is_not_rewritten_from_decoding(tmp_path,monkeypatch):
    _,source,_=make_document(tmp_path,'RMT-8-FÃ¤ng')
    monkeypatch.setattr(recovery.shutil,'which',lambda name:'/fake/tesseract')
    monkeypatch.setattr(recovery,'native_text_recovery_producer',lambda languages:{})
    monkeypatch.setattr(recovery.subprocess,'run',lambda *args,**kwargs:SimpleNamespace(
        returncode=0,stdout='RMT-8-FÃ¤ng'.encode()))
    receipt=recovery.inspect_native_text_regions(source,['deu'])
    assert receipt['confirmed_count']==0 and receipt['unresolved']


def test_regional_work_is_bounded_and_excess_evidence_is_not_silently_lost(tmp_path,monkeypatch):
    source=tmp_path/'several-lines.pdf'
    with fitz.open() as pdf:
        page=pdf.new_page()
        for y in [30,80,130]:
            page.insert_text((20,y),'RMT-8-FÃ¤ng')
        pdf.save(source)
    monkeypatch.setattr(recovery,'_MAX_REGIONS_PER_DOCUMENT',1)
    monkeypatch.setattr(recovery.shutil,'which',lambda name:'/fake/tesseract')
    monkeypatch.setattr(recovery,'native_text_recovery_producer',lambda languages:{})
    calls=[]
    def ocr(*args,**kwargs):
        calls.append(args)
        return SimpleNamespace(returncode=0,stdout='RMT-8-Fänge'.encode())
    monkeypatch.setattr(recovery.subprocess,'run',ocr)
    result=recovery.inspect_native_text_regions(source,['deu'])
    assert result['candidate_count']==3 and result['confirmed_count']==1
    assert len(calls)==2
    assert len(result['unresolved'])==2
    assert all(r['reason']=='region_budget_exceeded' for r in result['unresolved'])


def test_source_fragment_records_independent_native_and_prepared_failures():
    assert 'RMT-8-Finge' in SOURCE['prepared_item']['text']
    candidate=next(c for c in SOURCE['native_region']['candidates'] if c.get('confirmed_text')=='RMT-8-Fänge')
    assert candidate['original']=='RMT-8-FÃ¤ng' and candidate['decoded_hint']=='RMT-8-Fäng'
    assert SOURCE['repair']['replacement']=='RMT-8-Fänge'
    assert SOURCE['repair']['dpi']==[300,600]


def test_external_original_page_corroborates_complete_accented_word(tmp_path):
    library=os.environ.get('CORPUS_LIBRARY_DIR')
    if not library:
        pytest.skip('set CORPUS_LIBRARY_DIR for read-only original-source OCR replay')
    source=Path(library)/'B/Boysen-Ennen1987.pdf'
    if not source.exists():
        pytest.skip('original source PDF unavailable')
    producer=recovery.native_text_recovery_producer(['deu','eng'])
    if not producer['available']:
        pytest.skip('source replay requires Tesseract deu+eng')
    selected=tmp_path/'source-page-12.pdf'
    with fitz.open(source) as full,fitz.open() as page:
        page.insert_pdf(full,from_page=11,to_page=11)
        page.save(selected)
    result=recovery.inspect_native_text_regions(selected,['deu','eng'])
    candidate=next(c for r in result['regions'] for c in r['candidates'] if c['original']=='RMT-8-FÃ¤ng')
    assert candidate['confirmed_text']=='RMT-8-Fänge'
    assert result['unresolved']  # Other source damage remains explicitly uncertain.
