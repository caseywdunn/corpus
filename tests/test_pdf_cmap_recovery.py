"""Exact font decoding, guarded source glyphs, and preserved raw observations.

Chen rows are actual source/Docling captures. Synthetic controls are explicitly
separate; the driver double supplies captured source geometry, never guessed text.
"""
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from pipeline import pdf_cmap_recovery as cmap

FIXTURE = Path(__file__).parent / 'fixtures/pdf_cmap'
CAPTURE = json.loads((FIXTURE / 'chen_regions.json').read_text())


def fonts(case):
    return {name: {**font, 'mapping': {int(k): v for k, v in font['mapping'].items()},
                   'observed_codes': {int(k): v for k, v in font.get('observed_codes', {}).items()}}
            for name, font in case['fonts'].items()}


def captured_repair(case, **kwargs):
    return cmap.repair_cmap_region(case['item']['text'], case['glyphs'], fonts(case),
        verify_sign=lambda g, _: case['sign_decisions'][json.dumps(g['bbox'])]['verified'],
        verify_digit=lambda g: case['digit_decisions'][json.dumps(g['bbox'])], **kwargs)


def source_document():
    from docling_core.types.doc import DoclingDocument, DocItemLabel, ProvenanceItem, Size
    doc = DoclingDocument(name='Chen two captured regions')
    for case in CAPTURE['cases']:
        doc.add_page(page_no=case['page'], size=Size(**case['page_size']))
        doc.add_text(label=DocItemLabel(case['item']['label']), text=case['item']['text'],
                     orig=case['item']['orig'], prov=ProvenanceItem.model_validate(case['item']['prov'][0]))
    return doc


def test_actual_chen_scalar_decoding_retains_raw_source_and_exact_geometry():
    abstract, results = CAPTURE['cases']
    text, edits, unresolved = captured_repair(abstract)
    assert text == abstract['candidate_text'] and text.count('mg/m³') == 3
    assert 'HDYH' in edits[0]['original'] and edits[0]['matched_scalars'] > 1000
    assert len(edits[0]['formatting']) == 3
    assert all(e['raster']['verified'] and e['source_bbox'] for e in edits[0]['formatting'])
    assert all(e['source_bbox'] and e['cmap_sha256'] and e['encoding_sha256'] for e in edits[0]['mappings'])
    assert len(unresolved) == 5  # Unclean hyphen crops remain raw and explicit.
    text, edits, unresolved = captured_repair(results)
    assert text == results['candidate_text']
    assert 'Ｒ =0． 596' in text and 'P =0． 001' in text
    assert len(unresolved) == 1
    assert all(u['reason'] == 'cmap_scientific_glyph_not_source_confirmed' for u in unresolved)


def synthetic_region(native='abcde'):
    codes = '!#$%&'
    glyphs = [{'text': c, 'font': 'test', 'bbox': [i*10,0,i*10+9,10],
               'size': 10, 'origin': [i*10,10], 'raised': False} for i,c in enumerate(native)]
    font = {'test': {'mapping': dict(zip(map(ord, codes), native)), 'font_xref': 1,
                    'tounicode_xref': 2, 'cmap_sha256': 'synthetic-control'}}
    return codes, glyphs, font


@pytest.mark.parametrize('native', ['±abcd','⁻abcd','³abcd','＋abcd','＝abcd','×abcd','µabcd','−abcd'])
def test_cmap_never_bypasses_scientific_symbol_guards(native):
    raw, glyphs, font = synthetic_region(native)
    text, edits, unresolved = cmap.repair_cmap_region(raw, glyphs, font, verify_sign=lambda *_: False)
    assert text == '!abcd' and edits
    assert unresolved[0]['native'] == native[0]
    assert unresolved[0]['reason'] == 'cmap_scientific_glyph_not_source_confirmed'


def test_normal_decoded_text_is_unchanged_even_with_ambiguous_unused_map():
    raw, glyphs, font = synthetic_region()
    assert cmap.repair_cmap_region('abcde', glyphs, font) == ('abcde', [], [])
    font['test']['error'] = 'conflicting_mapping'
    assert cmap.repair_cmap_region('abcde', glyphs, font) == ('abcde', [], [])
    text, edits, unresolved = cmap.repair_cmap_region(raw, glyphs, font)
    assert text == raw and not edits and unresolved[0]['reason'] == 'ambiguous_or_unsupported_source_font_map'


def test_legitimate_decoded_text_omission_does_not_become_cmap_uncertainty():
    raw,glyphs,font=synthetic_region()
    for observed in ('abce','ab-cde','ab?de'):
        assert cmap.repair_cmap_region(observed,glyphs,font)==(observed,[],[])
    # Same source font, but actual positive raw-code evidence survives a
    # trailing missing scalar. Keep the uncertainty without inventing text.
    text,edits,unresolved=cmap.repair_cmap_region(raw[:-1],glyphs,font)
    assert text==raw[:-1] and not edits
    assert unresolved[0]['reason']=='cmap_region_not_one_to_one'
    # A real correctly decoded Type1 region with one omitted word is also
    # unrelated to this producer, despite containing a CMap font.
    case=CAPTURE['cases'][0]
    native=''.join(g['text'] for g in case['glyphs'])
    damaged=native.replace('species','',1)
    assert damaged!=native
    assert cmap.repair_cmap_region(damaged,case['glyphs'],fonts(case))==(damaged,[],[])


def test_one_to_one_and_unique_inverse_are_required_no_fuzzy_replacement():
    raw, glyphs, font = synthetic_region()
    for changed in [raw+'x', raw[:-1], raw[:2]+'x'+raw[3:]]:
        text, edits, unresolved = cmap.repair_cmap_region(changed, glyphs, font)
        assert text == changed and not edits and unresolved
    font['test']['mapping'][65] = 'a'
    assert cmap.repair_cmap_region(raw, glyphs, font)[2][0]['reason'] == 'ambiguous_reverse_cmap'


@pytest.mark.parametrize('source', [
    '2 beginbfchar <21> <0061> <21> <0062> endbfchar',
    '1 beginbfrange <21> <22> [<0061> <0062>] endbfrange',
    '1 beginbfchar <21> <D800> endbfchar',
    '1 beginbfchar <0121> <0061> endbfchar',
    '1 beginbfchar <21> <00610062> endbfchar',
])
def test_conflicting_non_scalar_and_unimplemented_cmaps_refuse(source):
    with pytest.raises(ValueError):
        cmap.parse_scalar_cmap(source)


def test_scalar_cmap_and_docling_exact_quote_normalization():
    assert cmap.parse_scalar_cmap('2 beginbfrange <21> <23> <0061> <60> <60> <007A> endbfrange') == {33:'a',34:'b',35:'c',96:'z'}
    raw, glyphs, font = synthetic_region('abcze')
    font['test']['mapping'].pop(ord('%'))
    font['test']['mapping'][96] = 'z'
    font['test']['observed_codes'] = {96: "'"}
    assert cmap.repair_cmap_region("!#$'&", glyphs, font)[0] == 'abcze'
    assert not cmap.repair_cmap_region('!#$%&', glyphs, font)[1]


def test_geometry_needs_raised_numeric_unit_not_flags_affiliation_or_baseline():
    case = CAPTURE['cases'][0]
    chars = deepcopy([g for g in case['glyphs'] if not g['text'].isspace()])
    indices = [i for i in range(len(chars)) if cmap._raised_unit_digit(chars,i)]
    assert len(indices) == 3
    i = indices[0]
    assert not chars[i]['raised']  # No PyMuPDF flags=1 dependency.
    baseline = deepcopy(chars)
    baseline[i]['origin'][1] = baseline[i-1]['origin'][1]
    baseline[i]['raised'] = True
    assert not cmap._raised_unit_digit(baseline,i)
    normal_size = deepcopy(chars)
    normal_size[i]['size'] = normal_size[i-1]['size']
    assert not cmap._raised_unit_digit(normal_size,i)
    prose = deepcopy(chars)
    prose[i-3]['text'] = 'a'  # mg/m -> ma/m is no supported unit.
    assert not cmap._raised_unit_digit(prose,i)
    footnote = deepcopy(chars)
    footnote[i-1]['text'] = 'n'
    assert not cmap._raised_unit_digit(footnote,i)


def test_failed_digit_evidence_keeps_decoded_scalar_and_explicit_uncertainty():
    case = CAPTURE['cases'][0]
    text, edits, unresolved = cmap.repair_cmap_region(case['item']['text'],case['glyphs'],fonts(case),
        verify_sign=lambda *_: False, verify_digit=lambda _: {'verified':False,'reason':'digit_ocr_disagreement'})
    assert 'mg/m 3' in text and 'mg/m³' not in text
    assert not edits[0]['formatting']
    assert sum(u['reason']=='cmap_unit_exponent_not_source_confirmed' for u in unresolved) == 3


def captured_driver(monkeypatch, tmp_path):
    """Driver instrumentation: source rows are captured; PDF identity is synthetic."""
    import fitz
    import pipeline.scientific_text
    cases = {c['page']-1:c for c in CAPTURE['cases']}
    class Pages:
        def __len__(self): return 8
        def __getitem__(self, i): return SimpleNamespace(number=i)
        def __enter__(self): return self
        def __exit__(self, *args): pass
    monkeypatch.setattr(fitz,'open',lambda _:Pages())
    monkeypatch.setattr(pipeline.scientific_text,'_source_lines',lambda p:[cases[p.number]['glyphs']])
    monkeypatch.setattr(cmap,'source_font_maps',lambda _,p:fonts(cases[p.number]))
    monkeypatch.setattr(cmap,'_verify_source_sign',lambda p,g,_:cases[p.number]['sign_decisions'][json.dumps(g['bbox'])]['verified'])
    monkeypatch.setattr(cmap._DigitVerifier,'verify',lambda _,p,g:cases[p.number]['digit_decisions'][json.dumps(g['bbox'])])
    path=tmp_path/'driver-identity-only.pdf'
    path.write_bytes(b'SYNTHETIC DRIVER IDENTITY; geometry is retained source capture')
    return path


def test_driver_partial_repair_is_idempotent_and_raw_orig_survives(monkeypatch,tmp_path):
    from docling_core.types.doc import DoclingDocument
    pdf=captured_driver(monkeypatch,tmp_path)
    doc=source_document()
    raw=[i.orig for i in doc.texts]
    report=cmap.recover_pdf_cmaps(doc,pdf)
    once=doc.model_dump(mode='json')
    assert len(report['repairs'])==2 and len(report['unresolved'])==6
    assert [i.orig for i in doc.texts]==raw
    assert cmap.recover_pdf_cmaps(doc,pdf)==report
    assert doc.model_dump(mode='json')==once
    assert all(i.meta.corpus__pdf_cmap[0]['producer']['digit_ocr'] for i in doc.texts)
    saved=tmp_path/'saved-docling.json'
    doc.save_as_json(saved)
    reloaded=DoclingDocument.load_from_json(saved)
    assert cmap.recover_pdf_cmaps(reloaded,pdf)==report
    assert reloaded.model_dump(mode='json')==once


@pytest.mark.parametrize('move',['bbox','page'])
def test_moved_item_does_not_replay_old_region_proof(monkeypatch,tmp_path,move):
    pdf=captured_driver(monkeypatch,tmp_path)
    doc=source_document()
    cmap.recover_pdf_cmaps(doc,pdf)
    if move=='bbox':
        doc.texts[0].prov[0].bbox.l+=30
    else:
        doc.texts[0].prov[0].page_no=8
    report=cmap.recover_pdf_cmaps(doc,pdf)
    assert all(r['item_ref']!=doc.texts[0].self_ref for r in report['repairs'])


def test_digit_ocr_budgets_cache_and_disagreement(monkeypatch):
    class Pix:
        def tobytes(self,_): return b'instrumented raster'
    rendered=[]
    page=SimpleNamespace(number=0,get_pixmap=lambda **kw:rendered.append(kw) or Pix())
    glyph={'text':'3','bbox':[10,10,14,16]}
    calls=[]
    monkeypatch.setattr(cmap.shutil,'which',lambda _: '/instrumented/tesseract')
    monkeypatch.setattr(cmap.subprocess,'run',lambda args,**kw:calls.append((args,kw)) or SimpleNamespace(stdout=b'3\n',returncode=0))
    verifier=cmap._DigitVerifier()
    first=verifier.verify(page,glyph)
    assert first['verified'] and len(calls)==2
    assert verifier.verify(page,glyph)==first and len(calls)==2 and len(rendered)==1
    assert all(0<kw['timeout']<=3 for _,kw in calls)
    assert all('tessedit_char_whitelist' not in ' '.join(args) for args,_ in calls)
    verifier.page_counts[0]=verifier.MAX_PAGE
    assert verifier.verify(page,{'text':'2','bbox':[20,20,24,26]})['reason']=='digit_ocr_budget_exhausted'
    assert len(calls)==2
    verifier=cmap._DigitVerifier()
    assert verifier.verify(page,{'text':'3','bbox':[0,0,1000,1000]})['reason']=='digit_ocr_pixel_budget_exhausted'
    assert verifier.verify(page,{'text':'3','bbox':[0,0,float('inf'),1]})['reason']=='digit_ocr_invalid_source_box'
    verifier=cmap._DigitVerifier()
    verifier.total_pixels=verifier.MAX_TOTAL_PIXELS
    assert verifier.verify(page,glyph)['reason']=='digit_ocr_pixel_budget_exhausted'
    verifier=cmap._DigitVerifier()
    verifier.cache={i: {} for i in range(verifier.MAX_DOCUMENT)}
    assert verifier.verify(page,glyph)['reason']=='digit_ocr_budget_exhausted'
    verifier=cmap._DigitVerifier()
    monkeypatch.setattr(cmap.subprocess,'run',lambda *a,**kw:SimpleNamespace(stdout=b'8\n',returncode=0))
    assert not verifier.verify(page,glyph)['verified']
    verifier=cmap._DigitVerifier()
    verifier.started-=verifier.DOCUMENT_SECONDS+1
    assert verifier.verify(page,glyph)['reason']=='digit_ocr_budget_exhausted'


def test_cmap_policy_and_models_invalidate_extraction_and_descendants(monkeypatch):
    from pipeline.build_inputs import config_fingerprints
    before=config_fingerprints({},panel_mode='ocr')
    monkeypatch.setattr(cmap,'PDF_CMAP_POLICY','changed-policy')
    after=config_fingerprints({},panel_mode='ocr')
    for stage in ('docling_extraction','text_chunking','figure_materialization','taxa_and_lexicon_extraction','figure_crossref'):
        assert before[stage]!=after[stage]
    import pipeline.source_spaces
    monkeypatch.setattr(pipeline.source_spaces,'source_spacing_producer',lambda:{'available':True,'traineddata_sha256':'different-model'})
    model=config_fingerprints({},panel_mode='ocr')
    assert model['docling_extraction']['extraction.pdf_cmap_producer']!=after['docling_extraction']['extraction.pdf_cmap_producer']


def test_original_pdf_replay_preserves_named_values_and_second_pass():
    root=os.environ.get('CORPUS_LIBRARY_ROOT')
    if not root:
        pytest.skip('optional original-PDF replay: set CORPUS_LIBRARY_ROOT to library repository')
    pdf=Path(root)/CAPTURE['source_pdf']
    assert hashlib.sha256(pdf.read_bytes()).hexdigest()==CAPTURE['source_sha256']
    producer=cmap.pdf_cmap_producer()
    if not producer['digit_ocr']['available']:
        pytest.skip('original digit replay needs installed Tesseract eng')
    doc=source_document()
    report=cmap.recover_pdf_cmaps(doc,pdf)
    assert doc.texts[0].text.count('mg/m³')==3
    assert 'Ｒ =0． 596' in doc.texts[1].text
    assert 'P =0． 001' in doc.texts[1].text
    assert len(report['unresolved'])==6  # Reviewed, independent raw hyphens.
    saved=doc.model_dump(mode='json')
    assert cmap.recover_pdf_cmaps(doc,pdf)==report
    assert doc.model_dump(mode='json')==saved


@pytest.mark.parametrize('failure',['null','nonstream','dangling','malformed'])
def test_malformed_optional_font_maps_do_not_fail_readable_extraction(failure):
    class PDF:
        def xref_get_key(self,*args): return ('xref','2 0 R')
        def xref_stream(self,_):
            if failure=='dangling': raise RuntimeError('xref missing')
            return {'null':None,'nonstream':'not bytes','malformed':b'invalid cmap'}[failure]
    page=SimpleNamespace(get_fonts=lambda **kw:[(1,'n/a','Type1','test','F0','',0)])
    mapped=cmap.source_font_maps(PDF(),page)
    assert mapped['test']['error'].startswith('unreadable_or_unsupported_source_map:')
    _,glyphs,_=synthetic_region()
    assert cmap.repair_cmap_region('abcde',glyphs,mapped)==('abcde',[],[])
    assert cmap.repair_cmap_region('abce',glyphs,mapped)==('abce',[],[])


def test_pinned_source_crops_keep_the_scientific_guard_decisions():
    import fitz
    manifest=json.loads((FIXTURE/'manifest.json').read_text())
    for name,expected in manifest['hashes'].items():
        assert hashlib.sha256((FIXTURE/name).read_bytes()).hexdigest()==expected
    for case in CAPTURE['cases']:
        for box,decision in case['sign_decisions'].items():
            pix=fitz.Pixmap(str(FIXTURE/decision['raster']))
            page=SimpleNamespace(get_pixmap=lambda **kw:pix)
            assert bool(cmap._verify_source_sign(page,{'bbox':json.loads(box)},decision['native']))==decision['verified']
