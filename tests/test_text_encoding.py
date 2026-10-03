"""Encoding recovery needs independent source agreement, never just plausibility."""
import copy
import json
import os
from pathlib import Path

import pytest

from pipeline.text_encoding import repair_region, recover_text_encoding

FIXTURE = Path(__file__).parent / 'fixtures/text_encoding/source_regions.json'
REGIONS = json.loads(FIXTURE.read_text())['regions']


def glyphs(region):
    if 'native_glyphs' in region:
        return [{'c': row[0], 'bbox': row[1:]} for row in region['native_glyphs']]
    return [{'c': char, 'bbox': [0, 0, 0, 0]} for char in region['native_text']]


@pytest.mark.parametrize('region', REGIONS, ids=lambda r: f"{Path(r['source']).stem}-p{r['physical_page']}-{r['item_ref']}")
def test_fresh_source_regions_preserve_chinese_caption_names_and_accents(region):
    repaired, edits, problems = repair_region(region['original'], glyphs(region))
    assert edits and not problems
    if region['source'].startswith('F/'):
        if region['item_ref'] == '#/texts/4':
            assert 'Maria Luz Fernández de Puelles' in repaired
            assert 'Javier Jansá' in repaired
        else:
            assert 'Oceanográfico' in repaired and 'Español' in repaired
            assert 'Oceanografía' in repaired
    else:
        assert repaired == ' '.join(region['native_text'].split())
        if region['physical_page'] == 30 and region['item_ref'] == '#/texts/1':
            assert '5301±8525' in repaired
            assert 'ind./100m3' in repaired
            assert 'Nanomia bijuga' in repaired
        if region['physical_page'] == 84:
            assert repaired == '圖1.台灣西南海域採樣測站之位置圖。'
    # Every edit retains the exact pre-repair substring and applies once.
    for edit in edits:
        a, b = edit['charspan']
        assert region['original'][a:b] == edit['original']
    again, further, problems = repair_region(repaired, glyphs(region))
    assert again == repaired and further == [] and problems == []


@pytest.mark.parametrize('text', [
    '管水母種類組成及豐度', 'Français, déjà vu, naïve façade; São Tomé.',
    'A literal Ã is unchanged.', 'A tilde ~ and an acute ´ are spacing symbols.',
    'Pacific fish and scientific figures; Alvarifio as quoted in the source.',
])
def test_correct_multilingual_text_and_spacing_symbols_are_unchanged(text):
    chars = [{'c': char, 'bbox': [i*10, 0, i*10+8, 10]} for i, char in enumerate(text)]
    assert repair_region(text, chars) == (text, [], [])


def test_similar_native_region_cannot_change_a_number_or_word():
    native = '測試中文資料海洋科學' * 8 + ' value1974'
    damaged = ('測試中文資料海洋科學' * 8 + ' value1965').encode('big5').decode('latin1')
    text, edits, problems = repair_region(damaged, [{'c': c, 'bbox': [0,0,0,0]} for c in native])
    assert text == damaged and not edits
    assert problems[0]['reason'] == 'possible_big5_bytes_without_matching_source_unicode'


@pytest.mark.parametrize('source', ['管水母種類組成及豐度', '謝辭'])
def test_missing_source_cannot_authorize_decoding(source):
    damaged = source.encode('big5').decode('latin1')
    text, edits, problems = repair_region(damaged, [])
    assert text == damaged and not edits and problems


def test_source_mojibake_is_flagged_without_guessing_missing_letters():
    text = 'RMT-8-FÃ¤ng'
    assert repair_region(text, [])[0:2] == (text, [])
    assert repair_region(text, [])[2][0]['reason'] == 'possible_utf8_mojibake_requires_source_review'


def test_mixed_repair_and_unresolved_evidence_does_not_create_new_accent_warning():
    chars = [{'c':'a','bbox':[0,0,10,10]}, {'c':'´','bbox':[2,0,6,10]}]
    repaired, edits, problems = repair_region('a´ FranÃ§ais', chars)
    assert repaired == 'á FranÃ§ais' and len(edits) == 1
    assert repair_region(repaired, chars) == (repaired, [], problems)


def test_repeated_unresolved_metadata_is_idempotent(tmp_path):
    fitz = pytest.importorskip('fitz')
    from docling_core.types.doc import DoclingDocument, DocItemLabel, Size, BoundingBox, ProvenanceItem
    pdf = fitz.open(); pdf.new_page(width=100,height=100); path=tmp_path/'source.pdf'; pdf.save(path); pdf.close()
    doc=DoclingDocument(name='uncertain')
    doc.add_page(page_no=1,size=Size(width=100,height=100))
    item=doc.add_text(label=DocItemLabel.TEXT,text='FranÃ§ais',
        prov=ProvenanceItem(page_no=1,bbox=BoundingBox(l=0,t=0,r=100,b=100,coord_origin='TOPLEFT'),charspan=(0,9)))
    recover_text_encoding(doc,path)
    before=copy.deepcopy(doc.export_to_dict())
    recover_text_encoding(doc,path)
    assert doc.export_to_dict()==before
    assert len(item.meta.corpus__text_encoding)==1


@pytest.mark.parametrize('region', REGIONS, ids=lambda r: f"source-{Path(r['source']).stem}-{r['physical_page']}-{r['item_ref']}")
def test_optional_original_pdf_replay_preserves_orig_and_provenance(region):
    library=os.environ.get('CORPUS_LIBRARY_DIR')
    if not library:
        pytest.skip('Set CORPUS_LIBRARY_DIR for original PDF geometry replay')
    from docling_core.types.doc import DoclingDocument, DocItemLabel, Size, BoundingBox, ProvenanceItem
    doc=DoclingDocument(name='source-replay')
    w,h=region['page_size'];doc.add_page(page_no=region['physical_page'],size=Size(width=w,height=h))
    item=doc.add_text(label=DocItemLabel(region['label']),text=region['original'],orig=region['original'],
        prov=ProvenanceItem(page_no=region['physical_page'],bbox=BoundingBox(**region['bbox']),charspan=(0,len(region['original']))))
    report=recover_text_encoding(doc,Path(library)/region['source'])
    assert report['repairs'] and not report['unresolved']
    assert item.orig==region['original']
    assert item.text==repair_region(region['original'],glyphs(region))[0]
    before=copy.deepcopy(doc.export_to_dict())
    recover_text_encoding(doc,Path(library)/region['source'])
    assert doc.export_to_dict()==before
