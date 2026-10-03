"""Source-confirmed table/image ownership and conservative abstention (#320)."""
import json
from pathlib import Path

import pytest
from docling_core.types.doc import DoclingDocument

from pipeline.source_layout import recover_table_caption_owners
from pipeline.table_structure import export_source_markdown


def source():
    path = Path(__file__).parent / 'fixtures/table_structure/Hosia2024-caption-owner.json'
    return DoclingDocument.model_validate(json.loads(path.read_text())['document'])


@pytest.mark.parametrize('top_left', [False, True])
def test_printed_table_caption_survives_save_and_chunk(tmp_path, monkeypatch, top_left):
    import docling.chunking
    from pipeline.chunking import chunk_text
    from tests.test_table_structure import WordTokenizer

    doc = source()
    if top_left:
        for item in [*doc.texts, *doc.tables, *doc.pictures]:
            for prov in item.prov:
                prov.bbox = prov.bbox.to_top_left_origin(doc.pages[prov.page_no].size.height)
    table, picture, caption = doc.tables[0], doc.pictures[0], doc.texts[0]
    assert not table.captions and caption.parent.cref == picture.self_ref
    original_text = [(t.text, t.orig, t.prov) for t in doc.texts]
    cells = table.data.model_dump()
    assert len(recover_table_caption_owners(doc)) == 1
    assert caption.parent.cref == table.self_ref
    assert [r.cref for r in table.captions] == [caption.self_ref]
    assert not picture.captions
    assert [r.cref for r in picture.children] == [t.self_ref for t in doc.texts[1:]]
    assert [(t.text, t.orig, t.prov) for t in doc.texts] == original_text
    assert table.data.model_dump() == cells
    assert recover_table_caption_owners(doc) == []
    path = tmp_path / 'docling_doc.json'
    doc.save_as_json(path)
    loaded = DoclingDocument.load_from_json(path)
    assert loaded.texts[0].parent.cref == loaded.tables[0].self_ref
    markdown = export_source_markdown(loaded)
    assert markdown.count(caption.text) == 1
    (tmp_path / 'text.json').write_text(json.dumps({'text': markdown}))
    real = docling.chunking.HybridChunker
    monkeypatch.setattr(docling.chunking, 'HybridChunker',
                        lambda **kw: real(tokenizer=WordTokenizer(limit=80), **kw))
    output = tmp_path / 'chunks.json'
    chunk_text(tmp_path / 'text.json', chunks_output=output)
    result = json.loads(output.read_text())
    assert result['chunker'] == 'hybrid_chunker'
    assert all(c['text'].strip() for c in result['chunks'])
    table_chunks = [c for c in result['chunks'] if c.get('tables')]
    assert table_chunks
    assert any(caption.text in c['text'] for c in table_chunks)


@pytest.mark.parametrize('case', [
    'figure_label', 'unnumbered', 'duplicate_table', 'existing_caption',
    'far_table', 'other_page', 'multipage_table', 'near_picture', 'intervening_text',
    'competing_caption', 'shifted_column',
])
def test_ambiguous_or_conflicting_ownership_is_unchanged(case):
    doc = source()
    table, picture, caption = doc.tables[0], doc.pictures[0], doc.texts[0]
    if case == 'figure_label':
        caption.text = caption.text.replace('TABLE 1', 'FIGURE 1')
    elif case == 'unnumbered':
        caption.text = 'Table of diagnostic features'
    elif case == 'duplicate_table':
        doc.add_table(data=table.data.model_copy(deep=True), prov=table.prov[0].model_copy(deep=True))
    elif case == 'existing_caption':
        table.captions.append(doc.texts[1].get_ref())
    elif case == 'far_table':
        table.prov[0].bbox.t -= 30
    elif case == 'other_page':
        table.prov[0].page_no = 6
    elif case == 'multipage_table':
        table.prov.append(table.prov[0].model_copy(deep=True))
    elif case == 'near_picture':
        picture.prov[0].bbox.t = caption.prov[0].bbox.b - 2
    elif case == 'intervening_text':
        prov = caption.prov[0].model_copy(deep=True)
        prov.bbox.t = caption.prov[0].bbox.b - 2
        prov.bbox.b = prov.bbox.t - 2
        doc.add_text(label='text', text='Intervening source text', prov=prov)
    elif case == 'competing_caption':
        other = doc.add_text(label='caption', text='TABLE 2 A competing caption',
                             prov=caption.prov[0].model_copy(deep=True), parent=picture)
        picture.captions.append(other.get_ref())
    elif case == 'shifted_column':
        table.prov[0].bbox.l += 100
    before = doc.model_dump()
    assert recover_table_caption_owners(doc) == []
    assert doc.model_dump() == before


def test_new_owner_policy_invalidates_extraction_and_consumers(monkeypatch):
    from pipeline import source_layout
    from pipeline.build_inputs import config_fingerprints

    current = config_fingerprints({}, panel_mode='off', surname_producer={})
    monkeypatch.setattr(source_layout, 'SOURCE_LAYOUT_POLICY', 'source_columns_v1')
    previous = config_fingerprints({}, panel_mode='off', surname_producer={})
    changed = {stage for stage in current if current[stage] != previous[stage]}
    assert changed == {'docling_extraction', 'text_chunking', 'taxa_and_lexicon_extraction',
                       'figure_materialization', 'figure_crossref'}
