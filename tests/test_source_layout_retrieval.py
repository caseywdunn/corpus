"""Named #304 source fragments through serialization, real chunking and serving.

The retained layout fixtures deliberately shorten some non-boundary prose.
These tests cover the captured boundary text, not new full-document extraction
or embedding quality. Only the token counter is local; Docling's real splitting,
merging and serialization code runs without downloading a model.
"""
import json
from types import SimpleNamespace

import pytest

from mcpsrv import app
from mcpsrv.tools.chunks import get_chunks, get_chunks_by_section
from pipeline.source_layout import recover_panel_caption_roles, repair_reading_order
from tests.test_source_text_integrity import CASES, document


def normalized(text):
    return ' '.join(text.split())


def refs(row):
    return {item['item_ref'] for item in row['source_items']}


def replay(tmp_path, monkeypatch, case, *, max_tokens=2000, repair=True, omit_heading=False):
    pytest.importorskip('docling')
    import docling.chunking
    from docling_core.transforms.chunker.tokenizer.base import BaseTokenizer
    from docling_core.types.doc import DoclingDocument
    from pipeline.chunking import chunk_text
    from pipeline.table_structure import export_source_markdown

    class LocalTokenizer(BaseTokenizer):
        def count_tokens(self, text):
            return len(text.split())

        def get_max_tokens(self):
            return max_tokens

        def get_tokenizer(self):
            # semchunk needs the callable when an actual item is split.
            return self.count_tokens

    real = docling.chunking.HybridChunker
    monkeypatch.setattr(docling.chunking, 'HybridChunker',
                        lambda **kw: real(tokenizer=LocalTokenizer(), **kw))
    doc = document(case)
    if omit_heading:
        doc.body.children = [r for r in doc.body.children
                             if not r.resolve(doc).text.startswith('Physalia minuta Church')]
    recover_panel_caption_roles(doc)
    order = repair_reading_order(doc) if repair else None
    source = tmp_path / 'docling_doc.json'
    doc.save_as_json(source)
    loaded = DoclingDocument.load_from_json(source)
    assert [r.cref for r in loaded.body.children] == [r.cref for r in doc.body.children]
    assert [(t.self_ref, t.text, t.label) for t in loaded.texts] == [
        (t.self_ref, t.text, t.label) for t in doc.texts]
    (tmp_path / 'text.json').write_text(json.dumps({'text': export_source_markdown(loaded)}))
    output = tmp_path / 'chunks.json'
    chunk_text(tmp_path / 'text.json', chunks_output=output, docling_doc_file=source)
    artifact = json.loads(output.read_text())
    assert artifact['chunker'] == 'hybrid_chunker'  # A fallback is not acceptance.
    sha = CASES[case]['sha256'][:12]
    monkeypatch.setattr(app, '_INDEX', SimpleNamespace(papers={sha: {'hash_dir': str(tmp_path)}}))
    # Real metadata-first discovery, bounded to this small source fragment.
    before = {p: p.read_bytes() for p in (source, output)}
    scan = get_chunks_by_section(sha, limit=50, with_text=False)
    assert len(scan) == artifact['total_chunks'] < 50
    assert all('text' not in row and 'len_chars' in row for row in scan)
    rows = get_chunks(sha, [row['chunk_id'] for row in scan])
    assert [row['text'] for row in rows] == [row['text'] for row in artifact['chunks']]
    for row, stored in zip(rows, artifact['chunks']):
        scope = row['context_projection']['records']['source_items']
        assert scope['available'] == len(stored['source_items'])
        assert scope['returned'] == len(row['source_items'])
        assert row['source_items'] == stored['source_items'][:scope['returned']]
        if scope['returned'] < scope['available']:
            assert scope['truncated'] and row['context_projection']['truncated']
    assert all(p.read_bytes() == value for p, value in before.items())
    return loaded, artifact, rows, sha, order


@pytest.mark.parametrize('max_tokens', [35, 60, 2000])
def test_church_continuations_retain_minuta_across_save_split_merge_and_serving(tmp_path, monkeypatch, max_tokens):
    doc, artifact, rows, sha, order = replay(tmp_path, monkeypatch, 'Church_etal2025-7-8', max_tokens=max_tokens)
    assert [page['page'] for page in order['reordered_pages']] == [8]
    continuation = next(t for t in doc.texts if t.text.startswith('principal tentacles'))
    continuation_rows = [row for row in rows if continuation.self_ref in refs(row)]
    assert continuation_rows
    if max_tokens == 35:
        assert len(continuation_rows) > 1  # Exercise actual item splitting.
    for row in continuation_rows:
        assert row['treatment_context']['name'] == 'Physalia minuta'
        assert row['treatment_context']['heading_page'] == 7
        assert row['section_type'] == 'diagnosis'
        assert row['section_class'] == 'description'
        assert any(p['page'] == 8 and p['item_ref'] == continuation.self_ref for p in row['source_items'])
        assert not any('Physalia physalis' in heading for heading in row['headings'])
    continuation_text = normalized(' '.join(row['text'] for row in continuation_rows))
    assert normalized(continuation.text) in continuation_text
    if max_tokens == 2000:
        assert 'Palpons and principal tentacles' in continuation_text

    diagnosis = get_chunks_by_section(sha, treatment_name='Physalia minuta', section_type='diagnosis', limit=50)
    assert {r['chunk_id'] for r in continuation_rows} <= {r['chunk_id'] for r in diagnosis}
    assert get_chunks_by_section(sha, treatment_name='Physalia physalis', section_type='diagnosis', limit=50) == []
    neighbor = get_chunks(sha, treatment_name='Physalia physalis')
    assert neighbor and any('Described as Holothuria physalis' in r['text'] for r in neighbor)
    assert all(continuation.self_ref not in refs(r) for r in neighbor)
    assert max(rows.index(row) for row in continuation_rows) < min(rows.index(row) for row in neighbor)

    captions = {t.self_ref for t in doc.texts if t.label.value == 'caption'}
    caption_rows = [row for row in rows if refs(row) & captions]
    assert caption_rows and set().union(*(refs(row) for row in caption_rows)) == captions
    assert all(refs(row) <= captions and row['treatment_context']['status'] == 'unknown'
               and row['section_type'] is None for row in caption_rows)
    assert not any(refs(row) & captions for row in diagnosis)
    assert artifact['treatment_context_policy']


def test_church_missing_heading_is_unknown_despite_literal_name_in_continuation(tmp_path, monkeypatch):
    doc, _, rows, sha, _ = replay(tmp_path, monkeypatch, 'Church_etal2025-7-8', omit_heading=True)
    continuation = next(t for t in doc.texts if t.text.startswith('principal tentacles'))
    selected = [row for row in rows if continuation.self_ref in refs(row)]
    assert selected and all(row['treatment_context'] == {'status': 'unknown', 'name': None} for row in selected)
    assert any('P. minuta' in row['text'] for row in selected)
    assert get_chunks(sha, treatment_name='Physalia minuta') == []
    assert get_chunks(sha, treatment_name='Physalia physalis')  # Adjacent heading still works.


def test_haddock_column_sentences_remain_contiguous_and_captions_stay_distinct(tmp_path, monkeypatch):
    doc, artifact, rows, _, order = replay(tmp_path, monkeypatch, 'Haddock_etal2005-1-1')
    assert [page['page'] for page in order['reordered_pages']] == [1]
    sentences = ['remnants of which were found inside two of our specimens.',
                 'spanning yellow to red (583, 620, and 680 nm)']
    for sentence in sentences:
        containing = [row for row in rows if sentence in normalized(row['text'])]
        assert len(containing) == 1
        assert all(p['page'] == 1 for p in containing[0]['source_items'])
    captions = {t.self_ref for t in doc.texts if t.label.value == 'caption'}
    caption_rows = [row for row in rows if refs(row) & captions]
    assert len(caption_rows) == 1
    assert refs(caption_rows[0]) == captions
    assert caption_rows[0]['text'].startswith('Fig. 1.')
    assert 'inside two of' not in caption_rows[0]['text']
    assert all('Scale bar' not in row['text'] for row in rows if not refs(row) & captions)
    first = next(row for row in rows if 'inside two of our specimens.' in normalized(row['text']))
    stored = next(row for row in artifact['chunks'] if row['chunk_id'] == first['chunk_id'])
    assert len(stored['source_items']) == 9
    assert first['context_projection']['records']['source_items'] == {
        'available': 9, 'returned': 6, 'truncated': True}
    # A mentioned genus is not evidence of an enclosing species treatment.
    assert all(row['treatment_context'] == {'status': 'unknown', 'name': None} for row in rows)


@pytest.mark.parametrize('case', ['Church_etal2025-7-8', 'Haddock_etal2005-1-1'])
def test_unrepaired_captured_order_reproduces_named_failure_through_same_chain(tmp_path, monkeypatch, case):
    doc, _, rows, _, _ = replay(tmp_path, monkeypatch, case, repair=False)
    if case.startswith('Church'):
        continuation = next(t for t in doc.texts if t.text.startswith('principal tentacles'))
        selected = [row for row in rows if continuation.self_ref in refs(row)]
        assert selected and all(row['treatment_context']['name'] == 'Physalia physalis' for row in selected)
    else:
        assert not any('inside two of our specimens.' in normalized(row['text']) for row in rows)
        assert not any('yellow to red (583, 620' in normalized(row['text']) for row in rows)
