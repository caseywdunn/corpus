"""A serialized table caption must not create a second blank passage."""
import json
from pathlib import Path

from docling_core.types.doc import DoclingDocument

from pipeline.chunking import chunk_text
from pipeline.table_structure import SourceTableSerializerProvider, export_source_markdown
from tests.test_table_structure import WordTokenizer


def test_source_table_caption_survives_without_blank_placeholder(tmp_path, monkeypatch):
    import docling.chunking
    fixture = Path(__file__).parent / 'fixtures/table_structure/Larimer1962-caption-placeholder.json'
    doc = DoclingDocument.model_validate(json.loads(fixture.read_text())['document'])
    real = docling.chunking.HybridChunker
    upstream = list(real(tokenizer=WordTokenizer(limit=30), merge_peers=False,
                         serializer_provider=SourceTableSerializerProvider()).chunk(doc))
    assert any(not c.text.strip() for c in upstream)  # Reproduce the producer defect.
    source = tmp_path / 'docling_doc.json'
    doc.save_as_json(source)
    (tmp_path / 'text.json').write_text(json.dumps({'text': export_source_markdown(doc)}))
    monkeypatch.setattr(docling.chunking, 'HybridChunker',
                        lambda **kw: real(tokenizer=WordTokenizer(limit=30), **kw))
    out = tmp_path / 'chunks.json'
    chunk_text(tmp_path / 'text.json', chunks_output=out)
    chunks = json.loads(out.read_text())['chunks']
    assert chunks and all(c['text'].strip() for c in chunks)
    text = '\n'.join(c['text'] for c in chunks)
    assert text.count(doc.texts[0].text) == 1
    # Removal changes neither the original populated passages nor their order.
    assert ''.join(c['text'] for c in chunks) == ''.join(c.text for c in upstream)
    assert [c['chunk_id'] for c in chunks] == [f'chunk_{i}' for i in range(len(chunks))]
    assert any(c.get('tables') for c in chunks)
