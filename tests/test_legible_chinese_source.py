"""A source-verified Chinese title must remain unchanged by encoding recovery.

This region is reconstructed from retained readable output and native geometry.
It is deliberately not presented as a fresh Docling capture or whole-paper test.
"""
import copy
import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from mcpsrv import app
from mcpsrv.tools.chunks import get_chunks, get_chunks_by_section
from pipeline.text_encoding import repair_region, recover_text_encoding


FIXTURE = Path(__file__).parent / 'fixtures/text_encoding/already_legible_chinese.json'
REGION = json.loads(FIXTURE.read_text())


def test_retained_legible_chinese_title_matches_native_source_without_repair():
    chars = [{'c': row[0], 'bbox': row[1:]} for row in REGION['native_glyphs']]
    assert ''.join(c['c'] for c in chars) == REGION['native_text'] == REGION['observed_text']
    assert repair_region(REGION['observed_text'], chars) == (REGION['observed_text'], [], [])


def test_optional_source_region_encoding_serialization_and_serving_remain_unchanged(tmp_path, monkeypatch):
    library = os.environ.get('CORPUS_LIBRARY_DIR')
    if not library:
        pytest.skip('Set CORPUS_LIBRARY_DIR for original Liu2012 source-region replay')
    import docling.chunking
    from docling_core.transforms.chunker.tokenizer.base import BaseTokenizer
    from docling_core.types.doc import DoclingDocument, DocItemLabel, Size, BoundingBox, ProvenanceItem
    from pipeline.chunking import chunk_text
    from pipeline.table_structure import export_source_markdown

    source = Path(library) / REGION['source']
    assert hashlib.sha256(source.read_bytes()).hexdigest() == REGION['source_sha256']
    doc = DoclingDocument(name='Liu2012-source-title-reconstruction')
    width, height = REGION['page_size']
    doc.add_page(page_no=REGION['physical_page'], size=Size(width=width, height=height))
    item = doc.add_text(label=DocItemLabel.TEXT, text=REGION['observed_text'], orig=REGION['observed_text'],
                        prov=ProvenanceItem(page_no=REGION['physical_page'], bbox=BoundingBox(**REGION['bbox']),
                                            charspan=(0, len(REGION['observed_text']))))
    before = copy.deepcopy(doc.export_to_dict())
    report = recover_text_encoding(doc, source)
    assert report['repairs'] == report['unresolved'] == []
    assert doc.export_to_dict() == before
    recover_text_encoding(doc, source)
    assert doc.export_to_dict() == before
    assert item.orig == item.text == REGION['observed_text']

    class LocalTokenizer(BaseTokenizer):
        def count_tokens(self, text):
            return len(text.split())

        def get_max_tokens(self):
            return 1000

        def get_tokenizer(self):
            return self.count_tokens

    real = docling.chunking.HybridChunker
    monkeypatch.setattr(docling.chunking, 'HybridChunker', lambda **kw: real(tokenizer=LocalTokenizer(), **kw))
    serialized = tmp_path / 'docling_doc.json'
    doc.save_as_json(serialized)
    (tmp_path / 'text.json').write_text(json.dumps({'text': export_source_markdown(doc)}))
    output = tmp_path / 'chunks.json'
    chunk_text(tmp_path / 'text.json', chunks_output=output, docling_doc_file=serialized)
    chunks = json.loads(output.read_text())
    assert chunks['chunker'] == 'hybrid_chunker'
    assert len(chunks['chunks']) == 1
    sha = REGION['source_sha256'][:12]
    monkeypatch.setattr(app, '_INDEX', SimpleNamespace(papers={sha: {'hash_dir': str(tmp_path)}}))
    saved = output.read_bytes()
    scan = get_chunks_by_section(sha, limit=1, with_text=False)
    assert len(scan) == 1 and 'text' not in scan[0]
    served = get_chunks(sha, [scan[0]['chunk_id']])[0]
    assert served['text'] == REGION['observed_text']
    assert served['source_items'][0]['page'] == 2
    assert 'text_integrity' not in served  # No fabricated repair receipt.
    assert output.read_bytes() == saved
