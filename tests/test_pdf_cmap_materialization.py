"""Captured Chen source evidence through unchanged build/serve consumers.

The source driver and OCR results replay actual pinned observations. Token
counting and numerical embedding are instrumentation; all production artifact
serialization, HybridChunker and embedding-input selection run unchanged.
"""
import json
from types import SimpleNamespace

from mcpsrv import app
from mcpsrv.tools.chunks import get_chunks
from tests.test_embedding_updates import build  # noqa: F401 -- pytest fixture
from tests.test_pdf_cmap_recovery import CAPTURE, captured_driver, source_document


def test_corrected_chen_saved_text_reaches_chunks_embedding_and_bounded_fetch(build,monkeypatch,tmp_path):  # noqa: F811
    import docling.chunking
    from docling_core.transforms.chunker.tokenizer.base import BaseTokenizer
    from pipeline.chunking import chunk_text
    from pipeline.pdf_cmap_recovery import recover_pdf_cmaps
    from pipeline.table_structure import export_source_markdown
    from pipeline.stages import _run_quality_gates

    class LocalTokenizer(BaseTokenizer):
        def count_tokens(self,text): return len(text.split())
        def get_max_tokens(self): return 2000
        def get_tokenizer(self): return self.count_tokens

    real=docling.chunking.HybridChunker
    monkeypatch.setattr(docling.chunking,'HybridChunker',lambda **kw:real(tokenizer=LocalTokenizer(),**kw))
    pdf=captured_driver(monkeypatch,tmp_path)
    doc=source_document()
    report=recover_pdf_cmaps(doc,pdf)
    folder=build.document(CAPTURE['source_sha256'][:12])
    doc.save_as_json(folder/'docling_doc.json')
    markdown=export_source_markdown(doc)
    (folder/'text.json').write_text(json.dumps({'text':markdown,'source_text_integrity':{'pdf_cmap':report}}))
    chunk_text(folder/'text.json',chunks_output=folder/'chunks.json')
    artifact=json.loads((folder/'chunks.json').read_text())
    assert artifact['chunker']=='hybrid_chunker'
    chunks=artifact['chunks']
    text=' '.join(c['text'] for c in chunks)
    assert text.count('mg/m³')==3 and 'Ｒ =0． 596' in text and 'P =0． 001' in text
    assert 'HDYH' not in text and 'corpus__pdf_cmap' not in text and 'same_region_exact_font_cmap' not in text
    assert any(c.get('text_integrity') for c in chunks)
    flags=_run_quality_gates(folder)
    assert any(f['gate']=='source_text_integrity' and f['detail'].startswith('pdf_cmap: 6') for f in flags)
    calls=[]
    def record(texts):
        calls.append(list(texts))
        return [[float(len(t)),1.] for t in texts]
    build.backend.embed=record
    assert build.run(folder)==len(chunks)
    assert calls[-1]==[c['text'] for c in chunks]
    assert build.problem(folder) is None
    stored={r['metadata']['chunk_id']:r['text'] for r in build.table.to_arrow().to_pylist()}
    monkeypatch.setattr(app,'_INDEX',SimpleNamespace(papers={folder.name:{'hash_dir':str(folder)}}))
    before={name:(folder/name).read_bytes() for name in ('docling_doc.json','text.json','chunks.json')}
    for chunk in chunks:
        served=get_chunks(folder.name,chunk_ids=[chunk['chunk_id']])
        assert len(served)==1 and served[0]['text']==stored[chunk['chunk_id']]==chunk['text']
        assert len(json.dumps(served).encode())<64*1024
    assert all((folder/name).read_bytes()==raw for name,raw in before.items())
