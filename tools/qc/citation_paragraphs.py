#!/usr/bin/env python3
"""Replay saved complete TEIs through citation materialization, without extraction.

Inputs are read-only build roots. A new output directory receives regenerated
paragraphs, actual reference authority, bounded tool replay and a measurement
receipt. Signature counts concern the supplied population, never all possible
text errors or fresh corpus extraction. No OCR, Grobid, embedding or model call.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import subprocess
from types import SimpleNamespace

from bib.authority import create_schema, phase1_corpus_papers, phase2_references
from mcpsrv import app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import get_excerpts_citing, get_intext_citations
from pipeline.citation_spans import CITATION_SPAN_POLICY
from pipeline.grobid_client import parse_tei_intext_citations, parse_tei_references
from pipeline.intext_citations import _matching_source_pdf

SIGNATURE = r'\(\d{4}[A-Za-z][A-Za-z.\- ]{1,30}\('
DETECTOR = re.compile(SIGNATURE)
REPO = Path(__file__).resolve().parents[2]


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def run(inputs, output):
    output = output.resolve()
    inputs = [root.resolve() for root in inputs]
    if any(output.is_relative_to(root) or root.is_relative_to(output) for root in inputs):
        raise ValueError('Output must be separate from every read-only input root')
    output.mkdir(parents=True, exist_ok=False)
    build = output / 'build'
    papers, cases = {}, []
    for root in inputs:
        for xml_path in sorted((root / 'documents').glob('*/grobid.tei.xml')):
            folder = xml_path.parent
            document = folder.name
            if document in papers:
                raise ValueError(f'Duplicate document: {document}; select one capture explicitly')
            xml = xml_path.read_text(encoding='utf-8')
            pdf = _matching_source_pdf(folder, xml)
            if pdf is None:
                raise ValueError(f'No matching prepared PDF/TEI provenance: {folder}')
            baseline = parse_tei_intext_citations(xml)
            rebuilt = parse_tei_intext_citations(xml, pdf_path=pdf)
            assert len(baseline['paragraphs']) == len(rebuilt['paragraphs'])
            references = parse_tei_references(xml)
            previous_refs = json.loads((folder / 'references.json').read_text())['references']
            assert references == previous_refs, f'Saved references differ from current parse: {folder}'
            metadata = json.loads((folder / 'metadata.json').read_text())
            destination = build / 'documents' / document
            destination.mkdir(parents=True)
            dump(destination / 'metadata.json', metadata)
            dump(destination / 'references.json', {'references': references})
            dump(destination / 'intext_citations.json', rebuilt)
            papers[document] = {'hash_dir': str(destination), 'title': metadata.get('title') or ''}
            before_hits = [i for i, text in enumerate(baseline['paragraphs']) if DETECTOR.search(text)]
            after_hits = [i for i, text in enumerate(rebuilt['paragraphs']) if DETECTOR.search(text)]
            changed = [i for i, (a, b) in enumerate(zip(baseline['paragraphs'], rebuilt['paragraphs'])) if a != b]
            record = {'document': document, 'input_root': str(root), 'title': metadata.get('title'),
                      'input_sha256': {name: sha(folder / name) for name in
                                       ('grobid.tei.xml', 'grobid.tei.xml.provenance.json', 'processed.pdf', 'metadata.json', 'references.json')},
                      'paragraphs': len(rebuilt['paragraphs']), 'reference_observations': len(references),
                      'before_hit_indices': before_hits, 'after_hit_indices': after_hits,
                      'changed_paragraphs': len(changed),
                      'text_sources': dict(Counter(row['text_source'] for row in rebuilt['citations'])),
                      'link_statuses': dict(Counter(row['validation_status'] for row in rebuilt['citations'])),
                      'detector_cases': [{'index': i, 'before': baseline['paragraphs'][i], 'after': rebuilt['paragraphs'][i]}
                                         for i in sorted(set(before_hits + after_hits))]}
            cases.append(record)
    if not cases:
        raise ValueError('No saved TEI inputs')
    db = build / 'biblio_authority.sqlite'
    with sqlite3.connect(db) as conn:
        create_schema(conn)
        phase1 = phase1_corpus_papers(conn, build)
        phase2 = phase2_references(conn, build)
        for case in cases:
            folder = build / 'documents' / case['document']
            refs = json.loads((folder / 'references.json').read_text())['references']
            observed = conn.execute('SELECT grobid_xml_id,raw_citation,title,year FROM reference_observations WHERE citing_corpus_hash=? ORDER BY ordinal', (case['document'],)).fetchall()
            assert observed == [(r['xml_id'], r['raw'], r['title'], r['year']) for r in refs]
        before = list(conn.iterdump())
        assert phase1_corpus_papers(conn, build) == 0
        assert phase2_references(conn, build) == (0, 0)
        assert list(conn.iterdump()) == before
        target_maps = {case['document']: dict(conn.execute('SELECT grobid_xml_id,cited_work_id FROM citations WHERE citing_corpus_hash=?', (case['document'],))) for case in cases}
    biblio = BiblioAuthority(db)
    old_index = app._INDEX
    app._INDEX = SimpleNamespace(biblio_db=biblio, papers=papers)
    excerpt_cache = {}
    try:
        for case in cases:
            folder = build / 'documents' / case['document']
            rebuilt = json.loads((folder / 'intext_citations.json').read_text())
            targets = target_maps[case['document']]
            for measured in case['detector_cases']:
                indices = [i for i, row in enumerate(rebuilt['citations']) if row['para_index'] == measured['index']]
                selected = indices[0]
                intext = get_intext_citations(case['document'], limit=1, offset=selected)
                assert intext['paragraphs'] == [measured['after']]
                assert intext['citations'] == [dict(rebuilt['citations'][selected], para_index=0)]
                measured['intext_replayed'] = True
                measured['excerpt_replayed'] = False
                for index in indices:
                    marker = rebuilt['citations'][index]
                    target = targets.get((marker.get('target_xml_id') or '').lstrip('#'))
                    if target is None:
                        continue
                    if target not in excerpt_cache:
                        rows, offset = [], 0
                        while True:
                            page = get_excerpts_citing(target, limit=200, offset=offset)
                            rows.extend(page['excerpts'])
                            if page['next_offset'] is None:
                                break
                            assert page['next_offset'] > offset
                            offset = page['next_offset']
                        excerpt_cache[target] = rows
                    found = [r for r in excerpt_cache[target] if r['citing_paper_hash'] == case['document'] and r['citation_index'] == index]
                    assert len(found) == 1 and found[0]['paragraph'] == measured['after']
                    measured['excerpt_replayed'] = True
                    measured['excerpt_work_id'] = target
                    measured['excerpt_citation_index'] = index
                    break
    finally:
        app._INDEX = old_index
        biblio.conn.close()
    source_files = ['pipeline/citation_spans.py', 'pipeline/grobid_client.py', 'pipeline/intext_citations.py', 'bib/authority.py', 'mcpsrv/tools/bibliography.py', 'tools/qc/citation_paragraphs.py']
    receipt = {'scope': 'Complete saved TEI/prepared PDF replay with matching producer receipt; current paragraph parser and actual authority materialization/tool functions. Not fresh extraction, all-library validation, or live MCP transport.',
               'detector': SIGNATURE, 'policy': CITATION_SPAN_POLICY,
               'code_base': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
               'consumed_code_sha256': {name: sha(REPO / name) for name in source_files},
               'input_roots': [str(p) for p in inputs], 'phase1': phase1, 'phase2': phase2,
               'unchanged_refresh': 'exact SQLite dump unchanged; zero phase1/phase2 work',
               'summary': {'documents': len(cases), 'paragraphs': sum(c['paragraphs'] for c in cases),
                           'reference_observations': sum(c['reference_observations'] for c in cases),
                           'before_hits': sum(len(c['before_hit_indices']) for c in cases),
                           'after_hits': sum(len(c['after_hit_indices']) for c in cases),
                           'intext_replays': sum(d['intext_replayed'] for c in cases for d in c['detector_cases']),
                           'excerpt_replays': sum(d['excerpt_replayed'] for c in cases for d in c['detector_cases'])},
               'cases': cases}
    dump(output / 'receipt.json', receipt)
    print(json.dumps(receipt['summary'], sort_keys=True))
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.input, args.output)


if __name__ == '__main__':
    main()
