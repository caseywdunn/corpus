"""Actual source runs and conservative boundary controls for #317."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re
import sqlite3
from types import SimpleNamespace

from lxml import etree
import pytest

from bib.authority import create_schema, phase1_corpus_papers, phase2_references
from mcpsrv import app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import get_excerpts_citing, get_intext_citations
from pipeline import citation_spans
from pipeline.grobid_client import parse_tei_intext_citations, parse_tei_references
from tests.test_citation_source_spans import CapturedSource

FIXTURES = Path(__file__).parent / 'fixtures/citation_spans/group_boundaries'
DETECTOR = re.compile(r'\(\d{4}[A-Za-z][A-Za-z.\- ]{1,30}\(')


def capture(name):
    proof = json.loads((FIXTURES / f'{name}.json').read_text())
    xml = (FIXTURES / f'{name}.tei.xml').read_text()
    assert hashlib.sha256(xml.encode()).hexdigest() == proof['fragment_sha256']
    return proof, xml


@pytest.mark.parametrize('name,expected', [
    ('Totton1965a', 'Leuckart (1853, 1854), Gegenbaur (1853), Vogt (1854)'),
    ('Tung2003', 'Margulis（1980, 1984）'),
])
def test_actual_source_fragments_materialize_and_serve_without_expansions(tmp_path, monkeypatch, name, expected):
    proof, xml = capture(name)
    baseline = parse_tei_intext_citations(xml)
    assert DETECTOR.search(baseline['paragraphs'][0])
    monkeypatch.setattr(citation_spans, 'PdfCitationSource', lambda _: CapturedSource(proof))
    data = parse_tei_intext_citations(xml, pdf_path=Path('captured.pdf'))
    assert expected in data['paragraphs'][0]
    assert not DETECTOR.search(data['paragraphs'][0])
    # Save/reparse is deterministic; valid repeated prose and raw expansions
    # remain in their own observations, not globally deleted.
    assert data == parse_tei_intext_citations(xml, pdf_path=Path('captured.pdf'))
    assert any(DETECTOR.search(''.join(o['surface'] for o in r['tei_observations']))
               for r in data['citations'])
    for row in data['citations']:
        if 'author_span' in row:
            assert data['paragraphs'][0][slice(*row['year_span'])] == row['citation_year']
    if name == 'Tung2003':
        for row in data['citations'][:2]:
            assert data['paragraphs'][0][slice(*row['author_span'])] == 'Margulis'
            assert row['source_span_evidence'] == {
                'method': 'unique_parenthesized_marker_at_box_edges',
                'source_interval': '。Margulis（1980, 1984）綜',
                'discarded_prefix': '。', 'discarded_suffix': '綜'}
            # Actual Grobid splits R. Ya. into a second author. Text repair
            # cannot establish the missing work link; keep it unresolved.
            assert row['target_xml_id'] is None
            assert row['validation_status'] == 'unresolved_author_year'
    sha = proof['source_sha256'][:12]
    doc = tmp_path / 'documents' / sha
    doc.mkdir(parents=True)
    refs = parse_tei_references(xml)
    (doc / 'metadata.json').write_text(json.dumps(proof['metadata']))
    (doc / 'references.json').write_text(json.dumps({'references': refs}))
    (doc / 'intext_citations.json').write_text(json.dumps(data))
    db = tmp_path / 'biblio.sqlite'
    with sqlite3.connect(db) as conn:
        create_schema(conn)
        assert phase1_corpus_papers(conn, tmp_path) == 1
        phase2_references(conn, tmp_path)
        observed = conn.execute('SELECT grobid_xml_id,raw_citation,title,year FROM reference_observations ORDER BY ordinal').fetchall()
        assert observed == [(r['xml_id'], r['raw'], r['title'], r['year']) for r in refs]
        targets = dict(conn.execute('SELECT grobid_xml_id,cited_work_id FROM citations'))
        before = list(conn.iterdump())
        assert phase2_references(conn, tmp_path) == (0, 0)
        assert list(conn.iterdump()) == before
    biblio = BiblioAuthority(db)
    monkeypatch.setattr(app, '_INDEX', SimpleNamespace(biblio_db=biblio, papers={sha: {'hash_dir': str(doc), 'title': name}}))
    try:
        served = get_intext_citations(sha, limit=200)
        assert served['paragraphs'] == data['paragraphs']
        assert served['citations'] == data['citations']
        # Use actual parsed/resolved targets only; unresolved Margulis markers
        # must not acquire an edge merely to make the excerpt test pass.
        found = 0
        for index, marker in enumerate(data['citations']):
            target = targets.get((marker.get('target_xml_id') or '').lstrip('#'))
            if target is None:
                continue
            response = get_excerpts_citing(target)
            row = next(r for r in response['excerpts'] if r['citing_paper_hash'] == sha and r['citation_index'] == index)
            assert row['paragraph'] == data['paragraphs'][0]
            assert expected in row['paragraph']
            found += 1
        assert found >= 1
    finally:
        biblio.conn.close()


def marker_source(proof):
    xml = etree.parse(str(FIXTURES / 'Tung2003.tei.xml'))
    refs = xml.findall('.//tei:body//tei:ref', citation_spans.NS)[:2]
    raw = ''.join(''.join(r.itertext()) for r in refs)
    return CapturedSource(proof).group_text(refs, raw)


@pytest.mark.parametrize('change', ['fully_contained', 'crossing_latin', 'crossing_digit', 'leading_letter', 'author_disagrees', 'year_disagrees'])
def test_crossing_glyph_is_not_permission_to_change_identity(change):
    proof, _ = capture('Tung2003')
    proof = deepcopy(proof)
    chars = proof['pdf_characters']['17']
    if change == 'fully_contained':
        chars[-1][1][2] = 485.0
    elif change == 'crossing_latin':
        chars[-1][0] = 'a'
    elif change == 'crossing_digit':
        chars[-1][0] = '7'
    elif change == 'leading_letter':
        chars[0][0] = 'B'
    elif change == 'author_disagrees':
        next(c for c in chars if c[0] == 'M')[0] = 'N'
    else:
        next(c for c in chars if c[0] == '4')[0] = '5'
    assert marker_source(proof) is None


def test_two_complete_markers_are_not_one_unique_parenthesized_marker():
    proof, _ = capture('Tung2003')
    chars = proof['pdf_characters']['17']
    first, last = chars[0], chars[-1]
    boxes = [(331.55, 659.07, 425.69, 673.04), (436.37, 660.71, 485.71, 673.29)]
    source = '。Margulis（1980, 1984） Margulis（1980, 1984）綜'
    assert citation_spans._bounded_marker(source, 'Margulis(1980Margulis( , 1984))', first, last, boxes) is None


@pytest.mark.parametrize('coords,expected', [
    ('1,10,112,20,10', True),  # adjacent line
    ('1,10,180,20,10', False),  # omitted intervening source prose
    ('2,10,10,20,10', True),  # next page continuation
    ('3,10,10,20,10', False),  # omitted physical page
    ('1,10,100,0,10', True),  # invalid geometry: preserve TEI fallback
    ('', True),
])
def test_source_group_boundaries_preserve_adjacent_and_unlocated_refs(coords, expected):
    left = etree.Element('ref', coords='1,10,100,20,10')
    right = etree.Element('ref', coords=coords)
    assert citation_spans._source_neighbors(left, right) == expected
