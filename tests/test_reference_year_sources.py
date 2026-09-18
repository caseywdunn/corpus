"""The four named #314 references, freshly parsed from their source pages."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from bib.authority import phase2_references
from bib.documents import find_work
from bib.reference_year import adjudicate, candidate_index
from mcpsrv import app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import get_missing_references
from pipeline.grobid_client import parse_tei_references
from tests.test_bibliographic_integrity import build


FIXTURE = Path(__file__).parent / 'fixtures/bibliographic_integrity/publication_year_sources'
DATA = json.loads((FIXTURE / 'cases.json').read_text())
CASES = DATA['cases']
TARGET = DATA['target']


def make_authority(path):
    return build(path, {c['paper_hash']: c['curated_entry'] for c in CASES + [TARGET]})


def write_refs(path, field):
    for case in CASES:
        destination = path / 'documents' / case['paper_hash'] / 'references.json'
        destination.write_text(json.dumps({'references': [case[field]]}))


def missing(db, monkeypatch):
    monkeypatch.setattr(app, '_INDEX', SimpleNamespace(biblio_db=BiblioAuthority(db)))
    return get_missing_references(min_citations=1)


def edges(conn):
    return list(conn.execute('SELECT citing_work_id,cited_work_id FROM citations ORDER BY 1,2'))


@pytest.mark.parametrize('case', CASES, ids=lambda c: c['historical_xml_id'])
def test_actual_source_page_capture_preserves_raw_and_historical_ids(case):
    capture = case['fresh_capture']
    fragment = (FIXTURE / capture['tei_fragment']).read_bytes()
    assert hashlib.sha256(fragment).hexdigest() == capture['fragment_sha256']
    assert parse_tei_references(fragment.decode()) == [case['fresh_reference']]
    assert case['retained_reference']['xml_id'] == case['historical_xml_id']
    assert case['retained_reference']['raw'] == ''
    assert case['fresh_reference']['raw']
    assert case['fresh_reference']['xml_id'] != case['historical_xml_id']
    assert capture['include_raw_citations'] is True
    assert capture['consolidate_citations'] == 0


def test_real_refreshed_observations_rebuild_counts_missing_and_preserve_old_rows(tmp_path, monkeypatch):
    incremental = tmp_path / 'incremental'
    conn, db = make_authority(incremental)
    target = find_work(conn, TARGET['paper_hash'])
    write_refs(incremental, 'retained_reference')
    phase2_references(conn, incremental)
    historical = list(conn.execute('SELECT * FROM reference_observations ORDER BY observation_id'))
    assert len(historical) == 4
    assert conn.execute('SELECT COUNT(*) FROM citations WHERE cited_work_id=?', (target,)).fetchone()[0] == 0
    assert sum(row['cited_by_count'] for row in missing(db, monkeypatch)) == 4
    # A producer-only rebuild cannot invent the missing historical raw text.
    conn.execute("UPDATE observation_work SET producer_version='old'")
    conn.commit()
    phase2_references(conn, incremental)
    assert list(conn.execute('SELECT * FROM reference_observations ORDER BY observation_id')) == historical
    assert conn.execute('SELECT COUNT(*) FROM citations WHERE cited_work_id=?', (target,)).fetchone()[0] == 0

    # This refresh uses actual new Grobid outputs, not hand-filled raw fields.
    write_refs(incremental, 'fresh_reference')
    phase2_references(conn, incremental)
    for row in historical:
        assert conn.execute('SELECT * FROM reference_observations WHERE observation_id=?', (row[0],)).fetchone() == row
    assert conn.execute('SELECT COUNT(*) FROM reference_observations').fetchone()[0] == 8
    assert conn.execute('SELECT COUNT(*) FROM citations WHERE cited_work_id=?', (target,)).fetchone()[0] == 4
    assert missing(db, monkeypatch) == []
    mappings = list(conn.execute('SELECT work_id,match_method FROM observation_work'))
    assert all(wid == target for wid, _ in mappings)
    assert [method for _, method in mappings].count('raw_publication_year_title_authors') == 3
    assert [method for _, method in mappings].count('title_year_authors_fuzzy') == 1
    decisions = {sha: json.loads(reasons) for sha, reasons in conn.execute('''
        SELECT o.citing_corpus_hash,q.reasons_json FROM observation_work ow
        JOIN reference_observations o USING(observation_id)
        JOIN reference_observation_quality q USING(observation_id)''')}
    assert decisions['d0894c24715f'][0]['line_hyphen_joined_for_comparison'] is True
    assert decisions['a16337443af7'][0]['raw_title_alignment'] == 'single_article_omission'

    clean = tmp_path / 'clean'
    clean_conn, clean_db = make_authority(clean)
    write_refs(clean, 'fresh_reference')
    phase2_references(clean_conn, clean)
    assert edges(clean_conn) == edges(conn)
    assert missing(clean_db, monkeypatch) == []
    changes = conn.total_changes
    assert phase2_references(conn, incremental) == (0, 0)
    assert conn.total_changes == changes


@pytest.mark.parametrize('change', [
    {'raw': ''},
    {'raw_replace': ('54, 25-90', '55, 25-90')},
    {'raw_replace': ('54, 25-90', '54, 25-91')},
    {'raw_replace': ('54, 25-90', '54, 25 90')},
    {'raw_replace': ('54, 25-90', '54, 25–90–105')},
    {'raw_replace': ('54, 25-90', '54, 25-90 - 105')},
    {'raw_replace': ('P.R. 1974', 'P.R. 1975')},
    {'raw_replace': ('P.R. 1974', 'P.R. 1974, 1975')},
    {'raw_replace': ('PUGH', 'P- UGH')},
    {'raw_replace': ('of siphonophores', 'of plankton')},
    {'authors': ['P Pugh', 'A Another']},
    {'doi': '10.9999/conflicting'},
    {'volume': '55'},
    {'pages': '25–91'},
])
def test_printed_article_omission_needs_all_independent_evidence(tmp_path, change):
    conn, _ = make_authority(tmp_path)
    ref = dict(CASES[3]['fresh_reference'])
    change = dict(change)
    if pair := change.pop('raw_replace', None):
        assert pair[0] in ref['raw']
        ref['raw'] = ref['raw'].replace(*pair)
    ref.update(change)
    assert adjudicate(ref, candidate_index(conn))[0] is None


def test_no_substantive_fuzzy_title_repair_or_multiple_article_omissions(tmp_path):
    conn, _ = make_authority(tmp_path)
    index = candidate_index(conn)
    ref = CASES[0]['fresh_reference']
    assert adjudicate(dict(ref, title=ref['title'].replace('cruise', 'crise')), index) == (None, [])
    omitted = CASES[3]['fresh_reference']
    assert adjudicate(dict(omitted, title=omitted['title'].removeprefix('The ')), index) == (None, [])


@pytest.mark.parametrize('title_change', [
    lambda title: title.replace('collected ', ''),
    lambda title: 'A reconsideration of ' + title,
    lambda title: title + ' and a comparison with other collections',
])
def test_substantive_omissions_and_title_prefix_or_suffix_are_not_candidates(tmp_path, title_change):
    conn, _ = make_authority(tmp_path)
    ref = CASES[3]['fresh_reference']
    assert adjudicate(dict(ref, title=title_change(ref['title'])), candidate_index(conn)) == (None, [])


def test_matching_explicit_parsed_locators_preserve_source_supported_omission(tmp_path):
    conn, _ = make_authority(tmp_path)
    ref = dict(CASES[3]['fresh_reference'], volume='54', pages='25–90')
    assert adjudicate(ref, candidate_index(conn))[0] == find_work(conn, TARGET['paper_hash'])


@pytest.mark.parametrize('same_doi', [False, True])
def test_article_omission_and_wrapping_do_not_disambiguate_editions(tmp_path, same_doi):
    entry = TARGET['curated_entry']
    second = dict(entry, _key='OtherEdition', edition='2',
                  doi=entry['doi'] if same_doi else '10.9999/other-edition')
    conn, _ = build(tmp_path, {'first': dict(entry, edition='1'), 'second': second})
    for case in (CASES[1], CASES[3]):
        work, reasons = adjudicate(case['fresh_reference'], candidate_index(conn))
        assert work is None
        assert reasons[0]['candidate_count'] == 2
