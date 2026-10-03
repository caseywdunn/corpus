"""Figure count and caption regressions from the September audit (#324/#332)."""
import json
from pathlib import Path

import pytest

from pipeline.figures import parse_panels_from_caption
from mcpsrv.bundle import _count_figures_and_chunks


@pytest.mark.parametrize('end', ['L', 'N', 'U', 'Z'])
def test_ranges_cover_the_full_uppercase_alphabet(end):
    panels = parse_panels_from_caption(
        f'FIGURE 1. A-C. live specimens. D-{end}. preserved specimens.'
    )
    assert [p['label'] for p in panels] == list(map(chr, range(65, ord(end) + 1)))


def test_individual_late_letters_follow_an_earlier_range():
    panels = parse_panels_from_caption('FIGURE 1. A-X. views. (Y) upper. (Z) lower.')
    assert [p['label'] for p in panels] == list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    assert panels[-1]['description'] == 'lower'


def test_source_sutherland_inline_genus_does_not_destroy_explicit_panel_inventory():
    fixture = json.loads((Path(__file__).parent / 'fixtures/figure_integrity/'
                          'sutherland_captions.json').read_text())
    caption = next(t['text'] for t in fixture['texts'] if t['self_ref'] == '#/texts/40')
    panels = parse_panels_from_caption(caption)
    assert [p['label'] for p in panels] == ['A', 'B']
    assert 'swimming by N. bijuga' in panels[0]['description']
    assert 'with N. bijuga' in panels[1]['description']


@pytest.mark.parametrize('caption', [
    'FIGURE 1. A-M. specimens. N. upper view.',
    'FIGURE 1. A-N. specimens viewed with N. upper view.',
    'FIGURE 1. A-M. specimens viewed with (N) upper view.',
])
def test_late_panel_definitions_and_independent_declarations_survive(caption):
    assert [p['label'] for p in parse_panels_from_caption(caption)] == list('ABCDEFGHIJKLMN')


def test_new_caption_labels_retire_old_completion_without_inventing_rois(tmp_path):
    from pipeline.figure_passes import _pass25_annotate_figures
    figures = tmp_path / 'figures.json'
    original_rois = [{'label': letter, 'roi_px': [0, 0, 10, 10]} for letter in 'ABCDEF']
    figures.write_text(json.dumps({'figures': [{
        'figure_id': 'fig12', 'figure_number': '12',
        'caption_text': 'FIGURE 12. (A-F) living. (G-N) fixed.',
        'pass3_status': 'completed', 'rois': original_rois,
    }]}))
    _pass25_annotate_figures(tmp_path / 'absent-text.json', figures)
    updated = json.loads(figures.read_text())['figures'][0]
    assert updated['panel_count_from_caption'] == 14
    assert updated['pass3_status'] == 'partial_caption_inventory'
    assert updated['pass3_inventory_reconciliation']['unlocated_labels'] == list('GHIJKLMN')
    assert updated['rois'] == original_rois
    before = figures.read_bytes()
    _pass25_annotate_figures(tmp_path / 'absent-text.json', figures)
    assert figures.read_bytes() == before


@pytest.mark.parametrize('target_kind', ['panel', 'figure', 'figure_discovery'])
def test_completion_reconciliation_preserves_full_coverage_and_numbered_targets(target_kind):
    from pipeline.figure_passes import _reconcile_panel_completion
    figure = {'pass3_status': 'completed', 'pass3_target_kind': target_kind,
              'panels_from_caption': [{'label': 'A'}, {'label': 'B'}],
              'rois': [{'label': 'A', 'roi_px': [0, 0, 5, 5]},
                       {'label': 'B', 'roi_px': [5, 0, 10, 5]}] if target_kind == 'panel' else []}
    _reconcile_panel_completion(figure)
    assert figure['pass3_status'] == 'completed'
    assert 'pass3_inventory_reconciliation' not in figure


@pytest.mark.parametrize('status', ['completed', 'vision_backend_failed'])
def test_new_vision_result_retires_inventory_note_only_after_success(tmp_path, monkeypatch, status):
    from pipeline import figure_passes
    figures = tmp_path / 'figures.json'
    image = tmp_path / 'figure.png'
    image.write_bytes(b'not opened by the controlled detection backend')
    note = {'previous_status': 'completed', 'unlocated_labels': ['B']}
    figures.write_text(json.dumps({'figures': [{
        'figure_id': 'fig1', 'figure_type': 'figure', 'file_path': str(image),
        'panels_from_caption': [{'label': 'A'}, {'label': 'B'}],
        'pass3_status': 'partial_caption_inventory', 'pass3_inventory_reconciliation': note,
    }]}))
    rois = [{'label': label, 'roi_px': box} for label, box in
            [('A', [0, 0, 5, 10]), ('B', [5, 0, 10, 10])]] if status == 'completed' else []
    monkeypatch.setattr(figure_passes, 'detect_figure_rois_via_vision',
                        lambda *args, **kwargs: {'pass3_status': status, 'rois': rois})
    figure_passes._pass3b_annotate_rois(figures, object())
    updated = json.loads(figures.read_text())['figures'][0]
    assert ('pass3_inventory_reconciliation' in updated) == (status == 'vision_backend_failed')


def test_specific_species_description_retains_shared_context():
    # Siebert et al. 2013, PDF page 5, Figure 2; source c9e7e8ae50a2.
    caption = ('FIGURE 2. (A–B) In situ photographs of holotype specimens. '
               '(A) Apolemia lanosa sp. nov., approximate length 2 m. '
               '(B) Apolemia rubriversa sp. nov., approximately 1.2 m in length. '
               '(A) Scale bar 2 cm.')
    panels = parse_panels_from_caption(caption)
    assert panels[0]['description'].startswith('Apolemia lanosa sp. nov.')
    assert panels[1]['description'].startswith('Apolemia rubriversa sp. nov.')
    assert 'In situ photographs of holotype specimens' in panels[0]['shared_descriptions'][0]
    assert 'Scale bar' not in panels[0]['description']


@pytest.mark.parametrize('caption,end', [
    # Siebert et al. 2013, PDF pages 7 and 16; source c9e7e8ae50a2.
    ('FIGURE 4. Apolemia lanosa sp. nov. Medium sized nectophore from the '
     'holotype specimen (living, A–C). Photographs and drawings of upper '
     '(A, D), lower (B, E), and lateral (C, F) views. Scale bar (5 mm) in C '
     'applies to figures. (G–S) Photographs of “upper” view of nectophores '
     'from the holotype and lateral view of the paratype specimen (T–U).', 'U'),
    ('FIGURE 12. Apolemia rubriversa sp. nov. (A–F) Medium sized nectophore '
     'of the holotype (live). Canals were colored red by prey pigment. '
     'Photographs and drawings of upper (A, D), lower (B, E) and lateral '
     '(C, F) views. (G–N) Available fixed nectophores of sample D195 shown '
     'in upper view. Scale bar (1 cm) in J applies to G–N.', 'N'),
])
def test_source_caption_late_ranges(caption, end):
    assert [p['label'] for p in parse_panels_from_caption(caption)] == list(
        map(chr, range(65, ord(end) + 1)),
    )


@pytest.mark.parametrize('number,expected', [('2', 'ABCD'), ('4', 'ABCDEFGHIJKLMNOPQRSTU'),
                                             ('12', 'ABCDEFGHIJKLMN')])
def test_complete_source_checked_siebert_captions(number, expected):
    fixture = json.loads((Path(__file__).parent / 'fixtures/figure_integrity/'
                          'siebert_panel_captions.json').read_text())
    caption = next(f['caption_text'] for f in fixture['figures'] if f['figure_number'] == number)
    panels = parse_panels_from_caption(caption)
    assert ''.join(p['label'] for p in panels) == expected
    if number == '2':
        assert 'Apolemia lanosa' in panels[0]['description']
        assert 'Apolemia rubriversa' in panels[1]['description']
        assert all('In situ photographs of holotype specimens' in p['shared_descriptions'][0]
                   for p in panels[:2])


@pytest.mark.parametrize('suffix', [
    'Compare with Figure 9 (M-Z).',
    'Photographs by N. Smith and Z. Jones.',
    'C.ped = pedicular canal; Z. rad = radial canal; H = hydroecium.',
])
def test_late_letters_do_not_invent_panels_from_other_context(suffix):
    panels = parse_panels_from_caption('FIGURE 1. A-B. specimens. ' + suffix)
    assert [p['label'] for p in panels] == ['A', 'B']


def test_bundle_rejects_a_stale_record_total(tmp_path):
    paper = tmp_path / 'abc123abc123'
    paper.mkdir()
    (paper / 'figures.json').write_text(json.dumps({
        'total_figures': 1,
        'figures': [{'filename': 'shared.png'}, {'filename': 'shared.png'}],
    }))
    with pytest.raises(ValueError, match='total_figures=1.*2 records'):
        _count_figures_and_chunks(tmp_path)
    # Historical artifacts that omit the field remain valid. Shared images
    # count twice because these are logical records, not a raster inventory.
    (paper / 'figures.json').write_text(json.dumps({
        'figures': [{'filename': 'shared.png'}, {'filename': 'shared.png'}],
    }))
    assert _count_figures_and_chunks(tmp_path) == (2, 0)


def test_legacy_index_uses_records_for_all_public_totals(tmp_path, monkeypatch):
    from mcpsrv.indexes import CorpusIndex
    from mcpsrv.tools import papers
    paper = tmp_path / 'documents' / 'abc123abc123'
    paper.mkdir(parents=True)
    (paper / 'metadata.json').write_text(json.dumps({'title': 'Grouped plate'}))
    (paper / 'figures.json').write_text(json.dumps({
        'total_figures': 1,
        'figures': [{'filename': 'shared.png'}, {'filename': 'shared.png'}],
    }))
    index = CorpusIndex(tmp_path)
    index.load()
    assert index.papers[paper.name]['n_figures'] == 2
    monkeypatch.setattr(papers, '_need_index', lambda: index)
    assert papers.corpus_summary()['n_figures_total'] == 2
    assert papers.list_papers()[0]['n_figures'] == 2
