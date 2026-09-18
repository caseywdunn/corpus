"""Figure count and caption regressions from the September audit (#324/#332)."""
import json

import pytest

from pipeline.figures import parse_panels_from_caption
from mcpsrv.bundle import _count_figures_and_chunks


@pytest.mark.parametrize('end', ['L', 'N', 'U', 'Z'])
def test_ranges_cover_the_full_uppercase_alphabet(end):
    panels = parse_panels_from_caption(
        f'FIGURE 1. A-C. live specimens. D-{end}. preserved specimens.'
    )
    assert [p['label'] for p in panels] == list(map(chr, range(65, ord(end) + 1)))


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
