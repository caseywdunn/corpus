"""Replay source geometry for Hosia and Sutherland caption failures (#322)."""
import json
from pathlib import Path
from types import SimpleNamespace as NS

import pytest

from pipeline.figures import (
    _text_fragments, expand_plate_figures, extract_caption_info,
    plate_legend_entries,
)


class Ref:
    def __init__(self, cref):
        self.cref = cref

    def resolve(self, document):
        return document.by_ref[self.cref]


def source_document(name):
    raw = json.loads((Path(__file__).parent / 'fixtures' / 'figure_integrity' /
                      f'{name}_captions.json').read_text())

    def item(raw_item):
        data = dict(raw_item)
        data['prov'] = [NS(**{**p, 'bbox': NS(**p['bbox'])}) for p in data['prov']]
        for key in ('parent',):
            if key in data:
                data[key] = Ref(data[key]['$ref'])
        for key in ('children', 'captions'):
            if key in data:
                data[key] = [Ref(r['$ref']) for r in data[key]]
        return NS(**data)

    doc = NS(texts=[item(t) for t in raw['texts']], pictures=[item(p) for p in raw['pictures']])
    doc.by_ref = {i.self_ref: i for i in doc.texts + doc.pictures}
    return doc


def test_hosia_figure4_retains_taxon_description_and_wrapped_scale_tail():
    doc = source_document('hosia')
    figure = extract_caption_info(doc.pictures[0], doc)
    body = doc.by_ref['#/texts/104'].text
    assert figure['caption_text'] == f'FIGURE 4. {body} bar 1 mm.'
    assert figure['figure_number'] == '4'
    assert figure['caption_status'] == 'bound'
    assert [part['text'] for part in figure['caption_fragments']] == [
        'FIGURE 4 bar 1 mm.', body,
    ]
    assert figure['caption_completeness'] == 'unverified'


def test_hosia_figure5_uses_owned_text_despite_noncaption_layout_labels():
    doc = source_document('hosia')
    body = doc.by_ref['#/texts/107']
    assert body.label == 'text'
    assert doc.by_ref['#/texts/106'].label == 'section_header'
    figure = extract_caption_info(doc.pictures[1], doc)
    assert figure['caption_text'] == f'FIGURE 5. {body.text}'
    assert figure['caption_text'].endswith('1 mm scale bar for bracts.')
    assert figure['caption_kind'] == 'prose_caption'
    assert len(figure['caption_fragments']) == 2


@pytest.mark.parametrize("legacy_link", [False, True])
def test_reconciliation_keeps_reconstructed_hosia_captions(legacy_link):
    doc = source_document('hosia')
    items = []
    for index, picture in enumerate(doc.pictures[:2]):
        prov = picture.prov[0]
        items.append({'docling_idx': index, 'page': 8,
                      'bbox': [prov.bbox.l, prov.bbox.b, prov.bbox.r, prov.bbox.t],
                      **extract_caption_info(picture, doc)})
    legends = plate_legend_entries([
        {'text': text, 'bbox': bbox}
        for item in doc.texts for text, bbox, page in _text_fragments(item) if page == 8
    ])
    if legacy_link:
        body = doc.by_ref['#/texts/104']
        bbox = body.prov[0].bbox
        items[0].update(figure_number=None, caption_text=body.text,
                        caption_bbox=[bbox.l, bbox.b, bbox.r, bbox.t])
    result = expand_plate_figures(items, {8: legends})
    assert len(result) == 2
    assert 'pneumatophore and upper nectosome' in result[0]['caption_text']
    assert result[0]['caption_text'].endswith('Scale bar 1 mm.')
    assert 'Composite illustration of mature' in result[1]['caption_text']


def test_nearby_unowned_prose_does_not_complete_bare_label():
    doc = source_document('hosia')
    doc.by_ref['#/texts/107'].parent = Ref('#/body')
    figure = extract_caption_info(doc.pictures[1], doc)
    assert figure['caption_text'] == 'FIGURE 5'
    assert figure['caption_kind'] == 'bare_label'
    assert figure['caption_completeness'] == 'unverified'


def test_sutherland_left_panel_uses_same_row_caption_not_results_prose():
    # Fresh one-page Docling extraction of PDF a65c78792a7e, page 5,
    # reproduces the audit's body-link failure before the correction.
    doc = source_document('sutherland')
    left = extract_caption_info(doc.pictures[0], doc)
    right = extract_caption_info(doc.pictures[1], doc)
    assert left['caption_text'] == right['caption_text'] == doc.by_ref['#/texts/40'].text
    assert left['figure_number'] == '5'
    assert 'result, N. bijuga traveled twice' not in left['caption_text']
    rejected = next(c for c in left['caption_candidates']
                    if c['caption_text'].startswith('result, N. bijuga'))
    assert rejected['rejection_reason'] == 'labelled_caption_on_same_picture_row'
    assert not rejected['chosen']


def test_a_different_row_does_not_override_an_unlabelled_caption():
    doc = source_document('sutherland')
    # A legitimate unlabelled caption remains usable; remove the actual Fig5
    # side caption so the only labelled alternative is Fig6 on the next row.
    doc.texts = [t for t in doc.texts if t.self_ref != '#/texts/40']
    doc.by_ref['#/texts/3'].text = 'Unlabelled but structurally linked specimen view.'
    doc.by_ref['#/texts/3'].prov[0].charspan = (0, len(doc.by_ref['#/texts/3'].text))
    result = extract_caption_info(doc.pictures[0], doc)
    assert result['caption_text'] == 'Unlabelled but structurally linked specimen view.'
    assert result['caption_kind'] == 'unlabelled_caption'


def test_reconstructed_hosia_figure1_keeps_the_reuse_exclusion():
    doc = source_document('hosia')
    result = extract_caption_info(doc.pictures[2], doc)
    assert result['figure_number'] == '1'
    assert 'Range of Nanomia nectophore shapes.' in result['caption_text']
    assert ('These images are not covered by the terms of the Creative Commons '
            'license of this publication.') in result['caption_text']
    assert result['caption_text'].endswith('please contact the relevant rights holder.')
