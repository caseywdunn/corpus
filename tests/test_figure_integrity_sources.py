"""Source replay uses the existing library; no extra fixture corpus is built.

CORPUS_LIBRARY_DIR=/path/to/library pytest tests/test_figure_integrity_sources.py
The expectations were checked against the source page, independently of the
layout detector. Automated geometry complements the recorded visual review.
"""
import hashlib
import json
import os
from pathlib import Path

import fitz
import pytest

from pipeline.figures import figure_rect_for_bbox, render_figures


CASES = json.loads((Path(__file__).parent / 'fixtures' / 'figure_integrity' /
                    'source_cases.json').read_text())['cases']


@pytest.mark.parametrize('case', CASES, ids=lambda c: f"issue-{c['issue']}")
def test_reviewed_source_extent_is_preserved(case, tmp_path):
    root = os.environ.get('CORPUS_LIBRARY_DIR')
    if not root:
        pytest.skip('set CORPUS_LIBRARY_DIR to replay existing source PDFs')
    source = Path(root) / case['library_path']
    assert hashlib.sha256(source.read_bytes()).hexdigest()[:12] == case['paper_hash']
    figure = {'figure_id': case['figure_id'], 'filename': 'figure.png',
              'extraction_method': 'docling', 'page': case['pdf_page'],
              'bbox': case['detected_bbox_bottom_left'],
              'bbox_coord_system': 'pdf_pts_bottom_left'}
    result = render_figures(source, [figure], tmp_path, native=True, pixel_cap=3000)
    assert result['rendered'] == 1
    with fitz.open(source) as doc:
        actual = figure_rect_for_bbox(figure['bbox'], figure['bbox_coord_system'],
                                      doc[case['pdf_page'] - 1].cropbox.height)
    assert actual.contains(fitz.Rect(case['source_image_bbox_top_left']))
    assert actual.y1 < case['caption_top']
    assert figure['bbox_boundary_evidence']['source'] == 'embedded_image_extent'
