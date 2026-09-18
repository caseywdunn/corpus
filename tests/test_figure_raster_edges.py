"""Preserve lettering clipped by layout detection at image boundaries (#329)."""
import io

import fitz
from PIL import Image, ImageDraw
import pytest

from pipeline.figures import complete_embedded_raster_bounds, render_figures


def _source_pdf(tmp_path, rotation=0, nearby_text=False):
    path = tmp_path / 'source.pdf'
    image = Image.new('RGB', (400, 200), 'white')
    draw = ImageDraw.Draw(image)
    draw.rectangle((385, 185, 398, 198), fill='red')  # identifying edge lettering
    buffer = io.BytesIO()
    image.save(buffer, format='PNG')
    doc = fitz.open()
    page = doc.new_page(width=400, height=500)
    page.insert_image(fitz.Rect(40, 80, 240, 180), stream=buffer.getvalue())
    if nearby_text:
        page.insert_text((235, 176), 'Adjacent prose', fontsize=5)
    else:
        page.insert_text((40, 199), 'FIGURE 1. Neighboring caption.', fontsize=10)
    page.set_rotation(rotation)
    doc.save(path)
    doc.close()
    return path


@pytest.mark.parametrize('rotation', [0, 90, 180, 270])
@pytest.mark.parametrize('coord', ['pdf_pts_top_left', 'pdf_pts_bottom_left'])
def test_render_recovers_source_edge_without_caption(tmp_path, rotation, coord):
    pdf = _source_pdf(tmp_path, rotation)
    clipped = [40, 80, 237, 177]
    bbox = clipped if coord.endswith('top_left') else [40, 323, 237, 420]
    figure = {'figure_id': 'figure_1', 'filename': 'fig_1.png', 'page': 1,
              'bbox': bbox, 'bbox_coord_system': coord, 'extraction_method': 'docling'}
    result = render_figures(pdf, [figure], tmp_path, native=False, fixed_scale=2)
    assert result['rendered'] == 1
    expected_bbox = [40, 80, 240, 180] if coord.endswith('top_left') else [40, 320, 240, 420]
    assert figure['bbox'] == expected_bbox
    assert figure['detected_bbox'] == bbox
    assert figure['bbox_boundary_evidence']['source'] == 'embedded_image_extent'
    with fitz.open(pdf) as doc:
        expected = doc[0].get_pixmap(matrix=fitz.Matrix(2, 2),
                                    clip=fitz.Rect(40, 80, 240, 180) * doc[0].rotation_matrix)
        with Image.open(tmp_path / 'fig_1.png') as actual:
            assert actual.tobytes() == expected.samples
            assert sorted(actual.size) == [200, 400]


def test_nearby_prose_vetoes_bound_expansion(tmp_path):
    with fitz.open(_source_pdf(tmp_path, nearby_text=True)) as doc:
        rect = fitz.Rect(40, 80, 237, 177)
        repaired, evidence = complete_embedded_raster_bounds(doc[0], rect)
        assert repaired == rect
        assert evidence is None


def test_partial_panel_is_not_expanded_to_whole_raster(tmp_path):
    with fitz.open(_source_pdf(tmp_path)) as doc:
        rect = fitz.Rect(40, 80, 140, 180)
        assert complete_embedded_raster_bounds(doc[0], rect) == (rect, None)


def test_fixed_repair_mode_leaves_correct_images_alone(tmp_path):
    pdf = _source_pdf(tmp_path)
    figure = {'filename': 'fig_1.png', 'page': 1, 'bbox': [40, 80, 240, 180],
              'bbox_coord_system': 'pdf_pts_top_left'}
    result = render_figures(pdf, [figure], tmp_path, native=False, repair_bounds_only=True)
    assert result['rendered'] == 0
    assert not (tmp_path / 'fig_1.png').exists()


def test_rerender_invalidates_boxes_in_the_old_image_frame(tmp_path):
    pdf = _source_pdf(tmp_path)
    figure = {'filename': 'fig_1.png', 'page': 1, 'bbox': [40, 80, 237, 177],
              'bbox_coord_system': 'pdf_pts_top_left', 'image_size_px': [394, 194],
              'rois': [{'type': 'panel', 'label': 'A', 'roi_px': [0, 0, 100, 100]}],
              'pass3_status': 'completed', 'panels_from_caption': [{'label': 'A'}]}
    render_figures(pdf, [figure], tmp_path, native=False)
    assert figure['rois'] == []
    assert figure['pass3_status'] == 'stale_image_geometry'
    assert figure['panels_from_caption'] == [{'label': 'A'}]
    assert figure['image_size_px'] == [400, 200]
    assert figure['roi_geometry_invalidated']['previous_image_size_px'] == [394, 194]


@pytest.mark.parametrize('rotation', [0, 90])
def test_page_crop_margins_use_local_unrotated_coordinates(tmp_path, rotation):
    original = _source_pdf(tmp_path)
    pdf = tmp_path / 'cropped.pdf'
    with fitz.open(original) as doc:
        doc[0].set_cropbox(fitz.Rect(20, 40, 350, 400))
        doc[0].set_rotation(rotation)
        doc.save(pdf)
    figure = {'filename': 'fig_1.png', 'page': 1, 'bbox': [20, 40, 217, 137],
              'bbox_coord_system': 'pdf_pts_top_left'}
    result = render_figures(pdf, [figure], tmp_path, native=False)
    assert result['rendered'] == 1
    assert figure['bbox'] == [20, 40, 220, 140]
    assert sorted(figure['image_size_px']) == [200, 400]
