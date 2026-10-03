"""Exercise both production backend paths across model/raster frames (#305)."""
import json
import sys
from types import SimpleNamespace

import pytest
from PIL import Image

from pipeline.figures import detect_figure_rois_via_vision
from pipeline.vision import ClaudeVisionBackend, LocalVLMBackend


def _response(size, pixels):
    w, h = size if pixels else (1, 1)
    return json.dumps({'panels': [
        {'label': label, 'panel_bbox_norm': [x * w, 0, (x + .5) * w, h],
         'label_bbox_norm': [x * w, 0, (x + .1) * w, .1 * h],
         'confidence': .9}
        for label, x in [('A', 0), ('B', .5)]
    ], 'embedded_figures': [
        {'figure_number': '1', 'panel_bbox_norm': [0, 0, w, h], 'confidence': .9},
    ]})


def _claude(response):
    backend = object.__new__(ClaudeVisionBackend)
    backend.model = 'stub'
    backend.max_tokens = 4096
    backend.client = SimpleNamespace(messages=SimpleNamespace(create=lambda **kw:
        SimpleNamespace(content=[SimpleNamespace(type='text', text=response)],
                        stop_reason='end_turn')))
    return backend


def _local(response, model_size, monkeypatch):
    # Stub inference only: the backend still prepares messages, reads the real
    # patch-grid contract, converts model output and normalizes public ROIs.
    import torch

    class Inputs(dict):
        def __getattr__(self, name):
            return self[name]

        def to(self, device):
            return self

    class Processor:
        image_processor = SimpleNamespace(patch_size=14)

        def apply_chat_template(self, messages, **kwargs):
            self.last_text = messages[1]['content'][1]['text']
            return self.last_text

        def __call__(self, **kwargs):
            return Inputs(input_ids=torch.tensor([[0]]), image_grid_thw=torch.tensor(
                [[1, model_size[1] // 14, model_size[0] // 14]],
            ))

        def batch_decode(self, *args, **kwargs):
            return [response]

    monkeypatch.setitem(sys.modules, 'qwen_vl_utils', SimpleNamespace(
        process_vision_info=lambda messages: ([messages[1]['content'][0]['image']], None),
    ))
    backend = object.__new__(LocalVLMBackend)
    backend._model_id = 'stub'
    backend._max_new_tokens = 4096
    backend._processor = Processor()
    backend._model = SimpleNamespace(device='cpu', generate=lambda **kwargs:
                                    torch.tensor([[0, 1]]))
    return backend


@pytest.mark.parametrize('backend_kind', ['claude', 'local'])
@pytest.mark.parametrize('raster_size', [(4000, 2400), (1400, 840)])
@pytest.mark.parametrize('pixels', [False, True])
def test_backend_boxes_use_the_served_raster_frame(
    tmp_path, monkeypatch, backend_kind, raster_size, pixels,
):
    w, h = raster_size
    path = tmp_path / 'figure.png'
    img = Image.new('RGB', raster_size, 'red')
    img.paste('blue', (w // 2, 0, w, h))
    img.save(path)
    if backend_kind == 'claude':
        model_size = (2000, 1200) if w > 2000 else raster_size
        backend = _claude(_response(model_size, pixels))
    else:
        # A processor frame independent of the source (Qwen's patch multiples).
        model_size = (840, 504)
        backend = _local(_response(model_size, pixels), model_size, monkeypatch)
    result = detect_figure_rois_via_vision(
        path, [{'label': 'A'}, {'label': 'B'}], backend,
    )
    assert result['image_size_px'] == list(raster_size)
    assert result['pass3_status'] == 'completed'
    a, b, embedded = result['rois']
    assert a['roi_px'] == [0, 0, w // 2, h]
    assert b['roi_px'] == [w // 2, 0, w, h]
    assert embedded['roi_px'] == [0, 0, w, h]
    assert a['label_bbox_px'] == [0, 0, w // 10, h // 10]
    assert b['label_bbox_px'] == [w // 2, 0, int(.6 * w), h // 10]
    assert img.crop(b['roi_px']).getextrema() == ((0, 0), (0, 0), (255, 255))
    for roi in result['rois']:
        provenance = roi['coordinate_provenance']['panel']
        assert provenance['model_input_size_px'] == list(model_size)
        assert provenance['raster_size_px'] == list(raster_size)
        assert provenance['model_coordinate_units'] == ('pixels' if pixels else 'normalized')


def test_expanded_caption_inventory_prevents_false_completion(tmp_path):
    from pipeline.figures import parse_panels_from_caption
    path = tmp_path / 'figure.png'
    Image.new('RGB', (1400, 840)).save(path)
    backend = _claude(_response((1400, 840), False))
    panels = parse_panels_from_caption('FIGURE 1. A-B. live specimens. C-N. preserved specimens.')
    result = detect_figure_rois_via_vision(path, panels, backend)
    assert len(panels) == 14
    assert result['pass3_status'] == 'partial_vision'
    assert {r.get('label') for r in result['rois'] if r['type'] == 'panel'} == {'A', 'B'}
