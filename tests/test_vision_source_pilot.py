"""The pilot must capture actual Qwen preprocessing, without loading weights."""
import json
from types import SimpleNamespace

import numpy as np
from PIL import Image
import pytest

from tools.qc.vision_source_pilot import ProcessorCapture, reconstruct_model_frame, run


def _image_and_processor():
    pytest.importorskip('torch')
    pytest.importorskip('transformers')
    from transformers.models.qwen2_vl.image_processing_qwen2_vl import Qwen2VLImageProcessor
    # Deliberately non-square with distinguishable axes and channels. A wrong
    # patch/merge/temporal inverse can retain the dimensions yet scramble it.
    y,x=np.indices((56,84))
    pixels=np.stack((x*3,y*4,(x+y)*2),axis=-1).astype('uint8')
    return Image.fromarray(pixels),Qwen2VLImageProcessor(do_resize=False)


def test_reconstruction_matches_real_stateless_qwen_processor_pixels():
    image,processor=_image_and_processor()
    result=processor(images=[image],return_tensors='pt')
    rebuilt=reconstruct_model_frame(result.pixel_values.numpy(),result.image_grid_thw[0].tolist(),processor)
    assert rebuilt.size==image.size
    assert np.array_equal(np.asarray(rebuilt),np.asarray(image))
    with pytest.raises(ValueError,match='Unsupported still-image'):
        reconstruct_model_frame(result.pixel_values.numpy(),[2,4,6],processor)


def test_capture_preserves_exact_patch_tensor_prompt_and_raw_response(tmp_path):
    import torch
    image,processor=_image_and_processor()
    class Wrapped:
        image_processor=processor
        def __call__(self,**kwargs):
            return processor(images=kwargs['images'],return_tensors='pt')
        def batch_decode(self,*args,**kwargs):
            return ['{"panels": [], "embedded_figures": []}']
    capture=ProcessorCapture(Wrapped(),tmp_path)
    result=capture(images=[image],text=['exact production prompt'])
    stored=torch.load(tmp_path/'processor-01-patches.pt',weights_only=True)
    assert torch.equal(stored['pixel_values'],result.pixel_values)
    assert torch.equal(stored['image_grid_thw'],result.image_grid_thw)
    assert (tmp_path/'processor-01-prompt.txt').read_text()=='exact production prompt'
    assert Image.open(tmp_path/'processor-01-frame.png').size==(84,56)
    assert capture.batch_decode(torch.tensor([[10,20]]),skip_special_tokens=True)==[
        '{"panels": [], "embedded_figures": []}']
    receipt=json.loads((tmp_path/'raw-responses.json').read_text())
    assert receipt[0]['generated_token_ids']==[[10,20]]
    assert receipt[0]['processor_call']==1


def test_run_refuses_cpu_fallback_before_loading_any_model(tmp_path,monkeypatch):
    import torch
    from pipeline import vision
    monkeypatch.setattr(torch.cuda,'is_available',lambda:False)
    monkeypatch.setattr(vision,'LocalVLMBackend',lambda **kw:pytest.fail('model must not load'))
    with pytest.raises(ValueError,match='CUDA unavailable'):
        run(SimpleNamespace(device='cuda',pilot=tmp_path))
