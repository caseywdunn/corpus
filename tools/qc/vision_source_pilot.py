#!/usr/bin/env python3
"""Capture a fresh local-VLM source pilot without changing corpus artifacts.

Existing panel_bbox_units.py inventories old ROI receipts; vlm_dtype_probe.py
compares converted boxes across dtypes. This companion captures the missing
evidence for a source review: exact prompts/responses, actual processor patch
tensors/grid, a reconstructed model-input raster, and production ROI crops.

Prepare on any host (no processor, OCR or model is loaded):
  python tools/qc/vision_source_pilot.py prepare --bundle BUNDLE --library LIBRARY \
      --out PILOT --case PAPER_HASH/FIGURE_ID/LABEL [--case ...]
Copy PILOT to the model host, then run with cached weights only:
  python tools/qc/vision_source_pilot.py run --pilot PILOT --device cuda \
      --dtype bfloat16 --expect-revision HUGGINGFACE_COMMIT

Keep the original pilot directory: run creates a separate capture/ directory
and refuses to overwrite an earlier capture. Successful numeric checks do NOT
establish scientific panel correctness. Complete each review.json against the
source-page image: target content, label, scale/context and neighboring panels.
Old receipts are comparison observations, never gold labels.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.qc.vlm_dtype_probe import caption_panel_labels, machine_report


def _write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def _hash(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _overlay(image, records, output, key="roi_px"):
    from PIL import ImageDraw
    image = image.convert("RGB")
    draw = ImageDraw.Draw(image)
    for i, record in enumerate(records):
        box = record.get(key)
        if box and len(box) == 4 and box[2] > box[0] and box[3] > box[1]:
            draw.rectangle(box, outline="red", width=3)
            draw.text((box[0]+3, box[1]+3), str(record.get("label") or i), fill="red")
    image.save(output)


def prepare(args):
    import fitz
    from PIL import Image
    from mcpsrv.figure_cache import figure_path

    bundle, library, output = args.bundle.resolve(), args.library.resolve(), args.out.resolve()
    if output.exists():
        raise ValueError("Use a new --out directory; prior evidence is never overwritten")
    if output.is_relative_to(bundle) or output.is_relative_to(library):
        raise ValueError("Pilot output must be outside the bundle and source library")
    output.mkdir(parents=True)
    cases = []
    for selector in args.case:
        parts = selector.split("/")
        if len(parts) != 3 or not all(re.fullmatch(r"[A-Za-z0-9_.-]+", p) for p in parts):
            raise ValueError("Each --case must be PAPER_HASH/FIGURE_ID/LABEL")
        paper_hash, figure_id, label = parts
        hd = bundle / "documents" / paper_hash
        summary = json.loads((hd / "summary.json").read_text())
        figures_path = hd / "figures.json"
        record = next(f for f in json.loads(figures_path.read_text())["figures"]
                      if f["figure_id"] == figure_id)
        source = next((library / p for p in summary["relative_paths"] if (library / p).is_file()), None)
        if source is None or not source.resolve().is_relative_to(library):
            raise ValueError(f"Original library PDF unavailable for {selector}")
        digest = _hash(source)
        if digest != summary["pdf_hash_full"]:
            raise ValueError(f"Original PDF differs from the selected receipt: {source}")
        page_number = record.get("source_page")
        if page_number is None:
            scan_path = hd / "scan_detection.json"
            scan = json.loads(scan_path.read_text()) if scan_path.is_file() else {}
            pages = scan.get("keeppages_selected") or []
            page_number = pages[record["page"]-1] if pages else record["page"]
        image_path = figure_path(hd, record)
        name = "-".join(parts)
        case_dir = output / name
        case_dir.mkdir()
        raster_name = "source-raster" + image_path.suffix
        shutil.copyfile(image_path, case_dir / raster_name)
        with Image.open(image_path) as raster:
            size = list(raster.size)
            _overlay(raster, record.get("rois", []), case_dir / "baseline-rois.png")
        with fitz.open(source) as pdf:
            pdf[page_number-1].get_pixmap(dpi=120).save(case_dir / "source-page.png")
        _write(case_dir / "baseline-figure.json", record)
        expected = caption_panel_labels(record)
        if label not in expected or len(expected) < 2:
            raise ValueError(f"Target label not in the multi-panel caption inventory: {selector}")
        case = {"case": name, "paper_hash": paper_hash, "figure_id": figure_id,
                "target_label": label, "expected_labels": expected,
                "raster": raster_name, "raster_size_px": size,
                "raster_sha256": _hash(image_path), "original_pdf": str(source),
                "original_pdf_sha256": digest, "physical_page": page_number,
                "figures_receipt_sha256": _hash(figures_path),
                "caption": record.get("caption_text") or ""}
        _write(case_dir / "case.json", case)
        cases.append(case)
    manifest = {"policy": "local-vlm-source-pilot-v1", "bundle": str(bundle),
                "bundle_manifest_sha256": _hash(bundle / "bundle_manifest.json"),
                "cases": cases, "model_run": False,
                "scope": "Fresh ROI inference on copied existing rasters; not a full extraction rebuild. "
                         "Baseline receipts may differ from the issue snapshot. Original seeded pilot "
                         "pass/failure controls must be supplied explicitly, never inferred from bounds."}
    _write(output / "pilot.json", manifest)
    return manifest


def reconstruct_model_frame(pixel_values, grid, image_processor):
    """Invert Qwen2/2.5-VL's documented patch layout, not a guessed resize.

    The tensor itself is also saved. Refuse an unfamiliar shape rather than
    labeling an arbitrary source resize as the actual model input.
    """
    import numpy as np
    from PIL import Image

    t, gh, gw = map(int, grid)
    p, m, temporal = (image_processor.patch_size, image_processor.merge_size,
                      image_processor.temporal_patch_size)
    values = np.asarray(pixel_values)
    if t != 1 or gh % m or gw % m or values.shape != (gh*gw, 3*temporal*p*p):
        raise ValueError(f"Unsupported still-image patch tensor/grid: {values.shape}, {grid}")
    frames = values.reshape(t,gh//m,gw//m,m,m,3,temporal,p,p)
    frames = frames.transpose(5,0,6,1,3,7,2,4,8).reshape(3,temporal,gh*p,gw*p)
    if not np.array_equal(frames[:,0], frames[:,-1]):
        raise ValueError("Unexpected distinct temporal frames for one source image")
    pixels = frames[:,0].astype(np.float32)
    if image_processor.do_normalize:
        pixels = (pixels * np.asarray(image_processor.image_std)[:,None,None]
                  + np.asarray(image_processor.image_mean)[:,None,None])
    if image_processor.do_rescale:
        pixels /= image_processor.rescale_factor
    return Image.fromarray(np.clip(np.rint(pixels),0,255).astype('uint8').transpose(1,2,0))


class ProcessorCapture:
    """Transparent proxy: capture production calls without changing inference."""
    def __init__(self, processor, directory):
        self.inner, self.directory = processor, directory
        self.calls, self.responses = [], []

    def __getattr__(self, name):
        return getattr(self.inner, name)

    def __call__(self, **kwargs):
        result = self.inner(**kwargs)
        call = len(self.calls)+1
        prefix = self.directory / f"processor-{call:02d}"
        import torch
        pixels = result.pixel_values.detach().cpu()
        grid = result.image_grid_thw.detach().cpu()
        tensors = {name:value.detach().cpu() for name,value in result.items()
                   if hasattr(value,'detach')}
        torch.save(tensors, str(prefix)+"-patches.pt")
        frame = reconstruct_model_frame(pixels.float().numpy(),grid[0].tolist(),self.inner.image_processor)
        frame.save(str(prefix)+"-frame.png")
        (Path(str(prefix)+"-prompt.txt")).write_text(kwargs['text'][0])
        entry = {"call":call, "grid_thw":grid.tolist(), "size_px":list(frame.size),
                 "patch_tensor_shape":list(pixels.shape), "patch_tensor_dtype":str(pixels.dtype),
                 "image_processor":self.inner.image_processor.to_dict()}
        _write(Path(str(prefix)+".json"),entry)
        self.calls.append(entry)
        return result

    def batch_decode(self, generated, **kwargs):
        decoded = self.inner.batch_decode(generated, **kwargs)
        self.responses.append({"processor_call":len(self.calls), "raw_text":decoded,
                               "generated_token_ids":generated.detach().cpu().tolist(),
                               "decode_options":kwargs})
        _write(self.directory / "raw-responses.json",self.responses)
        return decoded


class ModelCapture:
    def __init__(self, model, processor_capture):
        self.inner, self.capture = model, processor_capture
        self.calls = []

    def __getattr__(self, name):
        return getattr(self.inner, name)

    def generate(self, **kwargs):
        self.calls.append({"processor_call":len(self.capture.calls),
                           "max_new_tokens":kwargs['max_new_tokens'],"do_sample":kwargs['do_sample'],
                           "image_grid_thw":kwargs['image_grid_thw'].detach().cpu().tolist()})
        _write(self.capture.directory / "generation-calls.json",self.calls)
        return self.inner.generate(**kwargs)


def run(args):
    # Keep reproducibility separate from cache population. Downloads must be
    # an explicit prior operator action, never an accidental pilot side effect.
    os.environ['HF_HUB_OFFLINE'] = os.environ['TRANSFORMERS_OFFLINE'] = '1'
    import torch
    from PIL import Image
    from pipeline.figures import detect_figure_rois_via_vision
    from pipeline.model_provenance import vision_producer
    from pipeline.vision import LocalVLMBackend

    if args.device == 'cuda' and not torch.cuda.is_available():
        raise ValueError("CUDA unavailable; do not fall back to a CPU model run")
    if args.device == 'mps' and not torch.backends.mps.is_available():
        raise ValueError("MPS unavailable")
    identity = vision_producer('vision-local', args.model)
    if args.expect_revision and identity.get('revision') != args.expect_revision:
        raise ValueError(f"Cached model revision differs: {identity.get('revision')}")
    pilot = args.pilot.resolve()
    manifest = json.loads((pilot / 'pilot.json').read_text())
    for case in manifest['cases']:
        raster = pilot / case['case'] / case['raster']
        if _hash(raster) != case['raster_sha256']:
            raise ValueError(f"Copied raster changed: {raster}")
    output = pilot / 'capture'
    output.mkdir()  # Refuse overwriting a previous run.
    _write(output / 'environment.json',{'machine':machine_report(),'model_before_load':identity,
           'dtype':args.dtype,'device':args.device,'tool_sha256':_hash(Path(__file__))})
    backend = LocalVLMBackend(model=args.model,device=args.device,dtype=args.dtype)
    if args.expect_revision and backend.producer.get('revision') != args.expect_revision:
        raise ValueError("Loaded model revision differs from the requested source-pilot revision")
    _write(output / 'producer.json',backend.producer)
    processor, model = backend._processor, backend._model
    checks = []
    for case in manifest['cases']:
        source_dir = pilot / case['case']
        image_path = source_dir / case['raster']
        if _hash(image_path) != case['raster_sha256']:
            raise ValueError(f"Copied raster changed: {image_path}")
        destination = output / case['case']
        destination.mkdir()
        capture = ProcessorCapture(processor,destination)
        backend._processor = capture
        backend._model = ModelCapture(model,capture)
        result = detect_figure_rois_via_vision(image_path,
            [{'label':label} for label in case['expected_labels']],backend,caption_text=case['caption'])
        _write(destination / 'production-result.json',result)
        valid = bool(capture.calls) and result['pass3_status']=='completed'
        with Image.open(image_path) as raster:
            _overlay(raster,result.get('rois',[]),destination / 'fresh-rois.png')
            for index,roi in enumerate(result.get('rois',[])):
                box = roi['roi_px']
                bounded = 0 <= box[0] < box[2] <= raster.width and 0 <= box[1] < box[3] <= raster.height
                provenance = (roi.get('coordinate_provenance') or {}).get('panel',{})
                frame_ok = (provenance.get('raster_size_px')==list(raster.size)
                            and provenance.get('model_input_size_px')==capture.calls[-1]['size_px'])
                valid = valid and bounded and frame_ok
                if bounded:
                    raster.crop(box).save(destination / f'crop-{index:02d}.png')
        review = {'target_label':case['target_label'],'numeric_checks_passed':bool(valid),
                  'scientific_source_review':'pending','complete_target_content':None,
                  'label_preserved':None,'scale_and_context_preserved':None,
                  'neighboring_panel_substitution_absent':None,'known_good_control_preserved':None,
                  'reviewer':None,'notes':None}
        _write(destination / 'review.json',review)
        checks.append({'case':case['case'],**review})
    _write(output / 'summary.json',{'checks':checks,'gpu_acceptance':'source review still required'})
    return all(c['numeric_checks_passed'] for c in checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest='command',required=True)
    prep = commands.add_parser('prepare')
    prep.add_argument('--bundle',type=Path,required=True)
    prep.add_argument('--library',type=Path,required=True)
    prep.add_argument('--out',type=Path,required=True)
    prep.add_argument('--case',action='append',required=True)
    capture = commands.add_parser('run')
    capture.add_argument('--pilot',type=Path,required=True)
    capture.add_argument('--model',default='Qwen/Qwen2.5-VL-7B-Instruct')
    capture.add_argument('--expect-revision')
    capture.add_argument('--device',choices=['cuda','mps'],default='cuda')
    capture.add_argument('--dtype',choices=['float32','bfloat16','float16'],default='bfloat16')
    args = parser.parse_args()
    if args.command == 'prepare':
        prepare(args)
        print(f'Prepared {args.out}/pilot.json; no model loaded or run.')
        return 0
    return 0 if run(args) else 1


if __name__ == '__main__':
    raise SystemExit(main())
