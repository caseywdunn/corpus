"""Reproduce selection and reviewed-region rasters; never read candidate outputs."""
import argparse
import hashlib
import json
from pathlib import Path
import runpy


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('library',type=Path,help='Library repository containing library/ and transcriptions/')
    parser.add_argument('output',type=Path,help='New, empty output directory')
    args=parser.parse_args()
    here=Path(__file__).parent
    args.output.mkdir(parents=True,exist_ok=False)
    policy=json.loads((here/'selection-policy.json').read_text())
    manifest=args.library/'transcriptions/sources.json'
    assert hashlib.sha256(manifest.read_bytes()).hexdigest()==policy['source_manifest_sha256']
    inputs=json.loads((here/'transcription-inputs.json').read_text())
    for relative,expected in inputs.items():
        assert hashlib.sha256((args.library/relative).read_bytes()).hexdigest()==expected,relative
    namespace=runpy.run_path(str(here/'select.py'))
    namespace['main'].__globals__.update(ROOT=args.library,OUT=args.output)
    namespace['main']()
    assert (args.output/'selected.json').read_bytes()==(here/'selected.json').read_bytes()
    # Labels were frozen by image review, not generated from extraction or OCR.
    # Re-render their exact regions for a reviewer; never overwrite those labels.
    import fitz
    for label in json.loads((here/'labels.json').read_text()):
        source=args.library/label['source_pdf']
        with source.open('rb') as handle:
            assert hashlib.file_digest(handle,'sha256').hexdigest()==label['source_sha256']
        region=label['review']['rendered_source']
        with fitz.open(source) as pdf:
            page=pdf[label['physical_page']-1]
            page.get_pixmap(clip=fitz.Rect(region['bbox_top_left_pdf_points']),dpi=region['dpi']).save(
                args.output/f"review-{label['sample_number']:02}.png")


if __name__=='__main__':
    main()
