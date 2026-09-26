"""Render frozen review regions from hash-verified original PDFs; no OCR."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    import fitz

    review = Path(__file__).resolve().parent
    population = json.loads((review / "population.json").read_text())
    crops = json.loads((review / "crops.json").read_text())
    if sha256(review / "population.json") != crops["population_sha256"]:
        raise SystemExit("Frozen population does not match the crop manifest")
    rows = {row["id"]: row for row in population["rows"]}
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, required=True,
                        help="Library root containing the original B/, D_E/, H_J/, M/ PDF paths")
    parser.add_argument("--output", type=Path, required=True,
                        help="Directory for rendered PNGs; existing differing files are refused")
    parser.add_argument("--rows", nargs="+", choices=sorted(rows), default=list(rows))
    args = parser.parse_args()

    sources = {}
    for key in args.rows:
        row = rows[key]
        # Historical absolute paths remain immutable provenance; only the
        # original path relative to library/ is needed on another machine.
        relative = row["source"].split("/library/", 1)[1]
        source = args.library / relative
        if source not in sources:
            sources[source] = sha256(source)
        if sources[source] != row["source_sha256"]:
            raise SystemExit(f"Source hash mismatch: {source}")

    args.output.mkdir(parents=True, exist_ok=True)
    rendered = set()
    for key in args.rows:
        row, crop = rows[key], crops["crops"][key]
        target = args.output / crop["path"]
        if target.exists():
            if sha256(target) != crop["sha256"]:
                raise SystemExit(f"Refusing to overwrite differing crop: {target}")
        else:
            source = args.library / row["source"].split("/library/", 1)[1]
            with fitz.open(source) as document:
                page = document[crop["physical_page"] - 1]
                # Match the capture sequence: native-line inspection happened
                # before rendering and can populate MuPDF's font resources.
                page.get_text("dict")
                png = page.get_pixmap(clip=fitz.Rect(crop["bbox_top_left"]),
                                      dpi=crop["dpi"], alpha=False).tobytes("png")
            if hashlib.sha256(png).hexdigest() != crop["sha256"]:
                version = importlib.metadata.version("pymupdf")
                raise SystemExit(f"Render hash differs for {key}; using PyMuPDF {version}, "
                                 "original capture used 1.28.0. No crop written.")
            target.write_bytes(png)
        rendered.add(target)
    print(f"Verified {len(args.rows)} regions, {len(rendered)} PNGs, {len(sources)} source PDFs")


if __name__ == "__main__":
    main()
