#!/usr/bin/env python3
"""Reproduce only frozen source-review renders; no OCR, extraction or mutation."""
import argparse
import hashlib
import json
from pathlib import Path

import fitz


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    if args.out.resolve() == here or here in args.out.resolve().parents:
        parser.error("choose an output outside the frozen evidence directory")
    protocol = json.loads((here / "protocol.frozen.json").read_text())
    review = json.loads((here / "source-review.json").read_text())
    full = json.loads((here / "render_receipt.json").read_text())
    sources = {row["paper_hash"]: row for row in protocol["sources"]}
    for row in sources.values():
        path = args.library / row["source_relative_path"]
        if sha256(path) != row["source_sha256"]:
            raise ValueError(f"source hash mismatch: {path}")
    args.out.mkdir(parents=True, exist_ok=True)
    tasks = [
        (Path(row["crop_path"]).name, row, row["render_bbox_points"])
        for row in full["renders"]
    ]
    tasks.extend(
        (Path(row["crop"]).name, row, row["bbox_points"])
        for row in review["regions"].values()
    )
    results = []
    for name, row, box in tasks:
        source = sources[row["paper_hash"]]
        with fitz.open(args.library / source["source_relative_path"]) as pdf:
            page = pdf[row["physical_page"] - 1]
            image = page.get_pixmap(clip=fitz.Rect(box), dpi=row["dpi"], alpha=False)
            out = args.out / name
            image.save(out)
        actual = sha256(out)
        results.append({"file": name, "sha256": actual,
                        "expected_sha256": row["crop_sha256"],
                        "matches": actual == row["crop_sha256"]})
    receipt = {"pymupdf": fitz.VersionBind, "original_pymupdf": full["pymupdf"],
               "renders": results}
    (args.out / "reproduction.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"{sum(row['matches'] for row in results)}/{len(results)} rendered hashes match")
    return 0 if all(row["matches"] for row in results) else 2


if __name__ == "__main__":
    raise SystemExit(main())
