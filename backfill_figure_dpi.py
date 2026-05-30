#!/usr/bin/env python3
"""Re-render existing figure PNGs at higher DPI from stored bboxes (#121).

The pipeline previously rasterized docling figure crops at docling's
default ``images_scale=1.0`` — 72 dpi — and the saved PNG is never
resized downstream, so older bundles carry grainy figures that upscale
poorly in print (Quarto / Word). The extraction default is now
``figures.images_scale=2.0``, but that only affects *future* ingests.

This script lifts an *existing* bundle without re-running docling: every
figure record in ``figures.json`` already stores ``bbox`` +
``bbox_coord_system`` + ``page``, and the rendered source PDF is on disk
at ``<hash>/processed.pdf``. We re-render each figure's page region with
PyMuPDF at the target scale and overwrite the PNG in place. Pixel
dimensions in ``figures.json`` (only the PyMuPDF-path records carry
``width``/``height``) are refreshed when present.

Best-effort and idempotent-ish: a figure with no usable bbox/coord
system, or whose source PDF is missing, is skipped (counted, not
fatal). The re-rendered crop is bbox-exact and may frame slightly
differently from docling's original crop — acceptable for a quality
backfill; validate on one paper before a full run.

Usage:
    python backfill_figure_dpi.py <output_dir> [--scale 2.0] [--dry-run]
    python backfill_figure_dpi.py <output_dir> --scale 4.0 --only <HASH>
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Optional

logger = logging.getLogger("backfill_figure_dpi")

# Real figure crops worth re-rendering; skip journal furniture etc.
_RENDERABLE_TYPES = {"figure", "plate", "subpanel", None}


def _rect_for(bbox, coord_system: str, page_height: float):
    """Translate a stored bbox to a PyMuPDF (top-left origin) rect.

    docling stores ``[l, b, r, t]`` in bottom-left origin (``t > b``);
    PyMuPDF's path already stored top-left ``[x0, y0, x1, y1]``.
    """
    import fitz  # PyMuPDF

    l, a, r, b = (float(v) for v in bbox)
    if coord_system == "pdf_pts_bottom_left":
        # flip y: top-left y = page_height - (bottom-left y)
        return fitz.Rect(l, page_height - b, r, page_height - a)
    if coord_system == "pdf_pts_top_left":
        return fitz.Rect(l, a, r, b)
    return None  # unknown coord system → caller skips


def _backfill_paper(hash_dir: Path, scale: float, dry_run: bool) -> dict:
    import fitz  # PyMuPDF

    stats = {"rendered": 0, "skipped_no_bbox": 0, "skipped_no_src": 0, "errors": 0}
    figures_json = hash_dir / "figures.json"
    data = json.loads(figures_json.read_text(encoding="utf-8"))
    figures = data.get("figures", []) or []

    src_pdf = hash_dir / "processed.pdf"
    if not src_pdf.is_file():
        stats["skipped_no_src"] = len(figures)
        logger.warning("%s: no processed.pdf — skipping %d figures",
                       hash_dir.name, len(figures))
        return stats

    doc = fitz.open(src_pdf)
    try:
        changed = False
        for fig in figures:
            bbox = fig.get("bbox")
            coord = fig.get("bbox_coord_system")
            page = fig.get("page")
            fname = fig.get("filename")
            ftype = fig.get("figure_type")
            if ftype not in _RENDERABLE_TYPES:
                continue
            if not bbox or not coord or not page or not fname:
                stats["skipped_no_bbox"] += 1
                continue
            page_idx = int(page) - 1  # stored page is 1-based
            if page_idx < 0 or page_idx >= doc.page_count:
                stats["skipped_no_bbox"] += 1
                continue
            pg = doc[page_idx]
            rect = _rect_for(bbox, coord, pg.rect.height)
            if rect is None or rect.is_empty or rect.width <= 0 or rect.height <= 0:
                stats["skipped_no_bbox"] += 1
                continue
            png_path = hash_dir / "figures" / fname
            try:
                if dry_run:
                    new_w = round(rect.width * scale)
                    new_h = round(rect.height * scale)
                    logger.info("[dry-run] %s/%s → %d×%d px (scale %.1f)",
                                hash_dir.name, fname, new_w, new_h, scale)
                else:
                    pix = pg.get_pixmap(matrix=fitz.Matrix(scale, scale), clip=rect)
                    png_path.parent.mkdir(parents=True, exist_ok=True)
                    pix.save(png_path)
                    if "width" in fig or "height" in fig:
                        fig["width"], fig["height"] = pix.width, pix.height
                        changed = True
                    fig["images_scale"] = scale  # provenance of the re-render
                    changed = True
                stats["rendered"] += 1
            except Exception as e:  # one bad figure must not abort the paper
                logger.warning("%s/%s: render failed: %s", hash_dir.name, fname, e)
                stats["errors"] += 1
        if changed and not dry_run:
            figures_json.write_text(
                json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8",
            )
    finally:
        doc.close()
    return stats


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("output_dir", type=Path,
                        help="Corpus output dir (contains documents/<hash>/).")
    parser.add_argument("--scale", type=float, default=2.0,
                        help="Target images_scale (DPI = 72*scale). Default 2.0.")
    parser.add_argument("--only", default=None, metavar="HASH",
                        help="Re-render a single paper hash (validate before a "
                             "full run).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report intended renders without writing PNGs.")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    try:
        import fitz  # noqa: F401 — fail fast with a clear message
    except ImportError:
        logger.error("PyMuPDF (fitz) is required: pip install pymupdf")
        return 1

    docs = args.output_dir / "documents"
    if not docs.is_dir():
        logger.error("no documents/ under %s", args.output_dir)
        return 2

    if not 1.0 <= args.scale <= 8.0:
        logger.error("--scale must be in [1.0, 8.0] (got %s)", args.scale)
        return 2

    totals = {"rendered": 0, "skipped_no_bbox": 0, "skipped_no_src": 0, "errors": 0,
              "papers": 0}
    hash_dirs = sorted(d for d in docs.iterdir() if d.is_dir())
    if args.only:
        hash_dirs = [d for d in hash_dirs if d.name == args.only]
        if not hash_dirs:
            logger.error("no documents/%s", args.only)
            return 2

    for hash_dir in hash_dirs:
        if not (hash_dir / "figures.json").is_file():
            continue
        totals["papers"] += 1
        s = _backfill_paper(hash_dir, args.scale, args.dry_run)
        for k, v in s.items():
            totals[k] += v

    verb = "would re-render" if args.dry_run else "re-rendered"
    logger.info(
        "Done. %d papers; %s %d figures at scale %.1f (%d dpi); "
        "%d skipped (no bbox), %d skipped (no source), %d errors.",
        totals["papers"], verb, totals["rendered"], args.scale,
        round(72 * args.scale), totals["skipped_no_bbox"],
        totals["skipped_no_src"], totals["errors"],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
