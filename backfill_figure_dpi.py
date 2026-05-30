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

Two modes:
- ``--scale S`` (default): a single fixed scale for every figure
  (DPI = 72*S).
- ``--native``: preserve each figure's *native source resolution* — render
  at the densest overlapping embedded image's pixel density, which varies
  per figure and per scan. Vector figures have no native resolution, so
  they fall back to ``--vector-dpi`` (default 300). ``--max-dpi`` optionally
  caps pathological full-page scans.

Usage:
    python backfill_figure_dpi.py <output_dir> --native --only <HASH> --dry-run
    python backfill_figure_dpi.py <output_dir> --native               # full bundle
    python backfill_figure_dpi.py <output_dir> --native --vector-dpi 300 --max-dpi 600
    python backfill_figure_dpi.py <output_dir> --scale 4.0            # fixed mode
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Optional

from pipeline.figures import render_figures

logger = logging.getLogger("backfill_figure_dpi")


def _backfill_paper(hash_dir: Path, *, native: bool, fixed_scale: float,
                    vector_dpi: float, max_dpi, dry_run: bool) -> dict:
    """Re-render one paper's figures from its processed.pdf, via the shared
    pipeline.figures.render_figures — the exact logic the extraction stage's
    native default uses. Writes figures.json back when anything rendered."""
    figures_json = hash_dir / "figures.json"
    data = json.loads(figures_json.read_text(encoding="utf-8"))
    figures = data.get("figures", []) or []

    src_pdf = hash_dir / "processed.pdf"
    if not src_pdf.is_file():
        logger.warning("%s: no processed.pdf — skipping %d figures",
                       hash_dir.name, len(figures))
        return {"rendered": 0, "skipped_no_bbox": 0, "skipped_method": 0,
                "skipped_no_src": len(figures), "errors": 0}

    stats = render_figures(
        src_pdf, figures, hash_dir / "figures",
        native=native, fixed_scale=fixed_scale, vector_dpi=vector_dpi,
        max_dpi=max_dpi, dry_run=dry_run, label_prefix=f"{hash_dir.name}/",
    )
    stats["skipped_no_src"] = 0
    if not dry_run and stats["rendered"]:
        figures_json.write_text(
            json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8",
        )
    return stats

def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("output_dir", type=Path,
                        help="Corpus output dir (contains documents/<hash>/).")
    parser.add_argument("--native", action="store_true",
                        help="Preserve each figure's native source resolution: "
                             "render at the densest overlapping embedded image's "
                             "DPI (varies per figure/scan), falling back to "
                             "--vector-dpi for resolution-less vector figures. "
                             "Overrides --scale.")
    parser.add_argument("--scale", type=float, default=2.0,
                        help="Fixed mode: target images_scale (DPI = 72*scale). "
                             "Default 2.0. Ignored when --native is set.")
    parser.add_argument("--vector-dpi", type=float, default=300.0,
                        help="--native fallback DPI for vector figures (no native "
                             "resolution exists). Default 300.")
    parser.add_argument("--max-dpi", type=float, default=None,
                        help="--native ceiling, to bound pathological full-page "
                             "scans. Default: uncapped.")
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

    if not args.native and not 1.0 <= args.scale <= 8.0:
        logger.error("--scale must be in [1.0, 8.0] (got %s)", args.scale)
        return 2
    if args.native and args.vector_dpi <= 0:
        logger.error("--vector-dpi must be positive (got %s)", args.vector_dpi)
        return 2

    totals = {"rendered": 0, "skipped_no_bbox": 0, "skipped_method": 0,
              "skipped_no_src": 0, "errors": 0, "papers": 0}
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
        s = _backfill_paper(
            hash_dir, native=args.native, fixed_scale=args.scale,
            vector_dpi=args.vector_dpi, max_dpi=args.max_dpi, dry_run=args.dry_run,
        )
        for k, v in s.items():
            totals[k] += v

    verb = "would re-render" if args.dry_run else "re-rendered"
    mode_desc = (
        f"native (vector fallback {args.vector_dpi:.0f} dpi"
        + (f", capped {args.max_dpi:.0f} dpi)" if args.max_dpi else ")")
        if args.native else f"fixed scale {args.scale:.1f} ({round(72 * args.scale)} dpi)"
    )
    logger.info(
        "Done. %d papers; %s %d figures — mode: %s; "
        "%d skipped (no bbox), %d skipped (no source), %d errors.",
        totals["papers"], verb, totals["rendered"], mode_desc,
        totals["skipped_no_bbox"], totals["skipped_no_src"], totals["errors"],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
