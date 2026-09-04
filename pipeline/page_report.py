"""Generate an on-demand, self-contained page audit report (#274).

The report joins evidence that otherwise lives in four places: the rendered
``processed.pdf`` page, PDF word cells, Docling's selectable parsed text, and
``figures.json`` caption/figure/ROI decisions. It is operator-side output and
is deliberately never generated during a normal corpus build or copied into a
served bundle.

Usage::

    python -m pipeline.page_report OUTPUT/documents/<HASH>
    python -m pipeline.page_report OUTPUT/documents/<HASH> --pages 12-13,17

The output is one HTML file with bounded JPEG page renders embedded as data
URLs. No parallel ``visualizations/`` raster tree is created.
"""
from __future__ import annotations

import argparse
import base64
import html
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .figures import caption_evidence_summary


_DEFAULT_MAX_PAGES = 200
_DEFAULT_MAX_WIDTH = 1400
_DEFAULT_JPEG_QUALITY = 84


def _read_json(path: Path, default):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return default


def parse_page_spec(spec: Optional[str], total_pages: int) -> List[int]:
    """Parse ``1-3,8`` into validated, ordered one-based page numbers."""
    if not spec:
        return list(range(1, total_pages + 1))
    selected = set()
    for raw_part in spec.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if "-" in part:
            start_text, end_text = part.split("-", 1)
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as exc:
                raise ValueError(f"invalid page range {part!r}") from exc
            if end < start:
                raise ValueError(f"descending page range {part!r}")
            selected.update(range(start, end + 1))
        else:
            try:
                selected.add(int(part))
            except ValueError as exc:
                raise ValueError(f"invalid page number {part!r}") from exc
    outside = sorted(page for page in selected if page < 1 or page > total_pages)
    if outside:
        raise ValueError(
            f"page(s) outside processed.pdf's 1-{total_pages} range: {outside}"
        )
    if not selected:
        raise ValueError("page selection is empty")
    return sorted(selected)


def _docling_fragments(docling_data: Dict) -> Dict[int, List[Dict]]:
    """Reconstruct provenance-level text without importing Docling itself."""
    by_page: Dict[int, List[Dict]] = defaultdict(list)
    for item in docling_data.get("texts", []) or []:
        full_text = item.get("text") or ""
        for provenance in item.get("prov", []) or []:
            page = provenance.get("page_no")
            bbox = provenance.get("bbox") or {}
            if page is None:
                continue
            text = full_text
            span = provenance.get("charspan")
            if isinstance(span, list) and len(span) == 2:
                try:
                    start, end = int(span[0]), int(span[1])
                    if 0 <= start < end <= len(full_text):
                        text = full_text[start:end]
                except (TypeError, ValueError):
                    pass
            text = text.strip()
            if not text:
                continue
            by_page[int(page)].append({
                "label": str(item.get("label") or "text"),
                "text": text,
                "bbox": [
                    bbox.get("l"), bbox.get("b"), bbox.get("r"), bbox.get("t"),
                ],
                "coord_origin": bbox.get("coord_origin"),
            })
    return by_page


def _top_left_bbox(
    bbox: Optional[Sequence], page_height: float, coord_system: Optional[str],
) -> Optional[List[float]]:
    if not bbox or len(bbox) != 4:
        return None
    try:
        x0, y0, x1, y1 = (float(value) for value in bbox)
    except (TypeError, ValueError):
        return None
    x0, x1 = sorted((x0, x1))
    y0, y1 = sorted((y0, y1))
    if "bottom" in (coord_system or "").lower():
        y0, y1 = page_height - y1, page_height - y0
    return [x0, y0, x1, y1]


def _figure_page_bbox(figure: Dict, page_height: float) -> Optional[List[float]]:
    coord_system = figure.get("bbox_coord_system")
    if not coord_system:
        coord_system = (
            "pdf_pts_top_left"
            if figure.get("extraction_method") == "pymupdf"
            else "pdf_pts_bottom_left"
        )
    return _top_left_bbox(figure.get("bbox"), page_height, coord_system)


def _roi_page_bbox(
    figure: Dict, roi: Dict, page_height: float,
) -> Optional[List[float]]:
    figure_bbox = _figure_page_bbox(figure, page_height)
    size = figure.get("image_size_px")
    roi_px = roi.get("roi_px")
    if not figure_bbox or not size or len(size) != 2 or not roi_px:
        return None
    try:
        image_width, image_height = float(size[0]), float(size[1])
        rx0, ry0, rx1, ry1 = (float(value) for value in roi_px)
    except (TypeError, ValueError):
        return None
    if image_width <= 0 or image_height <= 0:
        return None
    fx0, fy0, fx1, fy1 = figure_bbox
    return [
        fx0 + min(rx0, rx1) / image_width * (fx1 - fx0),
        fy0 + min(ry0, ry1) / image_height * (fy1 - fy0),
        fx0 + max(rx0, rx1) / image_width * (fx1 - fx0),
        fy0 + max(ry0, ry1) / image_height * (fy1 - fy0),
    ]


def _safe_script_json(value) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":")).replace(
        "<", "\\u003c"
    )


def _render_page_data_url(page, max_width: int, jpeg_quality: int) -> str:
    import fitz

    scale = min(2.0, max_width / max(float(page.rect.width), 1.0))
    scale = max(scale, 0.5)
    pixmap = page.get_pixmap(
        matrix=fitz.Matrix(scale, scale), colorspace=fitz.csRGB, alpha=False,
    )
    encoded = base64.b64encode(
        pixmap.tobytes("jpeg", jpg_quality=jpeg_quality)
    ).decode("ascii")
    return f"data:image/jpeg;base64,{encoded}"


def _figure_caption_html(figure: Dict) -> str:
    evidence = caption_evidence_summary(figure)
    caption = figure.get("caption_text") or figure.get("caption") or ""
    title_bits = [str(figure.get("figure_id") or "unnamed figure")]
    if figure.get("figure_number"):
        title_bits.append(f"Fig. {figure['figure_number']}")
    summary = " · ".join(title_bits)
    parts = [
        "<details class='figure-evidence' open>",
        f"<summary>{html.escape(summary)}</summary>",
        "<div class='evidence-meta'>"
        f"status <strong>{html.escape(str(evidence['caption_status']))}</strong> · "
        f"confidence {html.escape(str(evidence['caption_confidence']))} · "
        f"kind {html.escape(str(evidence['caption_kind']))} · "
        f"source {html.escape(str(figure.get('caption_source')))} · "
        f"page distance {html.escape(str(evidence['caption_page_distance']))}"
        "</div>",
        f"<div class='caption'>{html.escape(caption) if caption else '<em>no caption</em>'}</div>",
    ]
    candidates = figure.get("caption_candidates") or []
    if candidates:
        parts.append("<ol class='candidates'>")
        for candidate in candidates:
            decision = "chosen" if candidate.get("chosen") else "rejected"
            reason = candidate.get("rejection_reason") or ""
            meta = ", ".join(
                str(value) for value in (
                    candidate.get("caption_source"),
                    candidate.get("confidence"),
                    f"page {candidate.get('caption_page')}",
                    f"gap {candidate.get('distance_pts')} pt",
                    reason,
                ) if value not in (None, "")
            )
            parts.append(
                f"<li><strong>{decision}</strong> ({html.escape(meta)}): "
                f"{html.escape(candidate.get('caption_text') or '')}</li>"
            )
        parts.append("</ol>")
    parts.append("</details>")
    return "".join(parts)


_CSS = """
:root { color-scheme: light; font-family: system-ui, sans-serif; color: #202124; }
body { margin: 0; background: #f4f5f7; }
header { position: sticky; top: 0; z-index: 5; padding: .8rem 1.2rem;
  background: rgba(255,255,255,.96); border-bottom: 1px solid #ccc; }
h1 { font-size: 1.2rem; margin: 0 0 .35rem; }
.meta, .stats, .evidence-meta { color: #5f6368; font-size: .86rem; }
.controls { display: flex; flex-wrap: wrap; gap: .8rem; margin-top: .55rem; }
.controls label { font-size: .85rem; }
main { max-width: 1800px; margin: 0 auto; padding: 1rem; }
.page-card { background: white; border: 1px solid #d9dce1; border-radius: 8px;
  margin: 0 0 1.2rem; padding: .9rem; scroll-margin-top: 7rem; }
.page-card h2 { font-size: 1.05rem; margin: 0 0 .3rem; }
.page-grid { display: grid; grid-template-columns: minmax(420px, 1.35fr) minmax(320px, 1fr);
  gap: 1rem; align-items: start; margin-top: .7rem; }
.page-image { position: relative; line-height: 0; border: 1px solid #bbb; background: #ddd; }
.page-image img, .page-image canvas { display: block; width: 100%; height: auto; }
.page-image canvas { position: absolute; inset: 0; pointer-events: none; }
.text-pane { min-width: 0; }
pre { white-space: pre-wrap; overflow-wrap: anywhere; font: .82rem/1.35 ui-monospace, monospace;
  max-height: 38rem; overflow: auto; background: #f7f7f8; padding: .7rem; }
.label { color: #777; user-select: none; }
.figure-evidence { margin: .5rem 0; border-top: 1px solid #eee; padding-top: .4rem; }
.figure-evidence summary { cursor: pointer; font-weight: 600; }
.caption { white-space: pre-wrap; margin: .3rem 0; }
.candidates { font-size: .85rem; padding-left: 1.4rem; }
.empty { color: #777; font-style: italic; }
@media (max-width: 900px) { .page-grid { grid-template-columns: 1fr; } }
"""


_JS = r"""
const layerStyles = {
  cells: {color:'#d93025', width:0.65, dash:[]},
  figures: {color:'#f29900', width:2.2, dash:[]},
  rois: {color:'#1a73e8', width:2.0, dash:[]},
  chosen: {color:'#188038', width:2.0, dash:[7,4]},
  rejected: {color:'#a142f4', width:1.6, dash:[4,4]}
};
function layerEnabled(name) {
  const box = document.querySelector(`[data-layer='${name}']`);
  return box && box.checked;
}
function drawCard(card) {
  const img = card.querySelector('img');
  const canvas = card.querySelector('canvas');
  const payload = JSON.parse(card.querySelector("script[type='application/json']").textContent);
  if (!img.complete || !img.naturalWidth) return;
  canvas.width = img.naturalWidth; canvas.height = img.naturalHeight;
  const ctx = canvas.getContext('2d');
  const sx = canvas.width / payload.page_width;
  const sy = canvas.height / payload.page_height;
  for (const layer of ['cells','figures','rois','chosen','rejected']) {
    if (!layerEnabled(layer)) continue;
    const style = layerStyles[layer];
    ctx.strokeStyle = style.color; ctx.fillStyle = style.color;
    ctx.lineWidth = style.width * Math.max(sx, sy); ctx.setLineDash(style.dash);
    for (const box of payload[layer]) {
      const [x0,y0,x1,y1,label] = box;
      ctx.strokeRect(x0*sx,y0*sy,(x1-x0)*sx,(y1-y0)*sy);
      if (label && layer !== 'cells') {
        ctx.font = `${Math.max(11, 10*sx)}px system-ui`;
        ctx.fillText(label, x0*sx + 3, Math.max(12, y0*sy - 3));
      }
    }
  }
  ctx.setLineDash([]);
}
function drawAll() { document.querySelectorAll('.page-image').forEach(drawCard); }
document.querySelectorAll("input[data-layer]").forEach(box => box.addEventListener('change', drawAll));
document.querySelectorAll('.page-image img').forEach(img => img.addEventListener('load', drawAll));
drawAll();
"""


def generate_page_report(
    hash_dir: Path,
    output_path: Optional[Path] = None,
    pages: Optional[Iterable[int]] = None,
    max_pages: int = _DEFAULT_MAX_PAGES,
    max_width: int = _DEFAULT_MAX_WIDTH,
    jpeg_quality: int = _DEFAULT_JPEG_QUALITY,
) -> Path:
    """Generate and return a self-contained audit report for one document."""
    import fitz

    hash_dir = Path(hash_dir).expanduser().resolve()
    pdf_path = hash_dir / "processed.pdf"
    docling_path = hash_dir / "docling_doc.json"
    missing = [str(path.name) for path in (pdf_path, docling_path) if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            f"{hash_dir} is missing required artifact(s): {', '.join(missing)}"
        )

    docling_data = _read_json(docling_path, {})
    figures_data = _read_json(hash_dir / "figures.json", {})
    metadata = _read_json(hash_dir / "metadata.json", {})
    scan = _read_json(hash_dir / "scan_detection.json", {})
    fragments_by_page = _docling_fragments(docling_data)
    figures = figures_data.get("figures", []) or []
    figures_by_page: Dict[int, List[Dict]] = defaultdict(list)
    candidates_by_page: Dict[int, List[Tuple[Dict, Dict]]] = defaultdict(list)
    for figure in figures:
        if figure.get("page") is not None:
            figures_by_page[int(figure["page"])].append(figure)
        for candidate in figure.get("caption_candidates", []) or []:
            if candidate.get("caption_page") is not None:
                candidates_by_page[int(candidate["caption_page"])].append(
                    (figure, candidate)
                )

    pdf = fitz.open(pdf_path)
    try:
        total_pages = pdf.page_count
        selected = (
            sorted(set(int(page) for page in pages))
            if pages is not None else list(range(1, total_pages + 1))
        )
        outside = [page for page in selected if page < 1 or page > total_pages]
        if outside:
            raise ValueError(f"page(s) outside 1-{total_pages}: {outside}")
        if not selected:
            raise ValueError("page selection is empty")
        if max_pages < 0:
            raise ValueError("max_pages must be zero or positive")
        if max_pages > 0 and len(selected) > max_pages:
            raise ValueError(
                f"report would render {len(selected)} pages; select pages or pass "
                f"--max-pages 0 to override the {max_pages}-page safety limit"
            )
        if not 400 <= max_width <= 4000:
            raise ValueError("max_width must be between 400 and 4000 pixels")
        if not 30 <= jpeg_quality <= 100:
            raise ValueError("jpeg_quality must be between 30 and 100")

        title = metadata.get("title") or metadata.get("filename") or hash_dir.name
        author_text = ", ".join(
            " ".join(filter(None, (author.get("forename"), author.get("surname"))))
            for author in metadata.get("authors", []) or []
        )
        common_provenance = " · ".join(filter(None, (
            f"scan: {scan.get('file_type')}" if scan.get("file_type") else None,
            f"OCR: {scan.get('ocr_mode')}" if scan.get("ocr_mode") else None,
            (
                "packs: " + ", ".join(scan.get("tesseract_packs") or [])
                if scan.get("tesseract_packs") else None
            ),
            (
                f"Docling schema: {docling_data.get('version')}"
                if docling_data.get("version") else None
            ),
            (
                f"pipeline: {figures_data.get('pipeline_version')}"
                if figures_data.get("pipeline_version") else None
            ),
        )))
        parts = [
            "<!doctype html><html><head><meta charset='utf-8'>",
            f"<title>Page audit — {html.escape(str(title))}</title>",
            f"<style>{_CSS}</style></head><body><header>",
            f"<h1>Page audit — {html.escape(str(title))}</h1>",
            "<div class='meta'>"
            + html.escape(" · ".join(filter(None, (
                author_text,
                str(metadata.get("year") or ""),
                f"hash {hash_dir.name}",
                f"{len(selected)} of {total_pages} processed pages",
            ))))
            + "</div>",
            f"<div class='meta'>{html.escape(common_provenance)}</div>",
            "<div class='controls'>",
            "<label><input type='checkbox' data-layer='cells'> PDF word cells (red)</label>",
            "<label><input type='checkbox' data-layer='figures' checked> figures (orange)</label>",
            "<label><input type='checkbox' data-layer='rois' checked> ROIs (blue)</label>",
            "<label><input type='checkbox' data-layer='chosen' checked> chosen captions (green)</label>",
            "<label><input type='checkbox' data-layer='rejected' checked> rejected captions (purple)</label>",
            "</div></header><main>",
        ]

        page_selection = scan.get("keeppages_selected") or []
        for page_number in selected:
            page = pdf[page_number - 1]
            page_width, page_height = float(page.rect.width), float(page.rect.height)
            words = page.get_text("words") or []
            cell_boxes = [
                [float(word[0]), float(word[1]), float(word[2]), float(word[3]), ""]
                for word in words
            ]
            page_figures = figures_by_page.get(page_number, [])
            figure_boxes, roi_boxes = [], []
            seen_figure_boxes = set()
            roi_count = rendered_roi_count = 0
            for figure in page_figures:
                box = _figure_page_bbox(figure, page_height)
                if box:
                    key = tuple(round(value, 2) for value in box)
                    label = str(figure.get("figure_number") or figure.get("figure_id") or "figure")
                    if key not in seen_figure_boxes:
                        figure_boxes.append([*box, label])
                        seen_figure_boxes.add(key)
                for roi in figure.get("rois", []) or []:
                    roi_count += 1
                    roi_box = _roi_page_bbox(figure, roi, page_height)
                    if roi_box:
                        rendered_roi_count += 1
                        roi_boxes.append([
                            *roi_box,
                            str(roi.get("label") or roi.get("figure_number") or "ROI"),
                        ])

            chosen_boxes, rejected_boxes = [], []
            for owner, candidate in candidates_by_page.get(page_number, []):
                box = _top_left_bbox(
                    candidate.get("caption_bbox"), page_height,
                    "pdf_pts_bottom_left",
                )
                if not box:
                    continue
                label = str(owner.get("figure_number") or owner.get("figure_id") or "caption")
                target = chosen_boxes if candidate.get("chosen") else rejected_boxes
                target.append([*box, label])

            fragments = fragments_by_page.get(page_number, [])
            parsed_text = "\n\n".join(
                f"[{fragment['label']}] {fragment['text']}" for fragment in fragments
            )
            source_page = (
                page_selection[page_number - 1]
                if len(page_selection) >= page_number else page_number
            )
            source_note = (
                f" · source page {source_page}" if source_page != page_number else ""
            )
            uncertain = sum(
                caption_evidence_summary(figure)["caption_status"] != "bound"
                for figure in page_figures
            )
            overlay = {
                "page_width": page_width,
                "page_height": page_height,
                "cells": cell_boxes,
                "figures": figure_boxes,
                "rois": roi_boxes,
                "chosen": chosen_boxes,
                "rejected": rejected_boxes,
            }
            parts.extend([
                f"<section class='page-card' id='page-{page_number}'>",
                f"<h2>Processed page {page_number}{html.escape(source_note)}</h2>",
                "<div class='stats'>"
                f"{len(words)} PDF words · {len(fragments)} Docling text items · "
                f"{len(parsed_text)} parsed characters · {len(page_figures)} figure records · "
                f"{rendered_roi_count}/{roi_count} ROIs mapped to page · "
                f"{len(chosen_boxes)} chosen and {len(rejected_boxes)} rejected caption boxes · "
                f"{uncertain} figures uncertain/unbound"
                "</div>",
                "<div class='page-grid'><div>",
                "<div class='page-image'>",
                f"<img alt='processed page {page_number}' src='{_render_page_data_url(page, max_width, jpeg_quality)}'>",
                "<canvas aria-hidden='true'></canvas>",
                f"<script type='application/json'>{_safe_script_json(overlay)}</script>",
                "</div></div><div class='text-pane'>",
                "<details open><summary>Selectable Docling text</summary>",
                f"<pre>{html.escape(parsed_text) if parsed_text else '<span class=\"empty\">no parsed text on this page</span>'}</pre>",
                "</details>",
                "<h3>Figure and caption evidence</h3>",
                "".join(_figure_caption_html(figure) for figure in page_figures)
                if page_figures else "<div class='empty'>no figure records on this page</div>",
                "</div></div></section>",
            ])
        parts.extend([f"</main><script>{_JS}</script></body></html>"])
    finally:
        pdf.close()

    output = Path(output_path) if output_path else hash_dir / "page_report.html"
    output = output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(parts), encoding="utf-8")
    return output


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("hash_dir", type=Path, help="One output/documents/<HASH> directory")
    parser.add_argument("--output", type=Path, default=None, help="Output HTML path")
    parser.add_argument("--pages", default=None, help="One-based selection, e.g. 1-3,8")
    parser.add_argument(
        "--max-pages", type=int, default=_DEFAULT_MAX_PAGES,
        help="Safety cap after selection; 0 disables (default: 200)",
    )
    parser.add_argument(
        "--max-width", type=int, default=_DEFAULT_MAX_WIDTH,
        help="Maximum embedded page-image width in pixels (default: 1400)",
    )
    parser.add_argument(
        "--jpeg-quality", type=int, default=_DEFAULT_JPEG_QUALITY,
        help="Embedded page JPEG quality, 30-100 (default: 84)",
    )
    args = parser.parse_args(argv)
    try:
        import fitz
        with fitz.open(args.hash_dir / "processed.pdf") as pdf:
            selected = parse_page_spec(args.pages, pdf.page_count)
        output = generate_page_report(
            args.hash_dir, args.output, selected, args.max_pages,
            args.max_width, args.jpeg_quality,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"page report: {exc}", file=sys.stderr)
        return 2
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
