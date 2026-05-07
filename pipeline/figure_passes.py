"""Figure passes 2.5 / 3a / 3b + cross-ref linker.

* :func:`_pass25_annotate_figures` — cheap always-on caption-panel
  parser + missing-figures scan over running text.
* :func:`_pass3a_annotate_rois` — Tesseract OCR-driven panel-label
  ROI detection on the figure PNG.
* :func:`_pass3b_annotate_rois` — vision-LLM ROI detection (Claude
  or local Qwen) over the same figure PNGs.
* :func:`_crossref_chunks_and_figures` — joins ``Fig. N`` mentions
  in chunk text to figure records by parsed figure_number.

All four read/write ``figures.json`` (and ``chunks.json`` for the
crossref). Pass 3a/3b are gated by CLI flags in run_pdf_processing_pipeline;
this module just exposes the work functions.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

from . import stamp_artifact
from .figures import (
    detect_figure_rois,
    detect_figure_rois_via_vision,
    detect_missing_figures,
    extract_caption_info,
    link_chunks_to_figures,
    parse_figure_number,
    parse_panels_from_caption,
)

logger = logging.getLogger(__name__)


def _pass25_annotate_figures(text_file: Path, figures_file: Path) -> None:
    """Pass 2.5 — caption-based panel parsing + missing-figures scan.

    Cheap, always-on. Reads ``text.json`` + ``figures.json``, and
    in-place-updates ``figures.json`` to add:

    * Per-figure:
        - ``panels_from_caption`` — ``[{label, description, kind}]`` parsed
          from the caption text (English/German/French/Russian figure
          prefixes + A./(A)/A-C patterns).
        - ``panel_count_from_caption`` — length of the list above.
    * Per-document (top-level of figures.json):
        - ``missing_figures`` — figure numbers cited in the running text
          that docling didn't extract. Each entry includes a
          best-effort ``caption_text_candidate`` pulled from the running
          text so the content isn't lost just because the image was.

    Does NOT change filenames or move images — that's Pass 3a's job.
    """
    if not figures_file.exists():
        return
    try:
        figures_data = json.load(figures_file.open(encoding="utf-8"))
    except Exception as e:
        logger.warning("Pass 2.5 skipped: couldn't read %s: %s", figures_file, e)
        return
    figures = figures_data.get("figures", []) or []

    for fig in figures:
        caption = fig.get("caption_text") or fig.get("caption") or ""
        panels = parse_panels_from_caption(caption)
        fig["panels_from_caption"] = panels
        fig["panel_count_from_caption"] = len(panels)

    # Missing-figures scan runs only if text.json is readable.
    running_text = ""
    if text_file.exists():
        try:
            td = json.load(text_file.open(encoding="utf-8"))
            running_text = td.get("text") or ""
        except Exception as e:
            logger.warning("Pass 2.5: could not read %s: %s", text_file, e)

    extracted_nums = {
        fig.get("figure_number") for fig in figures if fig.get("figure_number")
    }
    missing = detect_missing_figures(running_text, extracted_nums)
    figures_data["missing_figures"] = missing
    figures_data["total_missing_figures"] = len(missing)

    with figures_file.open("w", encoding="utf-8") as f:
        json.dump(stamp_artifact(figures_data), f, indent=2, ensure_ascii=False)

    # Small per-doc summary in the pipeline log — useful for spotting
    # when the corpus has extractions-vs-text gaps that need attention.
    n_panelled = sum(1 for f in figures if f.get("panel_count_from_caption", 0) > 1)
    logger.info(
        "Pass 2.5: %d/%d figures have multi-panel captions; %d missing-figure(s) inferred from text",
        n_panelled, len(figures), len(missing),
    )


def _pass3a_annotate_rois(figures_file: Path) -> None:
    """Pass 3a — OCR-driven panel/figure ROI detection on each real figure
    whose caption declares multiple panels.

    Opt-in. In-place modifies ``figures.json`` to add ``rois`` +
    ``pass3_status`` per figure, and ``image_size_px`` for figures whose
    images we opened. Skips figures with no multi-panel caption (no
    point paying OCR cost on single-panel figures).

    OCR reliability on line-art scientific figures varies a lot —
    dev_docs/PLAN.md §9 calls out vision-LLM fallback as Pass 3b for cases Pass 3a
    can't resolve. For now, figures where OCR finds no labels keep
    ``pass3_status = "no_labels_found"`` and no ROIs; downstream tools
    should fall back to whole-image retrieval + caption description.
    """
    if not figures_file.exists():
        return
    try:
        data = json.load(figures_file.open(encoding="utf-8"))
    except Exception as e:
        logger.warning("Pass 3a skipped: couldn't read %s: %s", figures_file, e)
        return
    figures = data.get("figures", []) or []
    n_ok = n_partial = n_none = n_skipped = 0
    for fig in figures:
        if fig.get("figure_type") not in ("figure", "plate", "subpanel"):
            continue
        panels = fig.get("panels_from_caption") or []
        if len(panels) <= 1:
            continue
        img_path = Path(fig.get("file_path", ""))
        if not img_path.exists():
            continue
        result = detect_figure_rois(img_path, panels)
        fig["rois"] = result.get("rois") or []
        fig["pass3_status"] = result.get("pass3_status")
        fig["ocr_token_count"] = result.get("ocr_token_count", 0)
        if result.get("image_size_px"):
            fig["image_size_px"] = result["image_size_px"]
        s = result.get("pass3_status")
        if s == "completed":
            n_ok += 1
        elif s == "partial_ocr":
            n_partial += 1
        elif s == "no_labels_found":
            n_none += 1
        else:
            n_skipped += 1
    with figures_file.open("w", encoding="utf-8") as f:
        json.dump(stamp_artifact(data), f, indent=2, ensure_ascii=False)
    logger.info(
        "Pass 3a: %d completed, %d partial, %d no-labels, %d skipped",
        n_ok, n_partial, n_none, n_skipped,
    )


def _pass3b_annotate_rois(figures_file: Path, vision_backend) -> None:
    """Pass 3b — vision-model-driven panel + compound-figure detection.

    Same contract as :func:`_pass3a_annotate_rois`: reads ``figures.json``,
    runs the backend on every real figure with a multi-panel caption,
    writes ROIs + ``pass3_status`` back in place.

    Runs INSTEAD of Pass 3a when both flags are set (Pass 3b supersedes
    Pass 3a because vision-LLM quality is consistently higher on
    scientific-figure panel detection). If Pass 3a had already populated
    ``rois`` for a figure, the vision result overwrites them.
    """
    if not figures_file.exists():
        return
    try:
        data = json.load(figures_file.open(encoding="utf-8"))
    except Exception as e:
        logger.warning("Pass 3b skipped: couldn't read %s: %s", figures_file, e)
        return
    figures = data.get("figures", []) or []

    n_ok = n_partial = n_none = n_compound = n_skipped = n_failed = 0
    for fig in figures:
        if fig.get("figure_type") not in ("figure", "plate", "subpanel"):
            continue
        panels = fig.get("panels_from_caption") or []
        if len(panels) <= 1:
            continue
        img_path = Path(fig.get("file_path", ""))
        if not img_path.exists():
            continue
        result = detect_figure_rois_via_vision(
            img_path, panels, vision_backend,
            caption_text=fig.get("caption_text") or fig.get("caption") or "",
        )
        fig["rois"] = result.get("rois") or []
        fig["pass3_status"] = result.get("pass3_status")
        fig["pass3_backend"] = result.get("pass3_backend")
        if result.get("image_size_px"):
            fig["image_size_px"] = result["image_size_px"]
        s = result.get("pass3_status") or ""
        if s.endswith("_compound"):
            n_compound += 1
        if s.startswith("completed"):
            n_ok += 1
        elif s.startswith("partial"):
            n_partial += 1
        elif s == "no_labels_found":
            n_none += 1
        elif s == "vision_backend_failed":
            n_failed += 1
        else:
            n_skipped += 1
    with figures_file.open("w", encoding="utf-8") as f:
        json.dump(stamp_artifact(data), f, indent=2, ensure_ascii=False)
    logger.info(
        "Pass 3b: %d completed, %d partial, %d no-labels, %d compound, "
        "%d skipped, %d backend-failed",
        n_ok, n_partial, n_none, n_compound, n_skipped, n_failed,
    )


def _crossref_chunks_and_figures(figures_file: Path, chunks_file: Path) -> None:
    """Run the bidirectional chunk ↔ figure cross-reference pass.

    Reads figures.json and chunks.json, calls
    :func:`figures.link_chunks_to_figures`, and writes both files back in
    place. No-op when either file is missing or unreadable.
    """
    if not figures_file.exists() or not chunks_file.exists():
        logger.debug("Skipping cross-ref: missing %s or %s", figures_file, chunks_file)
        return
    try:
        figures_data = json.load(figures_file.open(encoding="utf-8"))
        chunks_data = json.load(chunks_file.open(encoding="utf-8"))
    except Exception as e:
        logger.warning("Cross-ref skipped: could not load files: %s", e)
        return

    chunks = chunks_data.get("chunks", []) or []
    figures = figures_data.get("figures", []) or []
    if not chunks or not figures:
        logger.info(
            "Cross-ref: skipped (no chunks=%d or no figures=%d)",
            len(chunks), len(figures),
        )
        return

    link_chunks_to_figures(chunks, figures)

    # Write back — data was modified in place but be explicit about
    # re-serialization to keep JSON formatting consistent.
    chunks_data["chunks"] = chunks
    figures_data["figures"] = figures
    with figures_file.open("w", encoding="utf-8") as f:
        json.dump(stamp_artifact(figures_data), f, indent=2, ensure_ascii=False)
    with chunks_file.open("w", encoding="utf-8") as f:
        json.dump(stamp_artifact(chunks_data), f, indent=2, ensure_ascii=False)

    # Log a small summary so pipeline.log captures the link density.
    linked_figs = sum(1 for fig in figures if fig.get("referenced_in_chunks"))
    linked_chunks = sum(1 for ch in chunks if ch.get("figure_refs"))
    logger.info(
        "Cross-ref: %d/%d figures referenced, %d/%d chunks cite a figure",
        linked_figs, len(figures), linked_chunks, len(chunks),
    )
