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
    link_chunks_to_figures,
    parse_panels_from_caption,
    reconcile_plate_legend_numbers,
)

logger = logging.getLogger(__name__)

_VISION_PLATE_DISCOVERY_SOURCE = "vision_plate_discovery"


def _annotate_plate_figure_groups(figures, running_text: str) -> int:
    """Attach caption-derived numeric ROI targets to shared plates (#203).

    ``expand_plate_figures`` has already emitted one logical record per
    legend entry, all pointing at one image. This pass makes that relationship
    explicit without pretending the numbers are panel letters. The host alone
    is marked for one Pass 3 invocation; every sibling retains the same target
    list so the evidence is inspectable from any record.

    The running-text check deliberately excludes virtual siblings when it
    computes extracted numbers. That reconstructs the independent
    ``missing_figures`` signal that existed before plate expansion and records
    which proposed regions it corroborates.
    """
    for figure in figures:
        for field in (
            "plate_figures_from_caption",
            "plate_figure_count_from_caption",
            "plate_missing_figure_crosscheck",
            "plate_roi_host",
        ):
            figure.pop(field, None)

    pre_expansion_numbers = {
        str(figure.get("figure_number"))
        for figure in figures
        if figure.get("figure_number")
        and figure.get("shares_image_with") is None
    }
    independently_missing = {
        str(item.get("figure_number"))
        for item in detect_missing_figures(running_text, pre_expansion_numbers)
        if item.get("figure_number")
    }

    by_id = {
        figure.get("figure_id"): figure
        for figure in figures if figure.get("figure_id")
    }
    sibling_groups = {}
    for figure in figures:
        shared = figure.get("shares_image_with")
        if shared is not None:
            sibling_groups.setdefault(shared, []).append(figure)

    annotated = 0
    for shared, siblings in sibling_groups.items():
        host = by_id.get(f"docling_{shared}")
        if host is None:
            # Defensive fallback for an older/custom extractor that used a
            # different figure_id but preserved the shared filename.
            sibling_filename = siblings[0].get("filename")
            host = next((
                figure for figure in figures
                if figure.get("shares_image_with") is None
                and figure.get("filename") == sibling_filename
            ), None)
        if host is None:
            continue

        records = [host, *siblings]
        targets = []
        seen = set()
        for record in records:
            number = str(record.get("figure_number") or "").strip()
            if not number or number in seen:
                continue
            seen.add(number)
            targets.append({
                "label": number,
                "figure_number": number,
                "description": (
                    record.get("caption_text") or record.get("caption") or ""
                ),
                "kind": "figure",
                "figure_id": record.get("figure_id"),
                "missing_figure_crosscheck": number in independently_missing,
            })
        if len(targets) <= 1:
            continue

        corroborated = [
            target["figure_number"] for target in targets
            if target["missing_figure_crosscheck"]
        ]
        for record in records:
            record["plate_figures_from_caption"] = targets
            record["plate_figure_count_from_caption"] = len(targets)
            record["plate_missing_figure_crosscheck"] = corroborated
        host["plate_roi_host"] = True
        annotated += 1
    return annotated


def _roi_request(figure):
    """Return ``(targets, kind)`` for one Pass 3 invocation."""
    plate_targets = figure.get("plate_figures_from_caption") or []
    if figure.get("plate_roi_host") and len(plate_targets) > 1:
        return plate_targets, "figure"
    panels = figure.get("panels_from_caption") or []
    return panels, "panel"


def _plate_caption_text(targets) -> str:
    return "\n".join(
        f"Figure {target.get('figure_number')}: {target.get('description') or ''}"
        for target in targets
    )


def _apply_plate_roi_result(figures, host, targets, result) -> None:
    """Put each shared-plate ROI on its logical figure record."""
    by_id = {
        figure.get("figure_id"): figure
        for figure in figures if figure.get("figure_id")
    }
    rois = result.get("rois") or []
    for target in targets:
        record = by_id.get(target.get("figure_id"))
        if record is None:
            continue
        number = str(target.get("figure_number") or "")
        record["rois"] = [
            dict(roi) for roi in rois
            if str(roi.get("figure_number") or roi.get("label") or "") == number
        ]
        record["pass3_status"] = result.get("pass3_status")
        record["pass3_target_kind"] = "figure"
        record["plate_roi_source_figure_id"] = host.get("figure_id")
        if result.get("pass3_backend"):
            record["pass3_backend"] = result["pass3_backend"]
        if result.get("pass3_error"):
            record["pass3_error"] = result["pass3_error"]
        else:
            record.pop("pass3_error", None)
        if result.get("ocr_token_count") is not None:
            record["ocr_token_count"] = result["ocr_token_count"]
        if result.get("image_size_px"):
            record["image_size_px"] = result["image_size_px"]


def _is_bare_plate_for_discovery(figure) -> bool:
    """Whether Pass 3b may discover uncaptioned engravings on this image."""
    return (
        figure.get("figure_type") == "plate"
        and figure.get("caption_kind") == "bare_label"
        and figure.get("caption_status") == "bound"
        and figure.get("shares_image_with") is None
        and figure.get("image_shared_with") is None
        and not figure.get("plate_figures_from_caption")
    )


def _apply_plate_discovery_result(host, result):
    """Persist discovery evidence and materialize accepted logical figures."""
    candidates = result.get("figure_number_candidates") or []
    rois = result.get("rois") or []
    host["rois"] = rois
    host["pass3_status"] = result.get("pass3_status")
    host["pass3_target_kind"] = "figure_discovery"
    host["pass3_backend"] = result.get("pass3_backend")
    if result.get("image_size_px"):
        host["image_size_px"] = result["image_size_px"]
    if result.get("pass3_error"):
        host["pass3_error"] = result["pass3_error"]
    else:
        host.pop("pass3_error", None)
    host["plate_number_discovery"] = {
        "status": result.get("pass3_status"),
        "backend": result.get("pass3_backend"),
        "minimum_confidence": result.get("minimum_confidence"),
        "candidate_count": len(candidates),
        "accepted_count": sum(1 for candidate in candidates
                              if candidate.get("accepted")),
        "candidates": candidates,
    }

    accepted_by_number = {
        str(candidate.get("figure_number")): candidate
        for candidate in candidates if candidate.get("accepted")
    }
    records = []
    host_id = host.get("figure_id") or "plate"
    for roi in rois:
        number = str(roi.get("figure_number") or "").strip()
        evidence = accepted_by_number.get(number)
        if not number or evidence is None:
            continue
        records.append({
            "figure_id": f"{host_id}_vision_fig{number}",
            "filename": host.get("filename"),
            "file_path": host.get("file_path"),
            "extraction_method": host.get("extraction_method"),
            "figure_type": "figure",
            "figure_number": number,
            "figure_number_source": _VISION_PLATE_DISCOVERY_SOURCE,
            "figure_number_confidence": evidence.get("confidence"),
            # The model found a numbered region, not source caption prose.
            "caption_text": "",
            "caption_page": None,
            "caption_bbox": None,
            "caption_source": None,
            "caption_kind": None,
            "caption_status": "unbound",
            "caption_confidence": None,
            "caption_page_distance": None,
            "caption_candidates": [],
            "page": host.get("page"),
            "bbox": host.get("bbox"),
            "bbox_coord_system": host.get("bbox_coord_system"),
            "image_shared_with": host_id,
            "plate_source_figure_id": host_id,
            "plate_number": host.get("figure_number"),
            "plate_number_discovery_evidence": dict(evidence),
            "rois": [dict(roi)],
            "pass3_status": result.get("pass3_status"),
            "pass3_target_kind": "figure_discovery",
            "pass3_backend": result.get("pass3_backend"),
            "plate_roi_source_figure_id": host_id,
            "image_size_px": result.get("image_size_px"),
        })
    return records


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
    repaired_numbers = reconcile_plate_legend_numbers(figures, missing)
    if repaired_numbers:
        # Recompute from the corrected figure records so a resolved OCR loss
        # does not survive as a contradictory missing_figures entry.
        extracted_nums = {
            fig.get("figure_number") for fig in figures if fig.get("figure_number")
        }
        missing = detect_missing_figures(running_text, extracted_nums)
    plate_groups = _annotate_plate_figure_groups(figures, running_text)
    figures_data["missing_figures"] = missing
    figures_data["total_missing_figures"] = len(missing)

    with figures_file.open("w", encoding="utf-8") as f:
        json.dump(stamp_artifact(figures_data), f, indent=2, ensure_ascii=False)

    # Small per-doc summary in the pipeline log — useful for spotting
    # when the corpus has extractions-vs-text gaps that need attention.
    n_panelled = sum(1 for f in figures if f.get("panel_count_from_caption", 0) > 1)
    logger.info(
        "Pass 2.5: %d/%d figures have multi-panel captions; %d missing-figure(s) "
        "inferred from text; %d plate-legend OCR number(s) reconciled; "
        "%d shared plate(s) admitted to ROI detection",
        n_panelled, len(figures), len(missing), repaired_numbers, plate_groups,
    )


def _pass3a_annotate_rois(figures_file: Path) -> None:
    """Pass 3a — OCR-driven panel/figure ROI detection on each real figure
    whose caption declares multiple panels.

    Opt-in. In-place modifies ``figures.json`` to add ``rois`` +
    ``pass3_status`` per figure, and ``image_size_px`` for figures whose
    images we opened. Skips figures with neither a multi-panel caption nor a
    caption-enumerated shared plate (no point paying OCR cost otherwise).

    OCR reliability on line-art scientific figures varies a lot —
    vision-LLM fallback is Pass 3b, for cases Pass 3a can't resolve (see
    dev_docs/OVERVIEW.md). For now, figures where OCR finds no labels keep
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
        if (fig.get("figure_type") not in ("figure", "plate", "subpanel")
                and not fig.get("plate_roi_host")):
            continue
        # A virtual plate sibling points at the host's image. The host invokes
        # OCR once and distributes the numbered regions to all siblings.
        if (fig.get("shares_image_with") is not None
                and fig.get("plate_figures_from_caption")):
            continue
        targets, target_kind = _roi_request(fig)
        if len(targets) <= 1:
            continue
        img_path = Path(fig.get("file_path", ""))
        if not img_path.exists():
            continue
        result = detect_figure_rois(
            img_path, targets, target_kind=target_kind,
        )
        if target_kind == "figure":
            _apply_plate_roi_result(figures, fig, targets, result)
        else:
            fig["rois"] = result.get("rois") or []
            fig["pass3_status"] = result.get("pass3_status")
            fig["pass3_target_kind"] = "panel"
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
    """Pass 3b — vision-model-driven panel, compound and plate detection.

    Same contract as :func:`_pass3a_annotate_rois`: reads ``figures.json``,
    runs the backend on every real figure with a multi-panel caption or
    caption-enumerated shared plate, plus a tightly gated bare historical
    plate whose printed engraving numbers must be discovered from the image,
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
    # A refresh replaces derived discoveries rather than appending a second
    # generation. Their source record is immutable and remains in the list.
    figures = [
        figure for figure in figures
        if figure.get("figure_number_source") != _VISION_PLATE_DISCOVERY_SOURCE
    ]
    data["figures"] = figures
    for figure in figures:
        figure.pop("plate_number_discovery", None)
        if figure.get("pass3_target_kind") == "figure_discovery":
            for field in (
                "rois", "pass3_status", "pass3_target_kind", "pass3_backend",
                "pass3_error", "image_size_px",
            ):
                figure.pop(field, None)

    # Stamp the version into the Pass 3b log (#253). These logs recorded the
    # GPU and the config path but not the corpus version, so a run could not
    # be attributed to a build afterwards: the siphonophore_20260828 vision
    # pass could not be tied to a version at all, because `runs/` held only
    # records written after it and run.log attested to a different, later
    # run. One line makes the log self-describing.
    from . import PIPELINE_VERSION
    logger.info("Pass 3b starting | corpus %s | backend=%s | %d figure(s)",
                PIPELINE_VERSION, getattr(vision_backend, "name", "?"),
                len(figures))

    n_ok = n_partial = n_none = n_compound = n_discovery = n_skipped = n_failed = 0
    discovered_records = []
    for fig in figures:
        if (fig.get("figure_type") not in ("figure", "plate", "subpanel")
                and not fig.get("plate_roi_host")):
            continue
        if (fig.get("shares_image_with") is not None
                and fig.get("plate_figures_from_caption")):
            continue
        if _is_bare_plate_for_discovery(fig):
            targets, target_kind = [], "figure_discovery"
        else:
            targets, target_kind = _roi_request(fig)
            if len(targets) <= 1:
                continue
        img_path = Path(fig.get("file_path", ""))
        if not img_path.exists():
            continue
        result = detect_figure_rois_via_vision(
            img_path, targets, vision_backend,
            caption_text=(
                _plate_caption_text(targets) if target_kind == "figure"
                else fig.get("caption_text") or fig.get("caption") or ""
            ),
            target_kind=target_kind,
        )
        if target_kind == "figure_discovery":
            discovered_records.extend(_apply_plate_discovery_result(fig, result))
        elif target_kind == "figure":
            _apply_plate_roi_result(figures, fig, targets, result)
        else:
            fig["rois"] = result.get("rois") or []
            fig["pass3_status"] = result.get("pass3_status")
            fig["pass3_target_kind"] = "panel"
            fig["pass3_backend"] = result.get("pass3_backend")
            if result.get("pass3_error"):
                fig["pass3_error"] = result["pass3_error"]
            else:
                # A successful retry must not retain the prior failure reason.
                fig.pop("pass3_error", None)
            if result.get("image_size_px"):
                fig["image_size_px"] = result["image_size_px"]
        s = result.get("pass3_status") or ""
        if s.endswith("_compound"):
            n_compound += 1
        if s == "discovery_materialized":
            n_discovery += 1
        elif s.startswith("completed"):
            n_ok += 1
        elif s.startswith("partial"):
            n_partial += 1
        elif s == "no_labels_found":
            n_none += 1
        elif s == "vision_backend_failed":
            n_failed += 1
        else:
            n_skipped += 1
    figures.extend(discovered_records)
    with figures_file.open("w", encoding="utf-8") as f:
        json.dump(stamp_artifact(data), f, indent=2, ensure_ascii=False)
    logger.info(
        "Pass 3b: %d completed, %d partial, %d no-labels, %d compound, "
        "%d bare plates materialized, %d skipped, %d backend-failed",
        n_ok, n_partial, n_none, n_compound, n_discovery, n_skipped, n_failed,
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
