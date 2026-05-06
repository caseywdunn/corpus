#!/usr/bin/env python3
"""
Generalized corpus processing script that takes input_dir and output_dir arguments.
Recursively finds PDFs in input_dir, processes them with hash-based organization.
"""

import argparse
import hashlib
import json
import logging
import re
import shutil
import subprocess
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple
import multiprocessing
import tempfile
import os

from grobid_client import (
    GrobidClient,
    GrobidUnavailableError,
    parse_tei_header,
    parse_tei_intext_citations,
    parse_tei_references,
)
from bib import BibIndex, bib_entry_to_metadata
from figures import (
    extract_caption_info,
    parse_figure_number,
    parse_panels_from_caption,
    detect_missing_figures,
    detect_figure_rois,
    detect_figure_rois_via_vision,
    resolve_compound_figures,
    link_chunks_to_figures,
    generate_figures_report,
)
from taxa import (
    TaxonomyDB,
    extract_taxon_mentions,
    extract_lexicon_mentions,
    load_lexicon,
)


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pipeline foundations (#42 / #44 partial split)
# ---------------------------------------------------------------------------
#
# Cross-cutting infrastructure (config, logging, stage runner, failures,
# quality gates, resume helpers, hashing/orphan/summary I/O) lives in
# the ``pipeline/`` package. Imported here so existing references in
# the rest of this file keep working unchanged. Re-exported at the
# module level so callers that do ``from process_corpus import …``
# (tests, downstream tools) keep working during the deprecation
# window. Per-stage extraction is open as #45.
from pipeline.config import (
    _DEFAULT_CONFIG,
    _deep_merge,
    classify_section,
    load_config,
)
from pipeline import config as _pipeline_config
from pipeline.config import CONFIG  # mutable module-level singleton
from pipeline.log import per_pdf_file_log, setup_root_logging
from pipeline.stages import (
    _all_stage_artifacts_complete,
    _classify_exception,
    _file_sha256,
    _HugeDocumentError,
    _pdf_page_count,
    _REASON_CODES,
    _run_quality_gates,
    _safe_load_json,
    _should_run_stage,
    _stage,
    _stage_artifacts_present,
    _utcnow_iso,
)
from pipeline.io import (
    HASH_PREFIX_LEN,
    _verify_or_raise_collision,
    audit_orphans,
    calculate_pdf_hash,
    create_output_structure,
    create_summary_json,
    find_all_pdfs,
    get_relative_paths,
    short_hash,
)
from pipeline.annotate import _extract_taxa_and_anatomy
from pipeline.chunking import chunk_text, ingest_to_vector_db
from pipeline.metadata import _write_placeholder_metadata, _year_from_filename, extract_metadata
from pipeline.figure_passes import (
    _crossref_chunks_and_figures,
    _pass25_annotate_figures,
    _pass3a_annotate_rois,
    _pass3b_annotate_rois,
)
from pipeline.scan import (
    _annotate_pack_availability,
    _available_tesseract_langs,
    _compose_ocr_langs,
    _detect_language,
    _gibberish_score,
    _resolve_tesseract_packs,
    _text_layer_scripts,
    _visual_page_script,
    create_cell_visualizations,
    detect_scan_type,
    prepare_pdf,
)




def extract_docling_content(
    pdf_path: Path,
    text_output: Path,
    figures_output: Path,
    figures_dir: Path,
    visualizations_dir: Path,
    docling_doc_output: Optional[Path] = None,
    scan_file_type: Optional[str] = None,
):
    """Extract text and figures using docling, with PyMuPDF fallback for figures.

    Guarantees that only successfully saved figures are listed in figures.json.
    Also persists the parsed :class:`DoclingDocument` to
    ``docling_doc_output`` (defaults to ``<text_output.parent>/docling_doc.json``)
    so the chunking step can re-load the structured document and run
    HybridChunker without re-parsing the PDF.

    ``scan_file_type`` is the ``file_type`` from ``scan_detection.json``
    (``"scanned"`` / ``"born_digital"`` / ``"broken_text_layer"`` / None).
    Used to suppress the PyMuPDF figure fallback on rasterized scans, where
    each ``page.get_images()`` would just emit one full-page bitmap per page
    rather than a real figure.
    """
    document = None
    text_content = None

    # Ensure figure directory exists
    figures_dir.mkdir(parents=True, exist_ok=True)

    try:
        from docling.document_converter import DocumentConverter, PdfFormatOption
        from docling.datamodel.base_models import InputFormat
        from docling.datamodel.pipeline_options import PdfPipelineOptions
        # Configure converter to generate picture images explicitly
        pipeline_options = PdfPipelineOptions(
            do_ocr=False,
            do_table_structure=True,
            generate_picture_images=True,
            generate_page_images=False,
            do_picture_classification=True,
        )
        pdf_format_option = PdfFormatOption(pipeline_options=pipeline_options)
        converter = DocumentConverter(format_options={InputFormat.PDF: pdf_format_option})
        result = converter.convert(str(pdf_path))
        document = result.document

        # Extract text from docling if available
        text_content = {
            "title": document.name,
            "text": document.export_to_markdown(),
            "pages": len(document.pages) if hasattr(document, "pages") else None,
        }
    except ImportError as e:
        logger.warning("Docling not available (%s), will fall back for figures", e)
    except Exception as e:
        logger.warning("Docling extraction failed (%s), will fall back for figures", e)

    # If we couldn't get text via docling, write a placeholder text file
    if text_content is None:
        text_content = {
            "title": pdf_path.stem,
            "text": f"# {pdf_path.stem}\n\n[Text extraction unavailable — see logs]",
            "pages": None,
        }

    # Save text content (from docling or placeholder)
    with open(text_output, "w", encoding="utf-8") as f:
        json.dump(text_content, f, indent=2, ensure_ascii=False)

    # Persist the structured DoclingDocument so chunking can reload it
    # without re-parsing the PDF. This separation keeps each pipeline step
    # a pure function of the on-disk artifacts of the previous step.
    if document is not None:
        if docling_doc_output is None:
            docling_doc_output = text_output.parent / "docling_doc.json"
        try:
            document.save_as_json(docling_doc_output)
            logger.info("Saved DoclingDocument to %s", docling_doc_output)
        except Exception as e:
            logger.warning("Could not save DoclingDocument: %s", e)

    # Accumulate figure metadata only when file is saved
    figures_data = []

    # Helper to append a figure entry only if file exists
    def append_figure(figure_id: str, figure_path: Path, caption: str, meta: dict):
        if figure_path.exists() and figure_path.stat().st_size > 0:
            entry = {
                "figure_id": figure_id,
                "filename": figure_path.name,
                "file_path": str(figure_path),
                "caption_text": caption or "",  # canonical; see PLAN.md §3
            }
            entry.update(meta or {})
            figures_data.append(entry)
        else:
            logger.warning("Skipping missing/empty figure file: %s", figure_path)

    def _docling_prov_to_bbox_page(picture) -> dict:
        """Extract (page, bbox) from a docling picture's provenance, if any.

        Docling stores bboxes in bottom-left origin PDF points; we record the
        coordinate system alongside so downstream consumers don't have to
        guess.
        """
        meta: dict = {}
        prov = getattr(picture, "prov", None)
        if not prov:
            return meta
        first = prov[0]
        page_no = getattr(first, "page_no", None)
        bbox = getattr(first, "bbox", None)
        if page_no is not None:
            meta["page"] = page_no
        if bbox is not None:
            try:
                meta["bbox"] = [float(bbox.l), float(bbox.b), float(bbox.r), float(bbox.t)]
                meta["bbox_coord_system"] = "pdf_pts_bottom_left"
            except Exception as e:  # bbox shape varies between docling versions
                logger.warning("Could not serialize docling bbox: %s", e)
        return meta

    # Two-pass docling figure extraction (Phase D.2 — see PLAN.md).
    #
    # Pass 1: gather every docling Picture's image + bbox + caption in
    # memory. We don't write images to disk yet — the filename policy
    # depends on the figure type (real figure? plate? graphical furniture?
    # subpanel?) which we only know after dedup + classification.
    #
    # Pass 2: classify_figure() + dedupe_figures() prune redundant panel
    # crops and flag non-figures (tiny journal-furniture boxes, captionless
    # undersized elements).
    #
    # Pass 3: save surviving items under semantic filenames (fig_3.png,
    # plate_2.png, _graphic_5.png) so a directory listing is human-readable
    # and the MCP tools can default to real figures by filtering on
    # figure_type.
    if document is not None and hasattr(document, "pictures"):
        from figures import (
            classify_figure,
            dedupe_figures,
            compose_figure_filename,
            FIGURE_TYPE_FIGURE,
            FIGURE_TYPE_PLATE,
            FIGURE_TYPE_SUBPANEL,
        )

        raw_items: List[Dict] = []
        for idx, picture in enumerate(document.pictures or []):
            caption_info = extract_caption_info(picture, document)
            caption_text = caption_info.get("caption_text", "")
            figure_number = parse_figure_number(caption_text)
            bbox_meta = _docling_prov_to_bbox_page(picture)
            # Resolve image up-front (still in memory); we save after
            # classification once we know the intended filename.
            try:
                image = None
                if hasattr(picture, "get_image"):
                    image = picture.get_image(document) if callable(picture.get_image) else picture.get_image
            except Exception as e:
                logger.warning("Could not resolve docling figure %d image: %s", idx + 1, e)
                image = None

            raw_items.append({
                "docling_idx": idx + 1,
                "image": image,
                "caption_text": caption_text,
                "caption_page": caption_info.get("caption_page"),
                "caption_bbox": caption_info.get("caption_bbox"),
                "caption_source": caption_info.get("caption_source"),
                "figure_number": figure_number,
                "bbox": bbox_meta.get("bbox"),
                "page": bbox_meta.get("page"),
                "bbox_coord_system": bbox_meta.get("bbox_coord_system"),
            })

        # Pass 2: classify, then dedupe. Order matters — subpanels are
        # identified during dedupe and we don't want the initial
        # classify_figure() to call them "figure" and then the dedup pass
        # relabel them.
        for it in raw_items:
            if "figure_type" not in it:
                it["figure_type"] = classify_figure(it)
        items = dedupe_figures(raw_items)

        # Pass 3: save surviving images, record metadata.
        saved_count = 0
        for it in items:
            image = it.get("image")
            if image is None or not hasattr(image, "save"):
                logger.warning(
                    "No savable image for docling_%s (figure_type=%s); skipping",
                    it.get("docling_idx"), it.get("figure_type"),
                )
                continue
            filename = compose_figure_filename(it)
            # Extremely unlikely collision guard — if a filename is already
            # taken by a prior entry this pass, fall back to a docling-idx-
            # qualified name so nothing is overwritten silently.
            out_path = figures_dir / filename
            if out_path.exists():
                stem, ext = out_path.stem, out_path.suffix
                out_path = figures_dir / f"{stem}_docling_{it['docling_idx']}{ext}"
            try:
                image.save(str(out_path))
            except Exception as e:
                logger.warning("Could not save figure %s: %s", filename, e)
                continue
            saved_count += 1

            meta: Dict = {
                "extraction_method": "docling",
                "figure_type": it.get("figure_type"),
                "figure_number": it.get("figure_number"),
                "page": it.get("page"),
                "bbox": it.get("bbox"),
                "bbox_coord_system": it.get("bbox_coord_system"),
                "caption_page": it.get("caption_page"),
                "caption_bbox": it.get("caption_bbox"),
                "caption_source": it.get("caption_source"),
            }
            if it.get("figure_type") == FIGURE_TYPE_SUBPANEL:
                meta["primary_figure_docling_idx"] = it.get("primary_figure_docling_idx")
                meta["panel_letter"] = it.get("panel_letter")
            append_figure(
                figure_id=f"docling_{it['docling_idx']}",
                figure_path=out_path,
                caption=it.get("caption_text") or "",
                meta=meta,
            )

        # One-line summary of the classification outcome (useful in logs
        # + pipeline.log when reviewing extraction quality).
        type_counts: Dict[str, int] = {}
        for it in items:
            type_counts[it.get("figure_type", "?")] = type_counts.get(it.get("figure_type", "?"), 0) + 1
        raw_n = len(raw_items)
        dropped_dup = raw_n - len(items)
        logger.info(
            "Figure extraction: %d raw → %d kept (%d dup merged). types=%s",
            raw_n, len(items), dropped_dup, type_counts,
        )
        if saved_count == 0:
            logger.info("No figures saved via docling; will try PyMuPDF fallback")

    # Fallback: use PyMuPDF to extract embedded images from the PDF.
    # Skipped on scanned PDFs — every page is a rasterized bitmap, so
    # page.get_images() would just emit one full-page image per page rather
    # than rescuing real figures.
    if len(figures_data) == 0 and scan_file_type == "scanned":
        logger.info(
            "PyMuPDF figure fallback skipped: scan_detection.file_type=scanned "
            "(page bitmaps are not real figures)"
        )
    elif len(figures_data) == 0:
        try:
            from figures import classify_figure
            import fitz  # PyMuPDF
            doc = fitz.open(str(pdf_path))
            for page_num in range(len(doc)):
                page = doc.load_page(page_num)
                for img_index, img in enumerate(page.get_images()):
                    try:
                        xref = img[0]
                        # Attempt to recover the bbox of the image on the page.
                        # get_image_bbox accepts the image tuple directly; if
                        # that fails (PyMuPDF version differences), fall back
                        # to get_image_rects by xref.
                        bbox_list = None
                        try:
                            rect = page.get_image_bbox(img)
                            if rect and not rect.is_empty:
                                bbox_list = [float(rect.x0), float(rect.y0), float(rect.x1), float(rect.y1)]
                        except Exception:
                            try:
                                rects = page.get_image_rects(xref)
                                if rects:
                                    r = rects[0]
                                    bbox_list = [float(r.x0), float(r.y0), float(r.x1), float(r.y1)]
                            except Exception:
                                bbox_list = None

                        pix = fitz.Pixmap(doc, xref)
                        if pix.n - pix.alpha < 4:  # RGB or GRAY
                            figure_path = figures_dir / f"page{page_num + 1}_img{img_index + 1}.png"
                            pix.save(str(figure_path))
                            meta = {
                                "extraction_method": "pymupdf",
                                "page": page_num + 1,
                                "width": pix.width,
                                "height": pix.height,
                            }
                            if bbox_list is not None:
                                # PyMuPDF uses PDF-like coords but top-left origin
                                # (y grows downward). Record that so consumers
                                # don't confuse it with docling's bottom-left.
                                meta["bbox"] = bbox_list
                                meta["bbox_coord_system"] = "pdf_pts_top_left"
                            # Run the same size-threshold classifier docling
                            # figures get, so MCP _REAL_FIGURE_TYPES filtering
                            # works on PyMuPDF rescues too. No caption is
                            # available here, so size is the only signal —
                            # the classifier will return "graphical_element"
                            # for tiny page-corner thumbnails or "unclassified"
                            # for full-sized images without captions.
                            meta["figure_type"] = classify_figure({
                                "bbox": bbox_list,
                                "caption_text": "",
                                "figure_number": None,
                            })
                            append_figure(
                                figure_id=f"pymupdf_p{page_num + 1}_i{img_index + 1}",
                                figure_path=figure_path,
                                caption="",
                                meta=meta,
                            )
                        pix = None
                    except Exception as e:
                        logger.warning(
                            "Failed to save PyMuPDF image p%d i%d: %s",
                            page_num + 1, img_index + 1, e,
                        )
            doc.close()
        except ImportError:
            logger.warning("PyMuPDF not available; cannot run fallback image extraction")
        except Exception as e:
            logger.warning("PyMuPDF fallback failed: %s", e)

    # Write figures.json
    figures_info = {
        "figures": figures_data,
        "figures_directory": str(figures_dir),
        "total_figures": len(figures_data),
    }
    with open(figures_output, "w", encoding="utf-8") as f:
        json.dump(figures_info, f, indent=2)

    # Extract figure bboxes from the docling document before releasing it,
    # so the (large) document object can be garbage-collected before the
    # per-page image rendering loop.
    figure_bboxes_by_page: dict = {}
    if document is not None:
        if hasattr(document, "pictures") and document.pictures:
            for picture in document.pictures:
                if hasattr(picture, "prov") and picture.prov:
                    for prov_item in picture.prov:
                        if hasattr(prov_item, "page_no") and hasattr(prov_item, "bbox"):
                            page_no = prov_item.page_no
                            figure_bboxes_by_page.setdefault(page_no, []).append(
                                prov_item.bbox
                            )
    del document
    import gc; gc.collect()

    pdf_name = pdf_path.stem
    try:
        create_cell_visualizations(
            pdf_path, visualizations_dir, pdf_name,
            figure_bboxes_by_page=figure_bboxes_by_page,
        )
    except Exception as e:
        logger.warning("Creating visualizations failed: %s", e)









def main():
    parser = argparse.ArgumentParser(
        description="Process a corpus of PDFs with hash-based organization"
    )
    parser.add_argument("input_dir", type=Path, help="Input directory containing PDFs")
    parser.add_argument("output_dir", type=Path, help="Output directory for processed files")
    parser.add_argument("--resume", action="store_true", help="Skip documents whose summary.json already exists")
    parser.add_argument("--config", type=Path, default=None, help="Path to config.yaml (defaults to ./config.yaml)")
    parser.add_argument(
        "--grobid-url",
        default=os.environ.get("GROBID_URL", "http://localhost:8070"),
        help="Grobid service URL (default: $GROBID_URL or http://localhost:8070); "
        "pass --grobid-url='' or --no-grobid to skip metadata extraction",
    )
    parser.add_argument(
        "--no-grobid",
        action="store_true",
        help="Skip Grobid even if reachable (useful for dev iteration on non-metadata stages)",
    )
    parser.add_argument(
        "--bib",
        type=Path,
        default=None,
        help="Optional BibTeX file. Records whose 'file = {Foo.pdf}' field "
             "matches an input PDF supply the metadata header (title, authors, "
             "year, journal, DOI) instead of Grobid. References are still "
             "parsed from each PDF by Grobid.",
    )
    parser.add_argument(
        "--taxonomy-db",
        type=Path,
        default=None,
        help="Path to Darwin Core taxonomy SQLite "
             "(default: <output_dir>/taxonomy.sqlite). "
             "Build with: python ingest_taxonomy.py --source <dwc|dwca|worms> ...",
    )
    parser.add_argument(
        "--anatomy-lexicon",
        type=Path,
        default=None,
        help="Path to a domain-specific anatomy lexicon YAML. Optional and "
             "user-supplied — see demo/anatomy_lexicon.yaml for the format. "
             "Without this flag, anatomy extraction is skipped.",
    )
    parser.add_argument(
        "--lexicon",
        action="append",
        default=[],
        metavar="CATEGORY:PATH",
        help="Additional category-tagged lexicon YAML, e.g. "
             "--lexicon biogeography:/path/to/biogeo.yaml. Repeatable. "
             "Same YAML format as --anatomy-lexicon (see "
             "demo/anatomy_lexicon.yaml). Output lands in "
             "<hash>/<CATEGORY>.json; per-stage resume + corpus_status "
             "treat each category independently. (#24)",
    )
    parser.add_argument(
        "--no-taxa",
        action="store_true",
        help="Skip taxon + anatomy extraction",
    )
    parser.add_argument(
        "--content-aware-figures",
        action="store_true",
        help="Run Pass 3a — OCR-driven panel/figure ROI detection on multi-panel "
             "figures (adds ~1-2 s per figure with caption panels; opt-in because "
             "OCR reliability on line-art figures is mixed; see PLAN.md §9)",
    )
    parser.add_argument(
        "--vision-backend",
        choices=["claude", "local"],
        default=None,
        help="Run Pass 3b — vision-LLM-driven panel + compound-figure detection. "
             "'claude' uses the Anthropic API (needs ANTHROPIC_API_KEY); "
             "'local' uses an open-weights VLM on CUDA/MPS (Bouchet production). "
             "Supersedes --content-aware-figures when both are set.",
    )
    parser.add_argument(
        "--vision-model",
        default=None,
        help="Override the per-backend default vision model (e.g. "
             "claude-sonnet-4-6-20251001 for higher quality than Haiku).",
    )
    parser.add_argument(
        "--refresh-vision",
        action="store_true",
        help="With --resume and --vision-backend: instead of skipping hashes whose "
             "summary.json already exists, re-run ONLY Pass 3b on each hash's "
             "existing figures.json. No OCR/Docling/Grobid/chunking is re-done.",
    )
    parser.add_argument(
        "--audit-orphans",
        action="store_true",
        help="Read-only audit. List documents/<HASH>/ directories whose source "
             "PDF is no longer in input_dir, and LanceDB rows whose hash has no "
             "documents/ directory. Re-hashes input PDFs so it's path-independent. "
             "Does not delete anything.",
    )
    parser.add_argument(
        "--strict-network",
        action="store_true",
        help="Fail fast on the first transient external-service failure "
             "(Grobid 5xx, connect error, timeout) instead of retrying. "
             "Use for release-build runs where silent partial data is "
             "worse than aborting.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Discover input PDFs, apply --resume + batch filters, and report "
             "the work plan without processing anything. No external services "
             "are contacted; no files are written.",
    )
    parser.add_argument(
        "--batch-index",
        type=int,
        default=None,
        help="0-based index of the batch to process (use with --batch-size). "
             "Typically set to $SLURM_ARRAY_TASK_ID.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=256,
        help="Number of unique PDFs (hashes) per batch (default: 256).",
    )

    args = parser.parse_args()

    if args.refresh_vision and not args.vision_backend:
        parser.error("--refresh-vision requires --vision-backend to be set")

    setup_root_logging()

    if args.strict_network:
        from external import set_strict_network
        set_strict_network(True)

    # Mutate the pipeline.config singleton in place so all readers
    # (per-stage modules, quality_gates, etc.) see the loaded values.
    # Reassignment via ``global CONFIG = ...`` would only rebind the
    # local re-export — the canonical dict in pipeline.config would
    # stay at defaults.
    loaded = load_config(args.config)
    _pipeline_config.CONFIG.clear()
    _pipeline_config.CONFIG.update(loaded)

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.exists():
        logger.error("Input directory %s does not exist", input_dir)
        sys.exit(1)

    if args.audit_orphans:
        audit_orphans(input_dir, output_dir)
        return

    logger.info("Processing PDFs from: %s", input_dir)
    logger.info("Output directory: %s", output_dir)

    # One-time Grobid health check. If the service is unreachable at
    # startup we log and carry on with placeholder metadata for every
    # document rather than retrying (and logging) per PDF.
    grobid_client: Optional[GrobidClient] = None
    if args.no_grobid or not args.grobid_url:
        logger.info("Grobid skipped (--no-grobid or empty --grobid-url)")
    else:
        probe = GrobidClient(base_url=args.grobid_url)
        if probe.is_alive():
            logger.info("Grobid reachable at %s", args.grobid_url)
            grobid_client = probe
        else:
            logger.warning(
                "Grobid not reachable at %s — metadata will be placeholder. "
                "Start it with: docker compose up -d grobid",
                args.grobid_url,
            )

    # Optional BibTeX-driven metadata override. Loaded once and shared
    # across workers — entries are looked up by PDF basename.
    bib_index: Optional[BibIndex] = None
    if args.bib is not None:
        if not args.bib.exists():
            logger.error("Bib file %s does not exist", args.bib)
            sys.exit(1)
        try:
            bib_index = BibIndex.from_path(args.bib)
        except Exception as e:
            logger.error("Could not parse %s: %s", args.bib, e)
            sys.exit(1)

    # Open taxonomy snapshot and (if supplied) load anatomy lexicon. Both
    # are optional; missing inputs are logged and their output artifacts
    # are skipped.  The lexicon is opt-in via --anatomy-lexicon — there
    # is no default lookup because it's a domain-specific user input.
    taxonomy_db: Optional[TaxonomyDB] = None
    anatomy_lexicon: Optional[Dict[str, Dict]] = None
    taxonomy_fingerprint: Optional[Dict[str, Any]] = None
    anatomy_fingerprint: Optional[Dict[str, Any]] = None
    if not args.no_taxa:
        taxonomy_path = args.taxonomy_db or (args.output_dir / "taxonomy.sqlite")
        if taxonomy_path.exists():
            try:
                taxonomy_db = TaxonomyDB(taxonomy_path)
                # Stamp #29: hash once at startup so per-paper writes are
                # cheap. SHA-256 of taxonomy.sqlite is a stable identifier
                # that survives copy/move and changes any time the DB is
                # rebuilt.
                taxonomy_fingerprint = {
                    "path": str(taxonomy_path),
                    "sha256": _file_sha256(taxonomy_path),
                    "size": taxonomy_path.stat().st_size,
                }
                logger.info(
                    "Taxonomy snapshot loaded from %s (%d names, sha256=%s…)",
                    taxonomy_path, len(taxonomy_db.name_set()),
                    taxonomy_fingerprint["sha256"][:12],
                )
            except Exception as e:
                logger.warning(
                    "Could not open taxonomy snapshot %s: %s", taxonomy_path, e,
                )
        else:
            logger.warning(
                "Taxonomy snapshot %s not found — taxon extraction skipped. "
                "Build it with: python ingest_taxonomy.py %s --source <dwc|dwca|worms> ...",
                taxonomy_path, args.output_dir,
            )

        if args.anatomy_lexicon is not None:
            if args.anatomy_lexicon.exists():
                try:
                    anatomy_lexicon = load_lexicon(args.anatomy_lexicon)
                    anatomy_fingerprint = {
                        "path": str(args.anatomy_lexicon),
                        "sha256": _file_sha256(args.anatomy_lexicon),
                        "size": args.anatomy_lexicon.stat().st_size,
                    }
                    logger.info(
                        "Anatomy lexicon loaded from %s (%d terms, sha256=%s…)",
                        args.anatomy_lexicon, len(anatomy_lexicon),
                        anatomy_fingerprint["sha256"][:12],
                    )
                except Exception as e:
                    logger.warning(
                        "Could not load anatomy lexicon %s: %s",
                        args.anatomy_lexicon, e,
                    )
            else:
                logger.warning(
                    "Anatomy lexicon %s not found — anatomy extraction skipped",
                    args.anatomy_lexicon,
                )

    # Multi-category lexicons (#24). --lexicon CATEGORY:PATH, repeatable.
    extra_lexicons: Dict[str, Dict[str, Dict]] = {}
    extra_fingerprints: Dict[str, Dict[str, Any]] = {}
    for spec in (args.lexicon or []):
        if ":" not in spec:
            logger.error("--lexicon expects CATEGORY:PATH, got %r", spec)
            sys.exit(2)
        category, _, raw_path = spec.partition(":")
        category = category.strip().lower()
        path = Path(raw_path).expanduser()
        if not category or not category.replace("_", "").isalnum():
            logger.error(
                "--lexicon CATEGORY must be alphanumeric/underscore (got %r)",
                category,
            )
            sys.exit(2)
        if not path.exists():
            logger.warning(
                "--lexicon %s:%s not found — %s extraction skipped",
                category, path, category,
            )
            continue
        try:
            lex = load_lexicon(path)
            fp = {
                "path": str(path),
                "sha256": _file_sha256(path),
                "size": path.stat().st_size,
            }
            extra_lexicons[category] = lex
            extra_fingerprints[category] = fp
            logger.info(
                "Lexicon[%s] loaded from %s (%d terms, sha256=%s…)",
                category, path, len(lex), fp["sha256"][:12],
            )
        except Exception as e:
            logger.warning(
                "Could not load --lexicon %s:%s: %s", category, path, e,
            )

    # Vision backend for Pass 3b. Constructed once and reused so the
    # backend can keep long-lived state (API client, loaded model, etc.).
    vision_backend = None
    if args.vision_backend:
        try:
            from vision import get_vision_backend
            kwargs = {}
            if args.vision_model:
                kwargs["model"] = args.vision_model
            vision_backend = get_vision_backend(args.vision_backend, **kwargs)
            logger.info("Vision backend loaded: %s", vision_backend.name)
        except Exception as e:
            logger.error(
                "Could not load vision backend %r: %s — Pass 3b will be skipped",
                args.vision_backend, e,
            )

    # Create output directory structure
    documents_dir, vector_db_dir = create_output_structure(output_dir)

    # Find all PDFs and group by hash
    logger.info("Discovering PDFs...")
    pdf_map = find_all_pdfs(input_dir)

    logger.info("Found %d PDF file(s)", sum(len(paths) for paths in pdf_map.values()))
    logger.info("Unique PDFs (by hash): %d", len(pdf_map))

    # ── Pre-filter completed documents before batch slicing ─────────
    # When --resume is active and we're batching, remove fully-completed
    # hashes *before* dividing into batches. A hash is fully complete iff
    # summary.json exists AND every per-stage artifact is on disk (#28);
    # otherwise it has missing stages and the inner per-stage guards will
    # only re-run those.
    #
    # Exception: --refresh-vision (#27) explicitly targets the
    # already-completed population (it re-runs Pass 3b on existing
    # figures.json files). Skipping the pre-filter for that case
    # preserves the work-set so array tasks aren't silently emptied.
    if args.resume and args.batch_index is not None and not args.refresh_vision:
        expect_taxa = taxonomy_db is not None
        expect_anatomy = bool(anatomy_lexicon)
        before = len(pdf_map)
        kept = {}
        for h, paths in pdf_map.items():
            hd = documents_dir / short_hash(h)
            if (hd / "summary.json").exists() and _all_stage_artifacts_complete(
                hd, expect_taxa=expect_taxa, expect_anatomy=expect_anatomy,
            ):
                continue
            kept[h] = paths
        pdf_map = kept
        skipped = before - len(pdf_map)
        if skipped:
            logger.info(
                "Resume: filtered out %d fully-completed documents "
                "(%d remaining to process)", skipped, len(pdf_map),
            )

    # ── Batch slicing (for SLURM job arrays) ────────────────────────
    if args.batch_index is not None:
        all_hashes = sorted(pdf_map.keys())
        total = len(all_hashes)
        total_batches = -(-total // args.batch_size)  # ceiling division
        start = args.batch_index * args.batch_size
        end = start + args.batch_size
        batch_hashes = all_hashes[start:end]
        logger.info(
            "Batch %d/%d: processing hashes %d–%d (%d of %d total)",
            args.batch_index, total_batches - 1,
            start, min(end, total) - 1,
            len(batch_hashes), total,
        )
        if not batch_hashes:
            logger.warning(
                "Batch %d is empty (only %d hashes exist) — nothing to do",
                args.batch_index, total,
            )
            sys.exit(0)
        pdf_map = {h: pdf_map[h] for h in batch_hashes}

    if args.dry_run:
        n_would_full = n_would_partial = n_would_skip = 0
        expect_taxa = taxonomy_db is not None
        expect_anatomy = bool(anatomy_lexicon)
        for h in pdf_map:
            hd = documents_dir / short_hash(h)
            if not args.resume:
                n_would_full += 1
                continue
            if (hd / "summary.json").exists() and _all_stage_artifacts_complete(
                hd, expect_taxa=expect_taxa, expect_anatomy=expect_anatomy,
            ):
                n_would_skip += 1
            elif (hd / "summary.json").exists():
                # Partial — per-stage guards will run only the missing stages
                n_would_partial += 1
            else:
                n_would_full += 1
        logger.info(
            "Dry-run: %d unique PDF(s) in scope; would full-process %d, "
            "partial-process %d, skip %d (--resume = %s). Vision backend: %s. "
            "Grobid: %s. No files written.",
            len(pdf_map), n_would_full, n_would_partial, n_would_skip,
            "on" if args.resume else "off",
            args.vision_backend or "off",
            "off (--no-grobid or empty URL)" if (args.no_grobid or not args.grobid_url) else args.grobid_url,
        )
        return

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)

        for pdf_hash_full, pdf_paths in pdf_map.items():
            pdf_hash = short_hash(pdf_hash_full)
            hash_dir = documents_dir / pdf_hash

            # Detect prefix collision against any prior run in this dir.
            try:
                prior_matches = _verify_or_raise_collision(hash_dir, pdf_hash_full)
            except RuntimeError as e:
                logger.error(str(e))
                sys.exit(2)

            if args.resume and prior_matches:
                if args.refresh_vision and vision_backend is not None:
                    figures_file = hash_dir / "figures.json"
                    if not figures_file.exists():
                        logger.info(
                            "Skipping %s for vision refresh (no figures.json)", pdf_hash
                        )
                        continue
                    logger.info("Refreshing Pass 3b on %s", pdf_hash)
                    with per_pdf_file_log(hash_dir) as log_path:
                        logger.info("pipeline.log: %s (refresh-vision)", log_path)
                        try:
                            _pass3b_annotate_rois(figures_file, vision_backend)
                        except Exception as e:
                            logger.exception(
                                "Pass 3b refresh failed on %s: %s", pdf_hash, e
                            )
                    continue
                # Per-stage resume (#28): if every artifact is on disk,
                # skip the whole paper (fast path). Otherwise fall through
                # — run_pdf_processing_pipeline runs only the missing
                # stages.
                if _all_stage_artifacts_complete(
                    hash_dir,
                    expect_taxa=taxonomy_db is not None,
                    expect_anatomy=bool(anatomy_lexicon),
                ):
                    logger.info("Skipping %s (all stages complete)", pdf_hash)
                    continue
                logger.info(
                    "Resuming %s (re-running missing stages only)", pdf_hash
                )

            hash_dir.mkdir(exist_ok=True)

            # Use the first copy for processing (they're all identical by hash)
            primary_pdf = pdf_paths[0]

            logger.info(
                "Processing PDF hash %s (%d copies); primary file: %s",
                pdf_hash, len(pdf_paths), primary_pdf.relative_to(input_dir),
            )
            if len(pdf_paths) > 1:
                for path in pdf_paths[1:]:
                    logger.info("  additional copy: %s", path.relative_to(input_dir))

            # Run each document in a child process so that C-level crashes
            # (segfaults in docling/PyMuPDF) don't kill the whole batch.
            # We use the fork start method so the child inherits already-built
            # objects (grobid_client, taxonomy_db, vision_backend, _worker
            # closure) without pickling.
            #
            # macOS exception: PyTorch + Apple's Objective-C frameworks
            # (Metal/MPSGraph) are not fork-safe. The forked child crashes
            # the moment docling loads torch. There is no env-var workaround
            # that holds — OBJC_DISABLE_INITIALIZE_FORK_SAFETY just turns the
            # abort into a SIGSEGV. So on macOS we run inline and rely on
            # --resume to recover from the rare segfault.
            def _worker():
                with per_pdf_file_log(hash_dir) as log_path:
                    logger.info("pipeline.log: %s", log_path)
                    logger.info("pdf_hash_full: %s", pdf_hash_full)

                    processing_summary = run_pdf_processing_pipeline(
                        primary_pdf, hash_dir, temp_dir,
                        grobid_client=grobid_client,
                        taxonomy_db=taxonomy_db,
                        anatomy_lexicon=anatomy_lexicon,
                        content_aware_figures=args.content_aware_figures,
                        vision_backend=vision_backend,
                        bib_index=bib_index,
                        resume=args.resume,
                        taxonomy_fingerprint=taxonomy_fingerprint,
                        anatomy_fingerprint=anatomy_fingerprint,
                        extra_lexicons=extra_lexicons,
                        extra_fingerprints=extra_fingerprints,
                    )

                    summary_file = create_summary_json(
                        pdf_hash_full, pdf_paths, input_dir, hash_dir, processing_summary
                    )
                    logger.info("Created summary: %s", summary_file)

                    if processing_summary.get("status") == "success":
                        chunks_file = hash_dir / "chunks.json"
                        if chunks_file.exists():
                            logger.info("Writing vector-db ingestion marker...")
                            ingest_to_vector_db(chunks_file, vector_db_dir, pdf_hash)

            if sys.platform == "darwin":
                _worker()
            else:
                mp_ctx = multiprocessing.get_context("fork")
                proc = mp_ctx.Process(target=_worker)
                proc.start()
                proc.join()
                if proc.exitcode != 0:
                    if proc.exitcode < 0:
                        logger.error(
                            "Document %s killed by signal %d (segfault?) — skipping",
                            pdf_hash, -proc.exitcode,
                        )
                    else:
                        logger.error(
                            "Document %s worker exited with code %d — skipping",
                            pdf_hash, proc.exitcode,
                        )

    logger.info("Processing complete. Results saved to: %s", output_dir)
    logger.info("  Documents: %s", documents_dir)
    logger.info("  Vector DB: %s", vector_db_dir)


if __name__ == "__main__":
    main()
