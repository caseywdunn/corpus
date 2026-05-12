"""Docling extraction stage.

:func:`extract_docling_content` parses the (already-OCR'd-or-passed-
through) PDF into a structured :class:`DoclingDocument`, persists it
as ``docling_doc.json``, and emits ``text.json`` + ``figures.json``
plus per-figure PNGs. PyMuPDF fallback covers cases where docling
finds no Picture items, gated on the scan file_type so rasterized
plates don't get one-bitmap-per-page noise.
"""
from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

from . import stamp_artifact
from .figures import (
    detect_missing_figures,
    extract_caption_info,
    parse_figure_number,
)
from .scan import create_cell_visualizations

logger = logging.getLogger(__name__)


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
        logger.error("Docling not available (%s); cannot extract text", e)
        raise
    except Exception as e:
        logger.error("Docling extraction failed (%s); cannot extract text", e)
        raise

    # Save text content (docling succeeded — text_content is populated).
    # No placeholder fallback: a paper with no extracted text would be
    # embedded and served as fake content. The stage raises above so
    # _stage records a stage_failure and the outer pipeline marks
    # status=error; --resume retries cleanly because text.json isn't
    # written.
    with open(text_output, "w", encoding="utf-8") as f:
        json.dump(stamp_artifact(text_content), f, indent=2, ensure_ascii=False)

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
                "caption_text": caption or "",  # canonical; see dev_docs/PLAN.md §3
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

    # Two-pass docling figure extraction (Phase D.2 — see dev_docs/PLAN.md).
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
        from .figures import (
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
            from .figures import classify_figure
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
        json.dump(stamp_artifact(figures_info), f, indent=2)

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
