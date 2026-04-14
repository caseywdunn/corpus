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
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, List, Optional, Set
import tempfile
import os

from grobid_client import (
    GrobidClient,
    GrobidUnavailableError,
    parse_tei_header,
    parse_tei_references,
)


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Section-class normalizer (Phase B)
# ---------------------------------------------------------------------------

# Ordered list — first match wins. Patterns are case-insensitive and cover
# English, German, French, Russian, plus the taxonomy-paper conventions
# common in our siphonophore corpus (Description / Systematics).
_SECTION_PATTERNS = [
    ("abstract", r"\babstract\b|\bsummary\b|\bzusammenfassung\b|\br[ée]sum[ée]\b|резюме|сводка"),
    ("introduction", r"\bintroduction\b|\bвведение\b|\beinleitung\b"),
    ("methods", r"materials?\s+and\s+methods|\bmethods?\b|\bmethodology\b|study\s+area|материалы\s+и\s+методы|m[ée]thodes"),
    ("results", r"\bresults?\b|\bр[еé]зультаты\b|\bergebnisse\b|\br[ée]sultats\b"),
    ("description", r"\bdescription\b|\bsystematics?\b|\btaxonomy\b|\bsystematic\s+account\b|\bdiagnosis\b|описание"),
    ("discussion", r"\bdiscussion\b|обсуждение|\bdiskussion\b"),
    ("conclusion", r"\bconclusions?\b|выводы|\bschlussfolgerung"),
    ("acknowledgements", r"acknowledg(?:e?ment)s?|благодарности|danksagung|remerciements"),
    ("references", r"\breferences?\b|\bbibliograph|\bliterature\s+cited\b|литература"),
    ("appendix", r"\bappendix\b|приложение|anhang"),
]


def classify_section(headings: Optional[List[str]]) -> Optional[str]:
    """Map a chunk's heading trail to a canonical section class.

    ``headings`` is the list of headings provided by docling's
    HybridChunker — outermost first, nearest (most specific) last. We
    inspect the most specific heading first, then walk outward, and
    return the first canonical class that any heading matches. Returns
    None if no heading matches a known pattern.
    """
    if not headings:
        return None
    for h in reversed(headings):
        if not h:
            continue
        hlow = h.lower()
        for cls, pat in _SECTION_PATTERNS:
            if re.search(pat, hlow):
                return cls
    return None


# ---------------------------------------------------------------------------
# Config + logging helpers (Phase A)
# ---------------------------------------------------------------------------

# Defaults used if config.yaml is missing or a key is absent.
_DEFAULT_CONFIG = {
    "ocr": {
        "language_detection": True,
        "redo_ocr": True,
        "deskew": True,
        "rotate_pages": True,
        "optimize_level": 2,
    },
    "chunking": {
        "max_tokens": 8191,
        "model_name": "text-embedding-3-small",
        "merge_peers": True,
    },
    "embedding": {
        "model": "text-embedding-3-small",
        "dimensions": 1536,
        "batch_size": 100,
    },
    "vector_db": {
        "path": "output/lancedb",
        "table_name": "document_chunks",
    },
    "qc": {
        "enable_validation": True,
        "sample_size": 5,
    },
}


def _deep_merge(base: dict, overlay: dict) -> dict:
    """Recursively merge overlay into a copy of base. Overlay values win."""
    out = dict(base)
    for k, v in (overlay or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(config_path: Optional[Path] = None) -> dict:
    """Load config.yaml (if present) and merge it over the built-in defaults.

    Missing file is not an error — defaults are returned. Missing keys fall
    back to defaults. A malformed YAML is an error (loud failure).
    """
    if config_path is None:
        config_path = Path(__file__).parent / "config.yaml"
    if not config_path.exists():
        logger.info("No config.yaml at %s; using built-in defaults.", config_path)
        return dict(_DEFAULT_CONFIG)
    try:
        import yaml
    except ImportError as e:
        raise RuntimeError(
            "pyyaml is required to read config.yaml; pip install pyyaml"
        ) from e
    with open(config_path, "r", encoding="utf-8") as f:
        overlay = yaml.safe_load(f) or {}
    merged = _deep_merge(_DEFAULT_CONFIG, overlay)
    logger.info("Loaded config from %s", config_path)
    return merged


# Module-level config; populated by main(), also safe-loaded on import with
# defaults so helper functions can run outside main() (e.g., in tests).
CONFIG: dict = dict(_DEFAULT_CONFIG)


def setup_root_logging(level: int = logging.INFO) -> None:
    """Configure the root logger with a single stderr stream handler.

    Per-PDF file handlers are added/removed around each document by
    ``per_pdf_file_log``; they coexist with this stream handler so output
    appears both on the terminal and in the per-paper pipeline.log.
    """
    root = logging.getLogger()
    root.setLevel(level)
    # Avoid duplicate stream handlers if called twice.
    for h in root.handlers:
        if isinstance(h, logging.StreamHandler) and getattr(h, "_corpus_stream", False):
            return
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter("%(levelname)s %(name)s: %(message)s"))
    sh._corpus_stream = True  # marker so we don't re-add on repeated calls
    root.addHandler(sh)


@contextmanager
def per_pdf_file_log(hash_dir: Path):
    """Attach a FileHandler writing to ``<hash_dir>/pipeline.log`` for the
    duration of the context. All ``logging`` calls made inside the block are
    captured to that file in addition to the root stream handler.
    """
    hash_dir.mkdir(parents=True, exist_ok=True)
    log_path = hash_dir / "pipeline.log"
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
    )
    root = logging.getLogger()
    root.addHandler(fh)
    try:
        yield log_path
    finally:
        root.removeHandler(fh)
        fh.close()


def create_cell_visualizations(input_pdf: Path, output_dir: Path, pdf_name: str, document):
    """Create cell visualization PNGs using docling-parse with figure regions highlighted"""
    try:
        from docling_core.types.doc.page import TextCellUnit
        from docling_parse.pdf_parser import DoclingPdfParser, PdfDocument
        from PIL import ImageDraw
        
        logger.info("Creating cell visualizations...")
        
        pdf_parser = DoclingPdfParser()
        pdf_doc: PdfDocument = pdf_parser.load(path_or_stream=str(input_pdf))
        
        # Get figure bounding boxes from the docling document
        figure_bboxes_by_page = {}
        if hasattr(document, 'pictures') and document.pictures:
            for picture in document.pictures:
                if hasattr(picture, 'prov') and picture.prov:
                    for prov_item in picture.prov:
                        if hasattr(prov_item, 'page_no') and hasattr(prov_item, 'bbox'):
                            page_no = prov_item.page_no
                            if page_no not in figure_bboxes_by_page:
                                figure_bboxes_by_page[page_no] = []
                            figure_bboxes_by_page[page_no].append(prov_item.bbox)
        
        # Create visualization for word-level cells
        cell_unit = TextCellUnit.WORD
        
        for page_no, pred_page in pdf_doc.iterate_pages():
            img = pred_page.render_as_image(cell_unit=cell_unit)
            draw = ImageDraw.Draw(img)
            img_height = img.height
            
            # First, highlight figure regions in yellow
            if page_no in figure_bboxes_by_page:
                for bbox in figure_bboxes_by_page[page_no]:
                    # Convert bbox coordinates (assumes bottom-left origin) to PIL coordinates (top-left origin)
                    x0 = bbox.l
                    x1 = bbox.r
                    y0_raw = bbox.t  # top in document coordinates
                    y1_raw = bbox.b  # bottom in document coordinates
                    # Flip Y coordinates for PIL (top-left origin)
                    y0 = img_height - y0_raw
                    y1 = img_height - y1_raw
                    
                    # Ensure coordinates are in correct order
                    x0, x1 = min(x0, x1), max(x0, x1)
                    y0, y1 = min(y0, y1), max(y0, y1)
                    
                    rect = (x0, y0, x1, y1)
                    draw.rectangle(rect, fill="yellow", outline="orange", width=3)
            
            # Then, draw red rectangles around text cells (these will show over non-figure areas)
            for cell in pred_page.iterate_cells(unit_type=cell_unit):
                # Convert cell.rect to (x0, y0, x1, y1) and flip y for PIL
                x0 = min(getattr(cell.rect, "r_x0", 0), getattr(cell.rect, "r_x2", 0))
                x1 = max(getattr(cell.rect, "r_x0", 0), getattr(cell.rect, "r_x2", 0))
                y0_raw = min(getattr(cell.rect, "r_y0", 0), getattr(cell.rect, "r_y2", 0))
                y1_raw = max(getattr(cell.rect, "r_y0", 0), getattr(cell.rect, "r_y2", 0))
                y0 = img_height - y1_raw
                y1 = img_height - y0_raw
                rect = (x0, y0, x1, y1)
                
                # Check if this cell overlaps with any figure region
                is_in_figure = False
                if page_no in figure_bboxes_by_page:
                    for bbox in figure_bboxes_by_page[page_no]:
                        fig_x0, fig_x1 = bbox.l, bbox.r
                        fig_y0_raw, fig_y1_raw = bbox.t, bbox.b
                        fig_y0 = img_height - fig_y0_raw
                        fig_y1 = img_height - fig_y1_raw
                        fig_x0, fig_x1 = min(fig_x0, fig_x1), max(fig_x0, fig_x1)
                        fig_y0, fig_y1 = min(fig_y0, fig_y1), max(fig_y0, fig_y1)
                        
                        # Check for overlap
                        if (x0 < fig_x1 and x1 > fig_x0 and y0 < fig_y1 and y1 > fig_y0):
                            is_in_figure = True
                            break
                
                # Only draw cell boundaries for text regions (not in figures)
                if not is_in_figure:
                    draw.rectangle(rect, outline="red", width=2)
            
            # Save the visualization
            viz_filename = f"page_{page_no}_visualization.png"
            viz_path = output_dir / viz_filename
            img.save(str(viz_path))
            
        logger.info("Saved %d visualization PNGs", len(list(pdf_doc.iterate_pages())))

    except ImportError as e:
        logger.warning("Could not create visualizations — missing dependencies: %s", e)
    except Exception as e:
        logger.warning("Could not create visualizations: %s", e)


HASH_PREFIX_LEN = 12  # 48 bits; collision-safe up to ~1e7 documents


def calculate_pdf_hash(pdf_path: Path) -> str:
    """Calculate the full SHA-256 hex digest of a PDF file.

    The per-document directory uses a 12-char lowercase prefix of this digest
    (see :data:`HASH_PREFIX_LEN`); the full digest is recorded in each
    ``summary.json`` so we can verify on resume that a re-encountered prefix
    really identifies the same PDF (rather than a hash-prefix collision).
    """
    hasher = hashlib.sha256()
    with open(pdf_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def short_hash(full_hash: str) -> str:
    """Return the 12-char lowercase prefix used as the per-document dir name."""
    return full_hash[:HASH_PREFIX_LEN].lower()


def find_all_pdfs(input_dir: Path) -> Dict[str, List[Path]]:
    """Recursively find all PDFs under ``input_dir`` and group by full SHA-256.

    Returns a dict mapping full hex digest → list of paths that share it
    (duplicates). The caller derives the short directory name via
    :func:`short_hash`.
    """
    pdf_map: Dict[str, List[Path]] = {}
    for pdf_path in input_dir.rglob("*.pdf"):
        if not pdf_path.is_file():
            continue
        try:
            full_hash = calculate_pdf_hash(pdf_path)
            pdf_map.setdefault(full_hash, []).append(pdf_path)
        except Exception as e:
            logger.warning("Could not hash %s: %s", pdf_path, e)
    return pdf_map


def create_output_structure(output_dir: Path):
    """Create the output directory structure."""
    documents_dir = output_dir / "documents"
    vector_db_dir = output_dir / "vector_db"
    
    documents_dir.mkdir(parents=True, exist_ok=True)
    vector_db_dir.mkdir(parents=True, exist_ok=True)
    
    return documents_dir, vector_db_dir


def get_relative_paths(pdf_paths: List[Path], input_dir: Path) -> List[str]:
    """Get relative paths of PDFs from input directory."""
    return [str(path.relative_to(input_dir)) for path in pdf_paths]


def run_pdf_processing_pipeline(
    pdf_path: Path,
    hash_dir: Path,
    temp_dir: Path,
    grobid_client: Optional[GrobidClient] = None,
) -> Dict:
    """Run the per-PDF processing pipeline and return a summary dict.

    Called inside a :func:`per_pdf_file_log` block by the main loop, so all
    ``logging`` calls here are captured to ``<hash_dir>/pipeline.log``.
    """
    pdf_name = pdf_path.stem
    temp_pdf = temp_dir / f"{pdf_name}.pdf"

    # Copy PDF to temp directory for processing
    shutil.copy2(pdf_path, temp_pdf)

    # Create subdirectories in hash directory
    figures_dir = hash_dir / "figures"
    figures_dir.mkdir(exist_ok=True)
    visualizations_dir = hash_dir / "visualizations"
    visualizations_dir.mkdir(exist_ok=True)

    processing_summary = {
        "original_pdf": str(pdf_path),
        "processing_steps": [],
        "files_created": [],
        "errors": [],
    }

    try:
        logger.info("Detecting scan type...")
        detection_result = detect_scan_type(temp_pdf)
        detection_file = hash_dir / "scan_detection.json"
        with open(detection_file, "w") as f:
            json.dump(detection_result, f, indent=2)
        processing_summary["files_created"].append(str(detection_file))
        processing_summary["processing_steps"].append("scan_detection")

        logger.info("Preparing PDF...")
        processed_pdf = hash_dir / "processed.pdf"
        prepare_pdf(temp_pdf, detection_result, processed_pdf)
        processing_summary["files_created"].append(str(processed_pdf))
        processing_summary["processing_steps"].append("pdf_preparation")

        logger.info("Extracting text and figures...")
        text_file = hash_dir / "text.json"
        figures_file = hash_dir / "figures.json"
        docling_doc_file = hash_dir / "docling_doc.json"
        extract_docling_content(
            processed_pdf,
            text_file,
            figures_file,
            figures_dir,
            visualizations_dir,
            docling_doc_output=docling_doc_file,
        )
        processing_summary["files_created"].extend([str(text_file), str(figures_file)])
        if docling_doc_file.exists():
            processing_summary["files_created"].append(str(docling_doc_file))
        processing_summary["processing_steps"].append("docling_extraction")

        logger.info("Extracting metadata (Grobid)...")
        metadata_file = hash_dir / "metadata.json"
        references_file = hash_dir / "references.json"
        tei_file = hash_dir / "grobid.tei.xml"
        extract_metadata(
            processed_pdf,
            metadata_file,
            references_output=references_file,
            tei_output=tei_file,
            grobid_client=grobid_client,
        )
        processing_summary["files_created"].extend(
            [str(metadata_file), str(references_file)]
        )
        if tei_file.exists():
            processing_summary["files_created"].append(str(tei_file))
        processing_summary["processing_steps"].append("metadata_extraction")

        logger.info("Chunking text...")
        chunks_file = hash_dir / "chunks.json"
        chunk_text(text_file, metadata_file, chunks_file)
        processing_summary["files_created"].append(str(chunks_file))
        processing_summary["processing_steps"].append("text_chunking")

        processing_summary["status"] = "success"

    except Exception as e:
        processing_summary["status"] = "error"
        processing_summary["errors"].append(str(e))
        logger.exception("Error processing PDF: %s", e)

    return processing_summary


def detect_scan_type(pdf_path: Path) -> Dict:
    """Detect if PDF is scanned or born-digital."""
    try:
        import fitz  # PyMuPDF
        doc = fitz.open(pdf_path)
        total_text = ""
        total_chars = 0
        
        pages_to_check = min(3, len(doc))
        
        for page_num in range(pages_to_check):
            page = doc[page_num]
            page_text = page.get_text()
            total_text += page_text
            total_chars += len(page_text.strip())
        
        doc.close()
        
        meaningful_text = ''.join(total_text.split())
        avg_chars_per_page = total_chars / pages_to_check if pages_to_check > 0 else 0
        
        has_substantial_text = len(meaningful_text) > 500
        has_good_density = avg_chars_per_page > 100
        has_text = has_substantial_text and has_good_density
        
        return {
            "filename": pdf_path.name,
            "has_text": has_text,
            "needs_ocr": not has_text,
            "file_type": "born_digital" if has_text else "scanned"
        }
        
    except ImportError:
        # If PyMuPDF isn't available, avoid expensive OCR; proceed as born-digital.
        logger.warning("PyMuPDF not available; treating as born-digital (no OCR)")
        return {
            "filename": pdf_path.name,
            "has_text": True,
            "needs_ocr": False,
            "file_type": "born_digital",
        }
    except Exception as e:
        logger.warning("Error checking text content: %s", e)
        return {
            "filename": pdf_path.name,
            "has_text": False,
            "needs_ocr": True,
            "file_type": "scanned"
        }


def prepare_pdf(input_pdf: Path, detection_result: Dict, output_pdf: Path):
    """Prepare PDF by running OCR if needed."""
    needs_ocr = detection_result.get('needs_ocr', False)
    
    if needs_ocr:
        # Check if ocrmypdf is available
        if shutil.which("ocrmypdf") is None:
            logger.warning("ocrmypdf not found, copying original PDF")
            shutil.copy2(input_pdf, output_pdf)
            return

        optimize_level = str(CONFIG.get("ocr", {}).get("optimize_level", 2))
        logger.info("Running OCR on %s (detected as scanned)", input_pdf.name)
        cmd = [
            "ocrmypdf",
            "--force-ocr",
            "--optimize", optimize_level,
            "--color-conversion-strategy", "RGB",
            "--output-type", "pdf",
            str(input_pdf),
            str(output_pdf),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.warning("OCR failed, copying original PDF")
            logger.warning("OCR Error: %s", result.stderr)
            shutil.copy2(input_pdf, output_pdf)
        else:
            logger.info("OCR completed successfully")
    else:
        logger.info("Copying %s (detected as born-digital)", input_pdf.name)
        shutil.copy2(input_pdf, output_pdf)


def extract_docling_content(
    pdf_path: Path,
    text_output: Path,
    figures_output: Path,
    figures_dir: Path,
    visualizations_dir: Path,
    docling_doc_output: Optional[Path] = None,
):
    """Extract text and figures using docling, with PyMuPDF fallback for figures.

    Guarantees that only successfully saved figures are listed in figures.json.
    Also persists the parsed :class:`DoclingDocument` to
    ``docling_doc_output`` (defaults to ``<text_output.parent>/docling_doc.json``)
    so the chunking step can re-load the structured document and run
    HybridChunker without re-parsing the PDF.
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
                "caption": caption or "",
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

    # First try: use docling pictures if document is available
    if document is not None and hasattr(document, "pictures"):
        saved_count = 0
        for idx, picture in enumerate(document.pictures or []):
            figure_path = figures_dir / f"figure_{idx + 1}.png"

            # Caption (best-effort)
            caption = ""
            try:
                if hasattr(picture, "caption_text"):
                    caption = picture.caption_text(document) if callable(picture.caption_text) else (
                        str(picture.caption_text) if picture.caption_text else ""
                    )
            except Exception as e:
                logger.info("Caption extraction failed for figure %d: %s", idx + 1, e)

            # Save image if available
            try:
                image = None
                if hasattr(picture, "get_image"):
                    image = picture.get_image(document) if callable(picture.get_image) else picture.get_image
                if image is not None and hasattr(image, "save"):
                    image.save(str(figure_path))
                    saved_count += 1
                    meta = {"extraction_method": "docling"}
                    meta.update(_docling_prov_to_bbox_page(picture))
                    append_figure(
                        figure_id=f"docling_{idx + 1}",
                        figure_path=figure_path,
                        caption=caption,
                        meta=meta,
                    )
                else:
                    logger.warning("Docling did not return a savable image for figure %d", idx + 1)
            except Exception as e:
                logger.warning("Could not save docling figure %d: %s", idx + 1, e)

        if saved_count == 0:
            logger.info("No figures saved via docling; will try PyMuPDF fallback")

    # Fallback: use PyMuPDF to extract embedded images from the PDF
    if len(figures_data) == 0:
        try:
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

    # Create cell visualizations only if we have a document from docling
    if document is not None:
        pdf_name = pdf_path.stem
        try:
            create_cell_visualizations(pdf_path, visualizations_dir, pdf_name, document)
        except Exception as e:
            logger.warning("Creating visualizations failed: %s", e)


def _write_placeholder_metadata(pdf_path: Path, metadata_output: Path, references_output: Path):
    """Write empty-but-valid metadata.json and references.json.

    Used when Grobid is unavailable or its TEI can't be parsed; keeps
    downstream stages (chunking, embedding) functional and lets the
    summary show this document as metadata-less so we can triage later.
    """
    placeholder = {
        "filename": pdf_path.name,
        "title": "",
        "authors": [],
        "year": None,
        "journal": "",
        "doi": "",
        "abstract": "",
        "extraction_method": "placeholder",
    }
    with open(metadata_output, "w", encoding="utf-8") as f:
        json.dump(placeholder, f, indent=2, ensure_ascii=False)
    with open(references_output, "w", encoding="utf-8") as f:
        json.dump({"references": [], "total_references": 0}, f, indent=2)


def extract_metadata(
    pdf_path: Path,
    metadata_output: Path,
    references_output: Optional[Path] = None,
    tei_output: Optional[Path] = None,
    grobid_client: Optional[GrobidClient] = None,
):
    """Extract bibliographic metadata + references via Grobid.

    Runs ``/api/processFulltextDocument`` once, caches the raw TEI-XML at
    ``tei_output`` (if given), and parses it into ``metadata.json`` plus
    ``references.json``. If Grobid is unreachable or any step fails, falls
    back to placeholder files so downstream stages still run — the caller
    can inspect ``metadata["extraction_method"]`` to tell which path was
    taken.

    Parameters
    ----------
    pdf_path:
        The (possibly OCR'd) ``processed.pdf`` to send to Grobid.
    metadata_output:
        Output path for ``metadata.json`` (title, authors, year, …).
    references_output:
        Output path for ``references.json``. If None, defaults to
        ``<metadata_output.parent>/references.json``.
    tei_output:
        Cache path for the raw Grobid TEI. Existing non-empty TEI at this
        path is reused (skipping the Grobid call) — convenient for
        re-parsing without re-processing. If None, defaults to
        ``<metadata_output.parent>/grobid.tei.xml``.
    grobid_client:
        A live :class:`GrobidClient`, or None to skip Grobid entirely
        (placeholder-only mode).
    """
    hash_dir = metadata_output.parent
    if references_output is None:
        references_output = hash_dir / "references.json"
    if tei_output is None:
        tei_output = hash_dir / "grobid.tei.xml"

    tei_xml: Optional[str] = None

    # 1. Reuse cached TEI if present.
    if tei_output.exists() and tei_output.stat().st_size > 0:
        try:
            tei_xml = tei_output.read_text(encoding="utf-8")
            logger.info("Reusing cached Grobid TEI at %s", tei_output)
        except OSError as e:
            logger.warning("Could not read cached TEI %s: %s", tei_output, e)
            tei_xml = None

    # 2. Otherwise, call Grobid if a client is available.
    if tei_xml is None and grobid_client is not None:
        try:
            logger.info("Calling Grobid on %s ...", pdf_path.name)
            tei_xml = grobid_client.process_fulltext(pdf_path)
            tei_output.write_text(tei_xml, encoding="utf-8")
            logger.info("Wrote Grobid TEI to %s", tei_output)
        except GrobidUnavailableError as e:
            logger.warning("Grobid unavailable, using placeholder: %s", e)
        except Exception as e:
            logger.warning("Grobid call failed on %s: %s", pdf_path.name, e)

    # 3. Parse TEI if we have one; fall back on any parsing error.
    if tei_xml is not None:
        try:
            header = parse_tei_header(tei_xml)
            header["filename"] = pdf_path.name  # keep original filename for provenance
            with open(metadata_output, "w", encoding="utf-8") as f:
                json.dump(header, f, indent=2, ensure_ascii=False)

            refs = parse_tei_references(tei_xml)
            with open(references_output, "w", encoding="utf-8") as f:
                json.dump(
                    {"references": refs, "total_references": len(refs)},
                    f,
                    indent=2,
                    ensure_ascii=False,
                )
            logger.info(
                "Parsed Grobid TEI: %d authors, %d references, year=%s",
                len(header.get("authors", [])),
                len(refs),
                header.get("year"),
            )
            return
        except Exception as e:
            logger.warning(
                "Failed to parse Grobid TEI (%s); writing placeholder", e
            )

    # 4. Placeholder path.
    _write_placeholder_metadata(pdf_path, metadata_output, references_output)


def chunk_text(
    text_file: Path,
    metadata_file: Path,
    chunks_output: Path,
    docling_doc_file: Optional[Path] = None,
):
    """Chunk the document into structurally-aware pieces.

    When a serialized :class:`DoclingDocument` is available (``docling_doc.json``
    produced by :func:`extract_docling_content`), we drive docling's
    :class:`HybridChunker` — tokenizer-aware, respects headings, table and
    figure captions. Each chunk carries its heading trail and a derived
    ``section_class`` (see :func:`classify_section`).

    When no DoclingDocument is available (e.g., docling import failed
    upstream), we fall back to the prior naive character-window behavior so
    downstream stages still run.
    """
    with open(metadata_file, "r") as f:
        metadata = json.load(f)

    # Resolve default docling_doc_file relative to text_file's directory.
    if docling_doc_file is None:
        docling_doc_file = text_file.parent / "docling_doc.json"

    chunks: List[dict] = []
    chunker_name = "hybrid_chunker"

    if docling_doc_file.exists() and docling_doc_file.stat().st_size > 0:
        try:
            from docling.chunking import HybridChunker
            from docling_core.types.doc import DoclingDocument

            dl_doc = DoclingDocument.load_from_json(docling_doc_file)
            chunker = HybridChunker()
            for i, c in enumerate(chunker.chunk(dl_doc=dl_doc)):
                headings = list(getattr(c.meta, "headings", []) or [])
                captions = list(getattr(c.meta, "captions", []) or [])
                chunks.append(
                    {
                        "chunk_id": f"chunk_{i}",
                        "text": c.text,
                        "headings": headings,
                        "section_class": classify_section(headings),
                        "captions": captions,
                    }
                )
            logger.info("HybridChunker produced %d chunks", len(chunks))
        except Exception as e:
            logger.warning(
                "HybridChunker failed (%s); falling back to naive char chunker", e
            )
            chunks = []
            chunker_name = "naive_char_window"

    if not chunks:
        # Naive fallback.
        chunker_name = "naive_char_window"
        with open(text_file, "r", encoding="utf-8") as f:
            text_data = json.load(f)
        text = text_data.get("text", "")

        chunk_size = int(CONFIG.get("chunking", {}).get("max_tokens", 1000))
        if chunk_size <= 0:
            chunk_size = 1000

        for i in range(0, len(text), chunk_size):
            piece = text[i : i + chunk_size]
            chunks.append(
                {
                    "chunk_id": f"chunk_{len(chunks)}",
                    "text": piece,
                    "headings": [],
                    "section_class": None,
                    "captions": [],
                    "start_char": i,
                    "end_char": min(i + chunk_size, len(text)),
                }
            )
        logger.info("Naive chunker produced %d chunks", len(chunks))

    chunks_data = {
        "metadata": metadata,
        "chunker": chunker_name,
        "total_chunks": len(chunks),
        "chunks": chunks,
    }

    with open(chunks_output, "w", encoding="utf-8") as f:
        json.dump(chunks_data, f, indent=2, ensure_ascii=False)


def ingest_to_vector_db(chunks_file: Path, vector_db_dir: Path, pdf_hash: str):
    """Ingest chunks into vector database."""
    # Placeholder for vector database ingestion
    # This would typically involve embedding the text and storing in a vector database
    
    # Create a simple marker file for now
    ingestion_marker = vector_db_dir / f"{pdf_hash}_embedded.done"
    
    with open(ingestion_marker, 'w') as f:
        json.dump({
            "pdf_hash": pdf_hash,
            "chunks_file": str(chunks_file),
            "ingestion_timestamp": str(Path(chunks_file).stat().st_mtime),
            "status": "completed"
        }, f, indent=2)


def create_summary_json(
    pdf_hash_full: str,
    pdf_paths: List[Path],
    input_dir: Path,
    hash_dir: Path,
    processing_summary: Dict,
):
    """Write ``summary.json`` for one document.

    Records both the short directory prefix and the full SHA-256 so that
    ``--resume`` can verify that a re-encountered prefix refers to the same
    PDF (not a hash-prefix collision).
    """
    relative_paths = get_relative_paths(pdf_paths, input_dir)

    summary = {
        "pdf_hash": short_hash(pdf_hash_full),
        "pdf_hash_full": pdf_hash_full,
        "hash_algorithm": "sha256",
        "input_dir": str(input_dir),
        "relative_paths": relative_paths,
        "total_copies_found": len(pdf_paths),
        "processing_summary": processing_summary,
        "output_directory": str(hash_dir),
    }

    summary_file = hash_dir / "summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    return summary_file


def _verify_or_raise_collision(hash_dir: Path, pdf_hash_full: str) -> Optional[bool]:
    """If ``hash_dir/summary.json`` exists, verify its recorded full hash
    matches ``pdf_hash_full``. Returns True if it matches (resume-safe), False
    if no summary is present (fresh dir), and raises ``RuntimeError`` on a
    real hash-prefix collision.
    """
    summary_file = hash_dir / "summary.json"
    if not summary_file.exists():
        return False
    try:
        with open(summary_file, "r") as f:
            existing = json.load(f)
    except Exception as e:
        logger.warning("Could not read %s (%s); treating as incomplete", summary_file, e)
        return False
    existing_full = existing.get("pdf_hash_full")
    if existing_full is None:
        # Legacy summary from before we recorded full hashes. Trust it but warn.
        logger.warning(
            "Existing summary at %s has no pdf_hash_full; cannot verify "
            "against prefix collision. Treating as a match.",
            summary_file,
        )
        return True
    if existing_full != pdf_hash_full:
        raise RuntimeError(
            f"Hash-prefix collision detected at {hash_dir}: "
            f"existing summary records full hash {existing_full!r} but this "
            f"PDF hashes to {pdf_hash_full!r}. Increase HASH_PREFIX_LEN or "
            f"investigate duplicate inputs."
        )
    return True


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

    args = parser.parse_args()

    setup_root_logging()

    global CONFIG
    CONFIG = load_config(args.config)

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.exists():
        logger.error("Input directory %s does not exist", input_dir)
        sys.exit(1)

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

    # Create output directory structure
    documents_dir, vector_db_dir = create_output_structure(output_dir)

    # Find all PDFs and group by hash
    logger.info("Discovering PDFs...")
    pdf_map = find_all_pdfs(input_dir)

    logger.info("Found %d PDF file(s)", sum(len(paths) for paths in pdf_map.values()))
    logger.info("Unique PDFs (by hash): %d", len(pdf_map))

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
                logger.info("Skipping %s (already processed)", pdf_hash)
                continue

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

            # Everything below is captured to the per-PDF pipeline.log.
            with per_pdf_file_log(hash_dir) as log_path:
                logger.info("pipeline.log: %s", log_path)
                logger.info("pdf_hash_full: %s", pdf_hash_full)

                processing_summary = run_pdf_processing_pipeline(
                    primary_pdf, hash_dir, temp_dir, grobid_client=grobid_client
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

    logger.info("Processing complete. Results saved to: %s", output_dir)
    logger.info("  Documents: %s", documents_dir)
    logger.info("  Vector DB: %s", vector_db_dir)


if __name__ == "__main__":
    main()
