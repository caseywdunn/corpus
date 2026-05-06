"""Scan classification + OCR stage.

Largest single stage by line count. Bundles three concerns that
operate on the same input (PDF) and shared state (CONFIG, the
Tesseract pack map, the visual-script cross-check):

* Language + script detection on the existing text layer
  (``_detect_language``, ``_gibberish_score``, ``_text_layer_scripts``,
  ``_visual_page_script``).
* Tesseract pack resolution
  (``_available_tesseract_langs``, ``_resolve_tesseract_packs``,
  ``_compose_ocr_langs``).
* Three-class scan classification + ocrmypdf invocation
  (``detect_scan_type``, ``prepare_pdf``, ``_annotate_pack_availability``).
* Per-page QC visualizations (``create_cell_visualizations``).

Outputs: ``scan_detection.json`` (from detect_scan_type) and
``processed.pdf`` (from prepare_pdf, the OCR'd-or-passthrough result).
"""
from __future__ import annotations

import json
import logging
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from .config import CONFIG

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Language detection & OCR language selection (Phase C)
# ---------------------------------------------------------------------------

# ISO 639-1 → Tesseract traineddata code. Covers every language langdetect
# can return so we can check pack availability for any detected document. A
# detected language not in this map (or whose pack isn't installed) falls
# back to the CONFIG default union and triggers a warning.
_ISO_TO_TESSERACT = {
    # Latin-script European
    "en": "eng",
    "de": "deu",
    "fr": "fra",
    "la": "lat",
    "it": "ita",
    "es": "spa",
    "pt": "por",
    "nl": "nld",
    "pl": "pol",
    "sv": "swe",
    "no": "nor",
    "da": "dan",
    "fi": "fin",
    "ca": "cat",
    "cs": "ces",
    "hu": "hun",
    "ro": "ron",
    "sk": "slk",
    "sl": "slv",
    "hr": "hrv",
    "et": "est",
    "lv": "lav",
    "lt": "lit",
    "sq": "sqi",
    "cy": "cym",
    "af": "afr",
    "id": "ind",
    "sw": "swa",
    "tl": "tgl",
    "vi": "vie",
    "tr": "tur",
    # Cyrillic
    "ru": "rus",
    "uk": "ukr",
    "bg": "bul",
    "mk": "mkd",
    # Other scripts
    "el": "ell",
    "ar": "ara",
    "fa": "fas",
    "ur": "urd",
    "he": "heb",
    "hi": "hin",
    "mr": "mar",
    "ne": "nep",
    "bn": "ben",
    "gu": "guj",
    "pa": "pan",
    "ta": "tam",
    "te": "tel",
    "kn": "kan",
    "ml": "mal",
    "th": "tha",
    "ja": "jpn",
    "ko": "kor",
    "zh-cn": "chi_sim",
    "zh-tw": "chi_tra",
}


def _detect_language(text: str):
    """Return (iso_code, confidence) for the dominant language of ``text``.

    Returns (None, 0.0) on very short input or any detection failure. A seed
    is set on ``DetectorFactory`` for reproducibility (langdetect uses a
    randomized algorithm by default).
    """
    if not text or len(text.strip()) < 100:
        return None, 0.0
    try:
        from langdetect import detect_langs, DetectorFactory
    except ImportError:
        logger.warning("langdetect not installed; skipping language detection")
        return None, 0.0
    DetectorFactory.seed = 0
    try:
        # First ~10k chars is plenty; avoids slow detection on huge docs.
        results = detect_langs(text[:10000])
        if results:
            top = results[0]
            return top.lang, float(top.prob)
    except Exception as e:
        logger.debug("langdetect failed: %s", e)
    return None, 0.0


# Compiled once. Words of length <=2, or letters interleaved with digits, or
# all-digit tokens — all characteristic of text-layer corruption (Cyrillic
# mapped to Latin-1, hidden-layer OCR, etc.).
_DIGIT_IN_WORD_RE = re.compile(r"\d")


def _gibberish_score(text: str) -> float:
    """Fraction of tokens in ``text`` that look like text-layer garbage.

    A small score (<0.25) is normal even for clean English. Thresholds around
    0.5 reliably separate real prose from PDFs whose text layer maps Cyrillic
    glyphs to Latin-1 byte-by-byte (producing strings like "AKAllEMH5I HAYK").
    langdetect confidence alone does not catch this case because the mangled
    text is made of Latin bytes.
    """
    words = [w.strip(".,;:!?()[]{}\"'") for w in text.split()]
    words = [w for w in words if w]
    if not words:
        return 0.0
    sus = 0
    for w in words:
        if len(w) <= 2:
            sus += 1
            continue
        if w.isdigit():
            sus += 1
            continue
        if _DIGIT_IN_WORD_RE.search(w) and any(c.isalpha() for c in w):
            sus += 1
    return sus / len(words)


def _text_layer_scripts(text: str) -> Dict[str, float]:
    """Return the fraction of each major writing system present in ``text``.

    Letters outside the covered ranges are bucketed under ``"Other"``. Non-
    letter characters (digits, punctuation, whitespace) are ignored.

    Used as a cheap script fingerprint for cross-checking against Tesseract
    OSD. When a PDF's text layer is 95%+ Latin-family but the visual page
    image is Cyrillic (or Greek, etc.), the text layer is almost certainly
    a broken byte-for-byte mapping — the Stepanjants 1970 case.
    """
    if not text:
        return {}
    counts = {"Latin": 0, "Cyrillic": 0, "Greek": 0, "Arabic": 0, "Other": 0}
    for c in text:
        if not c.isalpha():
            continue
        cp = ord(c)
        if (0x0000 <= cp <= 0x024F) or (0x1E00 <= cp <= 0x1EFF) or (0x2C60 <= cp <= 0x2C7F):
            counts["Latin"] += 1
        elif 0x0400 <= cp <= 0x052F:
            counts["Cyrillic"] += 1
        elif 0x0370 <= cp <= 0x03FF:
            counts["Greek"] += 1
        elif 0x0600 <= cp <= 0x06FF:
            counts["Arabic"] += 1
        else:
            counts["Other"] += 1
    total = sum(counts.values())
    if total == 0:
        return {}
    return {k: v / total for k, v in counts.items() if v > 0}


def _visual_page_script(pdf_path: Path) -> Optional[str]:
    """Run Tesseract's script-detection OSD on sampled page images.

    Samples three pages spread through the document (≈25%, 50%, 75%) and
    runs Tesseract ``--psm 0`` OSD on each. Returns the first non-Latin
    script encountered, or ``"Latin"`` if every sampled page comes back
    Latin-family. ``None`` means OSD failed on every attempt.

    Why three pages spread through the document: historical papers often
    have a Latin-script cover/title page (journal metadata, publisher
    boilerplate) with content in a different script starting several pages
    in. The Stepanjants 1970 monograph is the canonical example — pages
    0–4 are Russian journal frontmatter in Latin-transliterated titles,
    pages 5+ are Cyrillic body text. Sampling a single early page misses
    the signal.

    Rendering at 200 DPI balances OSD accuracy against cost (~0.5–1.5 s
    per call per page). Called only for documents whose text layer is
    suspicious (dense + mostly Latin), not every PDF.
    """
    try:
        import fitz
    except ImportError:
        return None

    try:
        doc = fitz.open(pdf_path)
    except Exception as e:
        logger.debug("Could not open %s for OSD: %s", pdf_path, e)
        return None

    try:
        n = len(doc)
        if n == 0:
            return None
        if n >= 4:
            page_indices = sorted({n // 4, n // 2, (3 * n) // 4})
        else:
            page_indices = list(range(n))

        scripts_seen: List[str] = []
        for idx in page_indices:
            try:
                pix = doc[idx].get_pixmap(dpi=200)
                img_bytes = pix.tobytes("png")
            except Exception as e:
                logger.debug("Render failed on page %d: %s", idx, e)
                continue
            try:
                result = subprocess.run(
                    ["tesseract", "-", "-", "--psm", "0"],
                    input=img_bytes, capture_output=True, timeout=30,
                )
            except Exception as e:
                logger.debug("OSD subprocess failed on page %d: %s", idx, e)
                continue
            combined = (
                result.stdout.decode("utf-8", errors="replace")
                + result.stderr.decode("utf-8", errors="replace")
            )
            m = re.search(r"Script:\s*(\S+)", combined)
            if m:
                scripts_seen.append(m.group(1))

        if not scripts_seen:
            return None
        # Prefer any non-Latin-family script — these are the signal for
        # broken text layers. A single Cyrillic-page hit beats two Latin
        # hits, because Cyrillic body content is impossible for a text
        # layer that's 100% Latin-family to honestly represent.
        for s in scripts_seen:
            if s not in ("Latin", "Fraktur"):
                return s
        return scripts_seen[0]
    finally:
        doc.close()


# ISO 639-1 code → expected script (for cross-check). Latin-family languages
# also render as "Fraktur" for historical German — Tesseract reports that
# as a distinct script, so we treat Latin and Fraktur as compatible here.
# Script names match what Tesseract OSD emits ("Japanese", "Han", "Hangul",
# "Devanagari", etc.).
_LANG_EXPECTED_SCRIPT = {
    # Latin
    "en": "Latin", "de": "Latin", "fr": "Latin", "it": "Latin", "es": "Latin",
    "pt": "Latin", "nl": "Latin", "pl": "Latin", "sv": "Latin", "no": "Latin",
    "da": "Latin", "fi": "Latin", "la": "Latin", "ca": "Latin", "cs": "Latin",
    "hu": "Latin", "ro": "Latin", "sk": "Latin", "sl": "Latin", "hr": "Latin",
    "et": "Latin", "lv": "Latin", "lt": "Latin", "sq": "Latin", "cy": "Latin",
    "af": "Latin", "id": "Latin", "sw": "Latin", "tl": "Latin", "vi": "Latin",
    "tr": "Latin",
    # Cyrillic
    "ru": "Cyrillic", "uk": "Cyrillic", "bg": "Cyrillic", "sr": "Cyrillic",
    "mk": "Cyrillic", "be": "Cyrillic",
    # Other scripts
    "el": "Greek",
    "ar": "Arabic", "fa": "Arabic", "ur": "Arabic",
    "he": "Hebrew",
    "hi": "Devanagari", "mr": "Devanagari", "ne": "Devanagari",
    "bn": "Bengali",
    "gu": "Gujarati",
    "pa": "Gurmukhi",
    "ta": "Tamil", "te": "Telugu", "kn": "Kannada", "ml": "Malayalam",
    "th": "Thai",
    "ja": "Japanese",
    "ko": "Hangul",
    "zh-cn": "Han", "zh-tw": "Han",
}


_AVAILABLE_TESSERACT_LANGS_CACHE: Optional[frozenset] = None


def _available_tesseract_langs() -> frozenset:
    """Return the set of Tesseract language codes installed on this system.

    Shells out to ``tesseract --list-langs`` once and caches the result. If
    Tesseract isn't on PATH, returns an empty set (OCR will simply be a no-op
    with a warning).
    """
    global _AVAILABLE_TESSERACT_LANGS_CACHE
    if _AVAILABLE_TESSERACT_LANGS_CACHE is not None:
        return _AVAILABLE_TESSERACT_LANGS_CACHE
    if shutil.which("tesseract") is None:
        logger.warning("tesseract not found on PATH; OCR will be skipped")
        _AVAILABLE_TESSERACT_LANGS_CACHE = frozenset()
        return _AVAILABLE_TESSERACT_LANGS_CACHE
    try:
        out = subprocess.run(
            ["tesseract", "--list-langs"],
            capture_output=True, text=True, check=True,
        )
        # Output format: first line is a header, subsequent lines are codes.
        langs = {line.strip() for line in out.stdout.splitlines() if line.strip()}
        # Drop the header line (contains "List of available languages")
        langs = {l for l in langs if not l.lower().startswith("list of")}
        _AVAILABLE_TESSERACT_LANGS_CACHE = frozenset(langs)
    except Exception as e:
        logger.warning("Could not enumerate Tesseract languages: %s", e)
        _AVAILABLE_TESSERACT_LANGS_CACHE = frozenset()
    return _AVAILABLE_TESSERACT_LANGS_CACHE


# Tesseract OSD script name → ordered list of Tesseract packs to try when
# the visual page script disagrees with the text layer. The first installed
# pack from each list is used.
_SCRIPT_TO_TESSERACT = {
    "Cyrillic": ["rus"],
    "Greek": ["ell"],
    "Fraktur": ["deu_latf", "deu"],
    "Arabic": ["ara"],
    "Hebrew": ["heb"],
    "Devanagari": ["hin"],
    "Japanese": ["jpn"],
    "Hangul": ["kor"],
    "Han": ["chi_sim", "chi_tra"],
    "Thai": ["tha"],
    "Bengali": ["ben"],
    "Tamil": ["tam"],
    "Telugu": ["tel"],
    "Kannada": ["kan"],
    "Malayalam": ["mal"],
    "Gujarati": ["guj"],
    "Gurmukhi": ["pan"],
}


# Vendor watermark / wrapper signatures. When a PDF's text layer contains
# nothing but one of these markers (and below 5K chars across the sample),
# treat it as a scanned PDF whose actual content is image-only — even
# though docling/PyMuPDF can read the boilerplate banner. See
# detect_scan_type's vendor cross-check.
_VENDOR_BOILERPLATE = (
    "ProQuest ebrary",
    "biodiversitylibrary.org",
    "This page intentionally left blank",
)


def _resolve_tesseract_packs(
    detected_iso: Optional[str],
    visual_script: Optional[str] = None,
) -> List[str]:
    """Return targeted Tesseract packs for ``(detected_iso, visual_script)``.

    Pure precedence resolution — no fallback union, no automatic ``eng``
    suffix. Returns ``[]`` when no targeted pack is installed for the inputs,
    which is the signal callers use to decide whether to warn or to fall
    back to ``CONFIG["ocr"]["ocr_languages_default"]``.

    Precedence (same as :func:`_compose_ocr_langs`):

    1. ``visual_script`` (from Tesseract OSD) wins over ``detected_iso`` —
       langdetect's code is unreliable when the text layer is corrupt.
    2. Otherwise use ``_ISO_TO_TESSERACT[detected_iso]``. For German the
       Fraktur pack is added when installed (historical papers).
    """
    available = _available_tesseract_langs()
    if not available:
        return []
    chosen: List[str] = []
    if visual_script and visual_script in _SCRIPT_TO_TESSERACT:
        for tess in _SCRIPT_TO_TESSERACT[visual_script]:
            if tess in available and tess not in chosen:
                chosen.append(tess)
    elif detected_iso:
        tess = _ISO_TO_TESSERACT.get(detected_iso)
        if tess and tess in available:
            chosen.append(tess)
            if tess == "deu" and "deu_latf" in available:
                chosen.append("deu_latf")
    return chosen


def _compose_ocr_langs(
    detected_iso: Optional[str],
    visual_script: Optional[str] = None,
) -> List[str]:
    """Pick Tesseract language codes to pass to ocrmypdf ``-l``.

    Tries :func:`_resolve_tesseract_packs` first; if that returns nothing
    (unknown language, or pack not installed), falls back to
    ``CONFIG["ocr"]["ocr_languages_default"]`` filtered to installed packs.
    ``eng`` is always appended.
    """
    available = _available_tesseract_langs()
    if not available:
        return []

    chosen = _resolve_tesseract_packs(detected_iso, visual_script)
    if chosen:
        if "eng" in available and "eng" not in chosen:
            chosen.append("eng")
        return chosen

    cfg_default = CONFIG.get("ocr", {}).get(
        "ocr_languages_default",
        ["eng", "deu", "deu_latf", "fra", "rus", "lat"],
    )
    for lang in cfg_default:
        if lang in available and lang not in chosen:
            chosen.append(lang)
    if "eng" in available and "eng" not in chosen:
        chosen.append("eng")
    return chosen



_MAX_VIZ_PAGES = 200
_MAX_VIZ_WIDTH = 1600  # px — downscale wider renders to limit memory


def create_cell_visualizations(
    input_pdf: Path,
    output_dir: Path,
    pdf_name: str,
    figure_bboxes_by_page: dict,
):
    """Create cell visualization PNGs using docling-parse with figure regions highlighted.

    ``figure_bboxes_by_page`` maps page number to a list of docling bbox
    objects.  Extracted by the caller *before* releasing the docling
    document, so memory can be reclaimed before the (heavier) per-page
    rendering loop.
    """
    import gc

    try:
        from docling_core.types.doc.page import TextCellUnit
        from docling_parse.pdf_parser import DoclingPdfParser, PdfDocument
        from PIL import Image, ImageDraw

        logger.info("Creating cell visualizations...")

        pdf_parser = DoclingPdfParser()
        pdf_doc: PdfDocument = pdf_parser.load(path_or_stream=str(input_pdf))

        # Collect pages into a list (iterate_pages may be a one-shot generator).
        pages = list(pdf_doc.iterate_pages())
        if len(pages) > _MAX_VIZ_PAGES:
            logger.warning(
                "Skipping cell visualizations for %s (%d pages > %d max)",
                pdf_name, len(pages), _MAX_VIZ_PAGES,
            )
            return

        cell_unit = TextCellUnit.WORD
        n_saved = 0

        for page_no, pred_page in pages:
            img = pred_page.render_as_image(cell_unit=cell_unit)

            # Downscale wide renders to cap memory usage.
            if img.width > _MAX_VIZ_WIDTH:
                ratio = _MAX_VIZ_WIDTH / img.width
                img = img.resize(
                    (int(img.width * ratio), int(img.height * ratio)),
                    Image.LANCZOS,
                )

            draw = ImageDraw.Draw(img)
            img_height = img.height

            # First, highlight figure regions in yellow
            if page_no in figure_bboxes_by_page:
                for bbox in figure_bboxes_by_page[page_no]:
                    x0 = bbox.l
                    x1 = bbox.r
                    y0_raw = bbox.t
                    y1_raw = bbox.b
                    y0 = img_height - y0_raw
                    y1 = img_height - y1_raw
                    x0, x1 = min(x0, x1), max(x0, x1)
                    y0, y1 = min(y0, y1), max(y0, y1)
                    draw.rectangle((x0, y0, x1, y1), fill="yellow", outline="orange", width=3)

            # Then, draw red rectangles around text cells (not in figures)
            for cell in pred_page.iterate_cells(unit_type=cell_unit):
                x0 = min(getattr(cell.rect, "r_x0", 0), getattr(cell.rect, "r_x2", 0))
                x1 = max(getattr(cell.rect, "r_x0", 0), getattr(cell.rect, "r_x2", 0))
                y0_raw = min(getattr(cell.rect, "r_y0", 0), getattr(cell.rect, "r_y2", 0))
                y1_raw = max(getattr(cell.rect, "r_y0", 0), getattr(cell.rect, "r_y2", 0))
                y0 = img_height - y1_raw
                y1 = img_height - y0_raw
                rect = (x0, y0, x1, y1)

                is_in_figure = False
                if page_no in figure_bboxes_by_page:
                    for bbox in figure_bboxes_by_page[page_no]:
                        fig_x0, fig_x1 = bbox.l, bbox.r
                        fig_y0_raw, fig_y1_raw = bbox.t, bbox.b
                        fig_y0 = img_height - fig_y0_raw
                        fig_y1 = img_height - fig_y1_raw
                        fig_x0, fig_x1 = min(fig_x0, fig_x1), max(fig_x0, fig_x1)
                        fig_y0, fig_y1 = min(fig_y0, fig_y1), max(fig_y0, fig_y1)
                        if (x0 < fig_x1 and x1 > fig_x0 and y0 < fig_y1 and y1 > fig_y0):
                            is_in_figure = True
                            break

                if not is_in_figure:
                    draw.rectangle(rect, outline="red", width=2)

            viz_path = output_dir / f"page_{page_no}_visualization.png"
            img.save(str(viz_path))
            # Release per-page image memory immediately.
            del draw, img
            gc.collect()
            n_saved += 1

        logger.info("Saved %d visualization PNGs", n_saved)

    except ImportError as e:
        logger.warning("Could not create visualizations — missing dependencies: %s", e)
    except Exception as e:
        logger.warning("Could not create visualizations: %s", e)



def run_pdf_processing_pipeline(
    pdf_path: Path,
    hash_dir: Path,
    temp_dir: Path,
    grobid_client: Optional[GrobidClient] = None,
    taxonomy_db: Optional[TaxonomyDB] = None,
    lexicons: Optional[Dict[str, Dict[str, Dict]]] = None,
    content_aware_figures: bool = False,
    vision_backend=None,
    bib_index: Optional[BibIndex] = None,
    resume: bool = False,
    taxonomy_fingerprint: Optional[Dict[str, Any]] = None,
    lexicon_fingerprints: Optional[Dict[str, Dict[str, Any]]] = None,
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
        "started_at": _utcnow_iso(),
        "processing_steps": [],   # legacy: stage names, success only
        "files_created": [],
        "errors": [],             # legacy free-text; kept for backwards compat
        "stage_timings": [],      # #34: per-stage timing, success or fail
        "stage_failures": [],     # #34: structured failure records
        "skipped_stages": [],     # #28: stages whose artifact was already present
    }

    try:
        # Huge-document gate (#35) — skip-and-flag PDFs above max_pages
        # before any expensive stage runs.
        with _stage(processing_summary, "huge_document_check", hash_dir=hash_dir):
            max_pages = int(CONFIG.get("huge_document", {}).get("max_pages", 500))
            n_pages = _pdf_page_count(temp_pdf)
            if n_pages is not None:
                processing_summary["page_count"] = n_pages
                if n_pages > max_pages:
                    raise _HugeDocumentError(
                        f"PDF has {n_pages} pages (max: {max_pages}); skipping. "
                        "Pre-split this paper or raise huge_document.max_pages "
                        "to process it."
                    )
            processing_summary["processing_steps"].append("huge_document_check")

        # ── scan_detection ──────────────────────────────────────────
        detection_file = hash_dir / "scan_detection.json"
        if _should_run_stage("scan_detection", hash_dir=hash_dir,
                             resume=resume, processing_summary=processing_summary):
            with _stage(processing_summary, "scan_detection", hash_dir=hash_dir):
                logger.info("Detecting scan type...")
                detection_result = detect_scan_type(temp_pdf)
                with open(detection_file, "w") as f:
                    json.dump(detection_result, f, indent=2)
                processing_summary["files_created"].append(str(detection_file))
                processing_summary["processing_steps"].append("scan_detection")
        else:
            with open(detection_file) as f:
                detection_result = json.load(f)

        # ── pdf_preparation ─────────────────────────────────────────
        processed_pdf = hash_dir / "processed.pdf"
        if _should_run_stage("pdf_preparation", hash_dir=hash_dir,
                             resume=resume, processing_summary=processing_summary):
            with _stage(processing_summary, "pdf_preparation", hash_dir=hash_dir):
                logger.info("Preparing PDF...")
                prepare_pdf(temp_pdf, detection_result, processed_pdf)
                processing_summary["files_created"].append(str(processed_pdf))
                processing_summary["processing_steps"].append("pdf_preparation")

        # ── docling_extraction ─────────────────────────────────────
        text_file = hash_dir / "text.json"
        figures_file = hash_dir / "figures.json"
        docling_doc_file = hash_dir / "docling_doc.json"
        if _should_run_stage("docling_extraction", hash_dir=hash_dir,
                             resume=resume, processing_summary=processing_summary):
            with _stage(processing_summary, "docling_extraction", hash_dir=hash_dir):
                logger.info("Extracting text and figures...")
                extract_docling_content(
                    processed_pdf,
                    text_file,
                    figures_file,
                    figures_dir,
                    visualizations_dir,
                    docling_doc_output=docling_doc_file,
                    scan_file_type=detection_result.get("file_type"),
                )
                processing_summary["files_created"].extend([str(text_file), str(figures_file)])
                if docling_doc_file.exists():
                    processing_summary["files_created"].append(str(docling_doc_file))
                processing_summary["processing_steps"].append("docling_extraction")

        # ── metadata_extraction ─────────────────────────────────────
        metadata_file = hash_dir / "metadata.json"
        references_file = hash_dir / "references.json"
        tei_file = hash_dir / "grobid.tei.xml"
        if _should_run_stage("metadata_extraction", hash_dir=hash_dir,
                             resume=resume, processing_summary=processing_summary):
            with _stage(processing_summary, "metadata_extraction", hash_dir=hash_dir):
                logger.info("Extracting metadata (Grobid)...")
                bib_entry = bib_index.lookup(pdf_path.name) if bib_index is not None else None
                extract_metadata(
                    processed_pdf,
                    metadata_file,
                    references_output=references_file,
                    tei_output=tei_file,
                    grobid_client=grobid_client,
                    original_filename=pdf_path.name,
                    bib_entry=bib_entry,
                )
                processing_summary["files_created"].extend(
                    [str(metadata_file), str(references_file)]
                )
                if tei_file.exists():
                    processing_summary["files_created"].append(str(tei_file))
                processing_summary["processing_steps"].append("metadata_extraction")

        # ── text_chunking ──────────────────────────────────────────
        chunks_file = hash_dir / "chunks.json"
        if _should_run_stage("text_chunking", hash_dir=hash_dir,
                             resume=resume, processing_summary=processing_summary):
            with _stage(processing_summary, "text_chunking", hash_dir=hash_dir):
                logger.info("Chunking text...")
                chunk_text(text_file, chunks_output=chunks_file)
                processing_summary["files_created"].append(str(chunks_file))
                processing_summary["processing_steps"].append("text_chunking")

        with _stage(processing_summary, "figure_pass25_annotation", hash_dir=hash_dir):
            logger.info("Pass 2.5: annotating figures from captions + text...")
            _pass25_annotate_figures(text_file, figures_file)
            processing_summary["processing_steps"].append("figure_pass25_annotation")

        if vision_backend is not None:
            with _stage(processing_summary, "figure_pass3b_rois", hash_dir=hash_dir):
                logger.info("Pass 3b: vision-model-driven panel + compound detection...")
                _pass3b_annotate_rois(figures_file, vision_backend)
                processing_summary["processing_steps"].append("figure_pass3b_rois")
        elif content_aware_figures:
            with _stage(processing_summary, "figure_pass3a_rois", hash_dir=hash_dir):
                logger.info("Pass 3a: OCR-driven panel ROI detection...")
                _pass3a_annotate_rois(figures_file)
                processing_summary["processing_steps"].append("figure_pass3a_rois")

        # Pass 3c — resolve *_compound figures: split their ROIs, match
        # to missing_figures, rename the PNG to range notation. Cheap;
        # worth running whenever 3a or 3b has been run.
        if vision_backend is not None or content_aware_figures:
            with _stage(processing_summary, "figure_pass3c_resolve", hash_dir=hash_dir):
                logger.info("Pass 3c: compound figure resolution + file rename...")
                summary_3c = resolve_compound_figures(figures_file)
                logger.info(
                    "Pass 3c: %d resolved, %d renamed, %d unchanged, %d new records",
                    summary_3c.get("resolved", 0),
                    summary_3c.get("renamed", 0),
                    summary_3c.get("unchanged", 0),
                    summary_3c.get("new_records", 0),
                )
                processing_summary["processing_steps"].append("figure_pass3c_resolve")

        with _stage(processing_summary, "figure_crossref", hash_dir=hash_dir):
            logger.info("Linking chunks to figures...")
            _crossref_chunks_and_figures(figures_file, chunks_file)
            processing_summary["processing_steps"].append("figure_crossref")

        # ── taxa_anatomy_extraction ─────────────────────────────────
        # Stage runs only when a taxonomy DB or at least one lexicon
        # category is configured. The input_fingerprint captures the
        # taxonomy + per-category content hashes, so editing one
        # lexicon section forces this stage to re-run on --resume.
        run_taxa_anat = taxonomy_db is not None or bool(lexicons)
        if run_taxa_anat:
            taxa_anat_fingerprint: Dict[str, Any] = {}
            if taxonomy_db is not None and taxonomy_fingerprint is not None:
                taxa_anat_fingerprint["taxonomy"] = taxonomy_fingerprint
            if lexicon_fingerprints:
                taxa_anat_fingerprint["lexicons"] = lexicon_fingerprints
            if _should_run_stage(
                "taxa_anatomy_extraction",
                hash_dir=hash_dir,
                resume=resume,
                processing_summary=processing_summary,
                expected_fingerprint=taxa_anat_fingerprint,
            ):
                with _stage(processing_summary, "taxa_anatomy_extraction",
                            hash_dir=hash_dir, input_fingerprint=taxa_anat_fingerprint):
                    logger.info("Extracting taxa + lexicon mentions...")
                    taxa_anat_files = _extract_taxa_and_anatomy(
                        chunks_file,
                        hash_dir,
                        taxonomy_db,
                        lexicons,
                        taxonomy_fingerprint=taxonomy_fingerprint,
                        lexicon_fingerprints=lexicon_fingerprints,
                    )
                    processing_summary["files_created"].extend(str(p) for p in taxa_anat_files)
                    if taxa_anat_files:
                        processing_summary["processing_steps"].append("taxa_anatomy_extraction")

        with _stage(processing_summary, "figures_report", hash_dir=hash_dir):
            logger.info("Generating figures report...")
            report_path = generate_figures_report(hash_dir)
            if report_path:
                processing_summary["files_created"].append(str(report_path))
                processing_summary["processing_steps"].append("figures_report")

        # Quality gates (#36) — informational silent-failure detectors.
        # Run after success so artifacts are populated. A failed gate
        # records a quality_flag in summary.json but does not fail the
        # paper; corpus_status.py (#40) rolls these up for review.
        with _stage(processing_summary, "quality_gates", hash_dir=hash_dir):
            qgs = _run_quality_gates(hash_dir)
            processing_summary["quality_flags"] = qgs
            if qgs:
                logger.warning(
                    "Quality gates flagged %d issue(s): %s",
                    len(qgs), ", ".join(g["gate"] for g in qgs),
                )
            processing_summary["processing_steps"].append("quality_gates")

        processing_summary["status"] = "success"

    except Exception as e:
        processing_summary["status"] = "error"
        processing_summary["errors"].append(str(e))
        logger.exception("Error processing PDF: %s", e)

    processing_summary["ended_at"] = _utcnow_iso()
    return processing_summary


def _annotate_pack_availability(result: Dict) -> Dict:
    """Tag a :func:`detect_scan_type` result with Tesseract pack availability.

    Adds two fields:

    - ``tesseract_packs`` — the targeted Tesseract pack list resolved from
      ``(detected_language, visual_script)`` via
      :func:`_resolve_tesseract_packs`. Empty when no pack is installed for
      the detected language/script (the pipeline will then fall back to the
      configured default union).
    - ``tesseract_pack_available`` — ``True``/``False``/``None``. ``None``
      means no language was detected (e.g., a low-text scan), so we can't
      decide. Recording this for born-digital papers too is intentional —
      it lets you grep ``scan_detection.json`` across the corpus to find
      papers in unsupported languages.

    Emits a WARNING when ``needs_ocr=True`` but no targeted pack is
    installed; OCR will still run with the default union, but accuracy on
    the detected language/script will be poor. Install the appropriate
    pack (e.g. ``brew install tesseract-lang`` on macOS or
    ``apt-get install tesseract-ocr-<code>`` on Debian) for better results.
    """
    iso = result.get("detected_language")
    visual = result.get("visual_script")
    packs = _resolve_tesseract_packs(iso, visual)
    result["tesseract_packs"] = packs
    if iso is None and visual is None:
        result["tesseract_pack_available"] = None
    else:
        result["tesseract_pack_available"] = bool(packs)
    if result.get("needs_ocr") and not packs and (iso or visual):
        logger.warning(
            "No installed Tesseract pack for detected_language=%r visual_script=%r "
            "(%s); OCR will fall back to the default language union and may be "
            "low quality.",
            iso, visual, result.get("filename"),
        )
    return result


def detect_scan_type(pdf_path: Path) -> Dict:
    """Classify the text-layer state of a PDF into one of three buckets:

    - ``born_digital`` — dense, intelligible text layer in a detectable
      language. No OCR needed.
    - ``scanned`` — little or no text layer. Needs OCR (``--skip-text``
      default; pages with partial text are left alone).
    - ``broken_text_layer`` — dense text present but nonsense (e.g.,
      Cyrillic characters mapped 1:1 to Latin-1 bytes, as seen in the
      Stepanjants 1970 scan). Needs OCR with ``--force-ocr`` to replace
      the garbage.

    The returned dict is written to ``scan_detection.json`` and consumed
    by :func:`prepare_pdf`. Downstream stages can also read the
    ``detected_language`` field as a useful signal.
    """
    try:
        import fitz  # PyMuPDF
    except ImportError:
        logger.warning("PyMuPDF not available; treating as born-digital (no OCR)")
        return _annotate_pack_availability({
            "filename": pdf_path.name,
            "file_type": "born_digital",
            "has_text": True,
            "needs_ocr": False,
            "detected_language": None,
            "language_confidence": 0.0,
            "gibberish_score": 0.0,
            "ocr_mode": None,
        })

    # Read up to the first 5 pages (or all pages if the doc is shorter) —
    # enough text for language detection and gibberish scoring without
    # paying full-doc parsing cost just to triage.
    try:
        doc = fitz.open(pdf_path)
        pages_to_check = min(5, len(doc))
        total_text = ""
        for page_num in range(pages_to_check):
            total_text += doc[page_num].get_text()
        doc.close()
    except Exception as e:
        logger.warning("Error reading %s: %s; assuming scanned", pdf_path.name, e)
        return _annotate_pack_availability({
            "filename": pdf_path.name,
            "file_type": "scanned",
            "has_text": False,
            "needs_ocr": True,
            "detected_language": None,
            "language_confidence": 0.0,
            "gibberish_score": 0.0,
            "ocr_mode": "skip_text",
        })

    stripped = total_text.strip()
    total_chars = len(stripped)
    avg_chars_per_page = total_chars / max(pages_to_check, 1)
    has_substantial_text = total_chars > 500
    has_good_density = avg_chars_per_page > 100

    # Low-text path: treat as scanned and OCR it.
    if not (has_substantial_text and has_good_density):
        return _annotate_pack_availability({
            "filename": pdf_path.name,
            "file_type": "scanned",
            "has_text": False,
            "needs_ocr": True,
            "detected_language": None,
            "language_confidence": 0.0,
            "gibberish_score": 0.0,
            "total_chars_sampled": total_chars,
            "pages_checked": pages_to_check,
            "ocr_mode": "skip_text",
        })

    # Has enough text — now triage born_digital vs broken_text_layer.
    lang, conf = _detect_language(total_text)
    gib = _gibberish_score(total_text)
    scripts = _text_layer_scripts(total_text)
    latin_frac = scripts.get("Latin", 0.0)

    # --- Vendor-watermark cross-check ---
    # Some PDF exports (ProQuest eBooks, BHL scan exports) stamp a thin
    # per-page boilerplate banner on top of image-only content. Their text
    # layer clears the volume gate above but contains nothing but the
    # banner — the actual book content is in raster images. Re-route to
    # the OCR path so prepare_pdf can recover content from the images.
    if total_chars < 5000 and any(m in total_text for m in _VENDOR_BOILERPLATE):
        matched = next(m for m in _VENDOR_BOILERPLATE if m in total_text)
        logger.info(
            "Vendor-boilerplate-only text layer detected (%r) in %s; "
            "re-routing to OCR",
            matched, pdf_path.name,
        )
        return _annotate_pack_availability({
            "filename": pdf_path.name,
            "file_type": "scanned",
            "detection_reason": "vendor_boilerplate_only",
            "has_text": True,
            "needs_ocr": True,
            "detected_language": lang,
            "language_confidence": conf,
            "gibberish_score": gib,
            "text_layer_scripts": scripts,
            "visual_script": None,
            "total_chars_sampled": total_chars,
            "pages_checked": pages_to_check,
            "ocr_mode": "skip_text",
            "vendor_marker": matched,
        })

    threshold = float(CONFIG.get("ocr", {}).get("gibberish_threshold", 0.5))

    # --- Cheap gibberish path ---
    # High gibberish score is a direct signal the text layer is corrupt,
    # regardless of script. Catches the obvious cases.
    if gib > threshold:
        return _annotate_pack_availability({
            "filename": pdf_path.name,
            "file_type": "broken_text_layer",
            "detection_reason": "gibberish_score_above_threshold",
            "has_text": True,
            "needs_ocr": True,
            "detected_language": lang,
            "language_confidence": conf,
            "gibberish_score": gib,
            "text_layer_scripts": scripts,
            "visual_script": None,
            "total_chars_sampled": total_chars,
            "pages_checked": pages_to_check,
            "ocr_mode": "force_ocr",
        })

    # --- Visual-vs-text cross-check ---
    # If the text layer is almost entirely Latin-family characters *and*
    # gibberish is non-trivial, the Stepanjants case is possible: Cyrillic
    # glyphs mapped 1:1 to Latin-1 bytes. Confirm by running Tesseract
    # OSD on a rendered page and comparing the visual script to what the
    # text layer claims. Only invoked in the suspect zone so OSD cost is
    # bounded at corpus scale.
    visual = None
    if latin_frac > 0.90 and gib > 0.40:
        visual = _visual_page_script(pdf_path)
        # The expected script from langdetect — if it disagrees with what
        # OSD sees on the actual page image, the text layer is corrupt.
        expected = _LANG_EXPECTED_SCRIPT.get(lang or "", "Latin")
        # Latin/Fraktur are compatible (Fraktur is a Latin-family display
        # script for German). Any other mismatch is a red flag.
        compatible = {"Latin", "Fraktur", None}
        if (
            visual is not None
            and visual not in compatible
            and expected in compatible  # text layer "claims" Latin
        ):
            logger.info(
                "Text-layer/visual script mismatch: layer=Latin visual=%s — "
                "flagging as broken_text_layer",
                visual,
            )
            return _annotate_pack_availability({
                "filename": pdf_path.name,
                "file_type": "broken_text_layer",
                "detection_reason": "visual_script_mismatch",
                "has_text": True,
                "needs_ocr": True,
                "detected_language": lang,
                "language_confidence": conf,
                "gibberish_score": gib,
                "text_layer_scripts": scripts,
                "visual_script": visual,
                "total_chars_sampled": total_chars,
                "pages_checked": pages_to_check,
                "ocr_mode": "force_ocr",
            })

    return _annotate_pack_availability({
        "filename": pdf_path.name,
        "file_type": "born_digital",
        "detection_reason": "clean_text_layer",
        "has_text": True,
        "needs_ocr": False,
        "detected_language": lang,
        "language_confidence": conf,
        "gibberish_score": gib,
        "text_layer_scripts": scripts,
        "visual_script": visual,
        "total_chars_sampled": total_chars,
        "pages_checked": pages_to_check,
        "ocr_mode": None,
    })


def prepare_pdf(input_pdf: Path, detection_result: Dict, output_pdf: Path):
    """Run OCR if the detection result calls for it, else copy straight through.

    OCR behavior is driven by ``detection_result["ocr_mode"]``:

    - ``"skip_text"`` (default for scanned docs) — ``ocrmypdf --skip-text``
      so pages that already have a text layer are left untouched and only
      the blank pages get OCR'd. This is correct for mixed documents.
    - ``"force_ocr"`` (for broken text layers) — ``ocrmypdf --force-ocr`` to
      discard the garbage text layer and re-OCR everything. Used when
      :func:`detect_scan_type` flags a document as ``broken_text_layer``.

    Language selection uses ``detection_result["detected_language"]`` plus
    the config default union, filtered to what Tesseract actually has
    installed. See :func:`_compose_ocr_langs`.

    On any OCR failure we copy the original PDF through and log the error
    — downstream stages should still run on the untouched PDF rather than
    the pipeline aborting.
    """
    if not detection_result.get("needs_ocr"):
        logger.info("Copying %s (detected as %s)",
                    input_pdf.name, detection_result.get("file_type"))
        shutil.copy2(input_pdf, output_pdf)
        return

    if shutil.which("ocrmypdf") is None:
        logger.warning("ocrmypdf not found on PATH, copying original PDF")
        shutil.copy2(input_pdf, output_pdf)
        return

    ocr_mode = detection_result.get("ocr_mode", "skip_text")
    mode_flag = "--force-ocr" if ocr_mode == "force_ocr" else "--skip-text"

    langs = _compose_ocr_langs(
        detection_result.get("detected_language"),
        visual_script=detection_result.get("visual_script"),
    )
    if not langs:
        logger.warning(
            "No Tesseract languages available; copying original PDF (OCR skipped)"
        )
        shutil.copy2(input_pdf, output_pdf)
        return
    lang_arg = "+".join(langs)

    # Auto-degrade --optimize when pngquant isn't installed. ocrmypdf
    # requires pngquant for levels 2 and 3 and will exit 3 without it. We
    # still want OCR to proceed, just with smaller gains — drop to 1
    # (safe, internal-only optimizations) and log a nudge to install it.
    optimize_level = int(CONFIG.get("ocr", {}).get("optimize_level", 2))
    if optimize_level >= 2 and shutil.which("pngquant") is None:
        logger.warning(
            "pngquant not on PATH; reducing ocrmypdf --optimize from %d to 1. "
            "For tighter PDF output install pngquant (macOS: `brew install pngquant`).",
            optimize_level,
        )
        optimize_level = 1
    optimize_level = str(optimize_level)

    logger.info(
        "Running OCR on %s | file_type=%s mode=%s langs=%s",
        input_pdf.name,
        detection_result.get("file_type"),
        mode_flag,
        lang_arg,
    )
    cmd = [
        "ocrmypdf",
        mode_flag,
        "-l", lang_arg,
        "--optimize", optimize_level,
        "--color-conversion-strategy", "RGB",
        "--output-type", "pdf",
        str(input_pdf),
        str(output_pdf),
    ]

    ocr_timeout = float(CONFIG.get("stage_timeouts", {}).get("ocr", 1800))
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=ocr_timeout,
        )
    except subprocess.TimeoutExpired:
        # Re-raise so the wrapping _stage records reason_code=timeout.
        # The pipeline-level try/except still catches and continues to
        # the next paper, with the timeout structured into stage_failures[].
        logger.warning("OCR timed out after %.0fs on %s", ocr_timeout, input_pdf.name)
        raise

    if result.returncode != 0:
        logger.warning(
            "OCR failed (exit %d) on %s, copying original PDF",
            result.returncode, input_pdf.name,
        )
        # ocrmypdf writes human-readable errors to stderr; log a truncated
        # version so one bad doc doesn't flood the pipeline log.
        if result.stderr:
            logger.warning("ocrmypdf stderr (head): %s", result.stderr[:500])
        shutil.copy2(input_pdf, output_pdf)
    else:
        logger.info("OCR completed successfully (langs=%s)", lang_arg)
