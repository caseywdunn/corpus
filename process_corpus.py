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
from bib_metadata import BibIndex, bib_entry_to_metadata
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
    load_anatomy_lexicon,
    extract_taxon_mentions,
    extract_anatomy_mentions,
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
        "optimize_level": 2,
        # Used when language detection fails or the doc has a broken
        # text layer. Tesseract combines languages gracefully though
        # slowly; override in config.yaml to narrow for speed.
        "ocr_languages_default": ["eng", "deu", "deu_latf", "fra", "rus", "lat"],
        # Gibberish-score threshold for flagging a broken text layer.
        # 0.5 reliably catches the Cyrillic-as-Latin-1 case without
        # false-positiving on reference-heavy papers.
        "gibberish_threshold": 0.5,
    },
    "chunking": {
        "max_tokens": 8191,
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
    taxonomy_db: Optional[TaxonomyDB] = None,
    anatomy_lexicon: Optional[Dict[str, Dict]] = None,
    content_aware_figures: bool = False,
    vision_backend=None,
    bib_index: Optional[BibIndex] = None,
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
            scan_file_type=detection_result.get("file_type"),
        )
        processing_summary["files_created"].extend([str(text_file), str(figures_file)])
        if docling_doc_file.exists():
            processing_summary["files_created"].append(str(docling_doc_file))
        processing_summary["processing_steps"].append("docling_extraction")

        logger.info("Extracting metadata (Grobid)...")
        metadata_file = hash_dir / "metadata.json"
        references_file = hash_dir / "references.json"
        tei_file = hash_dir / "grobid.tei.xml"
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

        logger.info("Chunking text...")
        chunks_file = hash_dir / "chunks.json"
        chunk_text(text_file, metadata_file, chunks_file)
        processing_summary["files_created"].append(str(chunks_file))
        processing_summary["processing_steps"].append("text_chunking")

        logger.info("Pass 2.5: annotating figures from captions + text...")
        _pass25_annotate_figures(text_file, figures_file)
        processing_summary["processing_steps"].append("figure_pass25_annotation")

        if vision_backend is not None:
            logger.info("Pass 3b: vision-model-driven panel + compound detection...")
            _pass3b_annotate_rois(figures_file, vision_backend)
            processing_summary["processing_steps"].append("figure_pass3b_rois")
        elif content_aware_figures:
            logger.info("Pass 3a: OCR-driven panel ROI detection...")
            _pass3a_annotate_rois(figures_file)
            processing_summary["processing_steps"].append("figure_pass3a_rois")

        # Pass 3c — resolve *_compound figures: split their ROIs, match
        # to missing_figures, rename the PNG to range notation. Cheap;
        # worth running whenever 3a or 3b has been run.
        if vision_backend is not None or content_aware_figures:
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

        logger.info("Linking chunks to figures...")
        _crossref_chunks_and_figures(figures_file, chunks_file)
        processing_summary["processing_steps"].append("figure_crossref")

        logger.info("Extracting taxa + anatomy mentions...")
        taxa_anat_files = _extract_taxa_and_anatomy(
            chunks_file, hash_dir, taxonomy_db, anatomy_lexicon
        )
        processing_summary["files_created"].extend(str(p) for p in taxa_anat_files)
        if taxa_anat_files:
            processing_summary["processing_steps"].append("taxa_anatomy_extraction")

        logger.info("Generating figures report...")
        report_path = generate_figures_report(hash_dir)
        if report_path:
            processing_summary["files_created"].append(str(report_path))
            processing_summary["processing_steps"].append("figures_report")

        processing_summary["status"] = "success"

    except Exception as e:
        processing_summary["status"] = "error"
        processing_summary["errors"].append(str(e))
        logger.exception("Error processing PDF: %s", e)

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

    result = subprocess.run(cmd, capture_output=True, text=True)

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


# Years plausible for a siphonophore paper: 17xx (earliest Linnaean works) to
# the current decade. Tighter than a generic `\d{4}` so e.g. a trailing
# "20" in a specimen ID isn't picked up. Lookarounds exclude matches that
# are part of a longer digit run (e.g., "12005" shouldn't match 2005).
_FILENAME_YEAR_RE = re.compile(r"(?<!\d)(17\d{2}|18\d{2}|19\d{2}|20[0-4]\d)(?!\d)")


def _year_from_filename(name: str) -> Optional[int]:
    """Best-effort pub-year extraction from a PDF filename.

    Matches the first 4-digit year in the range 1700–2049, which covers the
    historical siphonophore literature without false-positiving on generic
    numbers (ISSNs, specimen counts, etc.). Returns None if no match. The
    Pugh-curated library uses an "Author(s)YYYY" convention almost
    universally, so this recovers years for the majority of the corpus where
    Grobid's header parser emits an empty ``<date/>``.
    """
    if not name:
        return None
    m = _FILENAME_YEAR_RE.search(name)
    if not m:
        return None
    return int(m.group(1))


def _write_placeholder_metadata(
    pdf_path: Path,
    metadata_output: Path,
    references_output: Path,
    original_filename: Optional[str] = None,
):
    """Write empty-but-valid metadata.json and references.json.

    Used when Grobid is unavailable or its TEI can't be parsed. Keeps
    downstream stages (chunking, embedding) functional and lets the
    summary show this document as metadata-less so we can triage later.
    Even in the placeholder path we attempt to recover a year from the
    filename — it's all the signal we have without Grobid.
    """
    effective_filename = original_filename or pdf_path.name
    fname_year = _year_from_filename(effective_filename)
    placeholder = {
        "filename": effective_filename,
        "title": "",
        "authors": [],
        "year": fname_year,
        "year_source": "filename" if fname_year is not None else None,
        "journal": "",
        "doi": "",
        "abstract": "",
        "extraction_method": "placeholder",
    }
    with open(metadata_output, "w", encoding="utf-8") as f:
        json.dump(placeholder, f, indent=2, ensure_ascii=False)
    with open(references_output, "w", encoding="utf-8") as f:
        json.dump({"references": [], "total_references": 0}, f, indent=2)
    with open(metadata_output.parent / "intext_citations.json", "w", encoding="utf-8") as f:
        json.dump({"paragraphs": [], "citations": []}, f, indent=2)


def extract_metadata(
    pdf_path: Path,
    metadata_output: Path,
    references_output: Optional[Path] = None,
    tei_output: Optional[Path] = None,
    grobid_client: Optional[GrobidClient] = None,
    original_filename: Optional[str] = None,
    bib_entry: Optional[Dict] = None,
):
    """Extract bibliographic metadata + references via Grobid.

    Runs ``/api/processFulltextDocument`` once, caches the raw TEI-XML at
    ``tei_output`` (if given), and parses it into ``metadata.json`` plus
    ``references.json``. If Grobid is unreachable or any step fails, falls
    back to placeholder files so downstream stages still run — the caller
    can inspect ``metadata["extraction_method"]`` to tell which path was
    taken.

    When ``bib_entry`` is supplied (a parsed BibTeX record matched to this
    PDF by filename), its title/authors/year/journal/DOI override the
    Grobid header — Grobid is still called for the reference list so we
    can build the citation graph, but the header parse is skipped. This
    lets a curated bibliography be the source of truth for header fields
    while keeping Grobid's reference extraction.

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
    bib_entry:
        Parsed BibTeX entry from :class:`bib_metadata.BibIndex`. When
        present, overrides Grobid's header parse for this document.
    """
    hash_dir = metadata_output.parent
    if references_output is None:
        references_output = hash_dir / "references.json"
    if tei_output is None:
        tei_output = hash_dir / "grobid.tei.xml"
    # Prefer the original PDF filename for both provenance and year-fallback
    # — by the time this runs, pdf_path is usually "processed.pdf" (post-OCR
    # or copied), whose name carries no year information.
    effective_filename = original_filename or pdf_path.name

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

    # 3. Parse references from TEI if available — independent of which
    # source we use for the header. Anything that fails leaves us with
    # an empty reference list (written below as a fallback).
    refs: List[dict] = []
    refs_parsed = False
    intext = {"paragraphs": [], "citations": []}
    if tei_xml is not None:
        try:
            refs = parse_tei_references(tei_xml)
            refs_parsed = True
        except Exception as e:
            logger.warning("Failed to parse Grobid TEI references (%s)", e)
        # In-text citation graph (issue #7).  Independent of refs parse —
        # we want partial recovery when one fails and the other doesn't.
        try:
            intext = parse_tei_intext_citations(tei_xml)
        except Exception as e:
            logger.warning("Failed to parse Grobid TEI in-text citations (%s)", e)

    # 4. Build the header. Bib record wins if present; otherwise Grobid's
    # header parse, with filename-year fallback. If both fail, placeholder.
    header: Optional[Dict] = None
    if bib_entry is not None:
        header = bib_entry_to_metadata(bib_entry, effective_filename)
        logger.info(
            "Using bib entry %r for header (overrides Grobid)",
            header.get("bib_key"),
        )
    elif tei_xml is not None:
        try:
            header = parse_tei_header(tei_xml)
            header["filename"] = effective_filename
            # Year fallback: Grobid emits <date/> (empty) on many papers
            # whose title page doesn't match its header template — especially
            # historical and non-English. The Pugh library's "AuthorYYYY"
            # filename convention lets us recover the year cheaply. Record
            # which source won in ``year_source`` so downstream consumers can
            # treat Grobid-extracted years (stronger) differently from
            # filename-derived years (weaker, depends on filename convention).
            if header.get("year") is None:
                fname_year = _year_from_filename(effective_filename)
                if fname_year is not None:
                    header["year"] = fname_year
                    header["year_source"] = "filename"
                    logger.info(
                        "Grobid returned year=None; recovered %d from filename %r",
                        fname_year, effective_filename,
                    )
                else:
                    header["year_source"] = None
            else:
                header["year_source"] = "grobid"
        except Exception as e:
            logger.warning(
                "Failed to parse Grobid TEI header (%s); writing placeholder", e
            )
            header = None

    # 5. Write outputs. If we got a header, write it + the reference list
    # (which may be empty if Grobid was unavailable). Otherwise placeholder.
    if header is not None:
        with open(metadata_output, "w", encoding="utf-8") as f:
            json.dump(header, f, indent=2, ensure_ascii=False)
        with open(references_output, "w", encoding="utf-8") as f:
            json.dump(
                {"references": refs, "total_references": len(refs)},
                f,
                indent=2,
                ensure_ascii=False,
            )
        with open(hash_dir / "intext_citations.json", "w", encoding="utf-8") as f:
            json.dump(intext, f, indent=2, ensure_ascii=False)
        logger.info(
            "Wrote metadata (method=%s, %d authors), %d references (parsed=%s), %d in-text citations",
            header.get("extraction_method"),
            len(header.get("authors", [])),
            len(refs),
            refs_parsed,
            len(intext["citations"]),
        )
        return

    # 6. Placeholder path.
    _write_placeholder_metadata(
        pdf_path, metadata_output, references_output,
        original_filename=effective_filename,
    )


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
        json.dump(figures_data, f, indent=2, ensure_ascii=False)

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

    OCR reliability on line-art scientific figures varies a lot — PLAN.md
    §9 calls out vision-LLM fallback as Pass 3b for cases Pass 3a
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
        json.dump(data, f, indent=2, ensure_ascii=False)
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
        json.dump(data, f, indent=2, ensure_ascii=False)
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
        json.dump(figures_data, f, indent=2, ensure_ascii=False)
    with chunks_file.open("w", encoding="utf-8") as f:
        json.dump(chunks_data, f, indent=2, ensure_ascii=False)

    # Log a small summary so pipeline.log captures the link density.
    linked_figs = sum(1 for fig in figures if fig.get("referenced_in_chunks"))
    linked_chunks = sum(1 for ch in chunks if ch.get("figure_refs"))
    logger.info(
        "Cross-ref: %d/%d figures referenced, %d/%d chunks cite a figure",
        linked_figs, len(figures), linked_chunks, len(chunks),
    )


def _extract_taxa_and_anatomy(
    chunks_file: Path,
    hash_dir: Path,
    taxonomy_db: Optional[TaxonomyDB],
    anatomy_lexicon: Optional[Dict[str, Dict]],
) -> List[Path]:
    """Run taxon + anatomy extraction on the per-PDF chunks and write
    ``taxa.json`` / ``anatomy.json``.

    Designed to degrade gracefully: if the taxonomy snapshot or anatomy
    lexicon isn't configured / available, the corresponding artifact is
    simply not emitted (not an error). This keeps the pipeline runnable
    on a dev machine that hasn't yet run ``ingest_taxonomy.py``.
    """
    out: List[Path] = []
    if not chunks_file.exists():
        return out
    try:
        chunks_data = json.load(chunks_file.open(encoding="utf-8"))
    except Exception as e:
        logger.warning("Skipping taxa/anatomy pass: couldn't read %s: %s", chunks_file, e)
        return out
    chunks = chunks_data.get("chunks") or []

    if taxonomy_db is not None:
        taxa_res = extract_taxon_mentions(chunks, taxonomy_db)
        taxa_file = hash_dir / "taxa.json"
        with taxa_file.open("w", encoding="utf-8") as f:
            json.dump(taxa_res, f, indent=2, ensure_ascii=False)
        out.append(taxa_file)
        logger.info(
            "Taxon mentions: %d (unique taxa: %d)",
            taxa_res["total_mentions"], taxa_res["unique_taxa"],
        )
    else:
        logger.info("Skipping taxon extraction (no taxonomy DB configured)")

    if anatomy_lexicon:
        anat_res = extract_anatomy_mentions(chunks, anatomy_lexicon)
        anat_file = hash_dir / "anatomy.json"
        with anat_file.open("w", encoding="utf-8") as f:
            json.dump(anat_res, f, indent=2, ensure_ascii=False)
        out.append(anat_file)
        logger.info(
            "Anatomy mentions: %d (unique terms: %d)",
            anat_res["total_mentions"], anat_res["unique_terms"],
        )
    else:
        logger.info("Skipping anatomy extraction (no lexicon configured)")
    return out


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
    if not args.no_taxa:
        taxonomy_path = args.taxonomy_db or (args.output_dir / "taxonomy.sqlite")
        if taxonomy_path.exists():
            try:
                taxonomy_db = TaxonomyDB(taxonomy_path)
                logger.info(
                    "Taxonomy snapshot loaded from %s (%d names)",
                    taxonomy_path, len(taxonomy_db.name_set()),
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
                    anatomy_lexicon = load_anatomy_lexicon(args.anatomy_lexicon)
                    logger.info(
                        "Anatomy lexicon loaded from %s (%d terms)",
                        args.anatomy_lexicon, len(anatomy_lexicon),
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
    # When --resume is active and we're batching, remove already-completed
    # hashes *before* dividing into batches. This ensures batches contain
    # only unprocessed documents, instead of some batches being entirely
    # skip-only while others carry all the remaining work.
    if args.resume and args.batch_index is not None:
        before = len(pdf_map)
        pdf_map = {
            h: paths for h, paths in pdf_map.items()
            if not (documents_dir / short_hash(h) / "summary.json").exists()
        }
        skipped = before - len(pdf_map)
        if skipped:
            logger.info(
                "Resume: filtered out %d already-completed documents "
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

            # Run each document in a child process so that C-level crashes
            # (segfaults in docling/PyMuPDF) don't kill the whole batch.
            # On Linux the default fork start method inherits all objects
            # (grobid_client, taxonomy_db, vision_backend) without pickling.
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

            proc = multiprocessing.Process(target=_worker)
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
