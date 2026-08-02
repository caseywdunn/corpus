"""Pipeline configuration + section-class normalizer.

Holds the built-in ``_DEFAULT_CONFIG`` and the loader that merges
``config.yaml`` over it. ``CONFIG`` is the live mutable singleton —
populated by :func:`pipeline.main.main` from the user's --config flag,
and also safe-loaded with defaults on import so helpers can read it
outside main() (tests, ad-hoc scripts).

Other pipeline modules read CONFIG via ``from .config import CONFIG``
or ``from . import config; config.CONFIG.get(...)``. The latter sees
post-load updates; the former captures the bound dict at import time.
Both work because CONFIG is mutated in place via :func:`load_config`'s
return value being assigned by main().
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Section-class normalizer
# ---------------------------------------------------------------------------
#
# Ordered list — first match wins. Patterns are case-insensitive and cover
# English, German, French, Russian, plus the taxonomy-paper conventions
# common in systematic literature (Description / Systematics).
_SECTION_PATTERNS = [
    ("abstract", r"\babstract\b|\bsummary\b|\bzusammenfassung\b|\br[ée]sum[ée]\b|резюме|сводка"),
    ("introduction", r"\bintroduction\b|\bвведение\b|\beinleitung\b"),
    # Before `methods`: a heading like "CONCLUSIONS ON THE METHODS OF
    # DEVELOPMENT" names its conclusions and only mentions methods, and
    # first-match-wins was labelling it `methods`.
    ("conclusion", r"\bconclusions?\b|выводы|\bschlussfolgerung"),
    # Russian headings in this literature run to "МАТЕРИАЛ И МЕТОДИКА
    # СБОРА И ОБРАБОТКИ" — singular материал, методика rather than
    # методы — which the stricter phrase missed.
    ("methods", r"materials?\s+and\s+methods|\bmethods?\b|\bmethodology\b|study\s+area|материал\w*\s+и\s+метод\w+|m[ée]thodes"),
    ("results", r"\bresults?\b|\bр[еé]зультаты\b|\bergebnisse\b|\br[ée]sultats\b"),
    # Morphology/anatomy headings are how descriptive taxonomic papers
    # label what IMRaD would call results — "GROSS MORPHOLOGY",
    # "Морфология дефинитивных колоний", "Bau der Siphonophoren".
    ("description", r"\bdescription\b|\bsystematics?\b|\btaxonomy\b|\bsystematic\s+account\b|\bdiagnosis\b|описание|\bmorpholog\w*|морфолог\w*|\banatom\w*|анатом\w*|\bbau\s+der\b"),
    ("discussion", r"\bdiscussion\b|обсуждение|\bdiskussion\b"),
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
# Config defaults + loader
# ---------------------------------------------------------------------------

# Defaults used if config.yaml is missing or a key is absent.
_DEFAULT_CONFIG = {
    "ocr": {
        "optimize_level": 2,
        # Used when language detection fails or the doc has a broken
        # text layer. Tesseract combines languages gracefully though
        # slowly; override in config.yaml to narrow for speed.
        "ocr_languages_default": [
            "eng",
            "deu",
            "deu_latf",
            "fra",
            "rus",
            "lat",
            # #46 — broader taxonomic-literature coverage
            "spa",
            "por",
            "chi_sim",
            "chi_tra",
            "jpn",
            "ell",
            "kor",
            # `grc` (Ancient Greek) is installed by environment.yaml but
            # left out of the default fallback union: low signal for
            # living taxonomy + adds noise to multi-pack OCR runs.
        ],
        # Gibberish-score threshold for flagging a broken text layer.
        # Raised from 0.5 → 0.65 to avoid the MilosMaley2005 false
        # positive (langdetect misfire as Swahili). Documents with
        # truly corrupt text (Cyrillic mojibake etc.) are still caught
        # by the visual_script_mismatch path before reaching this
        # threshold. Must stay in sync with config_schema.OcrConfig
        # default + the bundled config.template.yaml.
        "gibberish_threshold": 0.65,
        # Per-page --tesseract-timeout for ocrmypdf. A timeout blanks the
        # page and still exits 0, so this is deliberately generous.
        "tesseract_page_timeout": 900,
        # Re-OCR scanned pages that already carry a text layer of unknown
        # provenance. The gibberish/visual-script paths above only judge
        # what the text layer *says*; they cannot tell a typeset page
        # from a scan someone else OCR'd badly, which is most of a
        # historical corpus. `scan_page_fraction_min` is the fraction of
        # sampled pages that must be a single full-page image before the
        # document is treated as a scan — see scan._scanned_page_fraction.
        # Must stay in sync with config_schema.OcrConfig + the bundled
        # config.template.yaml.
        "reocr_scanned_text_layers": True,
        "scan_page_fraction_min": 0.40,
        # Identify a scan's language by OCRing a small sample, rather
        # than trusting a layer we've already rejected or falling back to
        # every installed pack at once. See scan._probe_language_by_ocr.
        "probe_language_by_ocr": True,
        # langdetect has no Latin profile, so Latin reads as a
        # low-confidence Romance language; below this floor we OCR with
        # the script union (which has `lat`) rather than a wrong pack.
        "probe_language_min_confidence": 0.85,
        # Bilingual original-plus-translation volumes are routine here,
        # so the probe may report more than one language per document.
        "probe_max_languages": 3,
        # Probe pages whose OCR output is this gibberish were read with
        # the wrong script's packs; discard their language vote.
        "probe_max_gibberish": 0.50,
        # 300, not 200 — see config.template.yaml. Lowering it broke
        # language detection on the hardest documents.
        "probe_dpi": 300,
        # Sampled from the middle 15-85% of pages — front matter and
        # references are both poor language-detection input.
        "probe_sample_pages": 5,
    },
    "chunking": {
        "max_tokens": 8191,
    },
    # Figure rasterization scale (#121). Saved figure DPI = 72 * scale,
    # fixed at extraction time. Must stay in sync with
    # config_schema.FiguresConfig.images_scale + config.template.yaml.
    # (figures.panel_detection lives in the pydantic FiguresConfig; the
    # pipeline reads that via the CLI flow, not this CONFIG global.)
    "figures": {
        # 'native' (default) renders each figure at its source's native
        # pixel density (#121); vector figures fall back to vector_dpi.
        # 'fixed' uses images_scale for every figure instead.
        "resolution_mode": "native",
        "vector_dpi": 300.0,
        "max_dpi": None,
        "images_scale": 2.0,
    },
    # Per-stage wallclock caps in seconds (#34). Hard-enforced where
    # technically feasible (subprocess.run on ocrmypdf; GrobidClient
    # request timeout). Stages without a hard-cap mechanism still record
    # timing in summary.json; the value is informational.
    "stage_timeouts": {
        "ocr": 1800,         # 30 min — long monographs can run scary long
        "grobid": 300,       # 5 min — enforced by GrobidClient default
        "docling": 600,      # 10 min — informational only (in-process)
        "vision_pass": 600,  # 10 min — informational only
    },
    # Huge-document gate (#35). PDFs above max_pages skip the pipeline
    # with a structured too_large failure rather than burning hours of
    # OCR / docling on something that may not even fit in memory. The
    # Haeckel 1888 *Challenger Siphonophorae* report (~600 pages of
    # plates) is the canary case. Chunked-OCR is an open follow-up;
    # the v0.2 implementation is skip-and-flag. Default raised from
    # 500 → 5000 in v0.3 — the 500 default was a v0.2 placeholder
    # that excluded too many real long-form monographs. Must stay in
    # sync with config_schema.HugeDocumentConfig + the bundled
    # config.template.yaml.
    "huge_document": {
        "max_pages": 5000,
    },
    # Silent-failure quality gates (#36). Each threshold is conservative
    # by default; tighten in config.yaml for stricter corpora. A gate
    # produces a flag in summary.json["quality_flags"]; nothing is
    # rejected automatically — operators decide what to do with flagged
    # papers using corpus_status.py (#40).
    "quality_gates": {
        "empty_text_min_chars": 500,
        "min_chars_per_page": 200,
        "max_gibberish_score": 0.5,
        "zero_refs_min_pages": 5,
        "min_median_chunk_chars": 50,
        "min_figure_mean_intensity": 10,
        "max_figures_sampled": 50,
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


# Default config.yaml lives at the repo root, one level up from this
# package. Resolved at call time so test fixtures can drop in a custom
# path.
_DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent.parent / "config.yaml"


def load_config(config_path: Optional[Path] = None) -> dict:
    """Load config.yaml (if present) and merge it over the built-in defaults.

    Missing file is not an error — defaults are returned. Missing keys fall
    back to defaults. A malformed YAML is an error (loud failure).
    """
    if config_path is None:
        config_path = _DEFAULT_CONFIG_PATH
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
