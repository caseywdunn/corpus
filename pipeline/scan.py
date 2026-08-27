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

import logging
import os
import re
import signal
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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


# Tags whose *script* subtag picks a different pack than the language alone
# would. This is the whole reason the public resolver takes BCP-47 rather than
# ISO 639: Tesseract ships separate packs for Fraktur German, the two Chinese
# script variants and Latin-script Serbian, and ISO has no way to name the
# distinction. `de` can only ever mean `deu`; `de-Latf` can mean `deu_latf`.
_BCP47_SCRIPT_PACKS = {
    # 19th-c. German set in Fraktur. `deu` follows as a fallback because a
    # Fraktur document's running heads, tables and citations are often set in
    # roman, and because `deu_latf` is a modern rename — a host carrying only
    # the older `frk` drops `deu_latf` at availability time and would
    # otherwise be left with nothing.
    ("de", "Latf"): ["deu_latf", "deu"],
    ("zh", "Hans"): ["chi_sim"],
    ("zh", "Hant"): ["chi_tra"],
    ("sr", "Latn"): ["srp_latn"],
    ("sr", "Cyrl"): ["srp"],
}

# langdetect emits `zh-cn` / `zh-tw`, which encode script as region and are
# not really ISO at all. They are still what `scan_detection.json` records, so
# they have to keep resolving; the region forms below exist for that
# compatibility and for tags copied out of a detection record by hand.
_BCP47_REGION_PACKS = {
    ("zh", "CN"): ["chi_sim"],
    ("zh", "SG"): ["chi_sim"],
    ("zh", "TW"): ["chi_tra"],
    ("zh", "HK"): ["chi_tra"],
    ("zh", "MO"): ["chi_tra"],
}

# Languages Tesseract has a pack for that `_ISO_TO_TESSERACT` cannot key,
# because they have no ISO 639-1 code. `grc` is the load-bearing one: a
# pre-1900 systematics paper quoting Aristotle needs it, and langdetect's
# vocabulary tops out at `el` (modern Greek).
_BCP47_LANG_PACKS = {
    "grc": ["grc"],
    "sr": ["srp"],
}

# No ISO 639-1 code, so a bare `zh` names a language with two scripts and no
# single pack. The union mirrors what the OSD branch already does for a `Han`
# verdict (`_SCRIPT_TO_TESSERACT`): when the script genuinely is unknown, both
# packs is the honest answer. Prefer `zh-Hans` / `zh-Hant` when it is known.
_BCP47_AMBIGUOUS = {
    "zh": ["chi_sim", "chi_tra"],
}


def _parse_bcp47(tag: str):
    """Return ``(language, script, region)`` from a BCP-47 tag, case-folded.

    Subtags are identified by shape, per BCP-47 itself: 4 alphabetic
    characters is a script, 2 alphabetic or 3 numeric is a region. Variants,
    extensions and private-use subtags are ignored — none of them bear on
    which Tesseract pack to load. Underscores are accepted because that is
    how the tag tends to come back from tooling that treats it as a locale.
    """
    parts = [p for p in re.split(r"[-_]", (tag or "").strip()) if p]
    if not parts:
        return "", None, None
    language = parts[0].lower()
    script = region = None
    for part in parts[1:]:
        if script is None and len(part) == 4 and part.isalpha():
            script = part.capitalize()
        elif region is None and (
            (len(part) == 2 and part.isalpha()) or (len(part) == 3 and part.isdigit())
        ):
            region = part.upper()
    return language, script, region


def bcp47_to_tesseract(tag: str) -> List[str]:
    """Map a BCP-47 language tag to the Tesseract packs that should read it.

    Public because the mapping is generic knowledge about PDFs and OCR, not
    about any one collection: a library's annotation pass resolves its
    per-document language into the `ocrlang` bib directive (#176) and needs
    this exact table. Duplicating it into each library repo is how the copies
    come to disagree with what `scan.py` actually loads, with nothing
    checking that they still agree.

    >>> bcp47_to_tesseract("ru")
    ['rus']
    >>> bcp47_to_tesseract("de-Latf")
    ['deu_latf', 'deu']
    >>> bcp47_to_tesseract("zh-Hant")
    ['chi_tra']
    >>> bcp47_to_tesseract("grc")
    ['grc']
    >>> bcp47_to_tesseract("xx")
    []

    Three things this deliberately does not do:

    **It does not check whether a pack is installed.** The caller is usually
    annotating a bib on a machine that will never run OCR. Availability is a
    property of the *build* host and is applied there, by
    :func:`_resolve_ocrlang_pin`, which drops missing packs with a warning.

    **It does not return vertical CJK companions.** BCP-47 describes language
    and script; vertical setting is typesetting, and `jpn_vert` is a different
    *model*, not a different language. It also must not be unioned with its
    horizontal sibling — measured against the gold set, `jpn_vert` alone
    scores 0.574 on vertical pages, `jpn` alone 0.246, and `jpn_vert+jpn`
    0.186, worse than either. See :func:`detect_vertical_cjk` (#196), which
    finds vertical pages geometrically instead.

    **Nothing calls it at run time.** The OCR language decision stays two
    tiers — an explicit `ocrlang` pin, else detection. Deriving packs from a
    bib field during a build would make *this table* an input to
    `processed.pdf` without being part of any fingerprint: improve the table
    and nothing invalidates, so documents keep their old OCR while the log
    reports the new `-l`. Resolving at annotation time instead leaves
    `ocrlang` a literal, directly-fingerprinted value, and a table change
    does nothing until the bib is rewritten — which invalidates visibly.
    """
    language, script, region = _parse_bcp47(tag)
    if not language:
        return []
    if script and (language, script) in _BCP47_SCRIPT_PACKS:
        return list(_BCP47_SCRIPT_PACKS[(language, script)])
    if region and (language, region) in _BCP47_REGION_PACKS:
        return list(_BCP47_REGION_PACKS[(language, region)])
    if language in _BCP47_LANG_PACKS:
        return list(_BCP47_LANG_PACKS[language])
    pack = _ISO_TO_TESSERACT.get(language)
    if pack:
        return [pack]
    if language in _BCP47_AMBIGUOUS:
        return list(_BCP47_AMBIGUOUS[language])
    return []


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


# Above this share of CJK characters, _gibberish_score returns 0.0
# rather than a meaningless number. 0.30 sits well clear of both
# populations measured on the reference corpus: Japanese and Chinese
# papers run 0.55-0.63, every Latin-script paper 0.00.
_CJK_SHARE_UNSCORABLE = 0.30


def _is_cjk(ch: str) -> bool:
    """True for Han, kana and Hangul characters."""
    cp = ord(ch)
    return (
        0x3040 <= cp <= 0x30FF      # hiragana + katakana
        or 0x3400 <= cp <= 0x4DBF   # CJK ext A
        or 0x4E00 <= cp <= 0x9FFF   # CJK unified
        or 0xF900 <= cp <= 0xFAFF   # CJK compatibility
        or 0xAC00 <= cp <= 0xD7AF   # Hangul syllables
        or 0x1100 <= cp <= 0x11FF   # Hangul jamo
    )


def _gibberish_score(text: str) -> float:
    """Fraction of tokens in ``text`` that look like text-layer garbage.

    A small score (<0.25) is normal even for clean English. Thresholds around
    0.5 reliably separate real prose from PDFs whose text layer maps Cyrillic
    glyphs to Latin-1 byte-by-byte (producing strings like "AKAllEMH5I HAYK").
    langdetect confidence alone does not catch this case because the mangled
    text is made of Latin bytes.

    CJK tokens are excluded from both numerator and denominator. The
    "token of <=2 characters is suspicious" rule encodes an assumption
    about alphabetic scripts that simply does not hold for Han, kana or
    Hangul, where one or two characters is an ordinary word and OCR
    output is legitimately full of short whitespace-separated runs.
    Counting them flagged Yamamori 2014 — which OCR'd correctly under
    `jpn` — as gibberish at 0.55. Scoring only the non-CJK remainder
    keeps the heuristic pointed at what it was built to detect while
    still judging the Latin-script parts of a mixed-script paper (that
    document is 88% Latin by character, so a document-level script test
    would not have caught this).
    """
    alpha_chars = [c for c in text if c.isalpha()]
    if alpha_chars:
        cjk_share = sum(_is_cjk(c) for c in alpha_chars) / len(alpha_chars)
        if cjk_share > _CJK_SHARE_UNSCORABLE:
            # On a CJK-dominant document this measure has no validity in
            # either direction, so refuse to produce a number rather than
            # produce a misleading one. Excluding CJK *tokens* is not
            # enough: what remains is page numbers, figure labels and
            # markup fragments — Yamamori 2014's non-CJK remainder is
            # ['##', 'Li', '\\_', "『'", '=', …] — which scores as
            # garbage no matter how good the OCR was. Yamamori OCR'd
            # correctly under `jpn` and still scored 0.82.
            #
            # Safe because the "is this a scan?" question is now answered
            # independently by _scanned_page_fraction; this heuristic is
            # no longer the only thing standing between a corrupt text
            # layer and the corpus.
            return 0.0

    words = [w.strip(".,;:!?()[]{}\"'") for w in text.split()]
    words = [w for w in words if w]
    scored = []
    for w in words:
        alpha = [c for c in w if c.isalpha()]
        if alpha and sum(_is_cjk(c) for c in alpha) / len(alpha) > 0.5:
            continue
        scored.append(w)
    if not scored:
        return 0.0
    sus = 0
    for w in scored:
        if len(w) <= 2:
            sus += 1
            continue
        if w.isdigit():
            sus += 1
            continue
        if _DIGIT_IN_WORD_RE.search(w) and any(c.isalpha() for c in w):
            sus += 1
    return sus / len(scored)


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


def _probe_sample_pages(n_pages: int, samples: int) -> List[int]:
    """Page indices to OCR-probe, spread through the document's body.

    Sampling is confined to the middle 15–85% of the page range. Both
    ends are furniture: front matter is covers, title pages, plates and
    publisher boilerplate — frequently in a different language from the
    body, and sometimes in no language at all — while the tail is
    references, which are a multilingual pile of proper nouns and the
    worst possible input to language detection. On a 314-page monograph
    a naive 75% sample lands squarely in the bibliography.

    Short documents can't afford the trim, so anything under ~7 pages
    just uses every page and lets per-page detection sort it out.
    """
    if n_pages <= 0:
        return []
    if n_pages < 7:
        return list(range(n_pages))[:samples]
    lo, hi = 0.15, 0.85
    span = hi - lo
    if samples == 1:
        fracs = [0.5]
    else:
        fracs = [lo + span * i / (samples - 1) for i in range(samples)]
    idxs = sorted({min(n_pages - 1, max(0, int(f * n_pages))) for f in fracs})
    return idxs


def _probe_language_by_ocr(
    pdf_path: Path, packs: List[str], pages: Optional[int] = None,
    dpi: Optional[int] = None,
) -> List[tuple]:
    """OCR sampled pages with ``packs`` and detect the language of each.

    Returns a list of ``(iso, confidence, n_pages)`` ordered by how many
    sampled pages voted for that language, most first. Empty on failure.

    Per-page rather than one verdict for the whole document, because
    bilingual PDFs are routine in this material: a Russian original with
    an English typescript translation stapled in from page 8, a French
    paper carrying Pugh's English translation. Concatenating the sample
    and detecting once would report whichever language happened to win
    the character count and silently drop the other — so the translation,
    or the original, would OCR under the wrong model.

    This exists because the cheap signals are unusable on exactly the
    documents that need OCR most. A scan's existing text layer is the
    thing we decided not to trust, so langdetect on it is worthless —
    it read Olfers 1824, a German Fraktur scan, as Catalan and routed it
    to the ``cat`` pack. Tesseract OSD only resolves *script*, and every
    Latin-script language looks alike to it.

    Doing a small amount of real OCR first breaks the circularity: OCR a
    handful of body pages with the broad script-appropriate union, and
    the output is clean enough for langdetect to be right. The full
    document can then be OCR'd once with the languages actually found,
    which is more accurate than running the whole union over every page —
    Tesseract's per-language models compete, and accuracy degrades as
    packs multiply.

    Cost is ``pages`` pages of OCR per scanned document, paid once at
    detection time and repaid by a better full-document pass.

    Callers must apply a confidence floor to the result. langdetect ships
    no Latin profile — ``la`` is simply not among its 55 languages — so
    Latin text cannot be identified, only mis-identified, and it comes
    back as a low-confidence Romance language (Linnaeus 1735 reads as
    Catalan at p=0.71 where every correct detection in the reference
    corpus scored 1.00). Rejecting anything under the floor sends those
    documents to the script-narrowed union, which does include ``lat``.
    That matters here specifically: Latin diagnoses and pre-Linnaean
    works are core taxonomic material, not an edge case.
    """
    if not packs:
        return []
    try:
        import fitz
    except ImportError:
        return []
    if shutil.which("tesseract") is None:
        return []
    try:
        doc = fitz.open(pdf_path)
    except Exception as e:
        logger.debug("Could not open %s for OCR language probe: %s", pdf_path, e)
        return []
    try:
        if pages is None:
            pages = int(CONFIG.get("ocr", {}).get("probe_sample_pages", 5))
        if dpi is None:
            # 300. Tempting to lower this — the probe is the slowest
            # part of detection — but 200 measurably degraded the result
            # on exactly the documents the probe exists for. Don't.
            dpi = int(CONFIG.get("ocr", {}).get("probe_dpi", 300))
        idxs = _probe_sample_pages(len(doc), pages)
        if not idxs:
            return []
        votes: Dict[str, List[float]] = {}
        for idx in idxs:
            try:
                img_bytes = doc[idx].get_pixmap(dpi=dpi).tobytes("png")
            except Exception as e:
                logger.debug("Probe render failed on page %d: %s", idx, e)
                continue
            # Ask OSD what script this page is in, then OCR it with only
            # that script's packs. Per page, not per document: that is
            # what lets a mixed-script volume be read correctly at all.
            page_packs = _probe_packs_for_script(
                _osd_script_for_page(img_bytes), packs
            )
            try:
                res = subprocess.run(
                    ["tesseract", "-", "-", "-l", "+".join(page_packs)],
                    input=img_bytes, capture_output=True, timeout=120,
                )
                page_text = res.stdout.decode("utf-8", errors="replace")
            except Exception as e:
                logger.debug("Probe OCR failed on page %d: %s", idx, e)
                continue
            # A plate, a blank verso or a half-title carries too little
            # text to detect from; counting it would let noise outvote
            # a real body page.
            if len(page_text.strip()) < 200:
                continue
            # Discard pages we evidently OCR'd with the wrong models.
            # OSD is not infallible — on Bernstein 1934 it read a
            # Cyrillic table as Latin, Tesseract dutifully transcribed it
            # with Latin packs ("Ta6auna 4 ... Bron. | NeNe cranguk"),
            # and langdetect called that Catalan at p=0.86. Confidence
            # cannot catch this; the gibberish score can, cleanly — the
            # garbage page scored 0.61 where that document's genuine
            # Russian and German pages scored 0.12-0.39.
            #
            # This is what separates a real second language from an
            # artifact: Bernstein 1934 turns out to be authentically
            # Russian *and* German (a substantial German section, normal
            # for pre-war Soviet zoology), and both survive the gate.
            # Only meaningful on Latin-dominant text. _gibberish_score
            # counts tokens of <=2 characters as suspicious, which is a
            # good signal for Latin-1 mojibake and a terrible one for
            # CJK, where Tesseract's output is legitimately full of short
            # whitespace-separated runs. Gating CJK on it threw away
            # Kawamura 1911a's Japanese pages — the exact content the
            # probe was added to find. Cyrillic sits in between (real
            # Russian body text scored 0.39), so restrict the gate to the
            # case it was built for: Latin characters pretending to be a
            # language they aren't.
            page_scripts = _text_layer_scripts(page_text)
            if page_scripts.get("Latin", 0.0) > 0.50:
                gib = _gibberish_score(page_text)
                if gib > _probe_gibberish_ceiling():
                    logger.debug(
                        "Probe page %d discarded: Latin-script OCR output "
                        "scores %.2f gibberish (likely OCR'd with the wrong "
                        "script's packs)", idx, gib,
                    )
                    continue
            iso, conf = _detect_language(page_text)
            if iso:
                votes.setdefault(iso, []).append(conf)
        # Most-voted first; confidence breaks ties. A language seen on one
        # page out of five is still reported — that is exactly the shape
        # of a stapled-in translation — and the caller decides whether the
        # evidence is strong enough to add its pack.
        return sorted(
            ((iso, max(cs), len(cs)) for iso, cs in votes.items()),
            key=lambda t: (t[2], t[1]),
            reverse=True,
        )
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
# Which Tesseract packs are written in which script. Used only to narrow
# the fallback language union when the script is known but the language
# isn't — see `restrict_to_script` in :func:`_compose_ocr_langs`. Not
# exhaustive: a code missing here simply never gets narrowed away.
_SCRIPT_PACK_FAMILIES = {
    "Latin": frozenset({
        "eng", "deu", "deu_latf", "fra", "lat", "ita", "spa", "por", "nld",
        "pol", "swe", "nor", "dan", "fin", "cat", "ces", "hun", "ron", "slk",
        "slv", "hrv", "est", "lav", "lit", "sqi", "cym", "afr", "ind", "swa",
        "tgl", "vie", "tur", "frk",
    }),
    "Cyrillic": frozenset({"rus", "ukr", "bul", "mkd", "srp", "srp_latn"}),
    "Greek": frozenset({"ell", "grc"}),
}

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

    No fallback union and no automatic ``eng`` suffix — those belong to
    :func:`_compose_ocr_langs`. Returns ``[]`` when no targeted pack is
    installed for the inputs, which is the signal callers use to decide
    whether to warn or to fall back to
    ``CONFIG["ocr"]["ocr_languages_default"]``.

    Both signals contribute (#176). ``visual_script`` leads, because OSD
    reads the page image and so stays right on exactly the documents
    langdetect gets wrong — a corrupt text layer. But it does not *replace*
    ``detected_iso``: OSD is a per-page guess from a single sampled page,
    and on this corpus it is wrong often enough that letting it win outright
    silently degraded 68 papers. In the 1,769-document siphonophore
    production corpus, 188 papers had a p>=0.99 Latin-script language
    detection overruled by a non-Latin OSD verdict — Bigelow 1914 read as
    Cyrillic, Alvarino 1976b as Greek, Broch 1928 as Japanese. All 188 were
    OCR'd. The 120 English ones were rescued by accident, because
    :func:`_compose_ocr_langs` appends ``eng`` anyway; the other 68 simply
    lost their correct pack (fra x35, spa x20, deu x7, ita x2, dan x2,
    hrv x1, nld x1) and OCR'd without it. None tripped the
    ``gibberish_after_ocr`` gate — they scored 0.26-0.45 against a 0.5
    threshold — so the damage was invisible.

    A union degrades gracefully where an exclusive choice does not:
    Tesseract accepts multi-pack ``-l`` fine, so a wrong OSD verdict now
    costs one surplus pack instead of the right one.

    Note ``_SCRIPT_TO_TESSERACT`` has no ``"Latin"`` key, so an OSD verdict
    of Latin never enters the first branch — the union only forms when OSD
    claims a *non-Latin* script, which is the disagreement case. For German
    the Fraktur pack is added when installed (historical papers).
    """
    available = _available_tesseract_langs()
    if not available:
        return []
    chosen: List[str] = []
    if visual_script and visual_script in _SCRIPT_TO_TESSERACT:
        for tess in _SCRIPT_TO_TESSERACT[visual_script]:
            if tess in available and tess not in chosen:
                chosen.append(tess)
    if detected_iso:
        tess = _ISO_TO_TESSERACT.get(detected_iso)
        if tess and tess in available and tess not in chosen:
            chosen.append(tess)
            if tess == "deu" and "deu_latf" in available and "deu_latf" not in chosen:
                chosen.append("deu_latf")
    return chosen


# Memory model for the OCR worker pool (#209), in GB. Each term was measured
# on a 12-core / 30 GB host building a 93-PDF corpuscle, reported with a
# process table in that issue:
#
#   ocrmypdf's own parent process      ~6.0    grows with document length
#   one Tesseract worker               ~1.9    on a dense 300-dpi scan
#   the Grobid JVM, for the whole run  ~3.4    -Xmx4g in docker-compose.yml
#
# plus headroom for docling's models, the page cache and the OS. A worker is
# budgeted at 2.5 rather than the measured 1.9 because the measurement is a
# single host and the failure mode is an OOM kill, not a slowdown.
_OCR_PARENT_GB = 6.0
_OCR_WORKER_GB = 2.5
_OCR_RESERVED_GB = 10.0          # Grobid JVM + docling + page cache + OS


def _host_ram_gb() -> Optional[float]:
    """Total physical RAM in GB, or ``None`` where it cannot be determined.

    ``os.sysconf`` rather than ``psutil``: psutil is present in the
    environment only as a transitive dependency of other packages, and this
    must not start failing because one of them drops it.
    """
    try:
        return (os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")) / 2 ** 30
    except (ValueError, OSError, AttributeError):
        return None


def _run_ocr(cmd: List[str], timeout: float):
    """Run ocrmypdf, and take its Tesseract workers down with it (#209).

    ``subprocess.run(..., timeout=)`` kills the direct child only. ocrmypdf's
    Tesseract workers are *grandchildren*, so on a timeout they were reparented
    to PID 1 and kept running — holding ~20 GB on the host that reported this,
    long after the pipeline had given up on the document.

    Running ocrmypdf in its own process group makes the whole tree killable in
    one call. The cost is that a terminal signal no longer reaches it by
    itself, so ``KeyboardInterrupt`` is forwarded explicitly; without that,
    Ctrl-C would leave exactly the orphans this is meant to prevent.

    A ``SIGKILL`` delivered to the pipeline from outside still cannot be
    handled here — nothing can — but that is now the only path that orphans
    the tree.
    """
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        start_new_session=True,
    )

    def _kill_tree():
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass                          # already gone

    try:
        stdout, stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        _kill_tree()
        proc.communicate()                # reap, so no zombie is left behind
        raise
    except BaseException:                 # KeyboardInterrupt, SystemExit, ...
        _kill_tree()
        proc.communicate()
        raise
    return subprocess.CompletedProcess(cmd, proc.returncode, stdout, stderr)


def _resolve_ocr_jobs() -> Optional[int]:
    """How many Tesseract workers ocrmypdf may run, or ``None`` for its default.

    ocrmypdf runs one worker per CPU and offers no way to reduce it from
    outside: it reads ``multiprocessing.cpu_count()``, which ignores CPU
    affinity, so neither ``taskset`` nor a cgroup CPU limit reaches it. A
    12-core host therefore reaches for ~20 GB of Tesseract on a dense scan and
    is OOM-killed. A cgroup memory limit contains the blast radius but does not
    prevent the overrun — the build still dies, just tidily.

    So the cap has to come from here. ``ocr.jobs`` sets it explicitly;
    otherwise it is derived from host RAM, and only ever *lowers* the worker
    count — on a machine with memory to spare this returns ``None`` and
    ocrmypdf's own default stands, so nothing gets slower.
    """
    configured = CONFIG.get("ocr", {}).get("jobs")
    if configured:
        return int(configured)

    ram = _host_ram_gb()
    cpus = os.cpu_count() or 1
    if ram is None:
        return None                       # unknown host; do not second-guess
    budget = ram - _OCR_RESERVED_GB - _OCR_PARENT_GB
    affordable = max(1, int(budget // _OCR_WORKER_GB))
    if affordable >= cpus:
        return None                       # RAM is not the binding constraint
    return affordable


def _vertical_cjk_hint(langs, detection_result) -> Optional[str]:
    """Tell the operator to consider a vertical CJK pack, or return ``None``.

    A hint rather than an automatic swap, because the measurement says no
    static rule is safe — see ``_VERTICAL_COMPANION``. Silent when the
    operator has already pinned packs (they have made this call themselves)
    and when the vertical model is not installed (advice to pin a pack that
    is not there sends them in circles).
    """
    if detection_result.get("ocrlang_honored"):
        return None
    available = _available_tesseract_langs()
    vertical = [
        _VERTICAL_COMPANION[lang] for lang in langs
        if lang in _VERTICAL_COMPANION and _VERTICAL_COMPANION[lang] in available
    ]
    if not vertical:
        return None
    return (
        "OCR'ing with a horizontal CJK model. If this document is set "
        "vertically, pin `ocrlang = {%s}` in the bib — worth about 2x the "
        "words recovered on vertical pages. Do not pin it alongside %s; the "
        "union scores worse than either alone."
        % ("+".join(vertical), "+".join(langs))
    )


# Below this many detected text lines a page does not get a vote on writing
# direction. A plate carrying two stray labels would otherwise swing the
# decision for a whole document.
_MIN_LINES_FOR_ORIENTATION_VOTE = 5

# A line this much wider than tall is horizontal; this much taller than wide is
# vertical. The gap between them is dead space, so a page of near-square blocks
# abstains rather than voting noise.
_HORIZONTAL_LINE_ASPECT = 1.6
_VERTICAL_LINE_ASPECT = 0.6


def _page_line_orientation(png_bytes: bytes, pack: str) -> tuple:
    """``(tall, wide)`` counts of Tesseract's detected text-line boxes.

    **Geometry, not recognition.** The obvious way to choose between `jpn` and
    `jpn_vert` is to OCR a page with each and keep the more confident result.
    Measured, that does not work: Tesseract reads a vertical column as a stack
    of single characters and is *confident* about each one, so on this corpus's
    vertical pages plain `jpn` scores 61.4 mean confidence against `jpn_vert`'s
    58.1 — it prefers the wrong model, and by enough to be decisive. (It does
    correctly reject `jpn_vert` on horizontal pages, 80.9 against 20.4, so the
    signal is not merely weak; it is asymmetric and misleading.)

    A line *box*, by contrast, is tall or wide whatever the glyphs inside it
    were read as. The separation is about three orders of magnitude — median
    width/height 0.04 on vertically-set pages against 21-53 on horizontal ones,
    including horizontal pages of the same document read with the same pack.
    """
    try:
        result = subprocess.run(
            ["tesseract", "-", "-", "-l", pack, "tsv"],
            input=png_bytes, capture_output=True, timeout=180,
        )
    except Exception as e:
        logger.debug("Orientation probe failed: %s", e)
        return 0, 0
    tall = wide = 0
    for line in result.stdout.decode("utf8", "replace").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 12 or parts[0] != "4":     # tsv level 4 == text line
            continue
        try:
            w, h = int(parts[8]), int(parts[9])
        except ValueError:
            continue
        if w < 20 or h < 20:                       # specks, not lines
            continue
        if w / h < _VERTICAL_LINE_ASPECT:
            tall += 1
        elif w / h > _HORIZONTAL_LINE_ASPECT:
            wide += 1
    return tall, wide


def _is_raster_page(page) -> bool:
    """True if the page is a full-bleed bitmap — i.e. a scan.

    Same test :func:`_scanned_page_fraction` applies, at page granularity.
    """
    area = abs(page.rect.width * page.rect.height)
    if not area:
        return False
    biggest = 0.0
    for block in page.get_text("dict").get("blocks", []):
        if block.get("type") == 1:
            x0, y0, x1, y1 = block.get("bbox", (0, 0, 0, 0))
            biggest = max(biggest, abs((x1 - x0) * (y1 - y0)))
    return biggest / area >= _FULL_PAGE_IMAGE_FRAC


def detect_vertical_cjk(pdf_path: Path, langs: List[str]) -> Optional[str]:
    """Return the ``_vert`` pack this document should use, or ``None``.

    Tesseract ships a separate model for vertically-set CJK and detection never
    reached it, so a vertically-set document was read by a horizontal model —
    worth about half the words on those pages (#196). The two packs cannot be
    combined: measured against a hand transcription, ``jpn_vert`` alone scores
    0.574 on vertical pages where ``jpn`` scores 0.246, while the union of both
    scores 0.186 — *worse than either* — because the models compete for the
    same glyphs. So this is an exclusive choice and it has to be made right.

    **The vote counts only raster pages.** Writing direction is a property of
    the page, not the document, and mixed volumes are routine here — an English
    typescript bound in front of a Japanese scan. But ``ocrmypdf`` takes one
    ``-l`` for the whole document, so one choice must serve. Born-digital pages
    resolve that: ``--redo-ocr`` preserves their existing text and the pack
    never touches them, so they should not get a vote on it. Counting only the
    pages OCR will actually rewrite turns an ambiguous 40% document-wide into
    an unambiguous majority.

    Both directions are equally costly to get wrong — the wrong pack scores
    ~0.21 either way — so this is a plain majority with no threshold to tune.
    """
    if shutil.which("tesseract") is None:
        return None
    candidates = [(lang, _VERTICAL_COMPANION[lang]) for lang in langs
                  if lang in _VERTICAL_COMPANION]
    if not candidates:
        return None                      # not a CJK document; nothing to decide
    available = _available_tesseract_langs()
    candidates = [(h, v) for h, v in candidates if v in available]
    if not candidates:
        return None
    horizontal_pack, vertical_pack = candidates[0]

    try:
        import fitz
    except ImportError:
        return None
    try:
        doc = fitz.open(pdf_path)
    except Exception as e:
        logger.debug("Could not open %s for orientation probe: %s", pdf_path, e)
        return None
    try:
        cfg = CONFIG.get("ocr", {})
        pages = int(cfg.get("probe_sample_pages", 5))
        dpi = int(cfg.get("probe_dpi", 300))
        vertical = horizontal = 0
        for idx in _probe_sample_pages(len(doc), pages):
            page = doc[idx]
            if not _is_raster_page(page):
                continue                 # born-digital: --redo-ocr keeps it
            try:
                img = page.get_pixmap(dpi=dpi).tobytes("png")
            except Exception:
                continue
            tall, wide = _page_line_orientation(img, horizontal_pack)
            if tall + wide < _MIN_LINES_FOR_ORIENTATION_VOTE:
                continue                 # too little on the page to judge
            if tall > wide:
                vertical += 1
            else:
                horizontal += 1
    finally:
        doc.close()

    if vertical > horizontal:
        logger.info(
            "%s: %d of %d sampled scan pages are set vertically; using %s "
            "instead of %s",
            pdf_path.name, vertical, vertical + horizontal,
            vertical_pack, horizontal_pack,
        )
        return vertical_pack
    return None


def _parse_ocrlang(raw: Optional[str]) -> List[str]:
    """Split an ``ocrlang`` bib value into Tesseract pack names.

    Written the way the operator already sees it in ocrmypdf's ``-l`` flag
    and in this module's own "Running OCR" log line — ``pol+eng``. Commas
    and bare whitespace are accepted too, because curators type both.
    """
    if not raw:
        return []
    out: List[str] = []
    for part in re.split(r"[+,\s]+", str(raw).strip().lower()):
        if part and part not in out:
            out.append(part)
    return out


def _resolve_ocrlang_pin(raw: Optional[str]) -> Tuple[List[str], List[str]]:
    """Return ``(honored_packs, dropped_packs)`` for an ``ocrlang`` value.

    A pin is honored only for packs Tesseract actually has installed.
    ocrmypdf exits non-zero on an unknown ``-l`` code, so passing one
    through would turn a typo into a failed OCR pass for that paper; and
    silently OCRing with a pack the operator did not ask for would be
    worse than ignoring the tag. Unknown names are therefore dropped and
    reported in ``scan_detection.json`` so the mistake is visible.
    """
    requested = _parse_ocrlang(raw)
    if not requested:
        return [], []
    available = _available_tesseract_langs()
    if not available:
        return [], requested
    honored = [p for p in requested if p in available]
    dropped = [p for p in requested if p not in available]
    return honored, dropped


def _osd_script_for_page(img_bytes: bytes) -> Optional[str]:
    """Tesseract OSD script name for one rendered page, or None."""
    try:
        res = subprocess.run(
            ["tesseract", "-", "-", "--psm", "0"],
            input=img_bytes, capture_output=True, timeout=30,
        )
    except Exception as e:
        logger.debug("OSD failed during probe: %s", e)
        return None
    combined = (
        res.stdout.decode("utf-8", errors="replace")
        + res.stderr.decode("utf-8", errors="replace")
    )
    m = re.search(r"Script:\s*(\S+)", combined)
    return m.group(1) if m else None


def _probe_packs_for_script(script: Optional[str], full: List[str]) -> List[str]:
    """Packs to OCR a probe page with, given its OSD script.

    Narrowing by script keeps the probe affordable. OCRing every sampled
    page with the whole 13-pack union costs ~85 s per document, because
    Tesseract runs each language model in turn; asking OSD first (a
    fraction of a second) and then OCRing with just that script's packs
    gets the same answer for a fraction of the work. A page OSD reports
    as Japanese needs `jpn`, not thirteen models competing.
    """
    if script and script not in ("Latin", "Fraktur"):
        packs = [p for p in _SCRIPT_TO_TESSERACT.get(script, []) if p in full]
        if packs:
            return packs
    latin = _SCRIPT_PACK_FAMILIES.get("Latin", frozenset())
    narrowed = [p for p in full if p in latin]
    return narrowed or full


def _probe_confidence_floor() -> float:
    """Minimum langdetect confidence for an OCR-probe result to be used."""
    return float(
        CONFIG.get("ocr", {}).get("probe_language_min_confidence", 0.85)
    )


def _probe_gibberish_ceiling() -> float:
    """Above this gibberish score, a probe page's OCR output is treated as
    an artifact of the wrong pack set and its language vote discarded."""
    return float(
        CONFIG.get("ocr", {}).get("probe_max_gibberish", 0.50)
    )


def _accept_probe_languages(votes: List[tuple]) -> List[tuple]:
    """Filter raw probe votes down to the languages we'll OCR with.

    Keeps every language that clears the confidence floor, capped at
    ``ocr.probe_max_languages`` (default 3). The cap exists because the
    point of probing was to *stop* handing Tesseract a big pack union;
    a bilingual original-plus-translation needs two, a trilingual volume
    three, and beyond that the detection is more likely noise than a
    genuinely polyglot document.
    """
    floor = _probe_confidence_floor()
    cap = int(CONFIG.get("ocr", {}).get("probe_max_languages", 3))
    return [v for v in votes if v[1] >= floor][:cap]


def _script_fallback_packs(script: Optional[str] = None) -> List[str]:
    """Installed packs from ``ocr_languages_default``, optionally narrowed
    to one script. The broad set used for the language probe and as the
    last-resort union when nothing better is known.
    """
    available = _available_tesseract_langs()
    cfg_default = CONFIG.get("ocr", {}).get(
        "ocr_languages_default",
        ["eng", "deu", "deu_latf", "fra", "rus", "lat"],
    )
    packs = [lang for lang in cfg_default if lang in available]
    if script:
        allowed = _SCRIPT_PACK_FAMILIES.get(script)
        if allowed:
            narrowed = [p for p in packs if p in allowed]
            if narrowed:
                return narrowed
    return packs


def _compose_ocr_langs(
    detected_iso: Optional[str],
    visual_script: Optional[str] = None,
    restrict_to_script: Optional[str] = None,
    extra_isos: Optional[List[str]] = None,
) -> List[str]:
    """Pick Tesseract language codes to pass to ocrmypdf ``-l``.

    Tries :func:`_resolve_tesseract_packs` first; if that returns nothing
    (unknown language, or pack not installed), falls back to
    ``CONFIG["ocr"]["ocr_languages_default"]`` filtered to installed packs.
    ``eng`` is always appended.

    ``restrict_to_script`` narrows that fallback union to packs written in
    one script. Tesseract slows roughly linearly in pack count and gets
    less accurate as they multiply, so when we know the *script* but not
    the language — a scan whose corrupt text layer is plainly Latin, where
    langdetect's specific guess is untrustworthy — there is no reason to
    also try Cyrillic, Greek and CJK. Cuts the default union from 13 packs
    to 7 on a Latin-script document.
    """
    available = _available_tesseract_langs()
    if not available:
        return []

    chosen = _resolve_tesseract_packs(detected_iso, visual_script)
    if chosen:
        # Additional languages the probe saw on other sampled pages — a
        # stapled-in translation, a bilingual volume. Appended rather
        # than replacing, so the dominant language still leads.
        for iso in extra_isos or []:
            for pack in _resolve_tesseract_packs(iso):
                if pack not in chosen:
                    chosen.append(pack)
        if "eng" in available and "eng" not in chosen:
            chosen.append("eng")
        return chosen

    cfg_default = CONFIG.get("ocr", {}).get(
        "ocr_languages_default",
        ["eng", "deu", "deu_latf", "fra", "rus", "lat"],
    )
    if restrict_to_script:
        allowed = _SCRIPT_PACK_FAMILIES.get(restrict_to_script)
        if allowed:
            narrowed = [lang for lang in cfg_default if lang in allowed]
            # Only narrow if something survives — a corpus whose configured
            # default set has no pack in this script keeps the full union
            # rather than falling through to nothing.
            if narrowed:
                cfg_default = narrowed
    for lang in cfg_default:
        if lang in available and lang not in chosen:
            chosen.append(lang)
    if "eng" in available and "eng" not in chosen:
        chosen.append("eng")
    return chosen



try:  # Pillow's stock decompression-bomb ceiling, captured before we raise it.
    from PIL import Image as _PILImage
    _PIL_DEFAULT_MAX_PIXELS = _PILImage.MAX_IMAGE_PIXELS
except Exception:  # pragma: no cover - Pillow is a hard dependency
    _PIL_DEFAULT_MAX_PIXELS = 178_956_970

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

        # Pillow's decompression-bomb guard trips on high-resolution
        # scan pages ("Decoded image exceeds size limit of 178956970
        # pixels"), which killed QC visualizations on 4 papers once
        # re-OCR started rendering scans we previously passed through.
        # The guard exists to stop a hostile *upload* exhausting memory;
        # here the image is one we just rendered from the operator's own
        # PDF, so the threat model doesn't apply. Raised rather than
        # disabled (None would remove the ceiling entirely), and scoped
        # to this function.
        Image.MAX_IMAGE_PIXELS = max(_PIL_DEFAULT_MAX_PIXELS or 0, 512_000_000)

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
    finally:
        # Restore Pillow's global guard whatever happened above — it is
        # process-wide state and later stages should keep the default.
        try:
            from PIL import Image as _Image
            _Image.MAX_IMAGE_PIXELS = _PIL_DEFAULT_MAX_PIXELS
        except Exception:
            pass



def _annotate_pack_availability(ocrlang: Optional[str], result: Dict,
                                pdf_path: Optional[Path] = None) -> Dict:
    """Tag a :func:`detect_scan_type` result with Tesseract pack availability.

    ``ocrlang`` is the operator's per-document override from the bib
    (#176), or None. When it names at least one installed pack it replaces
    the resolved ``tesseract_packs`` outright — see the ``ocrlang_*``
    fields below.

    Adds two fields:

    - ``tesseract_packs`` — the exact pack list :func:`prepare_pdf` will
      pass to ocrmypdf's ``-l``, produced by calling the same
      :func:`_compose_ocr_langs` with the same arguments. That includes the
      ``eng`` suffix and, when targeted resolution finds nothing, the
      configured fallback union.
    - ``tesseract_pack_available`` — ``True``/``False``/``None``. ``None``
      means no language was detected (e.g., a low-text scan), so we can't
      decide. This tracks *targeted* resolution only — whether this
      document's own language/script has a pack installed — so it stays
      meaningful even though ``tesseract_packs`` is never empty. Recording
      it for born-digital papers too is intentional: it lets you grep
      ``scan_detection.json`` across the corpus to find papers in
      unsupported languages.

    And, only when ``ocrlang`` is set:

    - ``ocrlang_requested`` — the raw string as the operator wrote it.
    - ``ocrlang_honored`` — whether it named any installed pack. False
      means the tag was ignored and detection drove the choice after all.
    - ``ocrlang_dropped`` — requested packs Tesseract doesn't have.

    ``tesseract_packs`` used to record targeted resolution alone while
    :func:`_compose_ocr_langs` appended ``eng`` on top, so a document OCR'd
    with ``rus+eng`` was filed as ``['rus']`` (#176). The field is read by
    operators grepping for bad OCR, and that gap is what made the
    OSD-precedence bug above look worse than it was — the diagnosis blamed
    a missing ``eng`` fallback that had been there all along. Deriving the
    value from the same function prepare_pdf calls keeps the two from
    drifting again.

    Emits a WARNING when ``needs_ocr=True`` but no targeted pack is
    installed; OCR will still run with the default union, but accuracy on
    the detected language/script will be poor. Install the appropriate
    pack (e.g. ``brew install tesseract-lang`` on macOS or
    ``apt-get install tesseract-ocr-<code>`` on Debian) for better results.
    """
    # A language inferred from a text layer we've decided not to trust is
    # not a usable pack hint — see the `language_trusted` note in
    # detect_scan_type's raster path. Drop it so pack resolution falls
    # through to the configured default union.
    iso = result.get("detected_language")
    if not result.get("language_trusted", True):
        iso = None
    visual = result.get("visual_script")
    extra_isos = (result.get("probe_languages") or [])[1:]

    # Targeted resolution drives the availability flag and the warning
    # below; it is the question "does this document's own language have a
    # pack?", which the fallback union would otherwise mask.
    targeted = _resolve_tesseract_packs(iso, visual)
    for extra in extra_isos:
        for pack in _resolve_tesseract_packs(extra):
            if pack not in targeted:
                targeted.append(pack)

    # Same call prepare_pdf makes, so the record cannot disagree with the
    # invocation.
    result["tesseract_packs"] = _compose_ocr_langs(
        iso,
        visual_script=visual,
        restrict_to_script=result.get("script_hint"),
        extra_isos=extra_isos,
    )
    # #196 — swap a horizontal CJK pack for its vertical counterpart when the
    # page geometry says the text is set vertically. Only reached when a CJK
    # pack was resolved, so it costs nothing on the rest of the corpus, and it
    # is a *replacement* rather than an addition: the union of both packs
    # measures worse than either alone.
    if pdf_path is not None:
        vertical_pack = detect_vertical_cjk(pdf_path, result["tesseract_packs"])
        if vertical_pack:
            # The vertical pack goes in ALONE — not merely in place of its
            # horizontal counterpart. Any other model competes for the same
            # glyphs and undoes it, and that includes `eng`. Measured on the
            # vertical pages of the document this was built for:
            #
            #     jpn+eng           0.246      (what detection used to pick)
            #     jpn_vert+jpn+eng  0.186
            #     jpn_vert+eng      0.176      (first attempt at this fix)
            #     jpn_vert          0.574
            #
            # Dropping `eng` is safe precisely because the vote counted only
            # raster pages: any born-digital Latin text in the volume is
            # preserved by --redo-ocr and never reaches Tesseract at all.
            result["tesseract_packs"] = [vertical_pack]
            result["vertical_cjk"] = True
    if iso is None and visual is None:
        result["tesseract_pack_available"] = None
    else:
        result["tesseract_pack_available"] = bool(targeted)

    # #176 — the operator's pin, applied last so it beats every inferred
    # signal. Detection still ran and its verdict stays recorded above:
    # the whole point is to be able to see what the pipeline believed and
    # what the operator overruled it with, side by side, in one file.
    honored, dropped = _resolve_ocrlang_pin(ocrlang)
    if ocrlang:
        result["ocrlang_requested"] = ocrlang
        result["ocrlang_honored"] = bool(honored)
        if dropped:
            result["ocrlang_dropped"] = dropped
    if honored:
        result["tesseract_packs"] = honored
        result["tesseract_pack_available"] = True
        if dropped:
            logger.warning(
                "ocrlang=%r on %s: no installed Tesseract pack for %s; "
                "OCRing with %s. Install the missing pack(s) (macOS: "
                "`brew install tesseract-lang`; Debian: "
                "`apt-get install tesseract-ocr-<code>`).",
                ocrlang, result.get("filename"), "+".join(dropped),
                "+".join(honored),
            )
        return result
    if ocrlang:
        logger.warning(
            "ocrlang=%r on %s names no installed Tesseract pack (%s); "
            "ignoring the override and falling back to detection. Check the "
            "spelling against `tesseract --list-langs` — these are Tesseract "
            "pack names (deu, fra), not ISO codes (de, fr).",
            ocrlang, result.get("filename"), "+".join(dropped),
        )
    if result.get("needs_ocr") and not targeted and (iso or visual):
        logger.warning(
            "No installed Tesseract pack for detected_language=%r visual_script=%r "
            "(%s); OCR will fall back to the default language union and may be "
            "low quality.",
            iso, visual, result.get("filename"),
        )
    return result


_FULL_PAGE_IMAGE_FRAC = 0.50

# At or above this scanned-page fraction the document is treated as a
# scan end to end and rasterizing it costs nothing. Below it there is
# real digital text worth preserving, so OCR uses --redo-ocr instead.
_MOSTLY_SCANNED_FRAC = 0.95

# Tesseract ships a separate model for vertically-set text in each CJK script.
# Detection never reaches them: `_ISO_TO_TESSERACT` maps "ja" to `jpn` and
# nothing downstream reconsiders it.
#
# This deliberately does NOT work like the Fraktur companion above, and the
# difference is measured rather than assumed. `deu`+`deu_latf` are added
# *together* because they degrade gracefully — a surplus pack costs little.
# The vertical models do not. Scored against a hand transcription of the same
# pages, as the fraction of printed words recovered:
#
#                        horizontal Japanese   vertical Japanese
#     jpn                       0.75                 0.25
#     jpn_vert                  0.21                 0.57
#     jpn_vert+jpn               --                  0.19
#
# The union is worse than *either* pack alone, because the two models compete
# for the same glyphs, and each pack is catastrophic on the other's direction.
# So there is no static rule to add here — the choice has to match the
# document, and `ocrmypdf` takes one pack list per document, so it cannot even
# be made per page. The supported answer is the `ocrlang` bib directive
# (#176), and this table exists only to tell the operator when to reach for it.
_VERTICAL_COMPANION = {
    "jpn": "jpn_vert",
    "chi_sim": "chi_sim_vert",
    "chi_tra": "chi_tra_vert",
    "kor": "kor_vert",
}


def _scanned_page_fraction(pdf_path: Path, pages_to_check: int = 8) -> Optional[float]:
    """Fraction of ``pages_to_check`` pages sampled across the document
    that are a scan —
    i.e. carry a *single* image covering at least half the page area.
    ``None`` if the PDF can't be read.

    This is the "is it a scan?" question, and it is independent of what
    the text layer says. A scanned page is one full-bleed bitmap; a
    born-digital page draws glyphs and vectors, with at most some inset
    figures. Measured over the 32-paper smoke corpus the two populations
    are perfectly separated with a wide margin — every scan scored
    0.50–1.00, every born-digital paper scored exactly 0.00 — which is
    why this is a page-count fraction rather than a mean area coverage.
    Mean coverage muddles the two: a two-page scan with one plate page
    lands at 0.50, close enough to a digital paper carrying a couple of
    large figures to be unsafe.

    Why this exists: the rest of :func:`detect_scan_type` only ever asks
    whether the *text layer reads sensibly*, so a scan carrying any
    plausible-looking OCR layer was classified ``born_digital`` and never
    re-OCR'd — 25 of 32 papers in that corpus, including every Fraktur
    one, whose layers are visibly corrupt yet score below the gibberish
    threshold.
    """
    try:
        import fitz  # PyMuPDF
    except ImportError:
        return None
    try:
        doc = fitz.open(pdf_path)
        # Sample across the whole document, not the first N pages. Mixed
        # volumes are common in this material and put the two halves in
        # page order: Kawamura 1911a is a born-digital English typescript
        # for pages 0-7 and a full-page Japanese scan from page 12, so
        # looking only at the front reported 0% raster and skipped OCR on
        # 13 pages of Japanese. Same reasoning as _probe_sample_pages,
        # but without trimming the ends — a scanned front cover or a
        # scanned plate section at the back is exactly what we want to
        # notice here.
        n_total = len(doc)
        if n_total == 0:
            doc.close()
            return None
        if n_total <= pages_to_check:
            indices = list(range(n_total))
        else:
            indices = sorted({
                min(n_total - 1, int(i * n_total / pages_to_check))
                for i in range(pages_to_check)
            })
        scanned_pages = 0
        counted = 0
        for i in indices:
            page = doc[i]
            page_area = abs(page.rect.width * page.rect.height)
            if not page_area:
                continue
            counted += 1
            biggest = 0.0
            # type == 1 is an image block; text blocks are type 0.
            for block in page.get_text("dict").get("blocks", []):
                if block.get("type") == 1:
                    x0, y0, x1, y1 = block.get("bbox", (0, 0, 0, 0))
                    biggest = max(biggest, abs((x1 - x0) * (y1 - y0)))
            if biggest / page_area >= _FULL_PAGE_IMAGE_FRAC:
                scanned_pages += 1
        doc.close()
        if not counted:
            return None
        return scanned_pages / counted
    except Exception as e:
        logger.warning(
            "Could not measure scanned-page fraction for %s: %s", pdf_path.name, e
        )
        return None


def detect_scan_type(pdf_path: Path, ocrlang: Optional[str] = None) -> Dict:
    """Classify the text-layer state of a PDF into one of three buckets:

    - ``born_digital`` — dense, intelligible text layer in a detectable
      language, on pages that are not raster scans. No OCR needed.
    - ``scanned`` — page images rather than digitally-set text. Needs OCR.
      Covers three cases: little or no text layer (``--skip-text``); a
      vendor-boilerplate-only layer; and a scan that already carries an
      OCR text layer of unknown provenance, which is re-OCR'd with
      ``--force-ocr`` so corpus's own OCR replaces it (#UX-P1.1).
    - ``broken_text_layer`` — dense text present but nonsense (e.g.,
      Cyrillic characters mapped 1:1 to Latin-1 bytes, as seen in the
      Stepanjants 1970 scan). Needs OCR with ``--force-ocr`` to replace
      the garbage.

    The policy is "OCR anything that isn't digitally native, even if it
    has been OCR'd before", because third-party text layers vary wildly
    in quality and a bad one is invisible downstream — it reads as clean
    text all the way into the embeddings. Set
    ``ocr.reocr_scanned_text_layers: false`` to trust existing layers
    instead.

    ``ocrlang`` is the optional per-document OCR language override from
    the paper's bib entry (#176) — a ``+``-joined list of Tesseract pack
    names, e.g. ``pol+eng``. When it names at least one installed pack it
    pins ``tesseract_packs`` outright, overruling both langdetect and
    Tesseract OSD. It selects *packs only*: it does not force OCR, so
    tagging a born-digital paper changes nothing. ``needs_ocr`` stays a
    detection decision.

    The returned dict is written to ``scan_detection.json`` and consumed
    by :func:`prepare_pdf`. Downstream stages can also read the
    ``detected_language`` field as a useful signal.
    """
    try:
        import fitz  # PyMuPDF
    except ImportError:
        logger.warning("PyMuPDF not available; treating as born-digital (no OCR)")
        return _annotate_pack_availability(ocrlang, {
            "filename": pdf_path.name,
            "file_type": "born_digital",
            "has_text": True,
            "needs_ocr": False,
            "detected_language": None,
            "language_confidence": 0.0,
            "gibberish_score": 0.0,
            "ocr_mode": None,
        }, pdf_path)

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
        return _annotate_pack_availability(ocrlang, {
            "filename": pdf_path.name,
            "file_type": "scanned",
            "has_text": False,
            "needs_ocr": True,
            "detected_language": None,
            "language_confidence": 0.0,
            "gibberish_score": 0.0,
            "ocr_mode": "skip_text",
        }, pdf_path)

    stripped = total_text.strip()
    total_chars = len(stripped)
    avg_chars_per_page = total_chars / max(pages_to_check, 1)
    has_substantial_text = total_chars > 500
    has_good_density = avg_chars_per_page > 100

    # Low-text path: treat as scanned and OCR it.
    if not (has_substantial_text and has_good_density):
        # There is no text layer to infer a language from, so without a
        # probe these documents OCR with the entire configured union —
        # every installed pack competing on every page, which is the
        # least accurate setting available. A short OCR sample buys a
        # real language and a targeted pack instead. Script is unknown
        # here (nothing to measure), so the probe uses the full union.
        accepted: List[tuple] = []
        if CONFIG.get("ocr", {}).get("probe_language_by_ocr", True):
            votes = _probe_language_by_ocr(pdf_path, _script_fallback_packs())
            accepted = _accept_probe_languages(votes)
            if accepted:
                logger.info(
                    "%s: no text layer to detect from; OCR probe found %s — "
                    "using for pack selection",
                    pdf_path.name,
                    ", ".join(f"{i} (p={c:.2f}, {n}pp)" for i, c, n in accepted),
                )
            elif votes:
                logger.info(
                    "%s: OCR probe returned only low-confidence languages "
                    "(%s) — OCRing with the full pack union instead",
                    pdf_path.name,
                    ", ".join(f"{i} p={c:.2f}" for i, c, _ in votes),
                )
        return _annotate_pack_availability(ocrlang, {
            "filename": pdf_path.name,
            "file_type": "scanned",
            "detection_reason": "no_text_layer",
            "has_text": False,
            "needs_ocr": True,
            "detected_language": accepted[0][0] if accepted else None,
            "language_confidence": accepted[0][1] if accepted else 0.0,
            "language_trusted": bool(accepted),
            "probe_languages": [i for i, _, _ in accepted],
            "gibberish_score": 0.0,
            "total_chars_sampled": total_chars,
            "pages_checked": pages_to_check,
            "ocr_mode": "skip_text",
        }, pdf_path)

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
        return _annotate_pack_availability(ocrlang, {
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
        }, pdf_path)

    # --- Raster-page cross-check -------------------------------------
    # Everything below this point reasons about the *content* of the text
    # layer. That cannot distinguish "digitally typeset" from "scanned,
    # then OCR'd by someone else at unknown quality" — and the latter is
    # most of a historical corpus. Ask the page geometry instead, and
    # re-OCR when the pages are bitmaps, whatever the layer says.
    _ocr_cfg = CONFIG.get("ocr", {})
    if _ocr_cfg.get("reocr_scanned_text_layers", True):
        frac_min = float(_ocr_cfg.get("scan_page_fraction_min", 0.40))
        coverage = _scanned_page_fraction(pdf_path, pages_to_check)
        if coverage is not None and coverage >= frac_min:
            logger.info(
                "%s: %.0f%% of sampled pages are full-page images — treating as "
                "a scan and re-OCRing over the existing text layer",
                pdf_path.name, 100 * coverage,
            )
            script_hint = "Latin" if latin_frac > 0.90 else None
            lang_trusted = latin_frac <= 0.90
            layer_lang = lang
            accepted: List[tuple] = []
            if _ocr_cfg.get("probe_language_by_ocr", True):
                # Probe every scan, not only the ones whose layer language
                # is untrustworthy. A layer containing Cyrillic or CJK does
                # pin *a* language reliably — but "reliably one language" is
                # wrong for a bilingual document, and those are common here:
                # Margulis 1976a is a Russian scan through p. 7 with an
                # English typescript from p. 8, Carré 1968 is French with an
                # English translation stapled in. Detecting from the layer
                # alone would OCR half of each under the wrong model.
                # Probe with the *full* pack union, never narrowed by
                # script_hint. That hint comes from the text layer, which
                # is the thing we distrust — and on a mixed-script volume
                # it is actively wrong: Kawamura 1911a's layer is the
                # English typescript in front, so a Latin-only probe OCR'd
                # its Japanese pages with Latin models and reported the
                # document as monolingual English. script_hint still
                # narrows the *fallback* union below, where it is only
                # ever consulted if the probe found nothing.
                votes = _probe_language_by_ocr(pdf_path, _script_fallback_packs())
                accepted = _accept_probe_languages(votes)
                if accepted:
                    logger.info(
                        "%s: OCR probe found %s (text layer claimed %r) — "
                        "using for pack selection",
                        pdf_path.name,
                        ", ".join(f"{i} (p={c:.2f}, {n}pp)" for i, c, n in accepted),
                        layer_lang,
                    )
                    lang, conf, lang_trusted = accepted[0][0], accepted[0][1], True
                elif votes and not lang_trusted:
                    logger.info(
                        "%s: OCR probe returned only low-confidence languages "
                        "(%s) — OCRing with the %s pack union instead "
                        "(langdetect cannot identify Latin, so low-confidence "
                        "Romance hits on historical material are suspect)",
                        pdf_path.name,
                        ", ".join(f"{i} p={c:.2f}" for i, c, _ in votes),
                        script_hint or "full",
                    )
            return _annotate_pack_availability(ocrlang, {
                "filename": pdf_path.name,
                "file_type": "scanned",
                "detection_reason": "raster_page_images",
                "has_text": True,
                "needs_ocr": True,
                "detected_language": lang,
                "language_confidence": conf,
                # We just declared this text layer untrustworthy, so a
                # language guessed *from* it inherits that. Script is the
                # discriminator: Cyrillic/CJK/Greek characters in the layer
                # are strong evidence langdetect got it right (they survive
                # bad OCR as themselves), but one Latin-script language is
                # indistinguishable from another once the glyphs are
                # mangled — which is how Olfers 1824, a German Fraktur
                # scan, was detected as Catalan and routed to the `cat`
                # pack. Untrusted ⇒ OCR with the configured default union
                # instead of one confidently-wrong pack.
                "language_trusted": lang_trusted,
                # Kept separate so the decision is auditable after the
                # fact: what the (untrusted) text layer claimed, versus
                # what OCRing a sample actually found. `probe_languages`
                # is ordered by how many sampled pages voted for each, so
                # a second entry is a real second language in the volume.
                "text_layer_language": layer_lang,
                "probe_languages": [i for i, _, _ in accepted],
                "probe_confidence": accepted[0][1] if accepted else None,
                # Script survives bad OCR even when the language guess
                # doesn't, so it's still a usable hint for narrowing the
                # fallback pack union when the probe comes back empty.
                "script_hint": script_hint,
                "gibberish_score": gib,
                "text_layer_scripts": scripts,
                "visual_script": None,
                "scanned_page_fraction": coverage,
                "total_chars_sampled": total_chars,
                "pages_checked": pages_to_check,
                # Not skip_text — that would preserve the very layer we
                # decided not to trust. Which of the two rewriting modes
                # depends on whether the whole volume is a scan: a mixed
                # one (Kawamura 1911a is a digital English typescript
                # bound in front of a Japanese scan) must not have its
                # digital half rasterized, so it gets --redo-ocr.
                "ocr_mode": (
                    "force_ocr" if coverage >= _MOSTLY_SCANNED_FRAC
                    else "redo_ocr"
                ),
            }, pdf_path)

    threshold = float(_ocr_cfg.get("gibberish_threshold", 0.65))

    # --- Cheap gibberish path ---
    # High gibberish score is a direct signal the text layer is corrupt,
    # regardless of script. Catches the obvious cases.
    if gib > threshold:
        return _annotate_pack_availability(ocrlang, {
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
        }, pdf_path)

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
            return _annotate_pack_availability(ocrlang, {
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
            }, pdf_path)

    return _annotate_pack_availability(ocrlang, {
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
        # Recorded on the no-OCR path too so the raster check is auditable
        # in both directions: grep for a high coverage here and you have
        # found a scan this classifier decided to trust.
        "scanned_page_fraction": _scanned_page_fraction(pdf_path, pages_to_check),
        "total_chars_sampled": total_chars,
        "pages_checked": pages_to_check,
        "ocr_mode": None,
    }, pdf_path)


# ocrmypdf reports per-page trouble on stderr and still exits 0. These
# are the messages that mean text was silently dropped or is suspect.
_OCR_WARNING_PATTERNS = (
    ("took too long to OCR", "page(s) hit the per-page OCR timeout and were "
                             "left blank — raise ocr.tesseract_page_timeout"),
    ("Suppressing OCR output", "page(s) had OCR text suppressed for improbable "
                               "aspect ratio"),
    ("possibly poor OCR", "page(s) flagged by Tesseract as poor quality"),
    ("Line cannot be recognized", "line(s) could not be recognized"),
)


def _ocr_timeout_for(pdf_path: Path, n_langs: int = 1) -> float:
    """Wall-clock cap for one document's ocrmypdf call, scaled by size.

    A flat cap cannot fit this corpus. Measured on a 12-core box with the
    v1.0 re-OCR path, throughput varies ~15x with how many Tesseract
    language packs are loaded, not just page count:

        Totton 1965a    314 pages, `eng` alone       1.3 s/page
        Linnaeus 1735    13 pages, 7-pack union     20.0 s/page

    Extrapolated to the largest document in the full siphonophore
    library (delle Chiaje 1830-31, 1,549 pages) that is 34 minutes at
    best and ~8.6 hours at worst — against a former flat default of 30
    minutes, which the *best* case already exceeded. A cluster array task
    with fewer cores than a workstation makes it worse.

    Raising the flat number instead would mean a genuinely stuck 3-page
    document also burns hours before failing, which is the thing a
    timeout exists to prevent. So: a floor for small documents plus a
    per-page allowance that scales with the real cost driver. On timeout
    the stage fails outright — ocrmypdf is killed mid-write — so the
    budget is deliberately generous, and on a cluster the SLURM walltime
    is the real backstop anyway.
    """
    cfg = CONFIG.get("stage_timeouts", {})
    base = float(cfg.get("ocr", 1800))
    # 30 s/page is already the *worst* measured rate (7-pack union at
    # 20 s/page, plus margin), so it is not scaled again by pack count —
    # doing that double-counts and produced a 51-hour budget for a
    # 1,549-page monograph. n_langs is accepted for logging and future
    # tuning, not applied.
    per_page = float(cfg.get("ocr_per_page", 30))
    try:
        import fitz
        with fitz.open(pdf_path) as doc:
            n_pages = len(doc)
    except Exception:
        return base
    return max(base, per_page * n_pages)


def _log_ocr_warnings(stderr: Optional[str], name: str) -> None:
    """Surface the ocrmypdf warnings that indicate lost or suspect text.

    These arrive on a *successful* exit, and the previous code only
    logged stderr when ocrmypdf failed — so "took too long to OCR -
    skipping" was discarded on every run it mattered.
    """
    if not stderr:
        return
    for needle, explanation in _OCR_WARNING_PATTERNS:
        n = stderr.count(needle)
        if n:
            level = logger.warning if "timeout" in explanation else logger.info
            level("[%s] ocrmypdf: %d %s", name, n, explanation)


def _warn_on_empty_ocr_pages(output_pdf: Path, name: str) -> None:
    """Warn when OCR produced no text at all on some pages.

    The end-state check, independent of whatever ocrmypdf said: a page
    that came out of OCR with zero characters is either genuinely blank,
    a plate, or a page we lost. Worth a line in the log either way,
    because the alternative is discovering it in the embeddings.
    """
    try:
        import fitz
        with fitz.open(output_pdf) as doc:
            empty = [i + 1 for i, pg in enumerate(doc) if not pg.get_text("text").strip()]
            total = len(doc)
    except Exception:
        return
    if empty and total:
        logger.warning(
            "[%s] %d/%d page(s) have no text after OCR (pages %s). Blank pages "
            "and plates are expected; a run of them is not.",
            name, len(empty), total,
            ", ".join(str(p) for p in empty[:12]) + ("…" if len(empty) > 12 else ""),
        )


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
    installed. See :func:`_compose_ocr_langs`. When the paper's bib entry
    carried an ``ocrlang`` override that named an installed pack (#176),
    that list is used verbatim instead.

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
    mode_flag = {
        "force_ocr": "--force-ocr",
        # Mixed volumes — a born-digital typescript bound with a scanned
        # original — must not be rasterized wholesale. --redo-ocr strips
        # and replaces existing OCR text while leaving genuine digital
        # text alone, so the good half survives intact.
        "redo_ocr": "--redo-ocr",
    }.get(ocr_mode, "--skip-text")

    if detection_result.get("vertical_cjk"):
        # #196 — the vertical/horizontal CJK choice was made from page
        # geometry during detection and is already applied to the recorded
        # packs. Take them verbatim rather than recomputing, which resolves
        # `ja` straight back to `jpn` and silently discards the decision.
        # Recomputing here is what made the run log disagree with
        # scan_detection.json on the first attempt: detection logged
        # "using jpn_vert" and OCR then ran `langs=jpn+eng`.
        langs = list(detection_result.get("tesseract_packs") or [])
    elif detection_result.get("ocrlang_honored"):
        # #176 — the operator pinned the packs in the bib. detect_scan_type
        # already validated them against `tesseract --list-langs` and wrote
        # the survivors here, so use them verbatim. No `eng` is appended:
        # being able to pin `jpn` alone, with no Latin pack competing for
        # the same glyphs, is most of the point of the override.
        langs = list(detection_result.get("tesseract_packs") or [])
    else:
        probe_langs = detection_result.get("probe_languages") or []
        langs = _compose_ocr_langs(
            detection_result.get("detected_language")
            if detection_result.get("language_trusted", True) else None,
            visual_script=detection_result.get("visual_script"),
            restrict_to_script=detection_result.get("script_hint"),
            extra_isos=probe_langs[1:],
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

    hint = _vertical_cjk_hint(langs, detection_result)
    if hint:
        logger.info("%s: %s", input_pdf.name, hint)

    # #209 — ocrmypdf otherwise runs one Tesseract worker per CPU, ~1.9 GB
    # each on a dense scan, and there is no way to reduce it from outside: it
    # reads multiprocessing.cpu_count(), which ignores CPU affinity, so neither
    # taskset nor a cgroup CPU limit reaches it.
    ocr_jobs = _resolve_ocr_jobs()

    logger.info(
        "Running OCR on %s | file_type=%s mode=%s langs=%s timeout=%.0fs%s",
        input_pdf.name,
        detection_result.get("file_type"),
        mode_flag,
        lang_arg,
        _ocr_timeout_for(input_pdf, len(langs)),
        f" jobs={ocr_jobs}" if ocr_jobs else "",
    )
    # Per-page OCR timeout. ocrmypdf's default is far too tight for what
    # this pipeline asks of Tesseract: on a dense 300-dpi historical scan
    # with several language packs loaded, a single page legitimately runs
    # for minutes. When the timeout fires ocrmypdf "copies the
    # preprocessed page into the final output" — a blank page — and still
    # exits 0, so the loss is invisible. Observed on Linnaeus 1735:
    # identical command, two runs, and pages 3/4/9/10 went from thousands
    # of characters to zero (56,631 chars → 15,292) purely on machine
    # load, with `OCR completed successfully` logged both times.
    page_timeout = str(int(CONFIG.get("ocr", {}).get("tesseract_page_timeout", 900)))
    cmd = [
        "ocrmypdf",
        mode_flag,
        "-l", lang_arg,
        "--tesseract-timeout", page_timeout,
        "--optimize", optimize_level,
        "--color-conversion-strategy", "RGB",
        "--output-type", "pdf",
        str(input_pdf),
        str(output_pdf),
    ]

    if ocr_jobs:
        cmd += ["--jobs", str(ocr_jobs)]

    ocr_timeout = _ocr_timeout_for(input_pdf, len(langs))
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

    if result.returncode != 0 and mode_flag == "--redo-ocr":
        # --redo-ocr is the fussiest of the three modes: it refuses on
        # PDFs whose existing text it can't cleanly attribute. Falling
        # back to --force-ocr still gets the scanned half OCR'd, at the
        # cost of rasterizing the digital half — worse than --redo-ocr,
        # much better than shipping the original untouched.
        logger.warning(
            "--redo-ocr failed (exit %d) on %s; retrying with --force-ocr",
            result.returncode, input_pdf.name,
        )
        if result.stderr:
            logger.warning("ocrmypdf stderr (head): %s", result.stderr[:500])
        cmd[1] = mode_flag = "--force-ocr"
        try:
            result = _run_ocr(cmd, ocr_timeout)
        except subprocess.TimeoutExpired:
            logger.warning(
                "OCR timed out after %.0fs on %s", ocr_timeout, input_pdf.name
            )
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
        logger.info(
            "OCR completed successfully (mode=%s langs=%s)", mode_flag, lang_arg
        )
        _log_ocr_warnings(result.stderr, input_pdf.name)
        _warn_on_empty_ocr_pages(output_pdf, input_pdf.name)
