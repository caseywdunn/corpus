"""Recover raw font codes only when a PDF's exact regional CMap proves them.

Some PDF parsers expose a Type1 font's encoded bytes instead of its ToUnicode
mapping (#303). This is not a dictionary or a general native-text replacement:
every observed nonspace scalar must re-encode exactly from the same region.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time
import unicodedata

from .source_layout import item_bounds

PDF_CMAP_POLICY = "exact-regional-type1-scalar-cmap-v1"
_MAX_REGION_CHARS = 8192
_GUARDED = set("±−-+=µμ－＋＝⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻₀₁₂₃₄₅₆₇₈₉₊₋")


def _guarded(char):
    return char in _GUARDED or unicodedata.category(char) == "Sm"


def parse_scalar_cmap(data):
    """Read the bounded one-byte/scalar subset; never guess unsupported maps."""
    if len(data) > 65_536:
        raise ValueError("scalar_cmap_size_limit")
    source = data.decode("ascii", errors="strict") if isinstance(data, bytes) else data
    source = re.sub(r"%[^\r\n]*", "", source)
    mapping = {}
    blocks = re.findall(r"(\d+)\s+begin(bfchar|bfrange)\b(.*?)end\2", source, re.S)
    if not blocks:
        raise ValueError("unsupported_or_missing_scalar_cmap")
    for count, kind, body in blocks:
        pattern = r"<([0-9a-fA-F]+)>\s*<([0-9a-fA-F]+)>"
        if kind == "bfrange":
            pattern += r"\s*<([0-9a-fA-F]+)>"
        rows = re.findall(pattern, body)
        if len(rows) != int(count) or re.sub(pattern, "", body).strip():
            raise ValueError("unsupported_scalar_cmap_form")
        for row in rows:
            first = int(row[0], 16)
            last = int(row[1], 16) if kind == "bfrange" else first
            encoded_unicode = row[-1]
            if len(row[0]) > 2 or first > last or last > 255 or len(encoded_unicode) != 4:
                raise ValueError("unsupported_scalar_cmap_range")
            initial = int(encoded_unicode, 16)
            for code in range(first, last + 1):
                value = initial + code - first
                if value > 0xFFFF or 0xD800 <= value <= 0xDFFF:
                    raise ValueError("non_scalar_cmap_value")
                char = chr(value)
                if code in mapping and mapping[code] != char:
                    raise ValueError("conflicting_scalar_cmap")
                mapping[code] = char
    return mapping


def _comparison_char(char):
    # Docling exposes Kangxi radicals for some ordinary Han characters. Only
    # these one-scalar compatibility equivalents participate in alignment;
    # never NFKC-fold superscripts, full-width operators, or numeric evidence.
    if 0x2E80 <= ord(char) <= 0x2FFF or 0xF900 <= ord(char) <= 0xFAFF:
        normalized = unicodedata.normalize("NFKC", char)
        if len(normalized) == 1:
            return normalized
    return char


def _compact(text):
    return [(i, _comparison_char(c)) for i, c in enumerate(text) if not c.isspace()]


def repair_cmap_region(text, chars, fonts, *, verify_sign=None, verify_digit=None):
    """Return exact scalar edits or uncertainty; preserve all raw observations."""
    old = _compact(text)
    source = [c for c in chars if not c["text"].isspace()]
    if not old or len(old) > _MAX_REGION_CHARS:
        return text, [], []
    if [c for _, c in old] == [_comparison_char(c["text"]) for c in source]:
        return text, [], []
    active = {c["font"] for c in source if c["font"] in fonts}
    if not active:
        return text, [], []
    reverse = {}
    for name in active:
        values = defaultdict(list)
        for code, char in fonts[name].get("mapping", {}).items():
            values[char].append(code)
        reverse[name] = values
    # Merely using a mapped font is not evidence of this decoder failure.
    # Gate even uncertainty on predominantly raw-code observations; otherwise
    # an ordinary layout omission would generate a misleading CMap warning.
    support, eligible, pairs = 0, 0, set()
    for (_, observed), glyph in zip(old, source):
        native = glyph["text"]
        codes = reverse.get(glyph["font"], {}).get(native, [])
        if len(codes) != 1:
            continue
        expected = fonts[glyph["font"]].get("observed_codes", {}).get(codes[0], chr(codes[0]))
        if _comparison_char(expected) == _comparison_char(native):
            continue
        eligible += 1
        if _comparison_char(expected) == observed:
            support += 1
            pairs.add((glyph["font"], observed, native))
    if support < 4 or len(pairs) < 3 or support < .8*eligible:
        return text, [], []
    ambiguous = [name for name in active if fonts[name].get("error")]
    if ambiguous:
        return text, [], [{"reason": "ambiguous_or_unsupported_source_font_map", "fonts": sorted(ambiguous)}]
    if len(old) != len(source):
        return text, [], [{"reason": "cmap_region_not_one_to_one", "observed_scalars": len(old), "source_scalars": len(source)}]
    proposed = []
    for (offset, observed), glyph in zip(old, source, strict=True):
        native = glyph["text"]
        codes = reverse.get(glyph["font"], {}).get(native, [])
        if len(codes) > 1:
            return text, [], [{"reason": "ambiguous_reverse_cmap", "font": glyph["font"], "native": native}]
        font = fonts.get(glyph["font"], {})
        expected = font.get("observed_codes", {}).get(codes[0], chr(codes[0])) if codes else native
        if _comparison_char(expected) != observed:
            return text, [], [{"reason": "cmap_reencoding_does_not_match_observation", "char_offset": offset}]
        if codes and observed != _comparison_char(native):
            # A raw ASCII byte is the specifically supported decoder failure.
            if not 33 <= codes[0] <= 126 or unicodedata.category(native).startswith("C"):
                return text, [], [{"reason": "unsupported_cmap_scalar_recovery", "char_offset": offset}]
            proposed.append((offset, text[offset], native, glyph, fonts[glyph["font"]]))
    if len(proposed) < 4 or len({(a, b, g["font"]) for _, a, b, g, _ in proposed}) < 3:
        return text, [], [{"reason": "insufficient_distinct_cmap_evidence"}]
    # Exact re-encoding establishes the candidate positions, but a wrong
    # ToUnicode sign must not bypass source-raster verification. Retain that
    # raw scalar and flag it while preserving independently proved positions.
    unresolved = []
    allowed = []
    for proposal in proposed:
        offset, observed, native, glyph, _ = proposal
        if _guarded(native) and (verify_sign is None or not verify_sign(glyph, native)):
            unresolved.append({"reason": "cmap_scientific_glyph_not_source_confirmed", "char_offset": offset,
                               "before": observed, "native": native, "source_bbox": glyph["bbox"]})
        else:
            allowed.append(proposal)
    if not allowed:
        return text, [], unresolved
    result = list(text)
    proofs = {}
    for offset, observed, native, glyph, font in allowed:
        result[offset] = native
        key = (font["font_xref"], observed, native)
        if key not in proofs:
            proofs[key] = {"font": glyph["font"], "font_xref": font["font_xref"],
                           "tounicode_xref": font["tounicode_xref"], "cmap_sha256": font["cmap_sha256"],
                           "encoding_xref": font.get("encoding_xref"), "encoding_sha256": font.get("encoding_sha256"),
                           "observed": observed, "replacement": native, "source_bbox": glyph["bbox"],
                           "first_char_offset": offset, "occurrences": 0}
        proofs[key]["occurrences"] += 1
    formatting = []
    allowed_offsets = {p[0] for p in allowed}
    for i, ((offset, _), glyph) in enumerate(zip(old, source, strict=True)):
        # This producer only formats the recovered scalar of a numeric unit;
        # no generic flags=4 / small-font / affiliation-footnote inference.
        if offset not in allowed_offsets or not _raised_unit_digit(source, i):
            continue
        evidence = verify_digit(glyph) if verify_digit else {"verified": False, "reason": "digit_verifier_unavailable"}
        record = {"char_offset": offset, "scalar": glyph["text"], "source_bbox": glyph["bbox"],
                  "source_origin": glyph["origin"], "source_size": glyph["size"],
                  "base_glyph": source[i-1], "raster": evidence}
        if evidence.get("verified"):
            replacement_digit = glyph["text"].translate(str.maketrans("23", "²³"))
            result[offset] = replacement_digit
            preceding = offset-1
            while preceding >= 0 and text[preceding].isspace():
                result[preceding] = ""
                preceding -= 1
            formatting.append({**record, "replacement": replacement_digit,
                               "evidence": "raised_unit_geometry_and_dual_unhinted_digit_ocr"})
        else:
            unresolved.append({**record, "reason": "cmap_unit_exponent_not_source_confirmed"})
    replacement = "".join(result)
    proof_hash = hashlib.sha256(json.dumps(source, ensure_ascii=False, sort_keys=True).encode()).hexdigest()
    return replacement, [{"charspan": [0, len(text)], "original": text, "replacement": replacement,
                          "evidence": "same_region_exact_font_cmap_reencoding", "matched_scalars": len(old),
                          "changed_scalars": len(allowed), "source_glyphs_sha256": proof_hash,
                          "mappings": list(proofs.values()), "formatting": formatting}], unresolved


def _raised_unit_digit(chars, i):
    glyph = chars[i]
    if glyph["text"] not in {"2", "3"} or i == 0:
        return False
    base = chars[i-1]
    # Examples: mg/m³, 2m³, 2mm². A letter in prose (including an
    # affiliation marker following a name) supplies no numeric-unit context.
    prefix = "".join(c["text"] for c in chars[max(0, i-20):i])
    if not re.search(r"(?:[mkµμ]?g/|[0-9][.,0-9]*)(?:[cmkµμ]?m)$", prefix):
        return False
    size = base["size"]
    return (size > 0 and .35 <= glyph["size"]/size <= .7
            and -.1*size <= glyph["bbox"][0]-base["bbox"][2] <= .6*size
            and .25*size <= base["origin"][1]-glyph["origin"][1] <= .8*size
            and glyph["bbox"][3] <= base["origin"][1]+.1*size)


class _DigitVerifier:
    """Bound independent raster work per page, document, pixels and time."""
    MAX_PAGE = 8
    MAX_DOCUMENT = 32
    MAX_CROP_PIXELS = 40_000
    MAX_TOTAL_PIXELS = 1_000_000
    DOCUMENT_SECONDS = 30
    CALL_SECONDS = 3

    def __init__(self):
        self.cache = {}
        self.page_counts = defaultdict(int)
        self.total_pixels = 0
        self.started = time.monotonic()
        self.executable = shutil.which("tesseract")

    def verify(self, page, glyph):
        import fitz
        key = (page.number, tuple(glyph["bbox"]), glyph["text"])
        if key in self.cache:
            return self.cache[key]
        result = {"verified": False, "dpi": 600, "modes": [7, 10], "language": "eng", "observations": []}
        if (len(self.cache) >= self.MAX_DOCUMENT or self.page_counts[page.number] >= self.MAX_PAGE
                or time.monotonic()-self.started >= self.DOCUMENT_SECONDS):
            return {**result, "reason": "digit_ocr_budget_exhausted"}
        self.page_counts[page.number] += 1
        self.cache[key] = result
        if not self.executable:
            result["reason"] = "digit_ocr_unavailable"
            return result
        if not all(math.isfinite(v) for v in glyph["bbox"]):
            result["reason"] = "digit_ocr_invalid_source_box"
            return result
        box = fitz.Rect(glyph["bbox"]) + (-.5, -.5, .5, .5)
        if box.is_empty or box.is_infinite:
            result["reason"] = "digit_ocr_invalid_source_box"
            return result
        # Bounding before rendering also limits hostile or malformed boxes.
        pixels = (int(box.width*600/72)+2)*(int(box.height*600/72)+2)
        if (pixels > self.MAX_CROP_PIXELS
                or self.total_pixels+pixels > self.MAX_TOTAL_PIXELS):
            result["reason"] = "digit_ocr_pixel_budget_exhausted"
            return result
        self.total_pixels += pixels
        data = page.get_pixmap(clip=box, dpi=600, colorspace=fitz.csGRAY).tobytes("png")
        result["crop_sha256"] = hashlib.sha256(data).hexdigest()
        result["crop_bbox"] = list(box)
        for mode in (7, 10):
            remaining = self.DOCUMENT_SECONDS-(time.monotonic()-self.started)
            if remaining <= 0:
                result["reason"] = "digit_ocr_budget_exhausted"
                return result
            try:
                process = subprocess.run([self.executable, "stdin", "stdout", "-l", "eng", "--psm", str(mode)],
                    input=data, capture_output=True, timeout=min(self.CALL_SECONDS, remaining))
            except (OSError, subprocess.TimeoutExpired) as exc:
                result["reason"] = type(exc).__name__
                return result
            raw = process.stdout.decode("utf-8", errors="replace")
            result["observations"].append({"mode": mode, "text": raw, "returncode": process.returncode})
            if process.returncode or raw.strip() != glyph["text"]:
                result["reason"] = "digit_ocr_disagreement"
                return result
        result["verified"] = True
        return result


def pdf_cmap_producer():
    from .source_spaces import source_spacing_producer
    installed = source_spacing_producer()
    return {"policy": PDF_CMAP_POLICY,
            "digit_ocr": {**{key: installed[key] for key in
                ("available", "version", "traineddata_sha256", "pymupdf_version") if key in installed},
                "language": "eng", "dpi": 600, "modes": [7, 10],
                "max_page": _DigitVerifier.MAX_PAGE, "max_document": _DigitVerifier.MAX_DOCUMENT,
                "max_crop_pixels": _DigitVerifier.MAX_CROP_PIXELS, "max_total_pixels": _DigitVerifier.MAX_TOTAL_PIXELS,
                "document_seconds": _DigitVerifier.DOCUMENT_SECONDS, "call_seconds": _DigitVerifier.CALL_SECONDS}}


def _font_map_record(pdf, font_xref):
    record = {"font_xref": font_xref}
    try:
        kind, value = pdf.xref_get_key(font_xref, "ToUnicode")
        if kind != "xref":
            return None
        cmap_xref = int(value.split()[0])
        record["tounicode_xref"] = cmap_xref
        data = pdf.xref_stream(cmap_xref)
        if not isinstance(data, bytes):
            raise ValueError("tounicode_is_not_a_stream")
        record["cmap_sha256"] = hashlib.sha256(data).hexdigest()
        record["mapping"] = parse_scalar_cmap(data)
        encoding_kind, encoding_value = pdf.xref_get_key(font_xref, "Encoding")
        if encoding_kind == "xref":
            encoding_xref = int(encoding_value.split()[0])
            encoding = pdf.xref_object(encoding_xref)
            record["encoding_xref"] = encoding_xref
            record["encoding_sha256"] = hashlib.sha256(encoding.encode()).hexdigest()
            differences = re.search(r"/Differences\s*\[([^\]]*)\]", encoding)
            if differences:
                observed_codes = {}
                code = None
                for token in re.findall(r"\d+|/[A-Za-z0-9_.]+", differences[1]):
                    if token.isdigit():
                        code = int(token)
                    elif code is not None:
                        # PDF /Encoding supplies quote-left/right glyphs;
                        # Docling PageAssembleModel.sanitize_text collapses
                        # both curly quotes to the single ASCII apostrophe.
                        if token in {"/quoteleft", "/quoteright"}:
                            observed_codes[code] = "'"
                        code += 1
                record["observed_codes"] = observed_codes
    except (ValueError, TypeError, RuntimeError, AttributeError, IndexError) as exc:
        # An optional source map must not break readable extraction when a
        # PDF contains a dangling xref or a non-stream ToUnicode object.
        record["error"] = f"unreadable_or_unsupported_source_map:{type(exc).__name__}"
    return record


def source_font_maps(pdf, page):
    result = {}
    for font_xref, _, subtype, name, *_ in page.get_fonts(full=True):
        if subtype != "Type1":
            continue
        record = _font_map_record(pdf, font_xref)
        if record is None:
            continue
        if name in result and result[name] != record:
            result[name] = {"error": "ambiguous_source_font_name"}
        else:
            result[name] = record
    return result


def _verify_source_sign(page, glyph, value):
    from .scientific_text import _ink, raster_agrees_with_sign
    ink = _ink(page, glyph)
    if ink is None:
        return False
    pixels = ink[0]
    if value in {"±", "−", "-", "－"}:
        return raster_agrees_with_sign(pixels, "±" if value == "±" else "−")
    # Micro-glyph restoration needs the separate numeric-unit producer's
    # descender/OCR evidence; a CMap alone cannot authorize it.
    if value in {"µ", "μ"}:
        return False
    height, width = pixels.shape
    rows = [i for i, count in enumerate(pixels.sum(axis=1)) if count >= .65*width]
    runs = []
    for row in rows:
        if not runs or row > runs[-1][-1]+1:
            runs.append([row])
        else:
            runs[-1].append(row)
    if value == "=":
        return (len(runs) == 2 and runs[1][0] - runs[0][-1] >= 2
                and pixels[rows, :].sum() >= .9*pixels.sum())
    if value == "+":
        cols = [i for i, count in enumerate(pixels.sum(axis=0)) if count >= .65*height]
        return (len(runs) == 1 and bool(cols) and .5 <= width/max(height, 1) <= 2
                and height*.2 < runs[0][0] < height*.8
                and width*.2 < min(cols) <= max(cols) < width*.8)
    return False


def recover_pdf_cmaps(document, pdf_path):
    """Materialize guarded regional font decoding before text serialization."""
    import fitz
    from docling_core.types.doc.common.meta import BaseMeta, FloatingMeta
    from .scientific_text import _source_lines
    pdf_path = Path(pdf_path)
    report = {"method": PDF_CMAP_POLICY, "producer": pdf_cmap_producer(), "repairs": [], "unresolved": []}
    with pdf_path.open("rb") as source:
        source_hash = hashlib.file_digest(source, "sha256").hexdigest()
    report["source_pdf_sha256"] = source_hash
    verifier = _DigitVerifier()
    with fitz.open(pdf_path) as pdf:
        cache = {}
        candidates = [(item, item_bounds(item, document), item.self_ref, item) for item in document.texts]
        for table in document.tables:
            if len(table.prov) != 1:
                continue
            page_no = table.prov[0].page_no
            for index, cell in enumerate(table.data.table_cells):
                if cell.bbox is not None and cell.text:
                    b = cell.bbox.to_top_left_origin(document.pages[page_no].size.height)
                    candidates.append((cell, (page_no, b.l, b.t, b.r, b.b), f"{table.self_ref}/cells/{index}", table))
        for item, box, ref, owner in candidates:
            if not box or not 1 <= box[0] <= len(pdf) or (item is owner and len(item.prov) != 1):
                continue
            notes = list(getattr(owner.meta, "corpus__pdf_cmap", []) or [])
            current_hash = hashlib.sha256(item.text.encode()).hexdigest()
            prior = next((entry for entry in notes if entry.get("method") == PDF_CMAP_POLICY
                          and entry.get("source_pdf_sha256") == source_hash and entry.get("item_ref") == ref
                          and entry.get("page") == box[0] and entry.get("source_bbox") == list(box[1:])
                          and entry.get("output_text_sha256") == current_hash
                          and entry.get("producer") == report["producer"]), None)
            if prior is not None:
                common = {key: prior[key] for key in ("item_ref", "page", "source_bbox", "source_pdf_sha256")}
                report["repairs"].extend({**common, **edit} for edit in prior["repairs"])
                report["unresolved"].extend({**common, **problem} for problem in prior["unresolved"])
                continue
            page_no = box[0]
            page = pdf[page_no-1]
            if page_no not in cache:
                cache[page_no] = (_source_lines(page), source_font_maps(pdf, page))
            lines, fonts = cache[page_no]
            if not fonts:
                continue
            # Preserve the PDF's native stream order. Sorting small raised
            # fragments by their upper edge moves unit digits before the line.
            _, left, top, right, bottom = box
            chars = [c for line in lines for c in line
                     if left-1 <= (c["bbox"][0]+c["bbox"][2])/2 <= right+1
                     and top-3 <= (c["bbox"][1]+c["bbox"][3])/2 <= bottom+3]
            repaired, edits, unresolved = repair_cmap_region(item.text, chars, fonts,
                verify_sign=lambda glyph, value: _verify_source_sign(page, glyph, value),
                verify_digit=lambda glyph: verifier.verify(page, glyph))
            common = {"item_ref": ref, "page": page_no, "source_bbox": list(box[1:]), "source_pdf_sha256": source_hash}
            report["repairs"].extend({**common, **edit} for edit in edits)
            report["unresolved"].extend({**common, **problem} for problem in unresolved)
            if edits:
                item.text = repaired
            if edits or unresolved:
                if owner.meta is None:
                    owner.meta = BaseMeta() if owner is item else FloatingMeta()
                entry = {**common, "method": PDF_CMAP_POLICY, "producer": report["producer"],
                         "input_text_sha256": current_hash,
                         "output_text_sha256": hashlib.sha256(item.text.encode()).hexdigest(),
                         "repairs": edits, "unresolved": unresolved}
                if entry not in notes:
                    notes.append(entry)
                owner.meta.corpus__pdf_cmap = notes
    return report
