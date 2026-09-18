"""Recover scientific glyphs from aligned source geometry, with raster checks.

Native PDF Unicode is evidence, not truth: legacy font encodings can map a
printed dash to ±. Mathematical sign changes therefore require the rendered
source glyph to agree. Superscripts require an actual raised source span.
"""
from __future__ import annotations

from difflib import SequenceMatcher
import re
import shutil
import subprocess

from .source_layout import item_bounds

SCIENTIFIC_TEXT_POLICY = "source_glyph_alignment_v1"
_SUPERSCRIPT = str.maketrans("0123456789+-−", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁻")
_SCIENTIFIC = set("±−µμ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")
_OLD_GLYPHS = set("+-−±µμm0123456789⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")
_LIGATURES = str.maketrans({"ﬁ": "fi", "ﬂ": "fl", "ﬀ": "ff", "ﬃ": "ffi", "ﬄ": "ffl"})


def _compact(text):
    chars, indexes = [], []
    for i, char in enumerate(text):
        for normalized in char.translate(_LIGATURES):
            if not normalized.isspace():
                chars.append(normalized)
                indexes.append(i)
    return "".join(chars), indexes


def _ink(page, glyph):
    import fitz
    import numpy as np
    rect = fitz.Rect(glyph["bbox"])
    pix = page.get_pixmap(clip=rect, dpi=600, colorspace=fitz.csGRAY)
    arr = np.asarray(pix.pil_image()) < 150
    if not arr.any():
        return None
    ys, xs = np.where(arr)
    return arr[min(ys):max(ys)+1, min(xs):max(xs)+1], (pix.y+max(ys))*72/600


def raster_agrees_with_sign(ink, sign):
    """Check horizontal bars plus the intervening vertical stem, not Unicode."""
    if ink is None or not ink.size:
        return False
    height, width = ink.shape
    rows = [i for i, count in enumerate(ink.sum(axis=1)) if count >= width*.6]
    runs = []
    for row in rows:
        if not runs or row > runs[-1][-1]+1:
            runs.append([row])
        else:
            runs[-1].append(row)
    if sign == "−":
        if len(runs) != 1:
            return False
        bar = runs[0]
        # A generous PDF font box may include the tip of a descender from
        # the preceding line. Require one dominant horizontal stroke.
        return bool(width >= 3*len(bar) and ink[bar, :].sum() >= .8*ink.sum())
    if sign != "±" or len(runs) != 2 or not (.5 <= width/max(height,1) <= 2):
        return False
    upper, lower = runs
    # A dash/date range has one bar. An equals sign has two bars but no
    # vertical stroke extending above the upper one.
    return (upper[0] > height*.2 and lower[0] > height*.65
            and bool(ink[:upper[0], int(width*.25):int(width*.75)+1].any()))


def _source_lines(page):
    import fitz
    lines = []
    for block in page.get_text("rawdict").get("blocks", []):
        for line in block.get("lines", []):
            chars = []
            for span in line["spans"]:
                for char in span["chars"]:
                    box = fitz.Rect(char["bbox"])*page.rotation_matrix
                    origin = fitz.Point(char["origin"])*page.rotation_matrix
                    chars.append({"text": char["c"], "bbox": list(box),
                                  "origin": [origin.x,origin.y], "font": span["font"],
                                  "size": span["size"], "raised": bool(span["flags"] & 1)})
            if chars:
                lines.append(chars)
    return lines


def _characters_in_box(lines, box):
    _, left, top, right, bottom = box
    selected = []
    for line in lines:
        chars = [c for c in line if left-1 <= (c["bbox"][0]+c["bbox"][2])/2 <= right+1
                 and top-3 <= (c["bbox"][1]+c["bbox"][3])/2 <= bottom+3]
        if chars:
            selected.append(chars)
    selected.sort(key=lambda line: (round(min(c["bbox"][1] for c in line)/3), min(c["bbox"][0] for c in line)))
    return [char for line in selected for char in line]


def _micro_ocr(page, glyph, next_glyph):
    """Corroborate an encoded m-shaped unit with a descender and rendered OCR."""
    import fitz
    executable = shutil.which("tesseract")
    if not executable:
        return False
    box = fitz.Rect(glyph["bbox"]) | fitz.Rect(next_glyph["bbox"])
    box += (-2,-2,2,2)
    pix = page.get_pixmap(clip=box, dpi=600)
    try:
        result = subprocess.run([executable, "stdin", "stdout", "--psm", "7"],
                                input=pix.tobytes("png"), capture_output=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0 and bool(re.search(r"[uµμ]m", re.sub(r"\s+", "", result.stdout.decode("utf-8", errors="replace"))))


def _native_text(page, chars):
    text, evidence = [], []
    compact = "".join(c["text"] for c in chars)
    for i, char in enumerate(chars):
        value = char["text"]
        method = None
        verified = False
        if value in {"±", "−"}:
            ink = _ink(page, char)
            verified = raster_agrees_with_sign(ink[0] if ink else None, value)
            method = "aligned_native_and_rendered_sign"
        elif value == "-" and char["raised"]:
            ink = _ink(page, char)
            verified = raster_agrees_with_sign(ink[0] if ink else None, "−")
            method = "aligned_native_and_rendered_sign"
        elif value == "+":
            # OCR text layers can themselves reduce a printed ± to +. The
            # two bars and upper stem must still be present in the source.
            ink = _ink(page, char)
            if raster_agrees_with_sign(ink[0] if ink else None, "±"):
                value, method, verified = "±", "rendered_plusminus_from_encoded_plus", True
        if char["raised"] and value in "0123456789+-−":
            # The geometry is independent of the potentially broken encoding.
            # A raised minus still has to render as a horizontal stroke.
            if value not in "+-−" or value == "+" or verified:
                value = value.translate(_SUPERSCRIPT)
                method, verified = "raised_source_span", True
        if value == "m" and i+1 < len(chars) and chars[i+1]["text"] == "m" and chars[i+1]["font"] != char["font"]:
            # Restrict to an actual numeric unit, never prose like "mammal".
            if re.search(r"\d\s*$", compact[:i]):
                ink = _ink(page, char)
                other = _ink(page, chars[i+1])
                if (ink and other and ink[1]-char["origin"][1] > .1*char["size"]
                        and other[1]-chars[i+1]["origin"][1] < .07*char["size"]):
                    value, method = "µ", "source_descender_font_and_rendered_ocr"
                    verified = _micro_ocr(page, char, chars[i+1])
        for part in value.translate(_LIGATURES):
            if not part.isspace():
                text.append(part)
                evidence.append({**char, "method": method, "verified": verified})
    return "".join(text), evidence


def repair_aligned_text(text, native, evidence):
    """Apply only narrow, independently evidenced symbol differences."""
    old, indexes = _compact(text)
    matcher = SequenceMatcher(None, old, native, autojunk=False)
    if not old or sum(block.size for block in matcher.get_matching_blocks()) < .85*min(len(old), len(native)):
        return text, [], []
    edits, unresolved = [], []
    candidates = []
    # Local anchors also locate repeated caption/body representations inside
    # one extraction item. Each copy must retain the same source sign (#303).
    for match in re.finditer("["+re.escape("".join(_SCIENTIFIC))+"]+", native):
        c,d = match.span()
        for size in (12, 8, 5, 3):
            left, right = native[max(0,c-size):c], native[d:d+size]
            if not left or not right:
                continue
            pattern = re.compile(re.escape(left)+"(?P<glyph>["+re.escape("".join(_OLD_GLYPHS))+r"]{0,4})"+re.escape(right))
            found = list(pattern.finditer(old))
            if found:
                candidates.extend((*m.span("glyph"), c,d) for m in found if m["glyph"] != native[c:d])
    replacements_by_span = {}
    for a,b,c,d in candidates:
        replacements_by_span.setdefault((a,b), set()).add(native[c:d])
    seen = set()
    for a,b,c,d in candidates:
        if (a,b) in seen:
            continue
        seen.add((a,b))
        if len(replacements_by_span[(a,b)]) != 1:
            unresolved.append({"before": old[a:b], "reason": "ambiguous_source_alignment"})
            continue
        proof = evidence[c:d]
        record = {"before": old[a:b], "after": native[c:d],
                  "source_bbox": [min(e["bbox"][0] for e in proof), min(e["bbox"][1] for e in proof),
                                  max(e["bbox"][2] for e in proof), max(e["bbox"][3] for e in proof)],
                  "evidence": sorted({e["method"] for e in proof if e["method"]})}
        if not all(e["verified"] for e in proof):
            unresolved.append({**record, "reason": "native_glyph_not_confirmed_by_source_raster"})
            continue
        start = indexes[a] if a < len(indexes) else len(text)
        end = indexes[b-1]+1 if b > a else start
        replacement = native[c:d]
        if replacement[0] in "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻":
            while start > 0 and text[start-1].isspace():
                start -= 1
        if replacement == "µ" and d < len(native) and native[d] == "m":
            while end < len(text) and text[end].isspace():
                end += 1
        edits.append({**record, "charspan": [start,end], "original": text[start:end], "replacement": replacement})
    edits.sort(key=lambda edit: edit["charspan"])
    # Conflicting/overlapping alignments are not safe to apply.
    overlaps = {i for i in range(len(edits)-1)
                if edits[i]["charspan"][1] > edits[i+1]["charspan"][0]}
    if overlaps:
        bad = overlaps | {i+1 for i in overlaps}
        unresolved.extend({**edits[i], "reason": "overlapping_source_alignment"} for i in sorted(bad))
        edits = [edit for i,edit in enumerate(edits) if i not in bad]
    result = text
    for edit in reversed(edits):
        start,end = edit["charspan"]
        result = result[:start]+edit["replacement"]+result[end:]
    return result, edits, unresolved


def prepare_scientific_text(document, pdf_path):
    """Repair structured text and table cells before all exports/embeddings."""
    import fitz
    from docling_core.types.doc.common.meta import BaseMeta, FloatingMeta
    report = {"method": SCIENTIFIC_TEXT_POLICY, "repairs": [], "unresolved": []}
    cache = {}
    with fitz.open(pdf_path) as pdf:
        candidates = []
        for item in document.texts:
            box = item_bounds(item, document)
            if box and len(item.prov) == 1:
                candidates.append((item, box, item.self_ref, item))
        for table in document.tables:
            if len(table.prov) != 1:
                continue
            page_no = table.prov[0].page_no
            for i,cell in enumerate(table.data.table_cells):
                if cell.bbox is not None and cell.text:
                    b = cell.bbox.to_top_left_origin(document.pages[page_no].size.height)
                    candidates.append((cell, (page_no,b.l,b.t,b.r,b.b), f"{table.self_ref}/cells/{i}", table))
        for item, box, ref, owner in candidates:
            page_no = box[0]
            if not 1 <= page_no <= len(pdf):
                continue
            if page_no not in cache:
                cache[page_no] = _source_lines(pdf[page_no-1])
            chars = _characters_in_box(cache[page_no], box)
            if not chars:
                continue
            native, evidence = _native_text(pdf[page_no-1], chars)
            repaired, edits, unresolved = repair_aligned_text(item.text, native, evidence)
            for entry in edits:
                report["repairs"].append({"item_ref": ref, "page": page_no, **entry})
            for entry in unresolved:
                report["unresolved"].append({"item_ref": ref, "page": page_no, **entry})
            if edits:
                # Docling's orig remains untouched; receipts retain exact old
                # substrings with their pre-repair offsets for table cells too.
                item.text = repaired
            if edits or unresolved:
                if owner.meta is None:
                    owner.meta = FloatingMeta() if owner is not item else BaseMeta()
                notes = list(getattr(owner.meta, "corpus__scientific_text", []))
                notes.append({"item_ref": ref, "page": page_no, "repairs": edits,
                              "unresolved": unresolved, "method": SCIENTIFIC_TEXT_POLICY})
                owner.meta.corpus__scientific_text = notes
    return report
