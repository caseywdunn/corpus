"""Keep independently verified original exponent evidence across OCR (#303).

Numeric-unit syntax selects candidates; it never supplies a digit. Both native
geometry and rendered ink must be raised, and two unhinted OCR modes must read
the already encoded digit. Prepared text must match the same source location.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
import math
from pathlib import Path
import re
import time

SOURCE_EXPONENT_POLICY = "original-unit-digit-geometry-dual-raster-v1"
_UNIT = re.compile(r"[0-9]+(?:[.,][0-9]+)?(?:[-−–][0-9]+(?:[.,][0-9]+)?)?(?:[cmkµμ]?m|[смк]?м)$")
_SUPER = str.maketrans("23", "²³")


def source_exponent_producer():
    from .pdf_cmap_recovery import pdf_cmap_producer
    return {"policy": SOURCE_EXPONENT_POLICY, "digit_ocr": pdf_cmap_producer()["digit_ocr"],
            "ink_dpi": 600, "ink_threshold": 150, "max_ink_pixels": 40_000}


def source_candidates(chars):
    """Candidate context and size/baseline evidence, never an expected power."""
    chars = [c for c in chars if not c["text"].isspace()]
    for i, glyph in enumerate(chars):
        if glyph["text"] not in {"2", "3"} or not i:
            continue
        base = chars[i-1]
        prefix = "".join(c["text"] for c in chars[max(0, i-40):i])
        match = _UNIT.search(prefix)
        size = base["size"]
        if (not match or not size > 0 or not .35 <= glyph["size"]/size <= .7
                or not -.1*size <= glyph["bbox"][0]-base["bbox"][2] <= .6*size
                or not .25*size <= base["origin"][1]-glyph["origin"][1] <= .8*size
                or glyph["bbox"][3] > base["origin"][1]+.1*size):
            continue
        prefix_chars = chars[i-len(match[0]):i]
        if (any(abs(c["origin"][1]-base["origin"][1]) > .3*size for c in prefix_chars)
                or any(not -.2*size <= right["bbox"][0]-left["bbox"][2] <= .8*size
                       for left, right in zip(prefix_chars, prefix_chars[1:]))):
            continue
        yield {"glyph": glyph, "base": base, "prefix": match[0]}


def _ink_box(page, glyph, verifier):
    import fitz
    import numpy as np
    if not all(math.isfinite(v) for v in glyph["bbox"]):
        return {"reason": "invalid_source_box"}
    box = fitz.Rect(glyph["bbox"])
    pixels = (int(box.width*600/72)+2)*(int(box.height*600/72)+2)
    if box.is_empty or box.is_infinite or pixels <= 0 or pixels > 40_000:
        return {"reason": "invalid_or_oversized_source_box"}
    if verifier.total_pixels+pixels > verifier.MAX_TOTAL_PIXELS:
        return {"reason": "source_ink_pixel_budget_exhausted"}
    verifier.total_pixels += pixels
    pix = page.get_pixmap(clip=box, dpi=600, colorspace=fitz.csGRAY)
    data = pix.tobytes("png")
    ink = np.asarray(pix.pil_image()) < 150
    result = {"crop_sha256": hashlib.sha256(data).hexdigest(), "crop_bbox": list(box)}
    if not ink.any():
        return {**result, "reason": "source_ink_missing"}
    ys, xs = np.where(ink)
    result["ink_bbox"] = [float(v)*72/600 for v in
                          (pix.x+min(xs), pix.y+min(ys), pix.x+max(xs)+1, pix.y+max(ys)+1)]
    return result


def raised_ink_agrees(base, digit):
    """Require actual ink above its baseline neighbour, not OCR font flags."""
    if not base or not digit:
        return False
    height = base[3]-base[1]
    return (height > 0 and .25*height <= digit[3]-digit[1] <= 1.2*height
            and digit[1] <= base[1]-.25*height and digit[3] <= base[3]-.3*height
            and -.1*height <= digit[0]-base[2] <= .8*height)


def inspect_source_exponents(pdf):
    """Inspect original pages before OCR, sharing a finite raster/OCR budget."""
    from .pdf_cmap_recovery import _DigitVerifier
    from .scientific_text import _source_lines
    verifier = _DigitVerifier()
    counts = defaultdict(int)
    report = {"method": SOURCE_EXPONENT_POLICY, "candidate_count": 0,
              "confirmed_count": 0, "candidates": []}
    for page_no, page in enumerate(pdf, 1):
        chars = [c for line in _source_lines(page) for c in line]
        for candidate in source_candidates(chars):
            counts[page_no] += 1
            report["candidate_count"] += 1
            entry = {"page": page_no, "page_size": [page.rect.width, page.rect.height],
                     **candidate, "status": "unresolved"}
            if (counts[page_no] > verifier.MAX_PAGE or report["candidate_count"] > verifier.MAX_DOCUMENT
                    or time.monotonic()-verifier.started >= verifier.DOCUMENT_SECONDS):
                entry["reason"] = "source_exponent_budget_exhausted"
            else:
                base = _ink_box(page, candidate["base"], verifier)
                digit = _ink_box(page, candidate["glyph"], verifier)
                entry["source_ink"] = {"base": base, "digit": digit}
                if not raised_ink_agrees(base.get("ink_bbox"), digit.get("ink_bbox")):
                    entry["reason"] = "source_ink_not_raised"
                else:
                    entry["digit_ocr"] = verifier.verify(page, candidate["glyph"])
                    if entry["digit_ocr"]["verified"]:
                        entry["status"] = "confirmed"
                        report["confirmed_count"] += 1
                    else:
                        entry["reason"] = entry["digit_ocr"].get("reason", "digit_ocr_unconfirmed")
            report["candidates"].append(entry)
    if report["candidate_count"]:
        report["producer"] = source_exponent_producer()
    return report


def _overlap(a, b):
    import fitz
    a, b = fitz.Rect(a), fitz.Rect(b)
    return not (a & b).is_empty and (a & b).get_area() >= .5*min(a.get_area(), b.get_area())


def _owners(document, page_no, source_box):
    from .source_layout import item_bounds
    owners = []
    for item in document.texts:
        box = item_bounds(item, document)
        if box and box[0] == page_no and _overlap(box[1:], source_box):
            owners.append((item, item, item.self_ref))
    for table in document.tables:
        if len(table.prov) != 1 or table.prov[0].page_no != page_no:
            continue
        for i, cell in enumerate(table.data.table_cells):
            if cell.bbox:
                box = cell.bbox.to_top_left_origin(document.pages[page_no].size.height)
                if _overlap((box.l, box.t, box.r, box.b), source_box):
                    owners.append((cell, table, f"{table.self_ref}/cells/{i}"))
    return owners


def _prepared_alignment(chars, candidate):
    """Locate exactly one prepared glyph at the original position and prefix."""
    selected = [(i, c) for i, c in enumerate(chars)
                if _overlap(c["bbox"], candidate["glyph"]["bbox"])]
    if len(selected) != 1:
        return None
    index, glyph = selected[0]
    original = candidate["glyph"]["text"]
    if glyph["text"] not in {original, original.translate(_SUPER), "?", "'", '"', "‘", "’", "“", "”"}:
        return None
    prefix = candidate["prefix"]
    if "".join(c["text"] for c in chars[max(0, index-len(prefix)):index]) != prefix:
        return None
    suffix = "".join(c["text"] for c in chars[index+1:index+9])
    if len(suffix) < 3:
        return None
    return glyph, prefix, suffix


def apply_source_exponents(document, pdf_path, recovery):
    """Use retained source proof only at a unique prepared and structured span."""
    import fitz
    from docling_core.types.doc.common.meta import BaseMeta, FloatingMeta
    from .scientific_text import _compact, _source_lines
    report = {"repairs": [], "unresolved": []}
    evidence = recovery.get("source_exponents", {})
    if evidence.get("method") != SOURCE_EXPONENT_POLICY or not evidence.get("candidates"):
        return report
    with Path(pdf_path).open("rb") as stream:
        prepared_hash = hashlib.file_digest(stream, "sha256").hexdigest()
    cache = {}
    with fitz.open(pdf_path) as pdf:
        for candidate in evidence["candidates"]:
            page_no = candidate["page"]
            source_box = candidate["glyph"]["bbox"]
            note = {"method": SOURCE_EXPONENT_POLICY, "status": "unresolved", "page": page_no,
                    "source_bbox": source_box, "original_native_digit": candidate["glyph"]["text"],
                    "source_pdf_sha256": recovery.get("source_pdf_sha256"),
                    "prepared_pdf_sha256": prepared_hash, "producer": evidence.get("producer"),
                    "source_evidence": candidate}
            owners = _owners(document, page_no, source_box)
            if not 1 <= page_no <= len(pdf):
                note["reason"] = "prepared_page_missing"
            elif any(abs(a-b) > .5 for a, b in zip(candidate["page_size"], list(pdf[page_no-1].rect)[2:])):
                note["reason"] = "prepared_page_geometry_changed"
            elif candidate["status"] != "confirmed":
                note["reason"] = candidate.get("reason", "source_exponent_unconfirmed")
            elif len(owners) != 1:
                note["reason"] = "structured_exponent_owner_ambiguous"
            else:
                if page_no not in cache:
                    cache[page_no] = [c for line in _source_lines(pdf[page_no-1]) for c in line
                                      if not c["text"].isspace()]
                alignment = _prepared_alignment(cache[page_no], candidate)
                if alignment is None:
                    note["reason"] = "prepared_exponent_alignment_unconfirmed"
                else:
                    glyph, prefix, suffix = alignment
                    item, owner, ref = owners[0]
                    compact, offsets = _compact(item.text)
                    pattern = prefix+glyph["text"]+suffix
                    replacement = candidate["glyph"]["text"].translate(_SUPER)
                    matches = list(re.finditer(re.escape(pattern), compact))
                    if not matches and prefix+replacement+suffix in compact:
                        continue  # Already applied; keep the original receipt.
                    if len(matches) != 1:
                        note["reason"] = "structured_exponent_alignment_unconfirmed"
                    else:
                        offset = offsets[matches[0].start()+len(prefix)]
                        if item.text[offset] == replacement:
                            continue
                        note.update(status="repaired", item_ref=ref, charspan=[offset, offset+1],
                                    original=item.text[offset], replacement=replacement,
                                    prepared_bbox=glyph["bbox"])
                        item.text = item.text[:offset]+replacement+item.text[offset+1:]
            report["repairs" if note["status"] == "repaired" else "unresolved"].append(note)
            for _, owner, ref in owners:
                if owner.meta is None:
                    owner.meta = FloatingMeta() if owner in document.tables else BaseMeta()
                notes = list(getattr(owner.meta, "corpus__native_text_recovery", []) or [])
                attached = {**note, "item_ref": ref}
                if attached not in notes:
                    notes.append(attached)
                owner.meta.corpus__native_text_recovery = notes
    return report
