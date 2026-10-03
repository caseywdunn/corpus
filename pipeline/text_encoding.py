"""Source-corroborated PDF encoding and accent recovery (#306, #312, #316).

The PDF's character geometry is independent evidence, not an unconditional
replacement text layer. Big5-compatible output is replaced only when the
same region's Unicode text re-encodes back to almost all of the observed
bytes. Spacing accents require overlap with the particular base glyph.
Unconfirmed candidates remain unchanged and carry review evidence.
"""
from __future__ import annotations

import re
import unicodedata
from collections import OrderedDict
from difflib import SequenceMatcher

TEXT_ENCODING_POLICY = "source-region-encoding-and-accent-geometry-v1"
_ACCENTS = {"´": "\u0301", "`": "\u0300", "˜": "\u0303", "~": "\u0303",
            "¨": "\u0308", "˚": "\u030a", "ˆ": "\u0302", "ˇ": "\u030c",
            "¸": "\u0327"}
_UTF8_DAMAGE = re.compile(r"(?:Ã\s*[\x80-\xbf]|Â\s*[\x80-\xbf])")


def _cjk(char):
    return "\u3400" <= char <= "\u9fff" or "\U00020000" <= char <= "\U0002ffff"


def _big5_candidate(text):
    # Ordinary Latin-1 prose has occasional accents, not dense byte pairs.
    high = sum(0x80 <= ord(c) <= 0xff for c in text)
    return high >= 4 and high / max(len(text), 1) >= 0.30


def _source_lines(page):
    return [line["spans"] for block in page.get_text("rawdict").get("blocks", [])
            for line in block.get("lines", [])]


def _region_chars(lines, box):
    x0, y0, x1, y1 = box
    out = []
    for spans in lines:
        chars = [c for span in spans for c in span.get("chars", [])
                 if x0 - .8 <= (c["bbox"][0] + c["bbox"][2]) / 2 <= x1 + .8
                 and y0 - .8 <= (c["bbox"][1] + c["bbox"][3]) / 2 <= y1 + .8]
        if chars:
            if out:
                out.append({"c": "\n", "bbox": [0, 0, 0, 0]})
            out.extend(chars)
    return out


def _indexed_compact(chars):
    text, indexes = [], []
    for i, char in enumerate(chars):
        for normalized in unicodedata.normalize("NFKD", "i" if char == "ı" else char):
            if not normalized.isspace():
                text.append(normalized)
                indexes.append(i)
    return "".join(text), indexes


def _accent_edits(text, chars):
    """Map an overlapping accent/base pair to its exact extracted offsets."""
    native, source_indexes = _indexed_compact([c["c"] for c in chars])
    target, target_indexes = _indexed_compact(text)
    aligned = {}
    for block in SequenceMatcher(None, native, target, autojunk=False).get_matching_blocks():
        for offset in range(block.size):
            aligned[source_indexes[block.a + offset]] = target_indexes[block.b + offset]
    edits, unresolved = [], []
    for i, glyph in enumerate(chars):
        accent = glyph["c"]
        if accent not in _ACCENTS:
            continue
        candidates = []
        # Restrict to neighboring letters in this source line. Whitespace
        # alone is insufficient: a literal spacing symbol does not overlap.
        for j in (i - 1, i + 1):
            if not 0 <= j < len(chars) or not chars[j]["c"].isalpha():
                continue
            base = chars[j]
            ax0, ay0, ax1, ay1 = glyph["bbox"]
            bx0, by0, bx1, by1 = base["bbox"]
            overlap = max(0, min(ax1, bx1) - max(ax0, bx0))
            if overlap < .6 * max(ax1 - ax0, .01):
                continue
            if not bx0 <= (ax0 + ax1) / 2 <= bx1:
                continue
            if max(0, min(ay1, by1) - max(ay0, by0)) < .5 * min(ay1 - ay0, by1 - by0):
                continue
            letter = "i" if base["c"] == "ı" else base["c"]
            composed = unicodedata.normalize("NFC", letter + _ACCENTS[accent])
            if len(composed) == 1 and composed != letter:
                candidates.append((j, composed))
        if len(candidates) != 1:
            continue
        j, composed = candidates[0]
        if i in aligned and j in aligned and aligned[i] == aligned[j] and text[aligned[i]] == composed:
            continue  # The same source pair was already composed on a prior pass.
        if i not in aligned or j not in aligned:
            unresolved.append({"reason": "accent_source_alignment_unconfirmed", "source_bbox": glyph["bbox"]})
            continue
        start, end = sorted((aligned[i], aligned[j]))
        end += 1
        if "".join(c for c in text[start:end] if not c.isspace()) not in (chars[j]["c"] + accent, accent + chars[j]["c"]):
            unresolved.append({"reason": "accent_source_alignment_unconfirmed", "source_bbox": glyph["bbox"]})
            continue
        edits.append({"charspan": [start, end], "original": text[start:end],
                      "replacement": composed, "source_bbox": glyph["bbox"],
                      "base_bbox": chars[j]["bbox"], "evidence": "overlapping_native_accent_and_base_glyph"})
    edits.sort(key=lambda e: e["charspan"])
    overlap = {i for i in range(len(edits) - 1) if edits[i]["charspan"][1] > edits[i + 1]["charspan"][0]}
    blocked = overlap | {i + 1 for i in overlap}
    unresolved.extend({"reason": "overlapping_accent_repairs", **edits[i]} for i in sorted(blocked))
    return [edit for i, edit in enumerate(edits) if i not in blocked], unresolved


def repair_region(text, chars):
    """Return repaired text and evidence; never infer missing name letters."""
    native = "".join(c["c"] for c in chars)
    edits, unresolved = [], []
    if _big5_candidate(text):
        confirmed = False
        if sum(_cjk(c) for c in native) >= 2:
            for encoding in ("big5", "cp950"):
                try:
                    expected = native.encode(encoding).decode("latin-1")
                except UnicodeError:
                    continue
                # Do not NFKC-normalize bytes: compatibility mappings such as
                # superscript digits would destroy their encoding identity.
                old = re.sub(r"\s+", "", text)
                source = re.sub(r"\s+", "", expected)
                alignment = SequenceMatcher(None, old, source, autojunk=False)
                matches = sum(m.size for m in alignment.get_matching_blocks())
                # In the byte-decoding failure Docling also clips line tails
                # to fit the doubled glyph count. Recover these only when
                # every observed byte agrees in order with this exact region;
                # substitutions/deletions do not get the relaxed source gate.
                clipped_tails = all(op[0] in ("equal", "insert") for op in alignment.get_opcodes())
                if not clipped_tails or matches / max(len(source), 1) < .65:
                    continue
                replacement = " ".join(native.split())
                edits.append({"charspan": [0, len(text)], "original": text,
                              "replacement": replacement, "evidence": "same_region_native_unicode_reencodes_to_observed_bytes",
                              "encoding": encoding, "observed_coverage": matches / len(old),
                              "source_coverage": matches / len(source),
                              "source_only_insertions": clipped_tails and matches < len(source)})
                text = replacement
                confirmed = True
                break
        if not confirmed:
            # Dense Latin-1 is merely a candidate. Require successful CJK
            # decoding before labeling an uncorroborated region suspicious.
            for encoding in ("big5", "cp950"):
                try:
                    decoded = text.encode("latin-1").decode(encoding)
                except UnicodeError:
                    continue
                if sum(_cjk(c) for c in decoded) >= 2:
                    unresolved.append({"reason": "possible_big5_bytes_without_matching_source_unicode", "original": text})
                    break
    if not edits:
        edits, accent_unresolved = _accent_edits(text, chars)
        unresolved.extend(accent_unresolved)
        for edit in reversed(edits):
            start, end = edit["charspan"]
            text = text[:start] + edit["replacement"] + text[end:]
    if _UTF8_DAMAGE.search(text):
        unresolved.append({"reason": "possible_utf8_mojibake_requires_source_review", "original": text})
    return text, edits, unresolved


def recover_text_encoding(document, pdf_path):
    """Repair text/cell content before exports; retain original Docling text."""
    import fitz
    from docling_core.types.doc.common.meta import BaseMeta, FloatingMeta

    report = {"method": TEXT_ENCODING_POLICY, "repairs": [], "unresolved": []}
    with fitz.open(pdf_path) as pdf:
        cache = OrderedDict()
        candidates = []
        for item in document.texts:
            if len(item.prov) == 1:
                prov = item.prov[0]
                candidates.append((item, prov.page_no, prov.bbox, item.self_ref, item))
        for table in document.tables:
            if len(table.prov) == 1:
                for i, cell in enumerate(table.data.table_cells):
                    if cell.bbox is not None and cell.text:
                        candidates.append((cell, table.prov[0].page_no, cell.bbox, f"{table.self_ref}/cells/{i}", table))
        for item, page_no, bbox, ref, owner in candidates:
            if not (_big5_candidate(item.text) or _UTF8_DAMAGE.search(item.text) or any(c in item.text for c in _ACCENTS)):
                continue
            if not 1 <= page_no <= len(pdf):
                report["unresolved"].append({"item_ref": ref, "page": page_no, "reason": "source_page_unavailable"})
                continue
            page = pdf[page_no - 1]
            if page_no not in cache:
                cache[page_no] = _source_lines(page)
                if len(cache) > 4:
                    cache.popitem(last=False)
            cache.move_to_end(page_no)
            b = bbox.to_top_left_origin(page.rect.height)
            chars = _region_chars(cache[page_no], (b.l, b.t, b.r, b.b))
            repaired, edits, unresolved = repair_region(item.text, chars)
            for entry in edits:
                report["repairs"].append({"item_ref": ref, "page": page_no, **entry})
            for entry in unresolved:
                report["unresolved"].append({"item_ref": ref, "page": page_no, **entry})
            if edits:
                item.text = repaired
            if edits or unresolved:
                if owner.meta is None:
                    owner.meta = BaseMeta() if owner is item else FloatingMeta()
                notes = list(getattr(owner.meta, "corpus__text_encoding", []))
                note = {"item_ref": ref, "page": page_no, "repairs": edits,
                        "unresolved": unresolved, "method": TEXT_ENCODING_POLICY}
                already_recorded = not edits and all(
                    any(problem in prior.get("unresolved", []) for prior in notes
                        if prior.get("method") == TEXT_ENCODING_POLICY and prior.get("item_ref") == ref)
                    for problem in unresolved)
                if note not in notes and not already_recorded:
                    notes.append(note)
                owner.meta.corpus__text_encoding = notes
    return report
