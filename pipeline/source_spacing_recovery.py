"""Retain source spaces lost by full-page OCR (#334).

Only exact native letters, nearby page anchors and two unhinted source-crop
readings can authorize insertion. Neither OCR nor a dictionary supplies words.
The input is the same kept-page PDF passed to OCR, so page ordinals map one to
one; the preparation receipt separately identifies original library pages.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
from pathlib import Path
import re
import time

from .source_spaces import corroborate_source_gaps, source_spacing_producer

SOURCE_SPACING_RECOVERY_POLICY = "original-prepared-exact-spaces-v1"
_LETTERS = re.compile(r"[^\W\d_]+", re.UNICODE)


def source_spacing_recovery_producer():
    return {**source_spacing_producer(), "policy": SOURCE_SPACING_RECOVERY_POLICY,
            "max_candidates_per_document": 32, "document_seconds": 60,
            "min_letters": 10, "max_letters": 64, "nearby_anchors": 2}


def _digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _letter_lines(page):
    """Keep punctuation barriers and actual glyph bounds when finding words."""
    import fitz
    lines = []
    for block in page.get_text("rawdict", flags=fitz.TEXTFLAGS_RAWDICT & ~fitz.TEXT_PRESERVE_IMAGES)["blocks"]:
        for line in block.get("lines", []):
            if tuple(line.get("dir", (1, 0))) != (1, 0):
                continue
            chars = [char for span in line.get("spans", []) for char in span.get("chars", [])]
            text = "".join(char["c"] for char in chars)
            # MuPDF normally emits one Unicode codepoint per glyph. A ligature
            # expanded into several codepoints has no unambiguous offset map.
            if len(text) != len(chars):
                continue
            words = []
            for match in _LETTERS.finditer(text):
                boxes = [char["bbox"] for char in chars[match.start():match.end()]]
                words.append({"text": match[0], "span": match.span(), "bbox": _union(boxes)})
            lines.append((text, words))
    return lines


def _union(boxes):
    return [min(b[0] for b in boxes), min(b[1] for b in boxes),
            max(b[2] for b in boxes), max(b[3] for b in boxes)]


def _aligned(source, prepared):
    import fitz
    a, b = fitz.Rect(source), fitz.Rect(prepared)
    overlap = a & b
    return (a.get_area() > 0 and b.get_area() > 0
            and overlap.get_area() >= .8 * min(a.get_area(), b.get_area())
            and .7 <= a.width / max(b.width, .01) <= 1.4
            and abs((a.x0+a.x1-b.x0-b.x1)/2) <= 5
            and abs((a.y0+a.y1-b.y0-b.y1)/2) <= 5)


def _candidates(source_page, prepared_page, producer):
    source_lines = _letter_lines(source_page)
    prepared_lines = _letter_lines(prepared_page)
    source_words = [word for _, words in source_lines for word in words]
    prepared_words = [word for _, words in prepared_lines for word in words]
    groups = defaultdict(list)
    for text, words in source_lines:
        for start in range(len(words)):
            for size in range(2, 5):
                selected = words[start:start+size]
                if len(selected) != size:
                    continue
                if not all(text[a["span"][1]:b["span"][0]].isspace()
                           for a, b in zip(selected, selected[1:])):
                    continue
                original = "".join(w["text"] for w in selected)
                if producer["min_letters"] <= len(original) <= producer["max_letters"]:
                    groups[original].append({"text": " ".join(w["text"] for w in selected),
                                             "bbox": _union([w["bbox"] for w in selected])})
    result = []
    for word in prepared_words:
        matches = [group for group in groups.get(word["text"], [])
                   if _aligned(group["bbox"], word["bbox"])]
        if not matches:
            continue
        candidate = {"original": word["text"], "prepared_bbox": word["bbox"]}
        if len(matches) != 1:
            result.append({**candidate, "reason": "ambiguous_original_phrase"})
            continue
        group = matches[0]
        candidate.update(source_bbox=group["bbox"], replacement=group["text"])
        anchors = []
        for anchor in source_words:
            if len(anchor["text"]) < 4 or abs(anchor["bbox"][1]-group["bbox"][1]) > 30:
                continue
            # Anchors lie outside the proposed repair, on either side nearby.
            if not (0 <= group["bbox"][0]-anchor["bbox"][2] <= 100
                    or 0 <= anchor["bbox"][0]-group["bbox"][2] <= 100):
                continue
            located = [p for p in prepared_words if p["text"] == anchor["text"]
                       and _aligned(anchor["bbox"], p["bbox"])]
            if len(located) == 1:
                anchors.append({"text": anchor["text"], "source_bbox": anchor["bbox"],
                                "prepared_bbox": located[0]["bbox"]})
        if len({a["text"] for a in anchors}) < producer["nearby_anchors"]:
            candidate["reason"] = "insufficient_page_alignment_anchors"
        candidate["anchors"] = anchors
        result.append(candidate)
    return result


def inspect_source_spacing(source_pdf, prepared_pdf, original_pages=None):
    """Measure original/prepared differences without changing either PDF."""
    import fitz
    producer = source_spacing_recovery_producer()
    report = {"method": SOURCE_SPACING_RECOVERY_POLICY, "producer": producer,
              "source_pdf_sha256": _digest(source_pdf),
              "prepared_pdf_sha256": _digest(prepared_pdf),
              "candidate_count": 0, "repairs": [], "unresolved": []}
    deadline = time.monotonic() + producer["document_seconds"]
    attempted = 0
    with fitz.open(source_pdf) as source, fitz.open(prepared_pdf) as prepared:
        if len(source) != len(prepared) or (original_pages and len(original_pages) != len(source)):
            report["unresolved"].append({"reason": "page_mapping_changed"})
            return report
        for index, (page, output_page) in enumerate(zip(source, prepared)):
            if time.monotonic() >= deadline or attempted >= producer["max_candidates_per_document"]:
                report["unresolved"].append({"page": index+1, "reason": "unexamined_pages_budget_exceeded",
                                             "unexamined_page_count": len(source)-index})
                break
            mapping = {"page": index+1, "original_page": original_pages[index] if original_pages else index+1}
            # Raw text boxes are unrotated while pixmap clips are rendered-page
            # coordinates. Refuse rotated pages until that transform is modeled.
            if (page.rotation or output_page.rotation or page.cropbox != output_page.cropbox
                    or page.rect != output_page.rect):
                report["unresolved"].append({**mapping, "reason": "page_geometry_changed"})
                continue
            candidates = _candidates(page, output_page, producer)
            report["candidate_count"] += len(candidates)
            for number, candidate in enumerate(candidates):
                note = {**mapping, **candidate}
                if "reason" in note:
                    report["unresolved"].append(note)
                    continue
                remaining = deadline-time.monotonic()
                if (number >= producer["max_candidates_per_page"]
                        or attempted >= producer["max_candidates_per_document"] or remaining < 1):
                    report["unresolved"].append({**note, "reason": "candidate_budget_exceeded"})
                    continue
                attempted += 1
                original, replacement = candidate["original"], candidate["replacement"]
                boundaries = [len(replacement[:m.start()].replace(" ", ""))
                              for m in re.finditer(" ", replacement)]
                proposal = {"text": original, "bbox": candidate["source_bbox"], "choices": [
                    {"original": original, "candidate": replacement, "start": 0,
                     "end": len(original), "boundaries": boundaries}]}
                _, observations = corroborate_source_gaps(page, [proposal], {
                    **producer, "timeout_seconds": min(producer["timeout_seconds"], remaining/2)})
                evidence = observations[0]
                evidence["route"] = SOURCE_SPACING_RECOVERY_POLICY
                evidence["proposal_evidence"] = "original_native_whitespace"
                note["crop_evidence"] = evidence
                if evidence["status"] == "verified":
                    report["repairs"].append(note)
                else:
                    report["unresolved"].append({**note, "reason": evidence["status"]})
    return report


def apply_source_spacing(document, prepared_pdf, receipt):
    """Insert confirmed spaces in one geometrically unique structured owner."""
    import fitz
    from docling_core.types.doc.common.meta import BaseMeta, FloatingMeta
    from .source_layout import item_bounds

    report = {"method": SOURCE_SPACING_RECOVERY_POLICY, "repairs": [], "unresolved": []}
    if not receipt or receipt.get("method") != SOURCE_SPACING_RECOVERY_POLICY:
        return report
    report.update({key: receipt[key] for key in ("producer", "source_pdf_sha256", "prepared_pdf_sha256")})
    report["unresolved"] = list(receipt.get("unresolved", []))
    if _digest(prepared_pdf) != receipt["prepared_pdf_sha256"]:
        report["unresolved"].append({"reason": "prepared_pdf_identity_changed"})
        return report
    for candidate in receipt["repairs"]:
        page_no, target = candidate["page"], fitz.Rect(candidate["prepared_bbox"])
        owners = []
        for item in document.texts:
            bounds = item_bounds(item, document)
            if bounds and bounds[0] == page_no and (fitz.Rect(bounds[1:]) & target).get_area() >= .8*target.get_area():
                owners.append((item, item, None))
        for table in document.tables:
            if len(table.prov) != 1 or table.prov[0].page_no != page_no:
                continue
            for index, cell in enumerate(table.data.table_cells):
                if cell.bbox:
                    box = cell.bbox.to_top_left_origin(document.pages[page_no].size.height)
                    if (fitz.Rect(box.l, box.t, box.r, box.b) & target).get_area() >= .8*target.get_area():
                        owners.append((cell, table, index))
        pattern = re.compile(r"(?<!\w)"+re.escape(candidate["original"])+r"(?!\w)")
        matches = [(item, owner, cell, match) for item, owner, cell in owners
                   for match in pattern.finditer(item.text)]
        if len(matches) != 1:
            if not matches and any(candidate["replacement"] in item.text for item, _, _ in owners):
                continue  # Applying the same receipt twice is a no-op.
            report["unresolved"].append({**candidate, "reason": "structured_owner_ambiguous_or_missing"})
            continue
        item, owner, cell, match = matches[0]
        note = {**candidate, "item_ref": owner.self_ref, "cell_index": cell,
                "charspan": list(match.span()), "method": SOURCE_SPACING_RECOVERY_POLICY,
                "source_pdf_sha256": receipt["source_pdf_sha256"],
                "prepared_pdf_sha256": receipt["prepared_pdf_sha256"], "status": "repaired"}
        item.text = item.text[:match.start()]+candidate["replacement"]+item.text[match.end():]
        if owner.meta is None:
            owner.meta = BaseMeta() if owner is item else FloatingMeta()
        owner.meta.corpus__source_spaces = [*(getattr(owner.meta, "corpus__source_spaces", []) or []), note]
        report["repairs"].append(note)
    return report
