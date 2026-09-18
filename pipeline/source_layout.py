"""Conservative source-geometry repairs before text export and chunking.

PDF text order is not reading order. Reorder only pages with independently
supported, non-overlapping prose columns; keep a receipt of every moved item.
"""
from __future__ import annotations

from collections import defaultdict
from statistics import median
import re

SOURCE_LAYOUT_POLICY = "source_columns_v1"

def _label(item):
    value = getattr(item, "label", "")
    return getattr(value, "value", str(value))


def item_bounds(item, document):
    """Return (page, left, top, right, bottom) in source page points."""
    prov = getattr(item, "prov", [])
    if not prov:
        children = [item_bounds(ref.resolve(document), document)
                    for ref in getattr(item, "children", [])]
        children = [box for box in children if box]
        if not children:
            return None
        # A spanning group retains its internal order. Its first source page
        # locates the group in the outer walk; later provenance is untouched.
        first_page = min(box[0] for box in children)
        children = [box for box in children if box[0] == first_page]
        return (children[0][0], min(b[1] for b in children), min(b[2] for b in children),
                max(b[3] for b in children), max(b[4] for b in children))
    p = prov[0]
    page = document.pages.get(p.page_no)
    if page is None:
        return None
    box = p.bbox.to_top_left_origin(page.size.height)
    return p.page_no, box.l, box.t, box.r, box.b


def repair_reading_order(document):
    """Repair evidenced column inversions without guessing ambiguous layouts.

    Headings and paragraph fragments stay in their printed column, including
    narrow text wrapping around a figure. Figure/table records follow prose on
    their source page so their own captions cannot interrupt sentence flow.
    """
    report = {"method": SOURCE_LAYOUT_POLICY, "reordered_pages": [], "unresolved": []}
    by_page = defaultdict(list)
    for ref in document.body.children:
        item = ref.resolve(document)
        box = item_bounds(item, document)
        if box is None:
            # Unknown geometry must not be silently placed into a treatment.
            report["unresolved"].append({"item_ref": item.self_ref, "reason": "missing_or_multipage_geometry"})
            return report
        by_page[box[0]].append((ref, item, box))
    ordered = []
    for page_no, entries in sorted(by_page.items()):
        prose = [entry for entry in entries if _label(entry[1]) not in
                 {"picture", "table", "caption", "page_header", "page_footer", "footnote"}]
        candidates = [entry for entry in prose if
                      _label(entry[1]) == "text" and len(getattr(entry[1], "text", "")) >= 60
                      and entry[2][3] - entry[2][1] >= 65]
        clusters = []
        for entry in sorted(candidates, key=lambda e: e[2][1]):
            x = entry[2][1]
            if clusters and abs(x - median(e[2][1] for e in clusters[-1])) < 12:
                clusters[-1].append(entry)
            else:
                clusters.append([entry])
        supported = [c for c in clusters if len(c) >= 2]
        starts = [median(e[2][1] for e in c) for c in supported]
        if len(starts) < 2 or any(b-a < 90 for a,b in zip(starts, starts[1:])):
            ordered.extend(entry[0] for entry in entries)
            continue
        # Every candidate column needs prose confined before the next one.
        if any(sum(e[2][3] <= starts[i+1]-4 for e in c) < 2
               for i,c in enumerate(supported[:-1])):
            ordered.extend(entry[0] for entry in entries)
            report["unresolved"].append({"page": page_no, "reason": "overlapping_prose_columns"})
            continue
        column_tops = [min(e[2][2] for e in c if i == len(starts)-1 or e[2][3] <= starts[i+1]-4)
                       for i,c in enumerate(supported)]
        spanning = [e for e in prose if e[2][1] <= starts[0]+12 and
                    e[2][3] >= starts[-1]+50 and e[2][2] > max(column_tops)+5]
        if spanning:
            ordered.extend(entry[0] for entry in entries)
            report["unresolved"].append({"page": page_no, "reason": "midpage_spanning_prose"})
            continue
        def order_key(entry):
            box = entry[2]
            col = min(range(len(starts)), key=lambda i: abs(box[1]-starts[i]))
            return col, box[2], box[1]
        sorted_prose = sorted(prose, key=order_key)
        extras = [e for e in entries if e not in prose]
        result = sorted_prose + extras
        old_refs = [e[1].self_ref for e in entries]
        new_refs = [e[1].self_ref for e in result]
        if old_refs != new_refs:
            report["reordered_pages"].append({"page": page_no, "column_starts": starts,
                                               "before": old_refs, "after": new_refs})
        ordered.extend(e[0] for e in result)
    if report["reordered_pages"]:
        document.body.children = ordered
    return report


def recover_panel_caption_roles(document):
    """Recover explicit (A), (B), ... continuation blocks directly under a caption."""
    from docling_core.types.doc import DocItemLabel
    changed = []
    items = [(item, item_bounds(item, document)) for item in document.texts]
    items = [(item, box) for item, box in items if box]
    for caption, box in items:
        if _label(caption) != "caption" or not re.match(r"(?:Fig(?:ure)?\.?|Plate)\s*\d", caption.text, re.I):
            continue
        bottom = box[4]
        for item, child in sorted(items, key=lambda pair: pair[1][2]):
            if child[0] != box[0] or child[2] < bottom-1:
                continue
            if child[2]-bottom > 16:
                break
            if _label(item) == "text" and abs(child[1]-box[1]) < 3 and re.match(
                    r"(?:\([A-Z](?:[–-][A-Z])?\)|See (?:Fig(?:ure)?\.?|Table)\b)", item.text):
                item.label = DocItemLabel.CAPTION
                changed.append({"item_ref": item.self_ref, "caption_ref": caption.self_ref,
                                "page": box[0], "reason": "aligned_explicit_panel_continuation"})
                bottom = child[4]
    return changed
