"""Conservative source-geometry repairs before text export and chunking.

PDF text order is not reading order. Reorder only pages with independently
supported, non-overlapping prose columns; keep a receipt of every moved item.
"""
from __future__ import annotations

from collections import defaultdict
from statistics import median
import re

SOURCE_LAYOUT_POLICY = "source_columns_table_captions_v2"

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


def recover_table_caption_owners(document):
    """Move an explicit table caption from a distant picture to its adjacent table.

    Require a unique, uncaptioned table immediately below the printed caption,
    with matching left edge and no intervening text. Existing table bindings,
    competing captions and uncertain/multipage geometry are left untouched.
    The source text and boxes are never changed (#320).
    """
    proposals = defaultdict(list)
    for caption in document.texts:
        if (_label(caption) != "caption" or len(caption.prov) != 1
                or not re.match(r"^Table\s+\d+\b", caption.text, re.I)
                or caption.parent is None):
            continue
        owner = caption.parent.resolve(document)
        box = item_bounds(caption, document)
        old = item_bounds(owner, document)
        if (_label(owner) != "picture" or len(owner.prov) != 1
                or box is None or old is None or old[0] != box[0]
                or not any(r.cref == caption.self_ref for r in owner.captions)):
            continue
        height = box[4] - box[2]
        width = box[3] - box[1]
        if height <= 0 or width <= 0:
            continue
        gap_limit = min(20, 2 * height)
        # A caption inside/near its current picture is ambiguous, even if a
        # nearby table also seems plausible.
        if old[2] - box[4] <= gap_limit:
            continue
        candidates = []
        for table in document.tables:
            tb = item_bounds(table, document)
            if (len(table.prov) != 1 or tb is None or tb[0] != box[0]
                    or not 0 <= tb[2] - box[4] <= gap_limit
                    or abs(tb[1] - box[1]) > height
                    or tb[3] < box[3] - height or tb[4] > old[2]):
                continue
            blockers = [item for item in document.texts
                        if item.self_ref != caption.self_ref
                        and (ib := item_bounds(item, document)) is not None
                        and ib[0] == box[0] and box[4] <= ib[2] < tb[2]
                        and min(ib[3], tb[3]) > max(ib[1], tb[1])]
            if not blockers:
                candidates.append(table)
        if len(candidates) == 1 and not candidates[0].captions:
            proposals[candidates[0].self_ref].append((caption, owner, candidates[0]))
    changed = []
    for candidates in proposals.values():
        if len(candidates) != 1:
            continue
        caption, owner, table = candidates[0]
        ref = caption.get_ref()
        owner.captions = [r for r in owner.captions if r.cref != caption.self_ref]
        owner.children = [r for r in owner.children if r.cref != caption.self_ref]
        table.captions.append(ref)
        if not any(r.cref == caption.self_ref for r in table.children):
            table.children.append(ref)
        caption.parent = table.get_ref()
        changed.append({"caption_ref": caption.self_ref, "previous_owner": owner.self_ref,
                        "table_ref": table.self_ref, "page": caption.prov[0].page_no,
                        "reason": "explicit_table_caption_adjacent_unique_table"})
    return changed


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
