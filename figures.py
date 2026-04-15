"""Figure + caption joint-object utilities and the figures_report.html
generator.

Phase D (see PLAN.md §3 "Figure+caption as a first-class object"). Three
concerns live here:

1. :func:`extract_caption_info` — pull caption text, page, and bbox for a
   docling Picture, first via docling's own caption linker, then via a
   proximity-based heuristic fallback over the document's TextItems, so
   figures without an explicit docling caption link still get one.
2. :func:`parse_figure_number` — best-effort extraction of the numeric
   identifier from a caption string ("Figure 3. ..." → "3"). Handles the
   common English / German / French / Russian forms.
3. :func:`link_chunks_to_figures` — post-chunking cross-reference pass
   that scans chunk text for "Fig. N" mentions and records bidirectional
   links (``chunk["figure_refs"]``, ``figure["referenced_in_chunks"]``).

Plus a small HTML report generator that a human can open to QC the whole
set of figures+captions for one paper at a glance.
"""

from __future__ import annotations

import html
import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


_PLATE_CAPTION_RE = re.compile(
    r"^\s*(?:plate|pl\.?|tafel|planche)\s*\d+",
    re.IGNORECASE,
)


# Multilingual "figure" prefix. Plate/Abbildung/Рисунок etc. all map to the
# same reference namespace — figure_refs_in_chunks just tracks the number.
_FIGURE_REF_RE = re.compile(
    r"""
    \b
    (?:fig(?:ure|\.?)|abb(?:ildung|\.?)|pl(?:ate|\.?)|plate|рис(?:унок|\.?))   # prefix
    \s*
    (\d+)                                                                        # number
    """,
    re.IGNORECASE | re.VERBOSE,
)

# Tighter version used only on caption text (the *leading* figure label).
_FIGURE_NUMBER_IN_CAPTION_RE = re.compile(
    r"^\s*(?:fig(?:ure|\.?)|abb(?:ildung|\.?)|pl(?:ate|\.?)|plate|рис(?:унок|\.?))\s*(\d+)",
    re.IGNORECASE,
)


def parse_figure_number(caption_text: str) -> Optional[str]:
    """Extract the figure number from the start of a caption string.

    Returns ``"3"`` for ``"Figure 3. Nectophore of Nanomia cara ..."``.
    Returns None if the caption doesn't start with a recognizable prefix —
    some captions (e.g., "Plate 2.") still match because we cover plate/
    pl/Abb/Pl/Рис. Best-effort, intended as a join key for body-text
    references, not as a canonical identifier.
    """
    if not caption_text:
        return None
    m = _FIGURE_NUMBER_IN_CAPTION_RE.match(caption_text)
    return m.group(1) if m else None


def _prov_to_bbox_and_page(prov_item) -> Tuple[Optional[List[float]], Optional[int]]:
    """Normalize a docling ProvenanceItem to (bbox_list, page_no) or Nones."""
    page_no = getattr(prov_item, "page_no", None)
    bbox = getattr(prov_item, "bbox", None)
    bbox_list: Optional[List[float]] = None
    if bbox is not None:
        try:
            bbox_list = [float(bbox.l), float(bbox.b), float(bbox.r), float(bbox.t)]
        except Exception as e:
            logger.debug("Could not serialize docling bbox: %s", e)
    return bbox_list, page_no


# ---------------------------------------------------------------------------
# Figure classification and deduplication
# ---------------------------------------------------------------------------

# PDF point heuristics. A full US Letter page is ~612×792 pts ≈ 485k pts²;
# A4 is ~595×842 pts ≈ 501k pts². We default to ~500k so these thresholds
# describe reasonable fractions of a page even without knowing the exact
# size — override ``page_area_pts`` when you do.
_DEFAULT_PAGE_AREA_PTS = 500_000
_MIN_FIGURE_AREA_FRAC = 0.03   # below 3% of page area + no caption → graphical
_MIN_FIGURE_DIM_PTS = 50       # below 50pt in either dim → graphical regardless


def _bbox_area(bbox: Optional[List[float]]) -> float:
    if not bbox or len(bbox) != 4:
        return 0.0
    return max(0.0, bbox[2] - bbox[0]) * max(0.0, bbox[3] - bbox[1])


def _bbox_overlap_fraction(a: List[float], b: List[float]) -> float:
    """Intersection area / smaller bbox area. Returns 0 when either bbox is
    missing or has zero area. Used to decide whether two same-numbered
    figures are redundant crops of the same panel or genuinely distinct
    subpanels."""
    if not (a and b and len(a) == 4 and len(b) == 4):
        return 0.0
    ix0 = max(a[0], b[0])
    iy0 = max(a[1], b[1])
    ix1 = min(a[2], b[2])
    iy1 = min(a[3], b[3])
    if ix1 <= ix0 or iy1 <= iy0:
        return 0.0
    inter = (ix1 - ix0) * (iy1 - iy0)
    area_a = _bbox_area(a)
    area_b = _bbox_area(b)
    smaller = min(area_a, area_b) if area_a and area_b else 0.0
    if smaller == 0:
        return 0.0
    return inter / smaller


# Figure-type strings emitted in figures.json. Consumers (MCP tools, the
# HTML report, downstream analytic queries) filter on these.
FIGURE_TYPE_FIGURE = "figure"
FIGURE_TYPE_PLATE = "plate"
FIGURE_TYPE_SUBPANEL = "subpanel"
FIGURE_TYPE_GRAPHICAL = "graphical_element"
FIGURE_TYPE_UNCLASSIFIED = "unclassified"


def classify_figure(item: Dict, page_area_pts: float = _DEFAULT_PAGE_AREA_PTS) -> str:
    """Classify one raw figure item into a ``figure_type``.

    ``item`` is the in-memory dict the two-pass extractor builds per picture:
    it carries ``bbox``, ``caption_text``, and ``figure_number``. The
    classifier never touches the saved image — it can be re-run on
    ``figures.json`` after the fact.

    Rules (first match wins):

    * **graphical_element** — bbox smaller than 50×50 pts OR below 3% of
      page area with no caption. Matches the journal-header/footer/logo
      furniture docling loves to extract from modern PDFs.
    * **plate** — caption starts with "Plate N" / "Pl. N" / "Tafel N" /
      "Planche N". Historical papers' non-figure non-furniture content
      lives here.
    * **figure** — has a caption AND a parseable figure number AND
      reasonable size.
    * **unclassified** — reasonable size but no caption or no number.
      Kept but flagged for the HTML report's manual-review section.
    """
    bbox = item.get("bbox")
    w = (bbox[2] - bbox[0]) if bbox and len(bbox) == 4 else 0
    h = (bbox[3] - bbox[1]) if bbox and len(bbox) == 4 else 0
    area = _bbox_area(bbox)
    caption = item.get("caption_text") or ""
    fig_num = item.get("figure_number")

    if w and h and (w < _MIN_FIGURE_DIM_PTS or h < _MIN_FIGURE_DIM_PTS):
        return FIGURE_TYPE_GRAPHICAL
    if area and area < _MIN_FIGURE_AREA_FRAC * page_area_pts and not caption:
        return FIGURE_TYPE_GRAPHICAL

    if _PLATE_CAPTION_RE.match(caption):
        return FIGURE_TYPE_PLATE

    if caption and fig_num:
        return FIGURE_TYPE_FIGURE

    return FIGURE_TYPE_UNCLASSIFIED


def _sort_by_reading_order(items: List[Dict]) -> List[Dict]:
    """Sort figure panels by approximate reading order — top row first,
    left-to-right within row.

    docling bboxes are PDF bottom-left origin, so ``bbox[3]`` (``t``, top)
    has the highest y on the page. We bucket top-y to the nearest 30 pts
    so panels in the same visual row sort together despite sub-pixel
    jitter from caption-box alignment, then break ties by ``bbox[0]``
    (left edge). This matches how panels are labelled A/B/C/D in
    scientific figures — A is top-left, B is top-right (or below A if
    vertical), C follows in the next row, etc.
    """
    def key(p):
        bb = p.get("bbox") or [0, 0, 0, 0]
        l, b, r, t = bb
        row_bucket = -round(t / 30) * 30  # negative → higher y sorts first
        return (row_bucket, l)
    return sorted(items, key=key)


def dedupe_figures(items: List[Dict]) -> List[Dict]:
    """Merge redundant crops, label multi-panel figures by position.

    Two modes per figure-number group:

    * **Coequal panels** (no one bbox encompasses the others) — the common
      PLOS-style case where Figure 3 is really 4 separate panel images
      labelled A/B/C/D. After removing overlap-duplicate crops, we sort
      the survivors by :func:`_sort_by_reading_order` and assign
      ``panel_letter = "a", "b", "c", …`` accordingly. Each panel keeps
      ``figure_type = "figure"`` — they're all real figures. The filenames
      become ``fig_3_a.png``, ``fig_3_b.png``, … (see
      :func:`compose_figure_filename`).
    * **Whole-figure + subpanels** — one bbox encompasses the others
      (overlap ≥ 80% of each sibling's area). The encompassing item is
      the primary (``figure_type = "figure"``, filename ``fig_3.png``);
      the contained items are reclassified as
      ``figure_type = "subpanel"`` with position-derived panel_letters.

    Items classified upstream as ``graphical_element`` or
    ``unclassified`` are **not** grouped by figure-number, even if the
    caption-proximity heuristic accidentally attached a figure label to
    them. This prevents the Siebert 2011 case where a 19×19 PLOS badge
    near Figure 1's caption ended up filed as ``fig_1_panel_b.png``.

    Returns a new list. Input order is preserved for items that weren't
    grouped; grouped items appear in reading order.
    """
    if not items:
        return items
    from collections import defaultdict

    by_num: Dict[str, List[Dict]] = defaultdict(list)
    passthrough: List[Dict] = []
    for it in items:
        num = it.get("figure_number")
        ftype = it.get("figure_type")
        # Only real figures and plates participate in dedup + panel
        # assignment. Graphical elements that happened to latch onto a
        # figure_number via the caption-proximity heuristic keep their
        # downgraded classification.
        if num and str(num) and ftype in (FIGURE_TYPE_FIGURE, FIGURE_TYPE_PLATE):
            by_num[str(num)].append(it)
        else:
            passthrough.append(it)

    kept: List[Dict] = list(passthrough)

    for num, group in by_num.items():
        if len(group) == 1:
            kept.append(group[0])
            continue

        # Step 1: remove overlap-duplicate crops. Sort largest-area first
        # and keep an item only if it doesn't overlap > 50% with any
        # already-kept item. This eliminates the whole-figure + per-panel-
        # crop duplication docling sometimes emits.
        sorted_group = sorted(group, key=lambda x: -_bbox_area(x.get("bbox")))
        unique: List[Dict] = []
        for cand in sorted_group:
            if any(
                _bbox_overlap_fraction(cand.get("bbox"), u.get("bbox")) > 0.5
                for u in unique
            ):
                logger.debug(
                    "Figure dedup: dropping docling_%s (Fig %s) "
                    "as overlap-duplicate",
                    cand.get("docling_idx"), num,
                )
                continue
            unique.append(cand)

        if len(unique) == 1:
            kept.append(unique[0])
            continue

        # Step 2: detect the whole-figure overview case. One bbox
        # encompasses every sibling (overlap ≥ 80% measured against the
        # sibling's bbox, so the overview's extra padding doesn't matter).
        largest = max(unique, key=lambda x: _bbox_area(x.get("bbox")))
        non_largest = [u for u in unique if u is not largest]
        encompassed = sum(
            1 for u in non_largest
            if _bbox_overlap_fraction(largest.get("bbox"), u.get("bbox")) > 0.8
        )
        has_whole = non_largest and encompassed >= len(non_largest)

        if has_whole:
            # Primary = largest (no panel letter); others are subpanels
            # with position-derived letters.
            ordered_panels = _sort_by_reading_order(non_largest)
            kept.append(largest)
            for idx, panel in enumerate(ordered_panels):
                panel["figure_type"] = FIGURE_TYPE_SUBPANEL
                panel["primary_figure_docling_idx"] = largest.get("docling_idx")
                panel["panel_letter"] = chr(ord("a") + idx)
                kept.append(panel)
        else:
            # Coequal panels. Sort by reading order, assign a, b, c, d, …
            # All keep figure_type = "figure".
            ordered = _sort_by_reading_order(unique)
            for idx, panel in enumerate(ordered):
                panel["panel_letter"] = chr(ord("a") + idx)
            kept.extend(ordered)

    return kept


# ---------------------------------------------------------------------------
# File-naming policy
# ---------------------------------------------------------------------------

_NAME_SAFE_RE = re.compile(r"[^a-z0-9]+")


def _slug(s: str) -> str:
    return _NAME_SAFE_RE.sub("_", (s or "").lower()).strip("_")


def compose_figure_filename(item: Dict) -> str:
    """Return a semantic filename for a classified figure item.

    Filename convention matches how the paper labels its figures:

    * Single-item figure → ``fig_{number}.png``
    * Multi-panel coequal figure (Fig 3 A/B/C/D) → ``fig_{number}_{letter}.png``
      (each panel is its own image, ``figure_type = "figure"``)
    * Whole-figure + contained subpanels:
      - whole → ``fig_{number}.png`` (``figure_type = "figure"``, no letter)
      - subpanels → ``fig_{number}_{letter}.png`` (``figure_type = "subpanel"``)
    * ``plate`` with number → ``plate_{number}.png`` (+ ``_{letter}`` if panelled)
    * ``figure`` / ``plate`` missing number → ``fig_docling_{idx}.png``
    * ``graphical_element`` → ``_graphic_{idx}.png`` (leading underscore
      sorts these separately from real content in a directory listing)
    * ``unclassified`` → ``_unclassified_{idx}.png``

    A collision-guard in the caller appends ``_docling_{idx}`` when a
    composed filename is already taken.
    """
    ftype = item.get("figure_type") or FIGURE_TYPE_UNCLASSIFIED
    idx = item.get("docling_idx")
    num = item.get("figure_number")
    letter = item.get("panel_letter")

    if ftype == FIGURE_TYPE_FIGURE:
        if num and letter:
            return f"fig_{num}_{letter}.png"
        return f"fig_{num}.png" if num else f"fig_docling_{idx}.png"
    if ftype == FIGURE_TYPE_PLATE:
        if num and letter:
            return f"plate_{num}_{letter}.png"
        return f"plate_{num}.png" if num else f"plate_docling_{idx}.png"
    if ftype == FIGURE_TYPE_SUBPANEL:
        if num and letter:
            return f"fig_{num}_{letter}.png"
        return f"_subpanel_{idx}.png"
    if ftype == FIGURE_TYPE_GRAPHICAL:
        return f"_graphic_{idx}.png"
    return f"_unclassified_{idx}.png"


def extract_caption_info(picture, document) -> Dict:
    """Return a dict with caption_text, caption_page, caption_bbox,
    bbox_coord_system, caption_source for a docling Picture.

    Tries in order:

    1. ``picture.captions`` — docling's own caption linker; when it
       succeeds, the TextItem it points to has the exact caption string and
       its own provenance (page + bbox). This is the good path.
    2. Proximity heuristic — scan all TextItems on the same page as the
       picture for one whose text starts with a figure-label prefix
       ("Figure", "Fig.", "Abb.", "Pl.", "Рис."), then pick the one whose
       bbox is vertically closest to the picture bbox. This catches the
       common case where docling's layout model didn't associate the
       caption with the picture object (plate-heavy monographs, end-matter
       figure sections, etc.).
    3. Nothing — return empty caption fields so downstream code can tell
       the caption genuinely wasn't found rather than silently blank.

    ``caption_source`` is always one of ``"docling_caption_link"``,
    ``"heuristic_proximity"``, or ``None``.
    """
    info: Dict = {
        "caption_text": "",
        "caption_page": None,
        "caption_bbox": None,
        "bbox_coord_system": None,
        "caption_source": None,
    }

    # --- Path 1: docling's own caption linker ---
    caps = getattr(picture, "captions", None)
    if caps:
        try:
            target = caps[0].resolve(document)
            text = (target.text or "").strip()
            if text:
                info["caption_text"] = text
                if getattr(target, "prov", None):
                    bbox_list, page_no = _prov_to_bbox_and_page(target.prov[0])
                    info["caption_page"] = page_no
                    info["caption_bbox"] = bbox_list
                    info["bbox_coord_system"] = "pdf_pts_bottom_left"
                info["caption_source"] = "docling_caption_link"
                return info
        except Exception as e:
            logger.debug("docling caption resolve failed: %s", e)

    # --- Path 2: proximity heuristic over same-page TextItems ---
    pic_prov = getattr(picture, "prov", None)
    if not pic_prov:
        return info
    pic_bbox_list, pic_page = _prov_to_bbox_and_page(pic_prov[0])
    if pic_page is None or pic_bbox_list is None:
        return info

    candidates: List[Tuple[float, object, List[float], int]] = []
    for text_item in getattr(document, "texts", []) or []:
        text = getattr(text_item, "text", "") or ""
        if not text:
            continue
        if not _FIGURE_NUMBER_IN_CAPTION_RE.match(text):
            continue
        prov = getattr(text_item, "prov", None)
        if not prov:
            continue
        bbox_list, page_no = _prov_to_bbox_and_page(prov[0])
        if page_no != pic_page or bbox_list is None:
            continue
        # Vertical distance: figure bottom to caption top (captions are
        # usually below the figure). In bottom-left coords, figure bottom is
        # bbox[1] (b) and caption top is bbox[3] (t). Use absolute gap.
        pic_top, pic_bottom = pic_bbox_list[3], pic_bbox_list[1]
        cap_top, cap_bottom = bbox_list[3], bbox_list[1]
        vertical_gap = min(
            abs(pic_bottom - cap_top),  # caption below figure
            abs(cap_bottom - pic_top),  # caption above figure (rare)
        )
        candidates.append((vertical_gap, text_item, bbox_list, page_no))

    if not candidates:
        return info
    candidates.sort(key=lambda c: c[0])
    _, best, best_bbox, best_page = candidates[0]
    info["caption_text"] = (best.text or "").strip()
    info["caption_page"] = best_page
    info["caption_bbox"] = best_bbox
    info["bbox_coord_system"] = "pdf_pts_bottom_left"
    info["caption_source"] = "heuristic_proximity"
    return info


# ---------------------------------------------------------------------------
# Chunk ↔ figure cross-reference pass
# ---------------------------------------------------------------------------


def link_chunks_to_figures(
    chunks: List[Dict],
    figures: List[Dict],
) -> Tuple[List[Dict], List[Dict]]:
    """Post-chunking bidirectional link between body-text figure mentions
    and figure objects.

    Scans each chunk's text for "Fig. N" / "Figure N" / "Pl. N" / etc.
    mentions, looks up figures whose parsed figure number matches, and
    records the link in both directions:

    * ``chunks[i]["figure_refs"]`` — list of ``figure_id`` referenced in
      that chunk.
    * ``figures[j]["referenced_in_chunks"]`` — list of ``chunk_id`` that
      mention that figure.

    Multiple figures can share a figure_number in the same doc (e.g.,
    "Fig. 1A" vs "Fig. 1B" or plate reuse); in that case the link is
    many-to-many within the matching number. Figures without a parseable
    number simply never match and come out with an empty
    ``referenced_in_chunks`` — accurate, and distinguishable from zero
    references to a parseable-number figure.

    Returns the mutated lists (they are also modified in place).
    """
    # Build figure_number -> [figure_id, ...]
    number_to_figure_ids: Dict[str, List[str]] = {}
    for fig in figures:
        num = fig.get("figure_number")
        if not num:
            continue
        number_to_figure_ids.setdefault(str(num), []).append(fig["figure_id"])
        # Initialize bucket if missing
        fig.setdefault("referenced_in_chunks", [])

    if not number_to_figure_ids:
        # Nothing to link; still initialize the fields on everything so
        # output schema is stable.
        for fig in figures:
            fig.setdefault("referenced_in_chunks", [])
        for ch in chunks:
            ch.setdefault("figure_refs", [])
        return chunks, figures

    for ch in chunks:
        text = ch.get("text", "") or ""
        seen_here: set = set()
        for m in _FIGURE_REF_RE.finditer(text):
            num = m.group(1)
            for fid in number_to_figure_ids.get(num, []):
                if fid in seen_here:
                    continue
                seen_here.add(fid)
        refs = sorted(seen_here)
        ch["figure_refs"] = refs
        for fid in refs:
            # Append chunk_id to figure's referenced_in_chunks
            for fig in figures:
                if fig["figure_id"] == fid:
                    lst = fig.setdefault("referenced_in_chunks", [])
                    cid = ch.get("chunk_id")
                    if cid and cid not in lst:
                        lst.append(cid)
                    break
    # Ensure every figure has the field (even if no refs)
    for fig in figures:
        fig.setdefault("referenced_in_chunks", [])
    return chunks, figures


# ---------------------------------------------------------------------------
# HTML QC report
# ---------------------------------------------------------------------------


_REPORT_CSS = """
body { font-family: -apple-system, system-ui, sans-serif; margin: 2em; color: #222; }
header { border-bottom: 1px solid #ddd; padding-bottom: 0.5em; margin-bottom: 1em; }
h1 { margin: 0 0 0.2em; font-size: 1.4em; }
h2 { margin: 2em 0 0.5em; font-size: 1.1em; color: #555; border-bottom: 1px solid #eee; padding-bottom: 0.3em; }
.meta { color: #666; font-size: 0.9em; }
.fig { display: flex; gap: 1em; margin: 1.5em 0; padding: 0.5em; border: 1px solid #eee; border-radius: 4px; }
.fig img { max-width: 320px; max-height: 300px; object-fit: contain; border: 1px solid #eee; background: #fafafa; }
.fig .body { flex: 1; min-width: 0; }
.fig .id { font-family: ui-monospace, monospace; font-size: 0.85em; color: #888; }
.fig .caption { margin-top: 0.25em; white-space: pre-wrap; }
.fig .tag { display: inline-block; font-size: 0.75em; padding: 0.1em 0.5em; border-radius: 4px; margin-right: 0.3em; }
.tag.docling { background: #e0f0ff; color: #03538b; }
.tag.heuristic { background: #fff3d6; color: #8a5a00; }
.tag.none { background: #ffe0e0; color: #8b0303; }
.tag.method { background: #eee; color: #555; }
.tag.type-figure { background: #e0ffe0; color: #035e03; }
.tag.type-plate  { background: #e8e0ff; color: #4a1f94; }
.tag.type-subpanel { background: #fff0d6; color: #6e4200; }
.tag.type-graphical_element { background: #eee; color: #555; }
.tag.type-unclassified { background: #ffe0e0; color: #8b0303; }
.refs { margin-top: 0.3em; color: #666; font-size: 0.85em; }
details > summary { cursor: pointer; color: #666; font-size: 0.95em; margin: 1em 0 0.5em; }
"""


# Ordered for rendering — real figures first, review bucket last.
_REPORT_SECTIONS = [
    ("Figures", {FIGURE_TYPE_FIGURE}),
    ("Plates", {FIGURE_TYPE_PLATE}),
    ("Subpanels", {FIGURE_TYPE_SUBPANEL}),
    ("Needs review (graphical elements / unclassified)",
     {FIGURE_TYPE_GRAPHICAL, FIGURE_TYPE_UNCLASSIFIED}),
]


def generate_figures_report(hash_dir: Path) -> Optional[Path]:
    """Write ``<hash_dir>/figures_report.html`` — a human-readable QC surface.

    Side-by-side thumbnails and captions for every figure, with source
    tags so it's easy to spot captions that came from the heuristic
    fallback (most likely to need manual review) and figures with no
    caption at all. Also shows which chunks reference each figure, when
    the cross-ref pass has run.

    Returns the path written, or None if figures.json was missing.
    """
    figures_file = hash_dir / "figures.json"
    if not figures_file.exists():
        return None
    try:
        figures_data = json.load(figures_file.open(encoding="utf-8"))
    except Exception as e:
        logger.warning("Could not read %s: %s", figures_file, e)
        return None

    figures = figures_data.get("figures", [])
    metadata = {}
    metadata_file = hash_dir / "metadata.json"
    if metadata_file.exists():
        try:
            metadata = json.load(metadata_file.open(encoding="utf-8"))
        except Exception:
            pass

    title = metadata.get("title") or metadata.get("filename") or hash_dir.name
    year = metadata.get("year")
    authors = ", ".join(
        f"{a.get('forename','').strip()} {a.get('surname','').strip()}".strip()
        for a in metadata.get("authors", []) or []
    )

    # Build HTML
    parts = [
        "<!doctype html>",
        "<html><head><meta charset='utf-8'>",
        f"<title>Figures — {html.escape(str(title))}</title>",
        f"<style>{_REPORT_CSS}</style>",
        "</head><body>",
        "<header>",
        f"<h1>{html.escape(str(title))}</h1>",
        f"<div class='meta'>{html.escape(authors)}"
        + (f" · {year}" if year else "")
        + f" · hash <code>{hash_dir.name}</code>"
        + f" · {len(figures)} figure(s)</div>",
        "</header>",
    ]

    # Group figures by type so real content comes first and the review
    # bucket (graphical elements + unclassified) can be collapsed.
    by_type: Dict[str, List[Dict]] = {}
    for fig in figures:
        by_type.setdefault(fig.get("figure_type") or FIGURE_TYPE_UNCLASSIFIED, []).append(fig)

    def _render_figure(fig: Dict) -> List[str]:
        img_rel = fig.get("filename") or ""
        img_path = f"figures/{img_rel}" if img_rel else ""
        src = fig.get("caption_source")
        src_tag = {
            "docling_caption_link": "<span class='tag docling'>docling-linked</span>",
            "heuristic_proximity": "<span class='tag heuristic'>heuristic</span>",
        }.get(src, "<span class='tag none'>no-caption</span>" if not fig.get("caption_text") else "")
        ftype = fig.get("figure_type") or FIGURE_TYPE_UNCLASSIFIED
        type_tag = f"<span class='tag type-{html.escape(ftype)}'>{html.escape(ftype)}</span>"

        fig_num = fig.get("figure_number") or ""
        page = fig.get("page")
        caption = fig.get("caption_text") or fig.get("caption") or ""
        refs = fig.get("referenced_in_chunks") or []

        out: List[str] = ["<div class='fig'>"]
        if img_path:
            out.append(f"<img src='{html.escape(img_path)}' alt='{html.escape(img_rel)}'>")
        else:
            out.append("<div style='width:320px;background:#fafafa;'></div>")
        out.append("<div class='body'>")
        # id line: figure_id · Fig. 3 · panel_b · page 5
        id_parts = [html.escape(fig.get("figure_id", ""))]
        if fig_num:
            id_parts.append(f"Fig. {html.escape(str(fig_num))}")
        if fig.get("panel_letter"):
            id_parts.append(f"panel {html.escape(fig['panel_letter'])}")
        if page is not None:
            id_parts.append(f"page {page}")
        out.append(f"<div class='id'>{' · '.join(id_parts)}</div>")
        out.append(f"<div>{type_tag}{src_tag}</div>")
        out.append(
            f"<div class='caption'>{html.escape(caption) if caption else '<em>(no caption)</em>'}</div>"
        )
        if refs:
            out.append(
                "<div class='refs'>referenced in "
                + ", ".join(f"<code>{html.escape(c)}</code>" for c in refs)
                + "</div>"
            )
        out.append("</div></div>")
        return out

    for heading, type_set in _REPORT_SECTIONS:
        figs = [f for t in type_set for f in by_type.get(t, [])]
        if not figs:
            continue
        # The "Needs review" section starts collapsed — it's for humans
        # to skim quickly, not central to the report's primary content.
        is_review = "Needs review" in heading
        if is_review:
            parts.append(f"<details><summary><h2 style='display:inline'>{html.escape(heading)} "
                         f"({len(figs)})</h2></summary>")
        else:
            parts.append(f"<h2>{html.escape(heading)} ({len(figs)})</h2>")
        for fig in figs:
            parts.extend(_render_figure(fig))
        if is_review:
            parts.append("</details>")

    parts.append("</body></html>")
    out = hash_dir / "figures_report.html"
    out.write_text("\n".join(parts), encoding="utf-8")
    return out
