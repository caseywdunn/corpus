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
.refs { margin-top: 0.3em; color: #666; font-size: 0.85em; }
"""


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

    for fig in figures:
        img_rel = fig.get("filename") or ""
        img_path = f"figures/{img_rel}" if img_rel else ""
        src = fig.get("caption_source")
        src_tag = {
            "docling_caption_link": "<span class='tag docling'>docling-linked</span>",
            "heuristic_proximity": "<span class='tag heuristic'>heuristic</span>",
        }.get(src, "<span class='tag none'>no-caption</span>" if not fig.get("caption_text") else "")
        method = fig.get("extraction_method", "")
        method_tag = f"<span class='tag method'>{html.escape(method)}</span>" if method else ""

        fig_num = fig.get("figure_number") or ""
        page = fig.get("page")
        caption = fig.get("caption_text") or fig.get("caption") or ""
        refs = fig.get("referenced_in_chunks") or []

        parts.append("<div class='fig'>")
        if img_path:
            parts.append(f"<img src='{html.escape(img_path)}' alt='{html.escape(img_rel)}'>")
        else:
            parts.append("<div style='width:320px;background:#fafafa;'></div>")
        parts.append("<div class='body'>")
        parts.append(
            f"<div class='id'>{html.escape(fig['figure_id'])}"
            + (f" · Fig. {html.escape(str(fig_num))}" if fig_num else "")
            + (f" · page {page}" if page is not None else "")
            + f"</div>"
        )
        parts.append(f"<div>{src_tag}{method_tag}</div>")
        parts.append(
            f"<div class='caption'>{html.escape(caption) if caption else '<em>(no caption)</em>'}</div>"
        )
        if refs:
            parts.append(
                "<div class='refs'>referenced in "
                + ", ".join(f"<code>{html.escape(c)}</code>" for c in refs)
                + "</div>"
            )
        parts.append("</div></div>")

    parts.append("</body></html>")
    out = hash_dir / "figures_report.html"
    out.write_text("\n".join(parts), encoding="utf-8")
    return out
