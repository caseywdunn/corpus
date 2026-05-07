"""Figure-keyed MCP tools.

Surfaces: get_figures_for_taxon, get_figures_for_lexicon_term,
get_figure, list_figure_rois, get_figure_roi_image, get_figure_image.
The image-returning tools wrap PIL crops in FastMCP's ``Image``
content type.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from mcp.server.fastmcp import Image

from ..app import _load_json, _need_index, mcp


# Restricted to the figure types that get returned by get_figures_for_*.
# Excludes graphical_element, plate_label, furniture, etc.
_REAL_FIGURE_TYPES = {"figure", "plate", "subpanel"}


@mcp.tool()
def get_figures_for_taxon(
    taxon_name: str,
    paper_hash: Optional[str] = None,
    limit: int = 50,
    include_all: bool = False,
) -> List[Dict]:
    """Figures from papers that mention the taxon, ranked by caption
    relevance.

    A figure whose caption names the taxon directly scores higher than
    a figure from a paper that merely mentions it elsewhere. Each
    result carries the on-disk figure image path.

    By default only returns items classified as ``figure`` or ``plate``
    (skipping journal furniture, subpanels of already-returned figures,
    and unclassifiable graphical elements). Pass ``include_all=True`` to
    see every extracted item including the review bucket.
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return [{"error": "no taxonomy snapshot configured"}]
    hit = idx.taxonomy_db.lookup(taxon_name)
    if not hit:
        return []
    aid = hit["accepted_taxon_id"]
    accepted_name_low = (hit["accepted_name"] or "").lower()
    matched_name_low = (hit["matched_name"] or "").lower()
    target_hashes = (
        [paper_hash] if paper_hash else idx.taxon_to_papers.get(aid, [])
    )

    rows: List[Dict] = []
    for h in target_hashes:
        p = idx.papers.get(h)
        if not p:
            continue
        figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
        for f in figs.get("figures", []) or []:
            ftype = f.get("figure_type")
            if not include_all and ftype not in _REAL_FIGURE_TYPES:
                continue
            caption = (f.get("caption_text") or f.get("caption") or "").lower()
            caption_hit = accepted_name_low in caption or (
                matched_name_low and matched_name_low in caption
            )
            rows.append({
                "paper_hash": h,
                "paper_title": p.get("title"),
                "figure_id": f.get("figure_id"),
                "figure_type": ftype,
                "page": f.get("page"),
                "caption_text": f.get("caption_text") or f.get("caption"),
                "figure_number": f.get("figure_number"),
                # Relative to the corpuscle's documents/ dir.
                # Call get_figure_image to fetch bytes.
                "image_path": f"{h}/figures/{f.get('filename') or ''}",
                "caption_has_taxon": caption_hit,
                "score": (100 if caption_hit else 0) + idx.taxon_mention_counts.get(aid, {}).get(h, 0),
            })
    rows.sort(key=lambda r: -r["score"])
    return rows[: int(limit)] if limit else rows


@mcp.tool()
def get_figures_for_lexicon_term(
    category: str,
    term: str,
    paper_hash: Optional[str] = None,
    limit: int = 50,
    include_all: bool = False,
) -> List[Dict]:
    """Figures whose captions mention a lexicon term, ranked by
    caption occurrence count.

    ``category`` is a top-level key in the corpus's lexicon
    (``anatomy``, ``biogeography``, …). ``term`` is matched
    case-insensitively as a substring; pass the canonical name or any
    declared synonym/translation. Returns real figures + plates by
    default; ``include_all=True`` includes the review bucket.
    """
    idx = _need_index()
    term_low = (term or "").strip().lower()
    if not term_low:
        return []
    # ``category`` is currently used for documentation + future filtering
    # by lexicon-membership. The match itself is caption-text search.
    _ = category
    target_hashes = (
        [paper_hash] if paper_hash else list(idx.papers.keys())
    )
    rows: List[Dict] = []
    for h in target_hashes:
        p = idx.papers.get(h)
        if not p:
            continue
        figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
        for f in figs.get("figures", []) or []:
            ftype = f.get("figure_type")
            if not include_all and ftype not in _REAL_FIGURE_TYPES:
                continue
            caption = (f.get("caption_text") or f.get("caption") or "")
            occ = caption.lower().count(term_low)
            if occ == 0:
                continue
            rows.append({
                "paper_hash": h,
                "paper_title": p.get("title"),
                "figure_id": f.get("figure_id"),
                "figure_type": ftype,
                "page": f.get("page"),
                "figure_number": f.get("figure_number"),
                "caption_text": caption,
                # Relative to the corpuscle's documents/ dir.
                # Call get_figure_image to fetch bytes.
                "image_path": f"{h}/figures/{f.get('filename') or ''}",
                "match_count": occ,
            })
    rows.sort(key=lambda r: -r["match_count"])
    return rows[: int(limit)] if limit else rows



@mcp.tool()
def get_figure(paper_hash: str, figure_id: str) -> Dict:
    """One figure's full record: caption, page, bbox, image path,
    cross-references."""
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return {"error": f"no such paper_hash: {paper_hash}"}
    figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
    for f in figs.get("figures", []) or []:
        if f.get("figure_id") == figure_id:
            return {
                **f,
                "paper_hash": paper_hash,
                "paper_title": p.get("title"),
                # Relative to the corpuscle's documents/ dir.
                "image_path": f"{paper_hash}/figures/{f.get('filename') or ''}",
            }
    return {"error": f"no such figure_id {figure_id!r} in paper {paper_hash}"}



@mcp.tool()
def list_figure_rois(paper_hash: str, figure_id: str) -> List[Dict]:
    """Return the per-panel / per-subfigure ROIs annotated on a figure.

    ROIs are populated by Pass 2.5 (caption-derived ``panels_from_caption``)
    and Pass 3a (OCR-derived ``rois`` with pixel bboxes). The caption-
    derived list gives labels + descriptions even when Pass 3a wasn't
    run or OCR didn't find the panel; ``rois`` gives pixel coordinates
    when available for :func:`get_figure_roi_image`.
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return [{"error": f"no such paper_hash: {paper_hash}"}]
    figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
    for f in figs.get("figures", []) or []:
        if f.get("figure_id") == figure_id:
            return [{
                "paper_hash": paper_hash,
                "figure_id": figure_id,
                "figure_number": f.get("figure_number"),
                "filename": f.get("filename"),
                "panels_from_caption": f.get("panels_from_caption") or [],
                "panel_count_from_caption": f.get("panel_count_from_caption", 0),
                "rois": f.get("rois") or [],
                "pass3_status": f.get("pass3_status"),
                "image_size_px": f.get("image_size_px"),
            }]
    return [{"error": f"no such figure_id {figure_id!r} in paper {paper_hash}"}]


@mcp.tool()
def get_figure_roi_image(
    paper_hash: str,
    figure_id: str,
    label: str,
) -> Dict:
    """Crop a panel ROI out of a figure image and return the crop's path.

    ``label`` is the panel letter (e.g. ``"A"``, ``"B"``) as stored in
    ``figures.json`` ``rois[*].label``. If Pass 3a found a pixel ROI
    for that label, we crop the figure PNG to that region and cache the
    result under ``<hash_dir>/figures/crops/{figure_id}__{label}.png``
    so repeated asks are free.

    Returns the crop path + caption description. If the label exists in
    ``panels_from_caption`` but Pass 3a couldn't find a pixel ROI for it
    (OCR missed the label), returns the whole figure's image path with
    ``crop: false`` — the LLM can still display the whole figure and
    reason about which region is which panel using the caption.
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return {"error": f"no such paper_hash: {paper_hash}"}
    hash_dir = Path(p["hash_dir"])
    figs = _load_json(hash_dir / "figures.json", default={}) or {}
    fig = next(
        (f for f in figs.get("figures", []) or [] if f.get("figure_id") == figure_id),
        None,
    )
    if fig is None:
        return {"error": f"no such figure_id {figure_id!r} in paper {paper_hash}"}

    whole_image = hash_dir / "figures" / (fig.get("filename") or "")
    caption_text = fig.get("caption_text") or fig.get("caption") or ""
    description_from_caption = next(
        (p["description"] for p in fig.get("panels_from_caption") or []
         if p.get("label") == label),
        None,
    )
    roi_entry = next(
        (r for r in fig.get("rois") or []
         if r.get("label") == label and r.get("roi_px")),
        None,
    )

    if roi_entry is None:
        # No pixel ROI — return whole figure with the caption info so the
        # LLM can display + caption-filter mentally.
        return {
            "paper_hash": paper_hash,
            "figure_id": figure_id,
            "label": label,
            "crop": False,
            # Relative to the corpuscle's documents/ dir.
            "image_path": f"{paper_hash}/figures/{whole_image.name}",
            "caption_text": caption_text,
            "description_from_caption": description_from_caption,
            "reason": "no_pixel_roi — Pass 3a didn't locate this label in the image",
        }

    # Cache crops so a second retrieval is free.
    crops_dir = hash_dir / "figures" / "crops"
    crops_dir.mkdir(parents=True, exist_ok=True)
    crop_path = crops_dir / f"{figure_id}__{label}.png"
    if not crop_path.exists():
        try:
            from PIL import Image
            with Image.open(whole_image) as im:
                x0, y0, x1, y1 = [int(v) for v in roi_entry["roi_px"]]
                # Clamp the ROI to the image bounds defensively — ROI
                # computation can exceed image dims at the edges.
                x0 = max(0, min(x0, im.width - 1))
                y0 = max(0, min(y0, im.height - 1))
                x1 = max(x0 + 1, min(x1, im.width))
                y1 = max(y0 + 1, min(y1, im.height))
                im.crop((x0, y0, x1, y1)).save(str(crop_path))
        except Exception as e:
            return {"error": f"could not crop figure: {e}"}

    return {
        "paper_hash": paper_hash,
        "figure_id": figure_id,
        "label": label,
        "crop": True,
        # Relative to the corpuscle's documents/ dir.
        "image_path": f"{paper_hash}/figures/crops/{crop_path.name}",
        "roi_px": roi_entry.get("roi_px"),
        "ocr_confidence": roi_entry.get("ocr_confidence"),
        "caption_text": caption_text,
        "description_from_caption": description_from_caption,
    }


@mcp.tool()
def get_figure_image(
    paper_hash: str,
    figure_id: str,
    label: Optional[str] = None,
) -> Image:
    """Return a figure (or panel crop) as inline PNG bytes.

    Use this when you need the image content itself; ``get_figure`` and
    ``get_figure_roi_image`` return paths plus metadata. Without
    ``label`` returns the whole figure. With ``label`` returns the
    panel crop, falling back to the whole figure when no pixel ROI was
    detected.
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        raise ValueError(f"no such paper_hash: {paper_hash}")

    hash_dir = Path(p["hash_dir"])
    figs = _load_json(hash_dir / "figures.json", default={}) or {}
    fig = next(
        (f for f in figs.get("figures", []) or [] if f.get("figure_id") == figure_id),
        None,
    )
    if fig is None:
        raise ValueError(f"no such figure_id {figure_id!r} in paper {paper_hash}")

    whole_image = hash_dir / "figures" / (fig.get("filename") or "")
    if not whole_image.exists():
        raise FileNotFoundError(f"figure file missing on disk: {whole_image}")

    if label is None:
        return Image(path=str(whole_image))

    roi_entry = next(
        (r for r in fig.get("rois") or []
         if r.get("label") == label and r.get("roi_px")),
        None,
    )
    if roi_entry is None:
        # No pixel ROI for this label — fall back to the whole figure.
        return Image(path=str(whole_image))

    crops_dir = hash_dir / "figures" / "crops"
    crops_dir.mkdir(parents=True, exist_ok=True)
    crop_path = crops_dir / f"{figure_id}__{label}.png"
    if not crop_path.exists():
        from PIL import Image as PILImage
        with PILImage.open(whole_image) as im:
            x0, y0, x1, y1 = [int(v) for v in roi_entry["roi_px"]]
            # Clamp defensively — ROI computation can exceed image dims.
            x0 = max(0, min(x0, im.width - 1))
            y0 = max(0, min(y0, im.height - 1))
            x1 = max(x0 + 1, min(x1, im.width))
            y1 = max(y0 + 1, min(y1, im.height))
            im.crop((x0, y0, x1, y1)).save(str(crop_path))
    return Image(path=str(crop_path))


