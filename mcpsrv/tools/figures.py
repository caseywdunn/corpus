"""Figure-keyed MCP tools.

Surfaces: get_figures_for_taxon, get_figures_for_lexicon_term,
get_figure, list_figure_rois, get_figure_roi_image, get_figure_image.
The image-returning tools wrap PIL crops in FastMCP's ``Image``
content type.

#51 — figure licensing. Each figure inherits license metadata from
its parent work (looked up at tool-call time from biblio_authority).
``get_figure_image`` refuses image bytes when ``publishable=false``,
returning a ValueError with the reason; clients with non-default
jurisdictions can read the raw license fields and re-derive.
``mcpsrv.main --allow-unpublishable`` bypasses for local rights-holder
cases.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from mcp.server.fastmcp import Image

from ..app import _load_json, _need_index, mcp


# Restricted to the figure types that get returned by get_figures_for_*.
# Excludes graphical_element, plate_label, furniture, etc.
_REAL_FIGURE_TYPES = {"figure", "plate", "subpanel"}


def _license_metadata_for_paper(paper_hash: str) -> Dict:
    """Look up license + publishable + attribution for a paper (#51).

    Returns ``{license, license_url, license_source, publishable,
    attribution}`` — all keys present even when the authority DB has
    no entry (everything falls back to ``unknown`` / ``False`` so
    clients see a consistent shape).
    """
    idx = _need_index()
    work = None
    if idx.biblio_db is not None:
        try:
            work = idx.biblio_db.get_work_by_corpus_hash(paper_hash)
        except Exception:
            work = None
    if not work:
        return {
            "license": None,
            "license_url": None,
            "license_source": "unknown",
            "publishable": False,
            "attribution": None,
        }
    return {
        "license": work.get("license"),
        "license_url": work.get("license_url"),
        "license_source": work.get("license_source") or "unknown",
        "publishable": bool(work.get("publishable")),
        "attribution": _attribution_string(work),
    }


def _attribution_string(work: Dict) -> Optional[str]:
    """Server-computed canonical attribution: ``Author Year. Title. doi:X``.

    Cheap to compute; clients can ignore it and re-derive from the raw
    fields if their citation style differs.
    """
    title = (work.get("title") or "").strip()
    year = work.get("year")
    doi = (work.get("doi") or "").strip()
    work_id = work.get("work_id")
    idx = _need_index()
    authors_str = ""
    if idx.biblio_db is not None and work_id:
        try:
            authors = idx.biblio_db.get_authors(work_id)
        except Exception:
            authors = []
        if authors:
            surnames = [a["surname"] for a in authors if a.get("surname")]
            if len(surnames) > 2:
                authors_str = f"{surnames[0]} et al."
            elif len(surnames) == 2:
                authors_str = f"{surnames[0]} & {surnames[1]}"
            elif len(surnames) == 1:
                authors_str = surnames[0]
    parts = []
    if authors_str:
        parts.append(authors_str)
    if year is not None:
        parts.append(f"({year})")
    if title:
        parts.append(title.rstrip(".") + ".")
    if doi:
        parts.append(f"doi:{doi}")
    return " ".join(parts) or None


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



# #76 — bounded defaults for the figure-dossier pair.
_FIGURE_DOSSIER_MAX_FIGURES_DEFAULT = 25
_FIGURE_DOSSIER_MAX_LINKED_CHUNKS_DEFAULT = 10
_FIGURE_CAPTION_PREVIEW_CHARS = 200


def _caption_preview(caption: str, n: int = _FIGURE_CAPTION_PREVIEW_CHARS) -> str:
    """Trim a caption to a preview length without breaking mid-word."""
    caption = (caption or "").strip()
    if len(caption) <= n:
        return caption
    cut = caption[: n].rsplit(" ", 1)[0]
    return cut + "…"


def _linked_chunks_for_figure(
    figure_id: str,
    chunks_by_id: Dict[str, Dict],
    max_chunks: int,
) -> List[Dict]:
    """Find chunks in this paper that reference ``figure_id`` via
    chunks.json ``figure_refs``. Returns lightweight entries (IDs +
    section + headings, no text — caller pairs with get_chunks)."""
    out: List[Dict] = []
    for cid, ch in chunks_by_id.items():
        if figure_id in (ch.get("figure_refs") or []):
            out.append({
                "chunk_id": cid,
                "section_class": ch.get("section_class"),
                "headings": ch.get("headings") or [],
            })
            if len(out) >= max_chunks:
                break
    return out


def _figure_dossier_entry(
    idx,
    paper_hash: str,
    figure_record: Dict,
    chunks_by_id: Dict[str, Dict],
    *,
    include_rois: bool,
    max_linked_chunks: int,
    extra_fields: Optional[Dict] = None,
) -> Dict:
    """Compose the per-figure dossier row. ROIs are summarised to
    counts + panel labels (caption-derived) rather than full ROI
    objects — those are still available via list_figure_rois."""
    p = idx.papers.get(paper_hash) or {}
    figure_id = figure_record.get("figure_id")
    entry: Dict[str, Any] = {
        "paper_hash": paper_hash,
        "paper_title": p.get("title"),
        "paper_year": p.get("year"),
        "figure_id": figure_id,
        "figure_type": figure_record.get("figure_type"),
        "page": figure_record.get("page"),
        "figure_number": figure_record.get("figure_number"),
        "caption_preview": _caption_preview(
            figure_record.get("caption_text") or figure_record.get("caption") or "",
        ),
        "image_path": f"{paper_hash}/figures/{figure_record.get('filename') or ''}",
        "linked_chunks": _linked_chunks_for_figure(
            figure_id, chunks_by_id, max_linked_chunks,
        ),
    }
    if include_rois:
        panels = figure_record.get("panels_from_caption") or []
        rois = figure_record.get("rois") or []
        entry["rois"] = {
            "panel_count_from_caption": figure_record.get(
                "panel_count_from_caption", 0,
            ),
            "panel_labels": [
                p.get("label") for p in panels if p.get("label")
            ],
            "n_rois_with_pixel_bbox": len(rois),
        }
    if extra_fields:
        entry.update(extra_fields)
    return entry


@mcp.tool()
def get_figure_dossier_for_taxon(
    taxon_name: str,
    max_figures: int = _FIGURE_DOSSIER_MAX_FIGURES_DEFAULT,
    max_linked_chunks: int = _FIGURE_DOSSIER_MAX_LINKED_CHUNKS_DEFAULT,
    include_rois: bool = True,
) -> Dict[str, Any]:
    """Figures linked to a taxon, each with its explanatory chunk IDs
    (#76). Supersedes ``get_figures_for_taxon`` + per-figure
    ``list_figure_rois`` + cross-ref against ``get_chunks_for_taxon``.
    Pair with ``get_chunks(paper_hash, chunk_ids=[...])`` to read the
    explanatory passages.

    Real figures + plates only (graphical_element / plate_label
    skipped). Ranked by caption-name match > mere paper-mention.

    Returns ``{taxon, n_papers_with_figures, n_figures, figures:
    [{paper_hash, paper_title, paper_year, figure_id, figure_type,
    page, figure_number, caption_preview, image_path,
    caption_has_taxon, linked_chunks: [{chunk_id, section_class,
    headings}], rois?: {panel_count_from_caption, panel_labels,
    n_rois_with_pixel_bbox}}]}``.
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return {"error": "no taxonomy snapshot configured"}
    hit = idx.taxonomy_db.lookup(taxon_name)
    if not hit:
        return {"not_found": True, "queried": taxon_name}

    aid = hit["accepted_taxon_id"]
    accepted_low = (hit.get("accepted_name") or "").lower()
    matched_low = (hit.get("matched_name") or "").lower()
    paper_hashes = list(idx.taxon_to_papers.get(aid, []))

    taxon_block: Dict[str, Any] = {
        "taxon_id": aid,
        "accepted_name": hit.get("accepted_name"),
    }
    if hit.get("rank"):
        taxon_block["rank"] = hit["rank"]

    scored: List[Dict] = []
    n_papers_with_figures = 0
    for h in paper_hashes:
        p = idx.papers.get(h)
        if not p:
            continue
        figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
        fig_list = figs.get("figures", []) or []
        chunks_data = _load_json(
            Path(p["hash_dir"]) / "chunks.json", default={},
        ) or {}
        chunks_by_id = {
            c.get("chunk_id"): c
            for c in chunks_data.get("chunks", []) or []
        }
        had_a_figure = False
        for f in fig_list:
            if f.get("figure_type") not in _REAL_FIGURE_TYPES:
                continue
            had_a_figure = True
            caption = (f.get("caption_text") or f.get("caption") or "").lower()
            caption_hit = (
                bool(accepted_low and accepted_low in caption)
                or bool(matched_low and matched_low in caption)
            )
            score = (100 if caption_hit else 0) + idx.taxon_mention_counts.get(
                aid, {},
            ).get(h, 0)
            entry = _figure_dossier_entry(
                idx, h, f, chunks_by_id,
                include_rois=include_rois,
                max_linked_chunks=max_linked_chunks,
                extra_fields={"caption_has_taxon": caption_hit},
            )
            scored.append((score, entry))
        if had_a_figure:
            n_papers_with_figures += 1

    scored.sort(key=lambda pair: -pair[0])
    figures_out = [entry for _, entry in scored[: max_figures]]
    return {
        "taxon": taxon_block,
        "n_papers_with_figures": n_papers_with_figures,
        "n_figures": len(figures_out),
        "figures": figures_out,
    }


@mcp.tool()
def get_figure_dossier_for_term(
    category: str,
    term: str,
    max_figures: int = _FIGURE_DOSSIER_MAX_FIGURES_DEFAULT,
    max_linked_chunks: int = _FIGURE_DOSSIER_MAX_LINKED_CHUNKS_DEFAULT,
    include_rois: bool = True,
) -> Dict[str, Any]:
    """Figures whose captions mention a lexicon term, each with its
    explanatory chunk IDs (#76).

    Replaces the chain ``get_figures_for_lexicon_term`` + per-figure
    ``list_figure_rois`` + cross-ref against ``get_chunks_for_topic``
    for the "show me figures depicting <term> and the passages that
    explain them" pattern.

    Category-agnostic. Match is caption-substring (case-insensitive)
    on ``term``. Returns the same shape as
    ``get_figure_dossier_for_taxon`` minus the taxon block, plus
    ``caption_match_count`` per figure (the substring hit count).
    """
    idx = _need_index()
    available = sorted(idx.lexicon_to_papers.keys())
    if category not in available:
        return {
            "error": "unknown_category",
            "queried_category": category,
            "available": available,
        }
    term_low = (term or "").strip().lower()
    if not term_low:
        return {"error": "empty_term"}

    scored: List[Dict] = []
    n_papers_with_figures = 0
    for h, p in idx.papers.items():
        figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
        fig_list = figs.get("figures", []) or []
        chunks_data = _load_json(
            Path(p["hash_dir"]) / "chunks.json", default={},
        ) or {}
        chunks_by_id = {
            c.get("chunk_id"): c
            for c in chunks_data.get("chunks", []) or []
        }
        had_a_figure = False
        for f in fig_list:
            if f.get("figure_type") not in _REAL_FIGURE_TYPES:
                continue
            caption = (f.get("caption_text") or f.get("caption") or "")
            occ = caption.lower().count(term_low)
            if occ == 0:
                continue
            had_a_figure = True
            entry = _figure_dossier_entry(
                idx, h, f, chunks_by_id,
                include_rois=include_rois,
                max_linked_chunks=max_linked_chunks,
                extra_fields={"caption_match_count": occ},
            )
            scored.append((occ, entry))
        if had_a_figure:
            n_papers_with_figures += 1

    scored.sort(key=lambda pair: -pair[0])
    figures_out = [entry for _, entry in scored[: max_figures]]
    return {
        "category": category,
        "term": term,
        "n_papers_with_figures": n_papers_with_figures,
        "n_figures": len(figures_out),
        "figures": figures_out,
    }


@mcp.tool()
def get_figure(paper_hash: str, figure_id: str) -> Dict:
    """One figure's full record: caption, page, bbox, image path,
    cross-references, plus license + publishable + attribution (#51)
    inherited from the parent work."""
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
                # #51 — license metadata inherited from the parent work.
                **_license_metadata_for_paper(paper_hash),
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

    #51: refuses to return image bytes when the parent work's
    ``publishable`` flag is false (license forbids reuse, or unknown
    license, or work is too recent for the configured PD cutoff).
    Override at server start with ``mcpsrv.main --allow-unpublishable``
    for local-only rights-holder cases. The raw license fields are
    always exposed via ``get_figure``, so clients with non-default
    jurisdictions can re-derive.
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        raise ValueError(f"no such paper_hash: {paper_hash}")

    # #51 — publishable gate. Server-level allow override carried on the
    # CorpusIndex. Refuses with a structured ValueError so clients
    # can branch on the message.
    if not getattr(idx, "allow_unpublishable", False):
        lic = _license_metadata_for_paper(paper_hash)
        if not lic["publishable"]:
            license_v = lic.get("license") or "unknown"
            raise ValueError(
                f"figure not publishable: license={license_v!r} "
                f"(source={lic['license_source']!r}). The image is not "
                f"returned to avoid downstream copyright issues. Read "
                f"get_figure({paper_hash!r}, {figure_id!r}) for the raw "
                f"license fields if your jurisdiction / use case differs, "
                f"or restart the server with --allow-unpublishable for "
                f"local rights-holder cases."
            )

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


@mcp.tool()
def get_figure_url(
    paper_hash: str,
    figure_id: str,
    label: Optional[str] = None,
) -> Dict:
    """Return a bearer-gated HTTP URL the caller can ``curl -o`` to
    land the figure PNG on disk *without* loading its bytes into the
    model's context window. Use instead of ``get_figure_image`` when
    file bytes must reach the filesystem (pandoc / LaTeX / PDF
    assembly) — the byte flow stays off the MCP JSON-RPC channel
    regardless of figure size.

    Fetch via ``curl -H "$auth_header" -o <path> "$url"``. Without
    ``label`` returns the whole figure; with ``label`` returns the
    panel crop if one exists (else falls back to the whole figure).
    Refuses unpublishable figures (#51) unless server started with
    ``--allow-unpublishable``.

    Returns ``{url, auth_header, mime_type, publishable, license,
    license_source}`` on success, ``{error: ...}`` on failure.
    """
    idx = _need_index()
    base = getattr(idx, "figure_url_base", None)
    if not base:
        return {
            "error": (
                "figure HTTP route is not available on this server. "
                "Possible causes: figure side-car failed to bind at "
                "startup (check server logs), or the server is running "
                "with an older mcpsrv that predates #69. "
                "Fall back to get_figure_image."
            ),
        }

    p = idx.papers.get(paper_hash)
    if not p:
        return {"error": f"no such paper_hash: {paper_hash}"}

    # Verify the figure record exists before handing out a URL — the
    # HTTP route would 404 anyway, but a JSON error here is cheaper
    # for the caller and surfaces the issue immediately.
    figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
    fig = next(
        (f for f in figs.get("figures", []) or [] if f.get("figure_id") == figure_id),
        None,
    )
    if fig is None:
        return {"error": f"no such figure_id {figure_id!r} in paper {paper_hash}"}

    lic = _license_metadata_for_paper(paper_hash)
    if not getattr(idx, "allow_unpublishable", False) and not lic.get("publishable"):
        return {
            "error": (
                f"figure not publishable: license={lic.get('license') or 'unknown'!r} "
                f"(source={lic.get('license_source')!r}). Server refuses to hand out "
                f"a URL the operator would only get a 403 from. Read "
                f"get_figure({paper_hash!r}, {figure_id!r}) for the raw license "
                f"fields, or restart the server with --allow-unpublishable for "
                f"local rights-holder cases."
            ),
            **lic,
        }

    url = f"{base}/figures/{paper_hash}/{figure_id}"
    if label:
        # Single key=value, no urlencoding gymnastics needed — figure_id
        # / label are constrained to [A-Za-z0-9._-] on the server side.
        url = f"{url}?label={label}"
    token = getattr(idx, "figure_auth_token", None)
    auth_header = f"Authorization: Bearer {token}" if token else None

    return {
        "url": url,
        "auth_header": auth_header,
        "mime_type": "image/png",
        "publishable": lic.get("publishable"),
        "license": lic.get("license"),
        "license_source": lic.get("license_source"),
        "fetch_hint": (
            "curl -fsSL -H \"$auth_header\" -o /tmp/fig.png \"$url\""
            if auth_header
            else "curl -fsSL -o /tmp/fig.png \"$url\"   # no auth configured"
        ),
    }


