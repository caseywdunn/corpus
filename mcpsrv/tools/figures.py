"""Figure-keyed MCP tools.

Surfaces: get_figures_for_taxon, get_figures_for_lexicon_term,
get_figure, list_figure_rois, get_figure_roi_image, get_figure_image.
The image-returning tools wrap PIL crops in MCPServer's ``Image``
content type.

#51 / #101 — figure licensing, keyed to the active output profile.
Each figure inherits license metadata from its parent work (looked up
at tool-call time from biblio_authority). The image-returning tools
gate on the per-call ``profile=`` (``report`` / ``manuscript`` /
``presentation``; see ``mcpsrv.profiles``): a ``strict`` profile refuses
bytes/URL when ``publishable=false``, the default permissive ``report``
allows them (in-chat display = fair use). The server ``--default-profile``
sets the fallback for calls that omit ``profile=``.

#154 — the *gate*, not the client, decides. Attribution fields
(``license`` / ``license_url`` / ``attribution``) go out in every profile
because captions need them; the clearance *determination*
(``publication_clearance`` / ``license_source``) is surfaced only under a
strict profile or on an explicit ``include_licensing=True``. Under the
permissive default the server has already authorized the figure, and
shipping a permission-shaped flag next to it caused figures to be
withheld in ordinary report use — on 86% of the served corpus that flag
meant only "we could not establish public domain", not "the
rightsholder refused".
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from mcp.server.mcpserver import Image
from pipeline.figures import EVIDENCE_FIGURE_TYPES, caption_evidence_summary

from ..app import _load_json, _need_index, _validated_limit, error, mcp
from ..profiles import get_profile, resolve_profile, unknown_profile_error


def _active_figure_profile(idx, profile: Optional[str]):
    """Resolve the active output profile for a figure call (#101):
    per-call ``profile=`` wins, else the server ``--default-profile``
    fallback, else built-in ``report``."""
    return resolve_profile(profile, getattr(idx, "default_profile", None))


def _figure_licensing_refusal(active, lic: Dict) -> Optional[str]:
    """Return a refusal reason when the active profile's figure-licensing
    policy forbids this figure, else None. Only the ``strict`` policy
    gates; ``permissive`` always allows.

    Deliberately keyed on the conservative ``publishable`` boolean rather
    than the five-state ``publication_clearance``: for *refusal* purposes,
    "could not establish public domain" and "explicitly restricted" do
    warrant the same answer. The distinction matters for what we *tell*
    the client, which is why the refusal message reports the state
    (#154 §2)."""
    if active.figure_licensing == "strict" and not lic.get("publishable"):
        state = lic.get("publication_clearance") or "no_record"
        because = {
            "restricted":
                "the recorded license forbids republication",
            "undetermined":
                "a license was recorded but not recognized (often a typo — "
                "check the .bib `license` field)",
            "no_record":
                "no license on file and the work is not old enough to be "
                "age-based public domain; this is an ABSENCE of evidence, "
                "not a refusal by the rightsholder",
        }.get(state, "clearance could not be established")
        return (
            f"figure withheld under profile {active.name!r} — "
            f"publication_clearance={state!r}: {because}. "
            f"license={lic.get('license') or 'none recorded'!r} "
            f"(source={lic.get('license_source')!r}). "
            "For in-chat display request profile='report'."
        )
    return None


# Restricted to the figure types that get returned by get_figures_for_*.
# Excludes graphical_element, plate_label, furniture, etc.
_REAL_FIGURE_TYPES = EVIDENCE_FIGURE_TYPES


def _resolve_lexicon_surfaces(idx, category: str, term: str) -> Dict:
    """Expand a lexicon query to every surface form of its concept (#143).

    Extraction is synonym-aware — each mention records the
    ``matched_text`` found and the ``canonical`` term it belongs to — but
    retrieval used to ignore that layer and substring-match the single
    string the caller passed. So querying ``wing`` missed captions saying
    ``ala``, querying ``ala`` missed ``wing``, and querying ``forewing``
    missed its sibling synonyms, even though the lexicon and the
    extraction layer both already knew they were the same concept.

    Returns ``{canonical, surfaces, resolved, queried}``:

    * ``resolved`` is False when the term isn't a known surface form in
      this category. We still search for the literal string then — a
      caller asking for something outside the lexicon should get a
      substring search rather than an empty list — but the caller can see
      that no expansion happened.
    * ``surfaces`` is lower-cased and always contains the query itself.

    Coverage note: surfaces come from forms that actually occurred in this
    corpus's scanned text, because a distilled bundle doesn't ship the
    lexicon YAML. A declared synonym nothing wrote won't be expanded.
    """
    queried = (term or "").strip().lower()
    surface_map = getattr(idx, "lexicon_surface_to_canonical", {}).get(category, {})
    canonical = surface_map.get(queried)
    if canonical is None:
        return {
            "canonical": None,
            "surfaces": [queried],
            "resolved": False,
            "queried": queried,
        }
    surfaces = set(
        getattr(idx, "lexicon_canonical_surfaces", {})
        .get(category, {})
        .get(canonical, set())
    )
    surfaces.add(queried)
    surfaces.add(canonical.lower())
    return {
        "canonical": canonical,
        "surfaces": sorted(surfaces),
        "resolved": True,
        "queried": queried,
    }


def _caption_surface_hits(caption_low: str, surfaces: List[str]) -> Dict:
    """Count caption occurrences across every surface form.

    Returns ``{occurrences, matched_surfaces}``. Occurrences sum over
    forms, so a caption using two synonyms of one concept ranks above one
    using a single form — which is the intent of ranking by occurrence.
    """
    total = 0
    matched: List[str] = []
    for s in surfaces:
        if not s:
            continue
        c = caption_low.count(s)
        if c:
            total += c
            matched.append(s)
    return {"occurrences": total, "matched_surfaces": matched}


def _caption_evidence_fields(figure: Dict) -> Dict:
    """Return the small, stable caption-binding summary for a figure.

    New builds persist these fields. The inference branch keeps older bundles
    truthful without rerunning caption association in the server: it only
    normalizes facts already present in the record (caption, source, and page).
    """
    return caption_evidence_summary(figure)


def _license_metadata_for_paper(paper_hash: str) -> Dict:
    """Look up license + clearance + attribution for a paper (#51).

    Returns ``{license, license_url, license_source, publishable,
    publication_clearance, attribution}`` — all keys present even when the
    authority DB has no entry, so callers see a consistent shape.

    **Internal view.** ``publishable`` is retained here because the
    strict-profile gate needs a single conservative boolean, but it must
    not go out on the wire unfiltered — see
    :func:`_license_fields_for_wire` and #154 §1/§2.
    """
    from bib.authority import CLEARANCE_NO_RECORD, clearance_state

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
            "publication_clearance": CLEARANCE_NO_RECORD,
            "attribution": None,
        }
    return {
        "license": work.get("license"),
        "license_url": work.get("license_url"),
        "license_source": work.get("license_source") or "unknown",
        "publishable": bool(work.get("publishable")),
        "publication_clearance": clearance_state(work),
        "attribution": _attribution_string(work),
    }


def _license_fields_for_wire(
    lic: Dict, active, include_licensing: bool = False,
) -> Dict:
    """The license fields that belong in a tool response (#154 §1).

    ``license`` / ``license_url`` / ``attribution`` go out in **every**
    profile — a caption needs them regardless of what the figure is for.

    The *clearance determination* does not. Under the permissive ``report``
    profile the server has already decided to serve the figure, so shipping
    ``publishable: false`` alongside it just invites the client to overrule
    a decision that was never theirs. That is not hypothetical: the default
    profile is ``report``, the server refuses nothing, and figures were
    still being withheld because a model reasonably reads itself as
    "publication-bound" and reads ``publishable`` as "may I show this".
    On 86% of the served corpus that flag was ``0`` meaning merely "we
    could not establish public domain".

    So the determination is included only when the active profile is
    actually gating on it, or when the caller explicitly asks via
    ``include_licensing=True``. When included it is reported as
    ``publication_clearance`` — a five-state description of what we know
    (see :func:`bib.authority.clearance_state`) — rather than the
    permission-shaped boolean.
    """
    out = {
        "license": lic.get("license"),
        "license_url": lic.get("license_url"),
        "attribution": lic.get("attribution"),
    }
    if include_licensing or active.figure_licensing == "strict":
        out["publication_clearance"] = lic.get("publication_clearance")
        out["license_source"] = lic.get("license_source")
    return out


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
    full_caption: bool = False,
    caption_only: bool = False,
) -> List[Dict]:
    """Figures from papers that mention the taxon, ranked by caption
    relevance.

    A figure whose caption names the taxon directly scores higher than
    a figure from a paper that merely mentions it elsewhere. **This means
    the list also includes figures whose caption does *not* name the
    taxon** — returned from any paper that mentions the taxon anywhere,
    with ``caption_has_taxon: false`` and a low ``score`` (the caption
    match contributes 100 to the score; a bare mention contributes only
    the mention count). For a precise "figures of this taxon" answer,
    filter on ``caption_has_taxon`` (or ``score``), or pass
    ``caption_only=True`` to return only caption-matched figures.

    By default only returns items classified as ``figure`` or ``plate``
    (skipping journal furniture, subpanels of already-returned figures,
    and unclassifiable graphical elements). Pass ``include_all=True`` to
    see every extracted item including the review bucket.

    ``caption_text`` is a preview (first ~200 chars) by default (#85);
    pass ``full_caption=True`` for the verbatim caption. Caption ownership is
    qualified by ``caption_status``, ``caption_confidence``,
    ``caption_page_distance``, and ``caption_kind``; do not treat an
    ``uncertain`` association as ordinary bound evidence.
    """
    try:
        n = _validated_limit(limit)
    except ValueError as e:
        return [error(str(e), "invalid_argument")]
    idx = _need_index()
    if idx.taxonomy_db is None:
        return [error("no taxonomy snapshot configured", "not_configured")]
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
            cap_full = f.get("caption_text") or f.get("caption") or ""
            caption = cap_full.lower()
            caption_hit = accepted_name_low in caption or (
                matched_name_low and matched_name_low in caption
            )
            if caption_only and not caption_hit:
                continue
            rows.append({
                "paper_hash": h,
                "paper_title": p.get("title"),
                "figure_id": f.get("figure_id"),
                "figure_type": ftype,
                "page": f.get("page"),
                "caption_text": cap_full if full_caption else _caption_preview(cap_full),
                "figure_number": f.get("figure_number"),
                # Relative to the corpuscle's documents/ dir.
                # Call get_figure_image to fetch bytes.
                "image_path": f"{h}/figures/{f.get('filename') or ''}",
                "caption_has_taxon": caption_hit,
                **_caption_evidence_fields(f),
                "score": (100 if caption_hit else 0) + idx.taxon_mention_counts.get(aid, {}).get(h, 0),
            })
    rows.sort(key=lambda r: -r["score"])
    return rows[:n]


@mcp.tool()
def get_figures_for_lexicon_term(
    category: str,
    term: str,
    paper_hash: Optional[str] = None,
    limit: int = 50,
    include_all: bool = False,
    full_caption: bool = False,
) -> Union[List[Dict], Dict]:
    """Figures whose captions mention a lexicon term, ranked by
    caption occurrence count.

    ``category`` is a top-level key in the corpus's lexicon
    (``anatomy``, ``biogeography``, …). The search is scoped to papers
    that registered at least one term in that category at extraction
    time — passing ``category="anatomy"`` will not surface figures
    whose only lexicon hit is a biogeography term. Unknown categories
    return ``{"error": "unknown_category", "available": [...]}``,
    matching ``get_figure_dossier_for_term``.

    ``term`` is **resolved through the lexicon** and matched against every
    surface form of its concept (#143), so a query for ``wing`` also finds
    captions that only ever say ``ala`` or ``forewing``, and vice versa.
    Pass the canonical name or any synonym — they are equivalent. Each row
    reports the ``matched_surfaces`` that actually hit, and the response is
    a list of rows plus a trailing note only when the term could *not* be
    resolved (then it degrades to a literal substring search, and
    ``resolved: false`` says so).

    Returns real figures + plates by default; ``include_all=True`` includes
    the review bucket.

    ``caption_text`` is a preview (first ~200 chars) by default (#85);
    pass ``full_caption=True`` for the verbatim caption. Caption ownership is
    qualified by ``caption_status``, ``caption_confidence``,
    ``caption_page_distance``, and ``caption_kind``.
    """
    try:
        n = _validated_limit(limit)
    except ValueError as e:
        return [error(str(e), "invalid_argument")]
    idx = _need_index()
    term_low = (term or "").strip().lower()
    if not term_low:
        return []
    available = sorted(idx.lexicon_to_papers.keys())
    if category not in available:
        return error("unknown lexicon category", "invalid_argument",
                     queried_category=category, available=available)
    # Restrict the candidate paper set to those tagged with at least
    # one term in this category at extraction time. Without this filter
    # the caption search runs across every paper and yields cross-
    # category false positives (e.g. an anatomy term that also occurs
    # in a biogeography-only paper).
    in_category: set = set()
    for paper_hashes in idx.lexicon_to_papers.get(category, {}).values():
        in_category.update(paper_hashes)
    target_hashes = (
        [paper_hash] if paper_hash else list(idx.papers.keys())
    )
    target_hashes = [h for h in target_hashes if h in in_category]
    # #143 — expand the query to every surface form of its concept.
    resolution = _resolve_lexicon_surfaces(idx, category, term_low)
    surfaces = resolution["surfaces"]
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
            hits = _caption_surface_hits(caption.lower(), surfaces)
            if hits["occurrences"] == 0:
                continue
            rows.append({
                "paper_hash": h,
                "paper_title": p.get("title"),
                "figure_id": f.get("figure_id"),
                "figure_type": ftype,
                "page": f.get("page"),
                "figure_number": f.get("figure_number"),
                "caption_text": caption if full_caption else _caption_preview(caption),
                # Relative to the corpuscle's documents/ dir.
                # Call get_figure_image to fetch bytes.
                "image_path": f"{h}/figures/{f.get('filename') or ''}",
                "match_count": hits["occurrences"],
                # Which synonym(s) actually hit — so a caller can tell a
                # canonical match from a synonym match (#143).
                "matched_surfaces": hits["matched_surfaces"],
                "canonical": resolution["canonical"],
                **_caption_evidence_fields(f),
            })
    rows.sort(key=lambda r: -r["match_count"])
    out = rows[:n]
    if not resolution["resolved"]:
        # Degraded to a literal substring search — say so rather than
        # letting an unrecognized term look like a synonym-aware result.
        out.append({
            "note": (
                f"{term!r} is not a known surface form in lexicon category "
                f"{category!r}; searched for the literal string only. Pass a "
                "canonical term or a synonym that occurs in the corpus for "
                "synonym-aware matching."
            ),
            "resolved": False,
        })
    return out



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
        **_caption_evidence_fields(figure_record),
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
    page, figure_number, caption_preview, caption_status,
    caption_confidence, caption_page_distance, caption_kind, image_path,
    caption_has_taxon, linked_chunks: [{chunk_id, section_class,
    headings}], rois?: {panel_count_from_caption, panel_labels,
    n_rois_with_pixel_bbox}}]}``.
    """
    idx = _need_index()
    if idx.taxonomy_db is None:
        return error("no taxonomy snapshot configured", "not_configured")
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

    ``term`` is resolved through the lexicon and matched against every
    surface form of its concept (#143) — a query for ``wing`` finds
    captions that only say ``ala``, and vice versa. Returns the same shape
    as ``get_figure_dossier_for_taxon`` minus the taxon block, plus
    ``caption_match_count`` per figure (summed across surface forms) and
    ``matched_surfaces`` naming the forms that hit. The response reports
    the resolved ``canonical`` and a ``resolved`` flag; when the term isn't
    a known surface form this degrades to a literal substring search with
    ``resolved: false``.
    """
    idx = _need_index()
    available = sorted(idx.lexicon_to_papers.keys())
    if category not in available:
        return error("unknown lexicon category", "invalid_argument",
                     queried_category=category, available=available)
    term_low = (term or "").strip().lower()
    if not term_low:
        return error("term must be non-empty", "invalid_argument")
    # #143 — same synonym expansion as get_figures_for_lexicon_term.
    resolution = _resolve_lexicon_surfaces(idx, category, term_low)
    surfaces = resolution["surfaces"]

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
            hits = _caption_surface_hits(caption.lower(), surfaces)
            occ = hits["occurrences"]
            if occ == 0:
                continue
            had_a_figure = True
            entry = _figure_dossier_entry(
                idx, h, f, chunks_by_id,
                include_rois=include_rois,
                max_linked_chunks=max_linked_chunks,
                extra_fields={
                    "caption_match_count": occ,
                    "matched_surfaces": hits["matched_surfaces"],
                },
            )
            scored.append((occ, entry))
        if had_a_figure:
            n_papers_with_figures += 1

    scored.sort(key=lambda pair: -pair[0])
    figures_out = [entry for _, entry in scored[: max_figures]]
    return {
        "category": category,
        "term": term,
        # #143 — what the term resolved to, and every form searched for.
        "canonical": resolution["canonical"],
        "resolved": resolution["resolved"],
        "surfaces_searched": surfaces,
        "n_papers_with_figures": n_papers_with_figures,
        "n_figures": len(figures_out),
        "figures": figures_out,
    }


@mcp.tool()
def get_figure(
    paper_hash: str,
    figure_id: str,
    profile: Optional[str] = None,
    include_licensing: bool = False,
) -> Dict:
    """One figure's full record: caption, page, bbox, image path,
    cross-references, plus the attribution fields (#51) inherited from the
    parent work.

    Caption ownership is explicit: ``caption_status`` is ``bound``,
    ``uncertain``, or ``unbound`` and is accompanied by confidence, kind, and
    page distance. New build artifacts also carry the full candidate trail.

    ``license`` / ``license_url`` / ``attribution`` are always present —
    you need them to caption the figure in any context.

    The publication-clearance *determination* is included only under a
    strict ``profile`` (``manuscript`` / ``presentation``), or when you
    pass ``include_licensing=True`` (#154). Under the default permissive
    ``report`` profile the server has already decided to serve the figure,
    so a clearance flag alongside it is just an invitation to withhold
    something that was cleared — display what you are given. When present,
    ``publication_clearance`` is one of ``public_domain`` /
    ``licensed_open`` / ``restricted`` / ``undetermined`` / ``no_record``;
    only ``restricted`` is positive evidence that republication was
    refused.
    """
    if profile is not None and get_profile(profile) is None:
        return unknown_profile_error(profile)
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return error(f"no such paper_hash: {paper_hash}", "not_found")
    active = _active_figure_profile(idx, profile)
    figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
    for f in figs.get("figures", []) or []:
        if f.get("figure_id") == figure_id:
            lic = _license_metadata_for_paper(paper_hash)
            return {
                **f,
                **_caption_evidence_fields(f),
                "paper_hash": paper_hash,
                "paper_title": p.get("title"),
                # Relative to the corpuscle's documents/ dir.
                "image_path": f"{paper_hash}/figures/{f.get('filename') or ''}",
                # #51 / #154 — attribution always; the clearance
                # determination only where it is actually being enforced.
                **_license_fields_for_wire(lic, active, include_licensing),
            }
    return error(f"no such figure_id {figure_id!r} in paper {paper_hash}", "not_found")



@mcp.tool()
def list_figure_rois(paper_hash: str, figure_id: str) -> List[Dict]:
    """Return the per-panel / per-subfigure ROIs annotated on a figure.

    ROIs are populated by Pass 2.5 (caption-derived ``panels_from_caption``)
    and Pass 3a (OCR-derived ``rois`` with pixel bboxes). The caption-
    derived list gives labels + descriptions even when Pass 3a wasn't
    run or OCR didn't find the panel; ``rois`` gives pixel coordinates
    when available for :func:`get_figure_roi_image`. The result also carries
    the caption status/confidence/kind/page-distance summary.
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return [error(f"no such paper_hash: {paper_hash}", "not_found")]
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
                **_caption_evidence_fields(f),
            }]
    return [error(f"no such figure_id {figure_id!r} in paper {paper_hash}", "not_found")]


@mcp.tool()
def get_figure_roi_image(
    paper_hash: str,
    figure_id: str,
    label: str,
    profile: Optional[str] = None,
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

    Caption status/confidence/kind/page distance qualify that description;
    callers should preserve the uncertainty when the binding is not strong.

    Honors the figure-licensing gate for the active ``profile``, like
    ``get_figure_image`` and ``get_figure_url`` (#154 §3). It previously
    did not, which made it a way around them: a client refused under
    ``manuscript`` could obtain the same pixels here — including, via the
    no-pixel-ROI fallback below, the whole uncropped figure.
    """
    if profile is not None and get_profile(profile) is None:
        return unknown_profile_error(profile)
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        return error(f"no such paper_hash: {paper_hash}", "not_found")
    # Gate before touching disk, matching the other enforcement points —
    # and before the `roi_entry is None` fallback, which returns the whole
    # figure and was the widest part of the bypass.
    active = _active_figure_profile(idx, profile)
    lic = _license_metadata_for_paper(paper_hash)
    refusal = _figure_licensing_refusal(active, lic)
    if refusal:
        return error(refusal, "forbidden")
    hash_dir = Path(p["hash_dir"])
    figs = _load_json(hash_dir / "figures.json", default={}) or {}
    fig = next(
        (f for f in figs.get("figures", []) or [] if f.get("figure_id") == figure_id),
        None,
    )
    if fig is None:
        return error(f"no such figure_id {figure_id!r} in paper {paper_hash}", "not_found")

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
            **_caption_evidence_fields(fig),
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
            return error(f"could not crop figure: {e}", "unavailable")

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
        **_caption_evidence_fields(fig),
        "description_from_caption": description_from_caption,
    }


@mcp.tool()
def get_figure_image(
    paper_hash: str,
    figure_id: str,
    label: Optional[str] = None,
    profile: Optional[str] = None,
) -> Image:
    """Return a figure (or panel crop) as inline PNG bytes.

    Use this when you need the image content itself; ``get_figure`` and
    ``get_figure_roi_image`` return paths plus metadata. Without
    ``label`` returns the whole figure. With ``label`` returns the
    panel crop, falling back to the whole figure when no pixel ROI was
    detected.

    #101: the figure-licensing gate is keyed to the active output
    ``profile`` (``report`` / ``manuscript`` / ``presentation``; pass it
    per call to reflect what you're producing). Under a ``strict``
    profile (manuscript / presentation) this refuses image bytes when
    the parent work's clearance could not be established. Under the
    default permissive ``report`` profile the image is returned (in-chat
    display is fair use).

    **Don't self-filter** (#154): pass the profile that matches what you
    are producing and let the server decide. If it returns bytes under
    ``report``, display them — a clearance determination is not included
    there precisely because it is not being enforced. ``get_figure(...,
    include_licensing=True)`` gives you the determination explicitly if
    you need to reason about it. Unknown
    profile names raise.
    """
    idx = _need_index()
    p = idx.papers.get(paper_hash)
    if not p:
        raise ValueError(f"no such paper_hash: {paper_hash}")

    if profile is not None and get_profile(profile) is None:
        raise ValueError(
            f"unknown profile {profile!r}; use list_output_profiles()"
        )
    active = _active_figure_profile(idx, profile)

    # #101 — figure-licensing gate, keyed to the active profile. Refuses
    # with a structured ValueError so clients can branch on the message.
    lic = _license_metadata_for_paper(paper_hash)
    refusal = _figure_licensing_refusal(active, lic)
    if refusal:
        raise ValueError(
            f"{refusal}. The image is not returned to avoid downstream "
            f"copyright issues. Read get_figure({paper_hash!r}, "
            f"{figure_id!r}) for the raw license fields, or pass "
            f"profile='report' for in-chat display."
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
    profile: Optional[str] = None,
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

    #101: the figure-licensing gate is keyed to the active ``profile``;
    a ``strict`` profile (manuscript / presentation) refuses an
    unpublishable figure, the default permissive ``report`` allows it.
    The resolved profile is **encoded into the returned URL** so the
    subsequent HTTP fetch enforces the same policy (otherwise a strict
    client could leak via the unprofiled path).

    Returns ``{url, auth_header, mime_type, profile, publishable,
    license, license_source, attribution}`` on success, ``{error: ...}``
    on failure (including ``unknown profile``).
    """
    idx = _need_index()
    if profile is not None and get_profile(profile) is None:
        return unknown_profile_error(profile)
    active = _active_figure_profile(idx, profile)
    base = getattr(idx, "figure_url_base", None)
    if not base:
        return error(
            "figure HTTP route is not available on this server. "
            "Possible causes: figure side-car failed to bind at "
            "startup (check server logs), or the server is running "
            "with an older mcpsrv that predates #69. "
            "Fall back to get_figure_image.",
            "unavailable",
        )

    p = idx.papers.get(paper_hash)
    if not p:
        return error(f"no such paper_hash: {paper_hash}", "not_found")

    # Verify the figure record exists before handing out a URL — the
    # HTTP route would 404 anyway, but a JSON error here is cheaper
    # for the caller and surfaces the issue immediately.
    figs = _load_json(Path(p["hash_dir"]) / "figures.json", default={}) or {}
    fig = next(
        (f for f in figs.get("figures", []) or [] if f.get("figure_id") == figure_id),
        None,
    )
    if fig is None:
        return error(f"no such figure_id {figure_id!r} in paper {paper_hash}", "not_found")

    lic = _license_metadata_for_paper(paper_hash)
    refusal = _figure_licensing_refusal(active, lic)
    if refusal:
        return error(
            f"{refusal}. Server refuses to hand out a URL the operator "
            f"would only get a 403 from. Read get_figure({paper_hash!r}, "
            f"{figure_id!r}) for the raw license fields, or pass "
            f"profile='report' for in-chat display.",
            "forbidden", profile=active.name, **lic,
        )

    # Encode the resolved profile into the URL so the HTTP route enforces
    # the same policy (the route defaults to the server fallback when no
    # profile is present — a strict client must not leak via that path).
    query = [f"profile={active.name}"]
    if label:
        # figure_id / label are constrained to [A-Za-z0-9._-] server-side,
        # so no urlencoding gymnastics are needed.
        query.append(f"label={label}")
    url = f"{base}/figures/{paper_hash}/{figure_id}?{'&'.join(query)}"
    token = getattr(idx, "figure_auth_token", None)
    auth_header = f"Authorization: Bearer {token}" if token else None

    return {
        "url": url,
        "auth_header": auth_header,
        "mime_type": "image/png",
        "profile": active.name,
        # #154 §1 — the server just authorized this URL under `active`, so
        # don't ship a clearance flag inviting the client to withhold it.
        # Attribution and license go out either way (a caption needs them);
        # the determination only under a strict profile.
        **_license_fields_for_wire(lic, active),
        "fetch_hint": (
            "curl -fsSL -H \"$auth_header\" -o /tmp/fig.png \"$url\""
            if auth_header
            else "curl -fsSL -o /tmp/fig.png \"$url\"   # no auth configured"
        ),
    }
