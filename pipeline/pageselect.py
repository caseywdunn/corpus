"""Per-document page selection: `keeppages` (#188).

A PDF in a scanned library is frequently not just the paper. Library cover
sheets, request slips and colour calibration targets get prepended by the
digitiser; a Russian original is bound in front of its English translation;
runs of blank versos pad the plates. Every one of those costs something, and
the costs compound: `detect_scan_type` samples pages to pick the OCR mode and
the language pack, so forty pages of scanner filler ahead of the body can
decide how the body is read. Then OCR pays full price for the filler, and a
calibration target becomes a figure.

`keeppages` is a physical page selection written on the bib entry:

    keeppages = {3--20}
    keeppages = {2,4,8--20,22--40,55}
    keeppages = {40--}

**Physical, 1-based positions in the file, never printed page numbers.** The
printed number is often absent, wrong, or restarts mid-volume — and on the
documents this targets it is precisely the front matter that has none. An
entry will routinely carry both `pages = {699--714}` (where the article sits
in its journal) and a `keeppages` that disagrees completely. That is correct,
not a data error.

The selection is applied *before* scan detection, so every later stage sees
only the selected pages, and page numbers in `figures.json` and `text.json`
are positions in the subset. :func:`selection_for` returns the resolved list,
which is itself the complete subset→source map — subset page *i* is
``selected[i - 1]`` — so `source_page` can be carried alongside `page`
without building a second structure to keep in sync.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import List, Optional, Tuple

logger = logging.getLogger(__name__)

# BibTeX page conventions: `--` for a range, `,` between terms, singles and
# ranges mixed freely. A single `-` is accepted too, because it is what an
# operator types half the time and there is no other meaning available here.
_TERM = re.compile(r"^\s*(\d+)\s*(?:-{1,3}\s*(\d*)\s*)?$")


class PageSelectionError(ValueError):
    """The `keeppages` value could not be parsed at all."""


def parse_keeppages(raw: Optional[str]) -> List[Tuple[int, Optional[int]]]:
    """Parse a `keeppages` value into ``[(start, end_or_None), ...]``.

    ``end`` is None for an open-ended range (`40--`, "page 40 to the end"),
    which is the common shape for trimming front matter off a document whose
    length the operator has not counted. A bare number is a one-page range.

    Raises :class:`PageSelectionError` on anything unparseable, because
    unlike `ocrlang` — where an unknown pack is dropped and detection carries
    on — there is no safe reading of a malformed page range. Silently keeping
    every page would look exactly like success.
    """
    if raw is None:
        return []
    text = raw.strip()
    if not text:
        return []
    terms: List[Tuple[int, Optional[int]]] = []
    for chunk in text.split(","):
        if not chunk.strip():
            continue
        m = _TERM.match(chunk)
        if not m:
            raise PageSelectionError(
                f"cannot parse page range {chunk.strip()!r} in keeppages={raw!r}; "
                "expected forms are 5, 3--20, or 40-- (1-based physical pages)"
            )
        start = int(m.group(1))
        if start < 1:
            raise PageSelectionError(
                f"page numbers are 1-based; got {start} in keeppages={raw!r}"
            )
        end_raw = m.group(2)
        if end_raw is None:                      # a bare number: one page
            terms.append((start, start))
        elif end_raw == "":                      # `40--`: open-ended
            terms.append((start, None))
        else:
            end = int(end_raw)
            if end < start:
                raise PageSelectionError(
                    f"range {start}--{end} runs backwards in keeppages={raw!r}; "
                    "a selection cannot reorder a document"
                )
            terms.append((start, end))
    return terms


def resolve_selection(terms, n_pages: int, name: str = "") -> Tuple[List[int], List[str]]:
    """Resolve parsed terms against a real page count.

    Returns ``(selected, warnings)`` where ``selected`` is a sorted,
    deduplicated list of 1-based page numbers.

    Sorted and deduplicated rather than taken in the order written: a
    selection says which pages *are* the document, and `40--50,10--20` must
    not become a way to reorder it. Overlapping terms collapse.

    A range running past the end is clamped with a warning rather than
    failing, matching how `ocrlang` drops packs that aren't installed — it is
    almost always a typo about where a document ends, and the pages that do
    exist are still the right ones. The warning is recorded in
    `scan_detection.json` so a clamp is visible rather than inferred.
    """
    warnings: List[str] = []
    if n_pages < 1:
        return [], warnings
    keep = set()
    for start, end in terms:
        stop = n_pages if end is None else end
        if start > n_pages:
            warnings.append(
                f"keeppages range {start}--{'' if end is None else end} starts "
                f"past the last page ({n_pages}); ignored"
            )
            continue
        if stop > n_pages:
            warnings.append(
                f"keeppages range {start}--{end} extends past the last page "
                f"({n_pages}); clamped"
            )
            stop = n_pages
        keep.update(range(start, stop + 1))
    selected = sorted(keep)
    for w in warnings:
        logger.warning("%s%s", f"{name}: " if name else "", w)
    return selected, warnings


def write_subset(src: Path, dst: Path, selected: List[int]) -> bool:
    """Write the selected pages of ``src`` to ``dst``. True if written.

    Returns False — leaving ``dst`` alone — when the selection is empty or
    covers the whole document, so an annotated document that keeps every page
    is byte-identical to an unannotated one and costs nothing.
    """
    if not selected:
        return False
    try:
        import fitz
    except ImportError:                          # pragma: no cover
        logger.warning("PyMuPDF unavailable; keeppages selection not applied")
        return False
    doc = fitz.open(src)
    try:
        if len(selected) == doc.page_count and selected == list(range(1, doc.page_count + 1)):
            return False
        out = fitz.open()
        try:
            for page_no in selected:
                out.insert_pdf(doc, from_page=page_no - 1, to_page=page_no - 1)
            out.save(dst)
        finally:
            out.close()
    finally:
        doc.close()
    return True


def selection_for(raw: Optional[str], pdf_path: Path,
                  name: str = "") -> Tuple[List[int], List[str]]:
    """Parse and resolve in one step against the real page count of ``pdf_path``.

    ``([], [])`` when there is no directive, which is the overwhelming case.
    """
    terms = parse_keeppages(raw)
    if not terms:
        return [], []
    try:
        import fitz
        with fitz.open(pdf_path) as doc:
            n_pages = doc.page_count
    except Exception as e:                       # pragma: no cover
        logger.warning("could not read page count of %s: %s", pdf_path, e)
        return [], []
    return resolve_selection(terms, n_pages, name=name)


# Keys in the persisted artifacts that hold a page number in the *subset*.
# `page` is on every figure and on chunk provenance; `caption_page` is
# separate because a caption can sit on a different leaf from its figure.
_PAGE_KEYS = ("page", "caption_page")


def _annotate(obj, selected: List[int]) -> int:
    """Recursively add ``source_page`` beside every subset page number."""
    n = 0
    if isinstance(obj, dict):
        for key in _PAGE_KEYS:
            value = obj.get(key)
            src_key = "source_page" if key == "page" else "caption_source_page"
            if isinstance(value, int) and 1 <= value <= len(selected):
                obj[src_key] = selected[value - 1]
                n += 1
        for v in obj.values():
            n += _annotate(v, selected)
    elif isinstance(obj, list):
        for v in obj:
            n += _annotate(v, selected)
    return n


def annotate_source_pages(paths, selected: List[int]) -> int:
    """Add ``source_page`` to every page number in ``paths``. Returns the count.

    Once pages are dropped, `page` in `figures.json` and `text.json` is a
    position in the *subset*, and that is the number a client sees when a
    figure is served. A figure reported as page 3 that is page 44 of the scan
    is a citation error, and a silent one — nothing downstream can detect it.

    So both numbers are carried: `page` stays subset-relative, because that is
    what indexes the artifacts actually on disk, and `source_page` says where
    it came from in the file the operator holds. Written only when a selection
    is active, so an unannotated document's artifacts are unchanged.

    Applied as a post-pass over the finished artifacts rather than threaded
    through each construction site: the map is a list lookup, the passes that
    rewrite `figures.json` would each need the same change, and one place to
    get wrong beats five.
    """
    import json

    total = 0
    for path in paths:
        path = Path(path)
        if not path.exists():
            continue
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError) as e:
            logger.warning("could not annotate source pages in %s: %s", path, e)
            continue
        n = _annotate(data, selected)
        if n:
            path.write_text(json.dumps(data, indent=2))
            total += n
    return total
