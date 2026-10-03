"""Materialize explicit figure-license exclusions during the build (#302).

The server consumes these facts; it never interprets caption prose. An
exclusion from the publication's license does not establish that every reuse
is forbidden. It means that the publication's license cannot clear the image.
See dev_docs/OVERVIEW.md, "Figure-specific rights evidence".
"""
from __future__ import annotations

from copy import deepcopy
import re

FIGURE_RIGHTS_VERSION = "figure-rights-v1"

# Deliberately require an explicit negation and license reference. A credit,
# "reproduced with permission", or a normal CC statement is not an exclusion.
_EXCLUSION = re.compile(
    r"\b(?:not\s+(?:covered|included)|excluded)\s+"
    r"(?:by|under|in|from)\s+(?:the\s+)?(?:terms\s+of\s+(?:the\s+)?)?"
    r"(?:creative\s+commons|cc[-\s]?by|(?:article|publication)(?:'s|’s)?)"
    r"[^.!?\n]{0,100}?\blicen[cs]e\b[^.!?\n]*",
    re.IGNORECASE,
)


def materialize_figure_rights(figures: list[dict], *, preserve_existing: bool = False) -> None:
    """Refresh caption evidence and propagate it across shared image records.

    All crops conservatively inherit the image's exclusions. Splitting a mixed
    image cannot grant clearance to a panel or whole-image fallback. A future
    reviewed panel-specific grant would need its own input/provenance and gates.
    Fresh extraction replaces obsolete derived evidence. Later figure passes
    preserve exclusions already observed before splitting/rebinding captions.
    """
    parents = list(range(len(figures)))

    def root(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i

    def join(a, b):
        parents[root(a)] = root(b)

    ids = {fig.get("figure_id"): i for i, fig in enumerate(figures)
           if fig.get("figure_id")}
    filenames = {}
    for i, fig in enumerate(figures):
        name = fig.get("filename")
        if name:
            if name in filenames:
                join(i, filenames[name])
            filenames[name] = i
        for key in ("image_shared_with", "plate_source_figure_id",
                    "plate_roi_source_figure_id"):
            parent = ids.get(fig.get(key))
            if parent is not None:
                join(i, parent)

    evidence = {}
    for i, fig in enumerate(figures):
        previous = fig.get("figure_rights") or {}
        if preserve_existing and previous.get("status") == "excluded_from_publication_license":
            evidence.setdefault(root(i), []).extend(deepcopy(previous.get("evidence") or []))
        caption = fig.get("caption_text") or fig.get("caption") or ""
        # Whitespace normalization also handles a line-wrapped license notice.
        normalized = " ".join(caption.split())
        for match in _EXCLUSION.finditer(normalized):
            start = max(normalized.rfind(".", 0, match.start()),
                        normalized.rfind("!", 0, match.start()),
                        normalized.rfind("?", 0, match.start())) + 1
            end = match.end()
            if end < len(normalized) and normalized[end] in ".!?":
                end += 1
            evidence.setdefault(root(i), []).append({
                "source": "caption_explicit_exclusion",
                "figure_id": fig.get("figure_id"),
                "page": fig.get("caption_page") or fig.get("page"),
                "caption_source": fig.get("caption_source"),
                "caption_status": fig.get("caption_status"),
                "text": normalized[start:end].strip(),
            })

    for i, fig in enumerate(figures):
        distinct = {tuple((k, str(v)) for k, v in sorted(item.items())): item
                    for item in evidence.get(root(i), [])}
        found = sorted(distinct.values(),
                       key=lambda e: (str(e["figure_id"]), e["text"]))
        fig["figure_rights"] = {
            "producer_version": FIGURE_RIGHTS_VERSION,
            "status": "excluded_from_publication_license" if found else "inherit",
            "scope": "whole_image_and_all_crops",
            "evidence": deepcopy(found),
        }
