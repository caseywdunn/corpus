"""Citation string formatters for the format_citation MCP tool (#79).

Hand-rolled, no new deps. One function per style. Default style is
``author-year`` — an APA-like reference-list + parenthetical-inline
shape, suitable for scientific writing and easy for downstream
LLM-mediated conversion to other styles.

Each formatter takes a ``fields`` dict with the shape:

    {
        "authors": [{"surname": "Totton", "forename": "A. K."}, ...],
        "year":    1965,
        "title":   "A synopsis of the Siphonophora",
        "journal": "Bulletin of the British Museum (Natural History)",
        "volume":  "7",   # optional
        "pages":   "1-200",  # optional
        "doi":     "10.5962/bhl.part.20269",  # optional
    }

and returns ``{"formatted": <reference-list entry>, "inline":
<parenthetical citation>}``. Fields beyond authors / year are best-
effort: any missing one is omitted from the output rather than
substituted with a placeholder, so partial records degrade
gracefully.
"""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

# Public set of supported styles. The MCP tool's ``style`` argument
# is validated against this; future styles (e.g. ``vancouver``,
# ``bibtex``) are TODOs that will extend this set + add a formatter.
SUPPORTED_STYLES = frozenset({"author-year"})


def format_citation(fields: Dict[str, Any], style: str = "author-year") -> Dict[str, str]:
    """Dispatch to the named style formatter. Raises ValueError on
    an unknown style — caller is expected to validate against
    :data:`SUPPORTED_STYLES` first if they care to fail soft."""
    if style == "author-year":
        return format_author_year(fields)
    raise ValueError(
        f"unknown citation style {style!r}; supported: "
        f"{sorted(SUPPORTED_STYLES)}"
    )


# ---------------------------------------------------------------------------
# author-year (APA-like)
# ---------------------------------------------------------------------------

def format_author_year(fields: Dict[str, Any]) -> Dict[str, str]:
    """APA-flavoured author-year style.

    Reference-list entry:
        Totton, A. K. (1965). A synopsis of the Siphonophora.
        *Bulletin of the British Museum (Natural History)*, *7*, 1-200.
        https://doi.org/10.5962/bhl.part.20269

    Parenthetical inline:
        (Totton, 1965)            — one author
        (Smith & Jones, 2010)     — two authors
        (Smith et al., 2010)      — three or more
    """
    authors = fields.get("authors") or []
    year = fields.get("year")
    title = _str(fields.get("title"))
    journal = _str(fields.get("journal"))
    volume = _str(fields.get("volume"))
    pages = _str(fields.get("pages"))
    doi = _str(fields.get("doi"))

    year_str = f"({year})" if year else "(n.d.)"
    author_str = _format_author_list(authors)

    parts: List[str] = []
    if author_str:
        parts.append(f"{author_str} {year_str}.")
    else:
        parts.append(f"{year_str}.")

    if title:
        parts.append(_ensure_terminal_period(title))

    if journal:
        journal_chunk = f"*{journal}*"
        if volume:
            journal_chunk += f", *{volume}*"
        if pages:
            journal_chunk += f", {pages}"
        parts.append(journal_chunk + ".")

    if doi:
        parts.append(f"https://doi.org/{doi}")

    formatted = " ".join(p for p in parts if p).strip()

    inline_year = str(year) if year else "n.d."
    inline_auth = _inline_author_phrase(authors)
    inline = f"({inline_auth}, {inline_year})" if inline_auth else f"({inline_year})"

    return {"formatted": formatted, "inline": inline}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _str(v: Any) -> str:
    """Coerce a possibly-None field to a stripped string."""
    if v is None:
        return ""
    return str(v).strip()


def _initials(forename: Optional[str]) -> str:
    """``"Anne Marie" → "A. M.", "A. K." → "A. K.", "" → ""``."""
    if not forename:
        return ""
    parts = [p for p in re.split(r"[\s.]+", forename.strip()) if p]
    return " ".join(f"{p[0].upper()}." for p in parts)


def _format_author_list(authors: List[Dict[str, Any]]) -> str:
    """``[{surname, forename}, ...] → "Smith, J., Jones, K., & Lee, Q."``.

    APA-style: surname-comma-initials, Oxford comma + ``&`` before
    the last entry. Authors with no surname are skipped. An author
    with no forename surfaces as just the surname (no trailing
    comma).
    """
    formatted: List[str] = []
    for a in authors or []:
        surname = _str(a.get("surname"))
        if not surname:
            continue
        init = _initials(a.get("forename"))
        formatted.append(f"{surname}, {init}".strip().rstrip(","))
    if not formatted:
        return ""
    if len(formatted) == 1:
        return formatted[0]
    if len(formatted) == 2:
        return f"{formatted[0]}, & {formatted[1]}"
    return ", ".join(formatted[:-1]) + ", & " + formatted[-1]


def _inline_author_phrase(authors: List[Dict[str, Any]]) -> str:
    """Parenthetical-citation author phrase: ``Smith``, ``Smith &
    Jones``, ``Smith et al.``"""
    surnames = [
        _str(a.get("surname"))
        for a in (authors or [])
        if _str(a.get("surname"))
    ]
    if not surnames:
        return ""
    if len(surnames) == 1:
        return surnames[0]
    if len(surnames) == 2:
        return f"{surnames[0]} & {surnames[1]}"
    return f"{surnames[0]} et al."


def _ensure_terminal_period(s: str) -> str:
    s = s.strip()
    if not s:
        return ""
    return s if s.endswith((".", "!", "?")) else s + "."
