"""Parse a BibTeX library and look up records by PDF filename.

Used by ``process_corpus.py`` when the user supplies ``--bib path/to/lib.bib``
to skip Grobid's header parser for any record whose ``file = {Foo.pdf}``
field matches an input PDF. References inside each PDF are still extracted
by Grobid; only the bibliographic *header* metadata (title, authors, year,
journal, DOI, abstract) is overridden.

The parser is intentionally minimal — it handles the subset of BibTeX that
``siphonophores.bib`` (machine-generated) emits: ``@type{key, field = {value},
...}`` with balanced-brace values. It is not a full bibtex parser; entries
that use ``"..."`` quoting or ``@string`` macros are best-effort.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional


logger = logging.getLogger(__name__)


def _strip_outer_braces(s: str) -> str:
    """Remove any braces used purely for bibtex grouping/escaping.

    Bibtex authors wrap protected casing or accented chars in ``{}``; the
    JSON metadata consumers don't want them, so we drop all ``{`` / ``}``
    after parsing. This is lossy but correct for the fields we care about.
    """
    return s.replace("{", "").replace("}", "").strip()


def _parse_fields(body: str) -> Dict[str, str]:
    """Parse the comma-separated field list inside an @entry{ ... } body.

    Expects ``body`` to be the substring *after* the entry key — i.e. the
    text between the first ``,`` (which separates the key from fields) and
    the closing ``}`` of the entry. Brace-balanced values are supported;
    quoted ``"..."`` values are too.
    """
    fields: Dict[str, str] = {}
    i, n = 0, len(body)
    while i < n:
        while i < n and body[i] in " \t\r\n,":
            i += 1
        if i >= n:
            break
        eq = body.find("=", i)
        if eq == -1:
            break
        name = body[i:eq].strip().lower()
        i = eq + 1
        while i < n and body[i] in " \t\r\n":
            i += 1
        if i >= n:
            break
        if body[i] == "{":
            depth, start, j = 1, i + 1, i + 1
            while j < n and depth > 0:
                c = body[j]
                if c == "{":
                    depth += 1
                elif c == "}":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            value = body[start:j]
            i = j + 1
        elif body[i] == '"':
            j = i + 1
            while j < n and body[j] != '"':
                if body[j] == "\\":
                    j += 1
                j += 1
            value = body[i + 1:j]
            i = j + 1
        else:
            j = i
            while j < n and body[j] != ",":
                j += 1
            value = body[i:j].strip()
            i = j
        if name:
            fields[name] = value.strip()
    return fields


def parse_bibtex(text: str) -> List[Dict]:
    """Parse a BibTeX file body into a list of entry dicts.

    Each entry has ``_type`` (e.g. ``"article"``), ``_key`` (the citation
    key), and one string-valued item per field. ``%`` comment lines and
    text outside of ``@type{...}`` blocks are skipped.
    """
    entries: List[Dict] = []
    i, n = 0, len(text)
    while i < n:
        at = text.find("@", i)
        if at == -1:
            break
        # Skip ``@`` inside obvious comment lines (``%`` to end-of-line).
        line_start = text.rfind("\n", 0, at) + 1
        if "%" in text[line_start:at]:
            i = text.find("\n", at)
            if i == -1:
                break
            continue
        brace = text.find("{", at)
        if brace == -1:
            break
        entry_type = text[at + 1:brace].strip().lower()
        # ``@comment{...}`` and ``@preamble{...}`` are not records.
        if entry_type in ("comment", "preamble", "string"):
            i = brace + 1
            continue
        depth, j = 1, brace + 1
        while j < n and depth > 0:
            c = text[j]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if depth != 0:
            logger.warning("Unbalanced braces near offset %d in bib; stopping parse", at)
            break
        body = text[brace + 1:j]
        comma = body.find(",")
        if comma == -1:
            i = j + 1
            continue
        key = body[:comma].strip()
        fields = _parse_fields(body[comma + 1:])
        fields["_type"] = entry_type
        fields["_key"] = key
        entries.append(fields)
        i = j + 1
    return entries


_AUTHOR_SPLIT_RE = re.compile(r"\s+and\s+", re.IGNORECASE)


def _split_authors(auth_str: str) -> List[Dict]:
    """Convert a bibtex author field into the dict-list format used in
    ``metadata.json`` (matching :func:`grobid_client.parse_tei_header`).

    Handles the common ``Last, First and Last, First`` form; falls back to
    treating a comma-less name as a surname-only entry (better than
    fabricating a forename).
    """
    if not auth_str:
        return []
    out: List[Dict] = []
    for raw in _AUTHOR_SPLIT_RE.split(auth_str):
        name = _strip_outer_braces(raw)
        if not name:
            continue
        if "," in name:
            surname, _, forename = name.partition(",")
            out.append({
                "forename": forename.strip(),
                "surname": surname.strip(),
                "affiliations": [],
            })
        else:
            out.append({
                "forename": "",
                "surname": name,
                "affiliations": [],
            })
    return out


_YEAR_RE = re.compile(r"\b(\d{4})\b")


def bib_entry_to_metadata(entry: Dict, filename: str) -> Dict:
    """Build a ``metadata.json`` dict from a parsed bibtex entry.

    Schema matches :func:`grobid_client.parse_tei_header` plus the
    ``filename`` / ``year_source`` fields the pipeline appends, so
    downstream consumers don't need to know which path produced the
    metadata. ``extraction_method`` is set to ``"bib"`` and the citation
    key is preserved as ``bib_key`` for traceability.
    """
    year: Optional[int] = None
    yr = (entry.get("year") or "").strip()
    if yr:
        m = _YEAR_RE.search(yr)
        if m:
            year = int(m.group(1))
    return {
        "filename": filename,
        "title": _strip_outer_braces(entry.get("title", "")),
        "authors": _split_authors(entry.get("author", "")),
        "year": year,
        "year_source": "bib" if year is not None else None,
        "journal": _strip_outer_braces(entry.get("journal", "")),
        "doi": _strip_outer_braces(entry.get("doi", "")),
        "abstract": _strip_outer_braces(entry.get("abstract", "")),
        "extraction_method": "bib",
        "bib_key": entry.get("_key", ""),
    }


class BibIndex:
    """Load a .bib file and index its entries by ``file = {...}`` value.

    Lookup is case-insensitive on the PDF filename (basename only). When
    multiple entries reference the same filename, the *last* one wins and
    a warning is logged — historically rare in our corpus, but possible
    with imperfect curation.
    """

    def __init__(self, entries: List[Dict]):
        self.entries = entries
        self._by_filename: Dict[str, Dict] = {}
        for e in entries:
            fname = (e.get("file") or "").strip()
            if not fname:
                continue
            # ``file`` can be a list like ``Foo.pdf;Bar.pdf`` in some bibtex
            # dialects; siphonophores.bib uses a single name, but split on
            # ``;`` and ``,`` defensively so we don't miss either side.
            for part in re.split(r"[;,]", fname):
                key = Path(part.strip()).name.lower()
                if not key:
                    continue
                if key in self._by_filename:
                    logger.warning(
                        "Bib filename collision on %r: keys %r and %r both claim it",
                        key, self._by_filename[key].get("_key"), e.get("_key"),
                    )
                self._by_filename[key] = e

    @classmethod
    def from_path(cls, path: Path) -> "BibIndex":
        text = Path(path).read_text(encoding="utf-8")
        entries = parse_bibtex(text)
        idx = cls(entries)
        logger.info(
            "Loaded %d bib entries (%d with file= field) from %s",
            len(entries), len(idx._by_filename), path,
        )
        return idx

    def lookup(self, filename: str) -> Optional[Dict]:
        """Return the entry whose ``file = {...}`` matches ``filename``.

        Matching is on the basename, lower-cased. Returns None if no entry
        is registered for that filename.
        """
        if not filename:
            return None
        return self._by_filename.get(Path(filename).name.lower())

    def __len__(self) -> int:
        return len(self.entries)

    def __contains__(self, filename: str) -> bool:
        return self.lookup(filename) is not None
