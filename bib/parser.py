"""Parse a BibTeX library and look up records by PDF filename.

Used by ``process_corpus.py`` when the user supplies ``--bib path/to/lib.bib``
to skip Grobid's header parser for any record whose ``file = {Foo.pdf}``
field matches an input PDF. References inside each PDF are still extracted
by Grobid; only the bibliographic *header* metadata (title, authors, year,
journal, DOI, abstract) is overridden.

The parser is intentionally minimal — it handles the subset of BibTeX
that machine-generated bibliographies emit: ``@type{key, field =
{value}, ...}`` with balanced-brace values. It is not a full bibtex
parser; entries that use ``"..."`` quoting or ``@string`` macros are
best-effort.
"""

from __future__ import annotations

import logging
import re
import unicodedata
from pathlib import Path
from typing import Dict, List, Optional


logger = logging.getLogger(__name__)


# Accent commands whose argument may abut them directly (``\'e``) as well
# as be braced (``\'{e}``). Mapped to Unicode combining marks and composed
# with NFC below, which is far smaller than an exhaustive letter table.
_LATEX_ACCENT_MARKS = {
    "`": "\u0300",   # grave
    "'": "\u0301",   # acute
    "^": "\u0302",   # circumflex
    '"': "\u0308",   # diaeresis
    "~": "\u0303",   # tilde
    "=": "\u0304",   # macron
    ".": "\u0307",   # dot above
}

# Accent commands spelled with a letter, which therefore *require* a brace
# or space before their argument (``\v{s}``, ``\c c``) — otherwise ``\cc``
# would be an unknown command rather than a cedilla.
_LATEX_LETTER_ACCENT_MARKS = {
    "u": "\u0306",   # breve
    "v": "\u030C",   # caron
    "H": "\u030B",   # double acute
    "c": "\u0327",   # cedilla
    "k": "\u0328",   # ogonek
    "r": "\u030A",   # ring above
}

# Standalone letter commands. Longest-first when substituted, so ``\oe``
# is not consumed as ``\o`` + a stray ``e``.
_LATEX_LETTERS = {
    r"\ss": "ß", r"\ae": "æ", r"\AE": "Æ", r"\oe": "œ", r"\OE": "Œ",
    r"\aa": "å", r"\AA": "Å", r"\o": "ø", r"\O": "Ø",
    r"\l": "ł", r"\L": "Ł", r"\i": "ı", r"\j": "ȷ",
    r"\dh": "ð", r"\DH": "Ð", r"\th": "þ", r"\TH": "Þ",
}

_ACCENT_RE = re.compile(
    r"\\([`'^\"~=.])\s*(?:\{([A-Za-z])\}|([A-Za-z]))"
)
_LETTER_ACCENT_RE = re.compile(
    r"\\([uvHckr])\s*(?:\{([A-Za-z])\}|\s([A-Za-z]))"
)
# Trailing whitespace is consumed because in LaTeX a space *terminates* a
# letter command rather than being content: ``\O rsted`` is "Ørsted", not
# "Ø rsted". A .bib that means to keep the space writes ``\ss{} and``,
# where the empty group ends the command and the space survives.
_LETTER_CMD_RE = re.compile(
    "(?:"
    + "|".join(re.escape(k) for k in sorted(_LATEX_LETTERS, key=len, reverse=True))
    + r")(?![A-Za-z])\s*"
)

# Characters a BibTeX writer must backslash-escape. Mirrored by
# ``bib/export.py``'s ``_BIBTEX_FIELD_ESCAPES`` so a round-trip is stable.
_ESCAPED_SPECIALS = ("&", "%", "_", "#", "$")

_WS_RE = re.compile(r"\s+")


def unescape_bib_value(s: str) -> str:
    """Convert one BibTeX field value to the plain text we serve.

    BibTeX values are LaTeX source, not text: ``&`` must be written
    ``\\&``, proper nouns are brace-protected (``{DNA}``), and accented
    letters are commands (``M{\\"u}ller``). Treating a value as plain text
    passes all of that through to ``metadata.json``, ``_serve/`` and every
    MCP tool that returns a title — so a client emitting a citation gets
    ``A.Braun \\& Vatke`` (#177). The ``.bib`` is not wrong in these cases;
    it is correctly escaped, and we were failing to decode it.

    Order matters. Accents are resolved first, because they are delimited
    by the very braces the last step drops; escaped specials next; braces
    only at the end.

    Brace removal is deliberately total and lossy: bibtex uses braces for
    case protection, which carries no meaning once the value is plain
    text. Escaped braces (``\\{`` / ``\\}``) unescape first so a value
    carrying one — Grobid OCR output does this, see #141 — leaves no stray
    backslash behind.
    """
    if not s:
        return ""
    out = s

    def _accent(match: re.Match, table: Dict[str, str]) -> str:
        letter = match.group(2) or match.group(3)
        return letter + table[match.group(1)]

    out = _ACCENT_RE.sub(lambda m: _accent(m, _LATEX_ACCENT_MARKS), out)
    out = _LETTER_ACCENT_RE.sub(
        lambda m: _accent(m, _LATEX_LETTER_ACCENT_MARKS), out
    )
    out = _LETTER_CMD_RE.sub(
        lambda m: _LATEX_LETTERS[m.group(0).rstrip()], out
    )

    for ch in _ESCAPED_SPECIALS:
        out = out.replace("\\" + ch, ch)

    # ``\\`` is a line break in bibtex, not a literal backslash.
    out = out.replace("\\\\", " ")
    out = out.replace("\\{", "{").replace("\\}", "}")
    out = out.replace("{", "").replace("}", "")

    # Compose the combining marks written above into precomposed
    # codepoints, so a served "Müller" compares equal to one typed
    # directly rather than differing by normalization form.
    out = unicodedata.normalize("NFC", out)
    return _WS_RE.sub(" ", out).strip()


# Retained under the old name because the behavior is a superset: every
# caller wanted "this bibtex value as text", which is what brace-stripping
# was approximating.
_strip_outer_braces = unescape_bib_value


def _scan_braced(text: str, open_idx: int) -> tuple:
    """Scan a brace-delimited run starting at ``text[open_idx] == "{"``.

    Returns ``(inner_start, close_idx, balanced)`` where ``close_idx`` is
    the index of the matching ``}`` (or ``len(text)`` when unbalanced) and
    ``balanced`` says whether a match was actually found.

    **A backslash escapes the next character** (#141). This is the root
    cause of the reported data loss: Grobid's OCR of a scanned reference
    list emits values like ``author = {Des, Ej\\{aims}``, and counting that
    ``\\{`` as an opening brace means the entry never closes. The parser
    then treated a *malformed value* as *the file ending here*. Note the
    quoted-string scanner below has always honored backslashes; only the
    brace path was missing it.
    """
    depth, j, n = 1, open_idx + 1, len(text)
    inner_start = open_idx + 1
    while j < n:
        c = text[j]
        if c == "\\":
            j += 2  # escaped char is literal — never counts as a brace
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return inner_start, j, True
        j += 1
    return inner_start, n, False


def _parse_fields(body: str, *, cite_key: str = "?") -> Dict[str, str]:
    """Parse the comma-separated field list inside an @entry{ ... } body.

    Expects ``body`` to be the substring *after* the entry key — i.e. the
    text between the first ``,`` (which separates the key from fields) and
    the closing ``}`` of the entry. Brace-balanced values are supported;
    quoted ``"..."`` values are too.

    An unbalanced brace inside a *value* used to swallow the rest of the
    entry with no warning at all — the quiet sibling of the entry-level
    bug in #141. It now logs, naming the entry and the field, so a
    half-parsed record is visible rather than merely wrong.
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
            start, j, balanced = _scan_braced(body, i)
            if not balanced:
                logger.warning(
                    "%s: unbalanced braces in field %r; the remaining fields "
                    "of this entry were not parsed",
                    cite_key, name or "<unnamed>",
                )
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


_TOP_LEVEL_AT_RE = re.compile(r"^[ \t]*@", re.MULTILINE)


def count_entry_markers(text: str) -> int:
    """Number of top-level ``@`` entry markers in a bib file body.

    Deliberately crude — it counts ``@`` at the start of a line, which is
    how every machine-generated bibliography we've seen formats entries.
    Used only to reconcile against the parsed count (#141), so
    over-counting a stray line-leading ``@`` inside a value produces a
    spurious warning at worst, never a dropped entry.
    """
    return len(_TOP_LEVEL_AT_RE.findall(text))


def _next_entry_start(text: str, from_idx: int) -> int:
    """Index of the next top-level ``@`` at or after ``from_idx``, or -1."""
    m = _TOP_LEVEL_AT_RE.search(text, from_idx)
    return m.end() - 1 if m else -1


def parse_bibtex(text: str) -> List[Dict]:
    """Parse a BibTeX file body into a list of entry dicts.

    Each entry has ``_type`` (e.g. ``"article"``), ``_key`` (the citation
    key), and one string-valued item per field. ``%`` comment lines and
    text outside of ``@type{...}`` blocks are skipped.

    **Malformed entries cost one entry, not the rest of the file** (#141).
    This used to ``break`` on the first unbalanced brace, so a single bad
    record silently discarded everything after it: on a real 19,834-entry
    export the parse stopped at ~1.75 MB and imported 2,258 entries, with
    one WARNING line and a summary that looked entirely plausible. Now the
    scan recovers at the next top-level ``@`` and, at the end, reconciles
    the parsed count against the number of ``@`` markers in the file so a
    shortfall is loud rather than invisible.
    """
    entries: List[Dict] = []
    skipped = 0
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
        _, j, balanced = _scan_braced(text, brace)
        if not balanced:
            # Recover rather than abandon the file. Resume at the next
            # top-level ``@`` *after this entry started*, so we make
            # forward progress even when the damage spans entries.
            skipped += 1
            resume = _next_entry_start(text, at + 1)
            logger.warning(
                "Unbalanced braces in entry near offset %d (type %r); "
                "skipping this entry and resuming at %s",
                at, entry_type,
                f"offset {resume}" if resume != -1 else "end of file",
            )
            if resume == -1:
                break
            i = resume
            continue
        body = text[brace + 1:j]
        comma = body.find(",")
        if comma == -1:
            i = j + 1
            continue
        key = body[:comma].strip()
        fields = _parse_fields(body[comma + 1:], cite_key=key or "?")
        fields["_type"] = entry_type
        fields["_key"] = key
        entries.append(fields)
        i = j + 1

    # Reconciliation (#141). The parsed count alone always looks
    # reasonable, which is what made the truncation invisible; comparing
    # it against the raw marker count is what makes a shortfall obvious.
    markers = count_entry_markers(text)
    if skipped:
        logger.warning(
            "bib parse: %d entries parsed, %d skipped as malformed "
            "(%d '@' markers in file)",
            len(entries), skipped, markers,
        )
    unaccounted = markers - len(entries) - skipped
    if unaccounted > 0:
        # Non-record markers (@comment/@preamble/@string) land here
        # legitimately, so this is a nudge rather than an alarm — but a
        # large gap means entries vanished without being counted.
        logger.warning(
            "bib parse: %d of %d '@' markers produced no entry and were not "
            "counted as malformed. @comment/@preamble/@string account for "
            "some; a large gap means entries were silently dropped.",
            unaccounted, markers,
        )
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
    # #51 + #54 — propagate license / skip metadata into metadata.json
    # so bib.authority can write it to the works.* columns.
    license_v = _strip_outer_braces(entry.get("license", ""))
    license_url = _strip_outer_braces(entry.get("licenseurl", ""))
    serve_raw = _strip_outer_braces(entry.get("serve", "")).strip().lower()
    serve_v: Optional[int] = None
    if serve_raw in {"0", "false", "no"}:
        serve_v = 0
    elif serve_raw in {"1", "true", "yes"}:
        serve_v = 1
    serve_reason = _strip_outer_braces(entry.get("servereason", ""))

    return {
        "filename": filename,
        "title": _strip_outer_braces(entry.get("title", "")),
        "authors": _split_authors(entry.get("author", "")),
        "year": year,
        "year_source": "bib" if year is not None else None,
        "journal": _strip_outer_braces(entry.get("journal", "")),
        "doi": _strip_outer_braces(entry.get("doi", "")),
        "abstract": _strip_outer_braces(entry.get("abstract", "")),
        "license": license_v or None,
        "license_url": license_url or None,
        "serve": serve_v,
        "serve_reason": serve_reason or None,
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
            # dialects; machine-generated bibs typically write a single
            # name, but split on ``;`` and ``,`` defensively so we don't
            # miss either side.
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
