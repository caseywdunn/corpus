"""Minimal Grobid client + TEI-XML parsers.

Grobid is a Java service that extracts bibliographic metadata and structure
from scientific PDFs. See ``docker-compose.yml`` for the local dev setup.

This module does two things:

1. :class:`GrobidClient` — a thin requests-based wrapper around the
   ``/api/processFulltextDocument`` endpoint. Falls back gracefully (raises
   :class:`GrobidUnavailableError`) when the service is not reachable so
   dev iteration without Docker running still produces usable output.

2. :func:`parse_tei_header` and :func:`parse_tei_references` — extract the
   structured bits we need (title, authors, year, journal, DOI, abstract;
   bibliography) from Grobid's TEI-XML, using lxml.

We intentionally do not pull in ``grobid-client-python`` — the endpoints
are simple enough to call directly and each added dependency is a thing to
maintain.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import requests
from lxml import etree

logger = logging.getLogger(__name__)

# TEI namespace used by Grobid's output
TEI_NS = "http://www.tei-c.org/ns/1.0"
NSMAP = {"tei": TEI_NS}


class GrobidUnavailableError(RuntimeError):
    """Raised when the Grobid service is unreachable.

    Distinct from per-request failures (timeouts, 5xx) so the caller can
    choose to fall back to a placeholder metadata step instead of aborting
    the whole run.
    """


class GrobidClient:
    """Thin wrapper over the subset of the Grobid REST API we use.

    Parameters
    ----------
    base_url:
        e.g. ``"http://localhost:8070"`` — the Grobid service URL.
    timeout:
        Per-request timeout in seconds. Grobid can be slow on long papers
        (minutes for a Haeckel monograph); set generously.
    """

    def __init__(self, base_url: str = "http://localhost:8070", timeout: float = 300.0):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def is_alive(self) -> bool:
        """Return True if the Grobid ``/api/isalive`` endpoint responds OK.

        Never raises — any connection error is caught and returned as
        False. Callers use this to decide whether to proceed with Grobid
        or to fall back.
        """
        try:
            r = requests.get(f"{self.base_url}/api/isalive", timeout=5)
            return r.status_code == 200 and r.text.strip().lower() == "true"
        except requests.RequestException as e:
            logger.debug("Grobid /isalive check failed: %s", e)
            return False

    def process_fulltext(
        self,
        pdf_path: Path,
        *,
        consolidate_header: int = 1,
        consolidate_citations: int = 0,
        include_raw_citations: bool = False,
    ) -> str:
        """Run ``/api/processFulltextDocument`` on a PDF and return TEI-XML.

        Parameters
        ----------
        consolidate_header, consolidate_citations:
            Grobid consolidation levels (0=off, 1=CrossRef lookup). Header
            consolidation meaningfully improves metadata quality for papers
            with recoverable DOIs; citation consolidation is slower and
            mostly useful at corpus scale — leave off by default.
        include_raw_citations:
            Add the ``includeRawCitations=1`` flag so the original citation
            string is preserved in the TEI.

        Raises
        ------
        GrobidUnavailableError
            Connection could not be established (service not running).
        RuntimeError
            Service returned a non-2xx response.
        """
        url = f"{self.base_url}/api/processFulltextDocument"
        params = {
            "consolidateHeader": str(consolidate_header),
            "consolidateCitations": str(consolidate_citations),
        }
        if include_raw_citations:
            params["includeRawCitations"] = "1"

        try:
            with open(pdf_path, "rb") as f:
                files = {"input": (pdf_path.name, f, "application/pdf")}
                r = requests.post(url, files=files, data=params, timeout=self.timeout)
        except requests.ConnectionError as e:
            raise GrobidUnavailableError(
                f"Cannot connect to Grobid at {self.base_url}: {e}"
            ) from e
        except requests.Timeout as e:
            raise RuntimeError(
                f"Grobid request timed out after {self.timeout}s on {pdf_path.name}"
            ) from e

        if r.status_code != 200:
            raise RuntimeError(
                f"Grobid returned HTTP {r.status_code} on {pdf_path.name}: "
                f"{r.text[:300]}"
            )
        return r.text


# ---------------------------------------------------------------------------
# TEI-XML parsing helpers
# ---------------------------------------------------------------------------


def _first_text(el) -> Optional[str]:
    """Collapse an lxml element's descendant text to a single whitespace-
    normalized string, or None if empty.
    """
    if el is None:
        return None
    text = " ".join(el.itertext())
    text = " ".join(text.split())
    return text or None


def _parse_tei(tei_xml: str):
    """Parse a TEI-XML string. Returns the root element.

    Wraps lxml errors in a clearer message.
    """
    try:
        return etree.fromstring(tei_xml.encode("utf-8"))
    except etree.XMLSyntaxError as e:
        raise RuntimeError(f"Could not parse TEI-XML: {e}") from e


def parse_tei_header(tei_xml: str) -> dict:
    """Extract bibliographic metadata from Grobid TEI-XML.

    Returns a dict with keys: title, authors (list of dicts), year, journal,
    doi, abstract. Missing fields are empty strings / lists / None; the
    caller decides how to treat absence.
    """
    root = _parse_tei(tei_xml)

    title_el = root.find(".//tei:titleStmt/tei:title", NSMAP)
    # Prefer biblStruct/analytic/title (the article-level title); fall back
    # to titleStmt/title which Grobid sometimes duplicates.
    analytic_title = root.find(".//tei:sourceDesc/tei:biblStruct/tei:analytic/tei:title", NSMAP)
    title = _first_text(analytic_title) or _first_text(title_el) or ""

    # Authors under analytic (article authors, not cited-work authors).
    authors: List[dict] = []
    for author_el in root.findall(".//tei:sourceDesc/tei:biblStruct/tei:analytic/tei:author", NSMAP):
        forename = _first_text(author_el.find("tei:persName/tei:forename", NSMAP))
        surname = _first_text(author_el.find("tei:persName/tei:surname", NSMAP))
        if not (forename or surname):
            continue
        affiliations = [
            _first_text(a) for a in author_el.findall("tei:affiliation", NSMAP)
        ]
        authors.append({
            "forename": forename or "",
            "surname": surname or "",
            "affiliations": [a for a in affiliations if a],
        })

    # Year: sourceDesc/biblStruct/monogr/imprint/date[@when]
    year: Optional[int] = None
    date_el = root.find(".//tei:sourceDesc/tei:biblStruct/tei:monogr/tei:imprint/tei:date", NSMAP)
    if date_el is not None:
        when = date_el.get("when") or _first_text(date_el) or ""
        # ``when`` looks like "1891", "1891-05", "1891-05-10"
        if when:
            try:
                year = int(when[:4])
            except ValueError:
                year = None

    # Journal / source title: monogr/title[@level='j']
    journal = _first_text(root.find(
        ".//tei:sourceDesc/tei:biblStruct/tei:monogr/tei:title[@level='j']", NSMAP
    )) or ""

    # DOI: analytic/idno[@type='DOI']
    doi = _first_text(root.find(
        ".//tei:sourceDesc/tei:biblStruct/tei:analytic/tei:idno[@type='DOI']", NSMAP
    )) or ""

    # Abstract: profileDesc/abstract
    abstract = _first_text(root.find(".//tei:profileDesc/tei:abstract", NSMAP)) or ""

    return {
        "title": title,
        "authors": authors,
        "year": year,
        "journal": journal,
        "doi": doi,
        "abstract": abstract,
        "extraction_method": "grobid",
    }


import re as _re

# --- Post-processing helpers for the reference list ---

# Titles that start with any of these are almost certainly not real citations —
# they're the author-address / ORCID / correspondence block at the end of the
# paper that Grobid's parser sometimes misclassifies as the next reference.
_ADDRESS_TITLE_RE = _re.compile(
    r"^\s*(?:addresses?|correspondence|corresponding author|funding|"
    r"acknowledg(?:e?ment)s?|competing interest|author contributions|"
    r"orcid|ethics statement|data availability|conflict of interest)[:.]?\s",
    _re.IGNORECASE,
)

# Author strings that are the continuation placeholder ("same as previous ref")
# scholarly papers often use em-dashes, hyphens, or underscores in the authors
# position for chained refs. Grobid captures them literally; we want to reuse
# the previous ref's authors instead.
_EM_DASH_AUTHORS_RE = _re.compile(r"^[\s_\-\u2012-\u2015\u2E3A\u2E3B]+$")

# Tokens that turn up as author "names" when Grobid fumbles a header / address
# block. Used to filter bogus refs whose authors list is dominated by these.
_BOGUS_AUTHOR_TOKENS = frozenset(
    s.lower() for s in ["Submitted", "Received", "Accepted", "Corresponding"]
)


def _looks_like_address_block(ref: dict) -> bool:
    """Return True if the ref is actually an author-address / end-matter block.

    The Marrus 2005 case in the demo set is the canonical example: Grobid
    turned the Sandholdt Road / NOC Southampton author addresses into a
    spurious reference. These entries share a signature — an ``ADDRESSES:``
    title, authors like ``Submitted``, or a journal value containing the
    institution block's tail ("P.R.P.) National Oceanography Centre" in the
    Marrus case).
    """
    title = ref.get("title") or ""
    if _ADDRESS_TITLE_RE.match(title):
        return True
    authors = ref.get("authors") or []
    if authors and all(a.lower() in _BOGUS_AUTHOR_TOKENS for a in authors):
        return True
    # A journal field that contains a parenthesized author initial + address
    # noise is another tell-tale for this kind of mis-parse.
    journal = ref.get("journal") or ""
    if _re.search(r"[A-Z]\.[A-Z]\.\)\s", journal):
        return True
    return False


def _resolve_em_dash_authors(refs: List[dict]) -> None:
    """In-place fix: when a ref's authors are literal em-dashes / underscores
    (the scholarly convention for 'same author as previous'), replace them
    with the previous ref's authors.

    If the previous ref is missing authors too, we walk back until we find
    one. If no prior ref has real authors, the field is left empty.
    """
    last_real_authors: List[str] = []
    for ref in refs:
        authors = ref.get("authors") or []
        # All tokens are dash/underscore placeholders?
        if authors and all(_EM_DASH_AUTHORS_RE.match(a) for a in authors):
            ref["authors"] = list(last_real_authors)
            ref["authors_continued_from_previous"] = True
        elif authors:
            last_real_authors = list(authors)


def parse_tei_references(tei_xml: str) -> List[dict]:
    """Extract the bibliography from Grobid TEI-XML.

    Returns a list of dicts, one per cited work, with best-effort fields:
    title, authors, year, journal, doi, raw (the original citation string
    if the TEI preserved it via includeRawCitations=1).

    Post-processing:
      * Address / end-matter blocks mistakenly parsed as references are
        dropped (see :func:`_looks_like_address_block`).
      * Em-dash / underscore author continuations are resolved to the
        previous reference's author list, matching how the paper reads.

    Dropped entries are counted in the caller's log but not emitted to the
    returned list.
    """
    root = _parse_tei(tei_xml)
    raw_refs: List[dict] = []

    # listBibl can appear under <back> or elsewhere; take all.
    for bs in root.findall(".//tei:listBibl/tei:biblStruct", NSMAP):
        # Title: prefer analytic/title (the article), fall back to
        # monogr/title (the book/journal).
        title = _first_text(bs.find("tei:analytic/tei:title", NSMAP)) or \
                _first_text(bs.find("tei:monogr/tei:title", NSMAP)) or ""

        authors: List[str] = []
        for ael in bs.findall(".//tei:author", NSMAP):
            forename = _first_text(ael.find("tei:persName/tei:forename", NSMAP)) or ""
            surname = _first_text(ael.find("tei:persName/tei:surname", NSMAP)) or ""
            name = (forename + " " + surname).strip()
            if name:
                authors.append(name)

        year = None
        date_el = bs.find(".//tei:imprint/tei:date", NSMAP)
        if date_el is not None:
            when = date_el.get("when") or _first_text(date_el) or ""
            if when:
                try:
                    year = int(when[:4])
                except ValueError:
                    year = None

        journal = _first_text(bs.find("tei:monogr/tei:title[@level='j']", NSMAP)) or ""
        doi = _first_text(bs.find(".//tei:idno[@type='DOI']", NSMAP)) or ""
        # Raw citation (only present if includeRawCitations=1 was requested)
        raw = _first_text(bs.find("tei:note[@type='raw_reference']", NSMAP)) or ""

        raw_refs.append({
            "xml_id": bs.get("{http://www.w3.org/XML/1998/namespace}id") or "",
            "title": title,
            "authors": authors,
            "year": year,
            "journal": journal,
            "doi": doi,
            "raw": raw,
        })

    # Em-dash author continuation: walk the list in-order BEFORE filtering.
    _resolve_em_dash_authors(raw_refs)

    # Drop address blocks and other non-citation residue.
    refs: List[dict] = []
    dropped = 0
    for ref in raw_refs:
        if _looks_like_address_block(ref):
            dropped += 1
            logger.info(
                "Dropping non-citation TEI entry %s: title=%r authors=%r",
                ref.get("xml_id"), (ref.get("title") or "")[:80],
                ref.get("authors"),
            )
            continue
        refs.append(ref)

    if dropped:
        logger.info("parse_tei_references: dropped %d non-citation entr%s",
                    dropped, "y" if dropped == 1 else "ies")

    return refs


def parse_tei_intext_citations(tei_xml: str) -> dict:
    """Extract the in-text citation graph (body → bibliography) from TEI.

    Walks every ``<ref type="bibr">`` in ``<text><body>`` and emits one
    record per citation, with the surface text, the enclosing section
    heading, and a pointer into a deduplicated paragraph list.

    The returned dict has two top-level fields:

    * ``"paragraphs"`` — list of paragraph text strings.  Multiple
      citations in the same paragraph share an entry rather than
      duplicating the (often long) text.
    * ``"citations"`` — list of dicts:

      * ``target_xml_id`` — ``"#b5"``-style ID into ``references.json``'s
        ``xml_id`` field.  ``None`` when Grobid couldn't resolve the
        surface text to a bibliography entry (~40% of cases in our
        sample); the surface is preserved for later fuzzy resolution.
      * ``surface`` — original text from the TEI ref element
        (e.g. ``"Furnestin, 1960"``).
      * ``section`` — text of the nearest ancestor ``<div>``'s heading,
        empty if none.
      * ``para_index`` — index into ``paragraphs``.

    The data is already on disk in ``grobid.tei.xml``; this is pure
    post-processing.  See issue #7.
    """
    root = _parse_tei(tei_xml)

    para_text_to_index: dict = {}
    paragraphs: List[str] = []
    citations: List[dict] = []

    # Restrict to body refs — header / listBibl have their own ref
    # elements that aren't in-text citations.
    for r in root.findall(".//tei:text/tei:body//tei:ref[@type='bibr']", NSMAP):
        target = r.get("target") or None
        surface = "".join(r.itertext()).strip()
        if not surface:
            continue

        # Section heading: nearest enclosing <div>'s <head>.  Grobid
        # nests divs (subsection within section); take the deepest div
        # whose head has text, otherwise the outermost.
        section = ""
        for div in r.iterancestors("{%s}div" % TEI_NS):
            head = div.find("tei:head", NSMAP)
            head_text = "".join(head.itertext()).strip() if head is not None else ""
            if head_text:
                section = head_text
                break

        # Enclosing paragraph text — the "excerpt" the user actually
        # wants when asking "show me passages citing paper X".
        para_el = next(r.iterancestors("{%s}p" % TEI_NS), None)
        if para_el is None:
            continue
        para_text = " ".join("".join(para_el.itertext()).split())
        if not para_text:
            continue

        para_index = para_text_to_index.get(para_text)
        if para_index is None:
            para_index = len(paragraphs)
            paragraphs.append(para_text)
            para_text_to_index[para_text] = para_index

        citations.append({
            "target_xml_id": target,
            "surface": surface,
            "section": section,
            "para_index": para_index,
        })

    return {"paragraphs": paragraphs, "citations": citations}
