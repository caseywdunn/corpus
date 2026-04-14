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


def parse_tei_references(tei_xml: str) -> List[dict]:
    """Extract the bibliography from Grobid TEI-XML.

    Returns a list of dicts, one per cited work, with best-effort fields:
    title, authors, year, journal, doi, raw (the original citation string
    if the TEI preserved it via includeRawCitations=1).
    """
    root = _parse_tei(tei_xml)
    refs: List[dict] = []

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

        refs.append({
            "xml_id": bs.get("{http://www.w3.org/XML/1998/namespace}id") or "",
            "title": title,
            "authors": authors,
            "year": year,
            "journal": journal,
            "doi": doi,
            "raw": raw,
        })

    return refs
