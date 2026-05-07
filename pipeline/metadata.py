"""Bibliographic metadata extraction stage (Grobid + bib fallback).

* :func:`extract_metadata` — calls Grobid on the processed PDF, falls
  back to a placeholder + a BibTeX-driven override when Grobid is
  unreachable. Writes ``metadata.json``, ``references.json``,
  ``intext_citations.json``, and (when available) ``grobid.tei.xml``.
* :func:`_year_from_filename` — extracts a 4-digit year from a PDF
  filename when Grobid emits an empty ``<date/>``.
* :func:`_write_placeholder_metadata` — degrades-gracefully fallback
  for runs that proceed without Grobid metadata.
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Optional

from bib import bib_entry_to_metadata

from . import stamp_artifact
from .grobid_client import (
    GrobidClient,
    GrobidUnavailableError,
    parse_tei_header,
    parse_tei_intext_citations,
    parse_tei_references,
)

logger = logging.getLogger(__name__)


# Years plausible for a published scientific paper: 17xx (earliest
# Linnaean works) to the current decade. Tighter than a generic `\d{4}`
# so e.g. a trailing "20" in a specimen ID isn't picked up. Lookarounds
# exclude matches that are part of a longer digit run (e.g., "12005"
# shouldn't match 2005).
_FILENAME_YEAR_RE = re.compile(r"(?<!\d)(17\d{2}|18\d{2}|19\d{2}|20[0-4]\d)(?!\d)")


def _year_from_filename(name: str) -> Optional[int]:
    """Best-effort pub-year extraction from a PDF filename.

    Matches the first 4-digit year in the range 1700–2049, covering
    historical scientific literature without false-positiving on
    generic numbers (ISSNs, specimen counts, etc.). Returns None if no
    match. Useful when the input library uses an "Author(s)YYYY"
    naming convention and Grobid's header parser emits an empty
    ``<date/>``.
    """
    if not name:
        return None
    m = _FILENAME_YEAR_RE.search(name)
    if not m:
        return None
    return int(m.group(1))


def _write_placeholder_metadata(
    pdf_path: Path,
    metadata_output: Path,
    references_output: Path,
    original_filename: Optional[str] = None,
):
    """Write empty-but-valid metadata.json and references.json.

    Used when Grobid is unavailable or its TEI can't be parsed. Keeps
    downstream stages (chunking, embedding) functional and lets the
    summary show this document as metadata-less so we can triage later.
    Even in the placeholder path we attempt to recover a year from the
    filename — it's all the signal we have without Grobid.
    """
    effective_filename = original_filename or pdf_path.name
    fname_year = _year_from_filename(effective_filename)
    placeholder = {
        "filename": effective_filename,
        "title": "",
        "authors": [],
        "year": fname_year,
        "year_source": "filename" if fname_year is not None else None,
        "journal": "",
        "doi": "",
        "abstract": "",
        "extraction_method": "placeholder",
    }
    with open(metadata_output, "w", encoding="utf-8") as f:
        json.dump(stamp_artifact(placeholder), f, indent=2, ensure_ascii=False)
    with open(references_output, "w", encoding="utf-8") as f:
        json.dump(stamp_artifact({"references": [], "total_references": 0}), f, indent=2)
    with open(metadata_output.parent / "intext_citations.json", "w", encoding="utf-8") as f:
        json.dump(stamp_artifact({"paragraphs": [], "citations": []}), f, indent=2)


def extract_metadata(
    pdf_path: Path,
    metadata_output: Path,
    references_output: Optional[Path] = None,
    tei_output: Optional[Path] = None,
    grobid_client: Optional[GrobidClient] = None,
    original_filename: Optional[str] = None,
    bib_entry: Optional[Dict] = None,
):
    """Extract bibliographic metadata + references via Grobid.

    Runs ``/api/processFulltextDocument`` once, caches the raw TEI-XML at
    ``tei_output`` (if given), and parses it into ``metadata.json`` plus
    ``references.json``. If Grobid is unreachable or any step fails, falls
    back to placeholder files so downstream stages still run — the caller
    can inspect ``metadata["extraction_method"]`` to tell which path was
    taken.

    When ``bib_entry`` is supplied (a parsed BibTeX record matched to this
    PDF by filename), its title/authors/year/journal/DOI override the
    Grobid header — Grobid is still called for the reference list so we
    can build the citation graph, but the header parse is skipped. This
    lets a curated bibliography be the source of truth for header fields
    while keeping Grobid's reference extraction.

    Parameters
    ----------
    pdf_path:
        The (possibly OCR'd) ``processed.pdf`` to send to Grobid.
    metadata_output:
        Output path for ``metadata.json`` (title, authors, year, …).
    references_output:
        Output path for ``references.json``. If None, defaults to
        ``<metadata_output.parent>/references.json``.
    tei_output:
        Cache path for the raw Grobid TEI. Existing non-empty TEI at this
        path is reused (skipping the Grobid call) — convenient for
        re-parsing without re-processing. If None, defaults to
        ``<metadata_output.parent>/grobid.tei.xml``.
    grobid_client:
        A live :class:`GrobidClient`, or None to skip Grobid entirely
        (placeholder-only mode).
    bib_entry:
        Parsed BibTeX entry from :class:`bib.BibIndex`. When
        present, overrides Grobid's header parse for this document.
    """
    hash_dir = metadata_output.parent
    if references_output is None:
        references_output = hash_dir / "references.json"
    if tei_output is None:
        tei_output = hash_dir / "grobid.tei.xml"
    # Prefer the original PDF filename for both provenance and year-fallback
    # — by the time this runs, pdf_path is usually "processed.pdf" (post-OCR
    # or copied), whose name carries no year information.
    effective_filename = original_filename or pdf_path.name

    tei_xml: Optional[str] = None

    # 1. Reuse cached TEI if present.
    if tei_output.exists() and tei_output.stat().st_size > 0:
        try:
            tei_xml = tei_output.read_text(encoding="utf-8")
            logger.info("Reusing cached Grobid TEI at %s", tei_output)
        except OSError as e:
            logger.warning("Could not read cached TEI %s: %s", tei_output, e)
            tei_xml = None

    # 2. Otherwise, call Grobid if a client is available.
    if tei_xml is None and grobid_client is not None:
        try:
            logger.info("Calling Grobid on %s ...", pdf_path.name)
            tei_xml = grobid_client.process_fulltext(pdf_path)
            tei_output.write_text(tei_xml, encoding="utf-8")
            logger.info("Wrote Grobid TEI to %s", tei_output)
        except GrobidUnavailableError as e:
            logger.warning("Grobid unavailable, using placeholder: %s", e)
        except Exception as e:
            logger.warning("Grobid call failed on %s: %s", pdf_path.name, e)

    # 3. Parse references from TEI if available — independent of which
    # source we use for the header. Anything that fails leaves us with
    # an empty reference list (written below as a fallback).
    refs: List[dict] = []
    refs_parsed = False
    intext = {"paragraphs": [], "citations": []}
    if tei_xml is not None:
        try:
            refs = parse_tei_references(tei_xml)
            refs_parsed = True
        except Exception as e:
            logger.warning("Failed to parse Grobid TEI references (%s)", e)
        # In-text citation graph (issue #7).  Independent of refs parse —
        # we want partial recovery when one fails and the other doesn't.
        try:
            intext = parse_tei_intext_citations(tei_xml)
        except Exception as e:
            logger.warning("Failed to parse Grobid TEI in-text citations (%s)", e)

    # 4. Build the header. Bib record wins if present; otherwise Grobid's
    # header parse, with filename-year fallback. If both fail, placeholder.
    header: Optional[Dict] = None
    if bib_entry is not None:
        header = bib_entry_to_metadata(bib_entry, effective_filename)
        logger.info(
            "Using bib entry %r for header (overrides Grobid)",
            header.get("bib_key"),
        )
    elif tei_xml is not None:
        try:
            header = parse_tei_header(tei_xml)
            header["filename"] = effective_filename
            # Year fallback: Grobid emits <date/> (empty) on many papers
            # whose title page doesn't match its header template — especially
            # historical and non-English. The Pugh library's "AuthorYYYY"
            # filename convention lets us recover the year cheaply. Record
            # which source won in ``year_source`` so downstream consumers can
            # treat Grobid-extracted years (stronger) differently from
            # filename-derived years (weaker, depends on filename convention).
            if header.get("year") is None:
                fname_year = _year_from_filename(effective_filename)
                if fname_year is not None:
                    header["year"] = fname_year
                    header["year_source"] = "filename"
                    logger.info(
                        "Grobid returned year=None; recovered %d from filename %r",
                        fname_year, effective_filename,
                    )
                else:
                    header["year_source"] = None
            else:
                header["year_source"] = "grobid"
        except Exception as e:
            logger.warning(
                "Failed to parse Grobid TEI header (%s); writing placeholder", e
            )
            header = None

    # 5. Write outputs. If we got a header, write it + the reference list
    # (which may be empty if Grobid was unavailable). Otherwise placeholder.
    if header is not None:
        with open(metadata_output, "w", encoding="utf-8") as f:
            json.dump(stamp_artifact(header), f, indent=2, ensure_ascii=False)
        with open(references_output, "w", encoding="utf-8") as f:
            json.dump(
                stamp_artifact(
                    {"references": refs, "total_references": len(refs)},
                ),
                f,
                indent=2,
                ensure_ascii=False,
            )
        with open(hash_dir / "intext_citations.json", "w", encoding="utf-8") as f:
            json.dump(stamp_artifact(intext), f, indent=2, ensure_ascii=False)
        logger.info(
            "Wrote metadata (method=%s, %d authors), %d references (parsed=%s), %d in-text citations",
            header.get("extraction_method"),
            len(header.get("authors", [])),
            len(refs),
            refs_parsed,
            len(intext["citations"]),
        )
        return

    # 6. Placeholder path.
    _write_placeholder_metadata(
        pdf_path, metadata_output, references_output,
        original_filename=effective_filename,
    )
