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
import hashlib
import logging
import re
import uuid
from pathlib import Path
from typing import Dict, List, Optional

from bib import bib_entry_to_metadata

from . import stamp_artifact
from . import external
from .config import CONFIG
from .stages import _file_sha256
from .version import __version__
from .grobid_client import (
    GrobidClient,
    GrobidUnavailableError,
    parse_tei_header,
    parse_tei_intext_citations,
    parse_tei_references,
    validate_fulltext_tei,
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
    grobid_result: Optional[Dict] = None,
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
    if grobid_result is not None:
        placeholder["grobid"] = grobid_result
    with open(metadata_output, "w", encoding="utf-8") as f:
        json.dump(stamp_artifact(placeholder), f, indent=2, ensure_ascii=False)
    with open(references_output, "w", encoding="utf-8") as f:
        json.dump(stamp_artifact({"references": [], "total_references": 0}), f, indent=2)
    with open(metadata_output.parent / "intext_citations.json", "w", encoding="utf-8") as f:
        json.dump(stamp_artifact({"paragraphs": [], "citations": []}), f, indent=2)


def _archive_tei(tei_output: Path, receipt_path: Path, hash_dir: Path):
    """Keep superseded evidence out of the active post-processing path."""
    if not tei_output.exists():
        return
    history = hash_dir / "metadata_cache_history"
    history.mkdir(exist_ok=True)
    archive_id = uuid.uuid4().hex
    tei_output.rename(history / f"{archive_id}-{tei_output.name}")
    if receipt_path.exists():
        receipt_path.rename(history / f"{archive_id}-{receipt_path.name}")
    logger.info("Archived disabled/unverified/stale Grobid cache under %s", history)


def extract_metadata(
    pdf_path: Path,
    metadata_output: Path,
    references_output: Optional[Path] = None,
    tei_output: Optional[Path] = None,
    grobid_client: Optional[GrobidClient] = None,
    original_filename: Optional[str] = None,
    bib_entry: Optional[Dict] = None,
    grobid_input: Optional[Dict] = None,
):
    """Extract bibliographic metadata + references via Grobid.

    Runs ``/api/processFulltextDocument`` once, caches the raw TEI-XML at
    ``tei_output`` (if given), and parses it into ``metadata.json`` plus
    ``references.json``. If Grobid is unreachable or any step fails, falls
    back to placeholder files so downstream stages still run. Returns a Grobid
    outcome and actual input receipt, also saved in ``metadata["grobid"]``.
    ``extraction_method=bib`` alone cannot prove that references were extracted.

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
        Cache path for the raw Grobid TEI. Reused only when its provenance
        receipt matches the prepared PDF, consolidation settings and TEI bytes.
        Unverified/stale files are archived out of the active path. If None, defaults to
        ``<metadata_output.parent>/grobid.tei.xml``.
    grobid_client:
        A live :class:`GrobidClient`, or None when disabled/unavailable.
    grobid_input:
        Expected capability/version from the shared resume resolver. With an
        enabled expectation and no client, verified TEI may still be reused.
        Direct calls without this argument treat no client as deliberate
        disablement, not as a transient outage.
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
    _g = CONFIG.get("grobid", {})
    if grobid_input is None:
        grobid_input = {"status": "complete" if grobid_client is not None else "disabled",
                        "service_version": None}
    enabled = grobid_input["status"] != "disabled"
    outcome = "unavailable" if enabled else "disabled"
    parse_failed = False
    tei_inputs = {
        "pdf_sha256": _file_sha256(pdf_path),
        "consolidate_header": int(_g.get("consolidate_header", 1)),
        "consolidate_citations": int(_g.get("consolidate_citations", 0)),
        "pipeline_version": __version__,
        "service_version": grobid_input.get("service_version"),
        "producer_id": _g.get("producer_id"),
    }
    receipt_path = tei_output.with_suffix(tei_output.suffix + ".provenance.json")

    # 1. Verify the cached payload as well as its inputs. Merely rerunning the
    # stage cannot fix an obsolete TEI cache after re-OCR or consolidation edits.
    if enabled and tei_output.exists() and tei_output.stat().st_size > 0:
        try:
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            cached = tei_output.read_text(encoding="utf-8")
            if (receipt.get("inputs") == tei_inputs
                    and receipt.get("tei_sha256") == hashlib.sha256(cached.encode("utf-8")).hexdigest()):
                validate_fulltext_tei(cached)
                tei_xml = cached
                outcome = "cached"
                logger.info("Reusing verified Grobid TEI at %s", tei_output)
        except (OSError, ValueError, AttributeError, RuntimeError):
            pass
    if tei_output.exists() and tei_xml is None:
        # Post-processing reads the active TEI directly; leaving a stale file
        # here would resurrect its citations even after placeholder fallback.
        _archive_tei(tei_output, receipt_path, hash_dir)

    # 2. Otherwise, call Grobid if a client is available.
    if enabled and tei_xml is None and grobid_client is not None:
        try:
            logger.info("Calling Grobid on %s ...", pdf_path.name)
            tei_xml = grobid_client.process_fulltext(
                pdf_path,
                consolidate_header=int(_g.get("consolidate_header", 1)),
                consolidate_citations=int(_g.get("consolidate_citations", 0)),
            )
            validate_fulltext_tei(tei_xml)
            tei_output.write_text(tei_xml, encoding="utf-8")
            receipt = {"inputs": tei_inputs,
                       "tei_sha256": hashlib.sha256(tei_xml.encode("utf-8")).hexdigest()}
            temporary = receipt_path.with_name(receipt_path.name + ".tmp-" + uuid.uuid4().hex)
            try:
                temporary.write_text(json.dumps(receipt, indent=2), encoding="utf-8")
                temporary.replace(receipt_path)
            finally:
                temporary.unlink(missing_ok=True)
            logger.info("Wrote Grobid TEI to %s", tei_output)
            outcome = "extracted"
        except Exception as e:
            if external.STRICT_NETWORK and (isinstance(e, GrobidUnavailableError)
                                            or external.is_transient(e.__cause__ or e)):
                raise
            tei_xml = None
            outcome = "request_failed"
            _archive_tei(tei_output, receipt_path, hash_dir)
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
            parse_failed = True
            logger.warning("Failed to parse Grobid TEI references (%s)", e)
        # In-text citation graph (issue #7).  Independent of refs parse —
        # we want partial recovery when one fails and the other doesn't.
        try:
            intext = parse_tei_intext_citations(tei_xml)
        except Exception as e:
            parse_failed = True
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
            parse_failed = True
            logger.warning(
                "Failed to parse Grobid TEI header (%s); writing placeholder", e
            )
            header = None

    if parse_failed:
        outcome = "parse_failed"
        _archive_tei(tei_output, receipt_path, hash_dir)
    complete = tei_xml is not None and not parse_failed
    result = {"outcome": outcome,
              "producer_id": _g.get("producer_id"),
              "fingerprint": {"status": "complete" if complete else ("unavailable" if enabled else "disabled"),
                              "service_version": grobid_input.get("service_version")}}

    # 5. Write outputs. If we got a header, write it + the reference list
    # (which may be empty if Grobid was unavailable). Otherwise placeholder.
    if header is not None:
        header["grobid"] = result
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
        return result

    # 6. Placeholder path.
    _write_placeholder_metadata(
        pdf_path, metadata_output, references_output,
        original_filename=effective_filename,
        grobid_result=result,
    )
    return result
