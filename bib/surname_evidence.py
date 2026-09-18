"""Apply source-verified surname evidence to derived reference views only."""

from __future__ import annotations

import json
from pathlib import Path


def load_source_reports(output_dir):
    """Read the versioned extraction report; caller fingerprints this mapping."""
    result = {}
    for path in sorted((Path(output_dir) / "documents").glob("*/text.json")):
        text = json.loads(path.read_text())
        report = (text.get("source_text_integrity") or {}).get("surnames")
        if report:
            result[path.parent.name] = report
    return result


def supported_reference(ref, report):
    """Require a verified image reading AND a complete curated title/year.

    The report never edits immutable reference observations. Missing titles,
    OCR-damaged titles and conflicting source candidates remain unresolved.
    """
    from .authority import extract_surname_from_ref_author, normalize_for_key

    if not report:
        return ref, []
    title = normalize_for_key(ref.get("title") or "")
    author_list = list(ref.get("authors") or [])
    reasons = []
    for position, author in enumerate(author_list):
        surname = extract_surname_from_ref_author(author)
        candidates = []
        leads = []
        for decision in report.get("decisions") or []:
            if (
                decision.get("year") != ref.get("year")
                or decision.get("original") != surname
            ):
                continue
            if decision.get("status") in (
                "source_supports_observed_spelling",
                "quoted_or_sic_context",
            ):
                continue
            leads.append(decision)
            if (
                decision.get("status") != "verified"
                or len(decision.get("candidates") or []) != 1
            ):
                continue
            curated = decision["candidates"][0]
            sources = [
                source
                for source in curated.get("sources", [])
                if len(title) >= 25 and normalize_for_key(source["title"]) == title
            ]
            if sources:
                candidates.append((decision, sources))
        replacements = {d["replacement"] for d, _ in candidates}
        if len(replacements) != 1:
            if leads:
                reasons.append(
                    {
                        "code": "possible_ocr_surname_conflict",
                        "author_position": position,
                        "observed_author": author,
                        "requires_source_review": True,
                        "candidate_surnames": sorted(
                            {
                                c["surname"]
                                for d in leads
                                for c in d.get("candidates", [])
                            }
                        )[:5],
                        "source_statuses": sorted({d["status"] for d in leads}),
                        "basis": "source_or_complete_reference_identity_not_uniquely_confirmed",
                    }
                )
            continue
        decision, sources = sorted(
            candidates, key=lambda c: (c[0].get("page", 0), c[0].get("item_ref", ""))
        )[0]
        start = author.rfind(surname)
        if start < 0:
            continue
        replacement = decision["replacement"]
        author_list[position] = (
            author[:start] + replacement + author[start + len(surname) :]
        )
        reasons.append(
            {
                "code": "source_ocr_surname_supported",
                "author_position": position,
                "observed_author": author,
                "derived_author": author_list[position],
                "year": ref["year"],
                "curated_bib_keys": [s["bib_key"] for s in sources],
                "basis": "source_image_ocr_consensus_and_exact_curated_title_year",
                "source_pdf_sha256": report.get("source_pdf_sha256"),
                "source_page": decision.get("page"),
                "item_ref": decision.get("item_ref"),
                "crop_sha256": decision.get("crop_sha256"),
                "catalog_sha256": report.get("catalog_sha256"),
                "producer_policy": report.get("policy"),
            }
        )
    return (
        (dict(ref, authors=author_list), reasons)
        if author_list != ref.get("authors")
        else (ref, reasons)
    )
