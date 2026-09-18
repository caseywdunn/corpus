#!/usr/bin/env python3
"""Backfill intext_citations.json for papers processed before issue #7.

The pipeline now emits an ``intext_citations.json`` next to each paper's
``references.json`` — but papers processed by an older pipeline don't
have one.  This tool walks ``<output_dir>/documents/*/grobid.tei.xml``
and parses each into ``intext_citations.json``. Uses the prepared PDF for
source citation text only when its provenance receipt matches the TEI and
PDF bytes; otherwise preserves TEI observations with explicit source status.
No Grobid call, no re-OCR. Old TEI without ref coordinates requires the
metadata stage to run against Grobid again for source-backed text repair.

Usage:
    python -m pipeline.intext_citations /path/to/output_dir

Idempotency (#30): by default each hash dir is processed only when
``intext_citations.json`` is missing. Re-running on an unchanged corpus
is a near-zero-work no-op (one stat() per hash dir). Adding a new paper
picks it up automatically. ``--force`` re-parses every paper —
deterministic re-write of the same content given the same TEI.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from pathlib import Path

from pipeline.grobid_client import parse_tei_intext_citations

logger = logging.getLogger("corpus.intext")


def _matching_source_pdf(hash_dir: Path, tei_xml: str):
    """Only use coordinates against the prepared PDF that produced the TEI."""
    pdf = hash_dir / "processed.pdf"
    receipt = hash_dir / "grobid.tei.xml.provenance.json"
    try:
        proof = json.loads(receipt.read_text(encoding="utf-8"))
        if proof.get("tei_sha256") != hashlib.sha256(tei_xml.encode("utf-8")).hexdigest():
            return None
        digest = hashlib.sha256()
        with pdf.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return pdf if digest.hexdigest() == proof.get("inputs", {}).get("pdf_sha256") else None
    except (OSError, ValueError, AttributeError, TypeError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output_dir", type=Path,
        help="Corpus output directory (contains documents/<hash>/ subdirs)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-parse even if intext_citations.json already exists",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Report which papers would be parsed without writing anything.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    docs = args.output_dir / "documents"
    if not docs.is_dir():
        if args.dry_run:
            # See pipeline/embed.py: on a corpuscle that has not been
            # extracted yet, "0 papers" is the plan, not a failure.
            logger.info(
                "Dry-run: %s does not exist yet; the extract step creates "
                "it on a real run. 0 papers to parse.", docs,
            )
            return 0
        logger.error("Not a corpus output dir: %s (no documents/)", args.output_dir)
        return 1

    n_total = n_parsed = n_skipped = n_failed = n_no_tei = 0
    n_citations = 0
    for hash_dir in sorted(docs.iterdir()):
        if not hash_dir.is_dir():
            continue
        n_total += 1
        tei_path = hash_dir / "grobid.tei.xml"
        out_path = hash_dir / "intext_citations.json"

        if not tei_path.exists():
            n_no_tei += 1
            continue  # placeholder / Grobid-failed papers — leave alone
        if out_path.exists() and not args.force:
            n_skipped += 1
            continue

        if args.dry_run:
            n_parsed += 1
            logger.debug("would parse %s", hash_dir.name)
            continue

        try:
            tei_xml = tei_path.read_text(encoding="utf-8")
            data = parse_tei_intext_citations(
                tei_xml, pdf_path=_matching_source_pdf(hash_dir, tei_xml),
            )
        except Exception as e:
            logger.warning("%s: parse failed: %s", hash_dir.name, e)
            n_failed += 1
            continue

        out_path.write_text(json.dumps(data, indent=2, ensure_ascii=False))
        n_parsed += 1
        n_citations += len(data["citations"])
        if n_parsed % 200 == 0:
            logger.info(
                "Parsed %d/%d papers, %d citations so far",
                n_parsed, n_total, n_citations,
            )

    if args.dry_run:
        logger.info(
            "Dry-run: %d papers (would parse %d; %d already up-to-date; "
            "%d have no grobid.tei.xml). No files written.",
            n_total, n_parsed, n_skipped, n_no_tei,
        )
    else:
        # `backfill_intext` is an idempotent top-up. In-text citations are
        # normally written during the main extraction stage, so on a healthy
        # re-run every paper is "already up-to-date" and nothing is parsed
        # here — that is success, not a no-op failure. Spell it out so
        # "0 parsed, N skipped" doesn't read like the feature did nothing.
        logger.info(
            "Done: %d papers — %d newly parsed, %d already up-to-date, "
            "%d without grobid.tei.xml, %d parse-failed; %d citations from "
            "newly parsed papers (papers already up-to-date kept the in-text "
            "citations written during extraction).",
            n_total, n_parsed, n_skipped, n_no_tei, n_failed, n_citations,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
