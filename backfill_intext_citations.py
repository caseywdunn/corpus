#!/usr/bin/env python3
"""Backfill intext_citations.json for papers processed before issue #7.

The pipeline now emits an ``intext_citations.json`` next to each paper's
``references.json`` — but papers processed by an older pipeline don't
have one.  This tool walks ``<output_dir>/documents/*/grobid.tei.xml``
and parses each into ``intext_citations.json``.  Pure post-processing of
the cached TEI; no Grobid call, no re-OCR.

Usage:
    python backfill_intext_citations.py /path/to/output_dir

By default skips papers whose ``intext_citations.json`` already exists.
``--force`` re-parses every paper.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from grobid_client import parse_tei_intext_citations

logger = logging.getLogger("backfill_intext")


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
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    docs = args.output_dir / "documents"
    if not docs.is_dir():
        logger.error("Not a corpus output dir: %s (no documents/)", args.output_dir)
        return 1

    n_total = n_parsed = n_skipped = n_failed = 0
    n_citations = 0
    for hash_dir in sorted(docs.iterdir()):
        if not hash_dir.is_dir():
            continue
        n_total += 1
        tei_path = hash_dir / "grobid.tei.xml"
        out_path = hash_dir / "intext_citations.json"

        if not tei_path.exists():
            continue  # placeholder / Grobid-failed papers — leave alone
        if out_path.exists() and not args.force:
            n_skipped += 1
            continue

        try:
            tei_xml = tei_path.read_text(encoding="utf-8")
            data = parse_tei_intext_citations(tei_xml)
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

    logger.info(
        "Done: %d papers (%d parsed, %d skipped, %d failed); %d citations total",
        n_total, n_parsed, n_skipped, n_failed, n_citations,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
