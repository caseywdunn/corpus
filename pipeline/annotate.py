"""Per-paper taxon + lexicon annotation stage.

Walks the chunks of one paper and emits ``taxa.json`` (when a taxonomy
snapshot is configured) plus one ``<category>.json`` per category in
the configured ``lexicons`` mapping.

Each artifact carries an ``input_fingerprint`` (#29) so corpus_status
can detect stale annotations after the input file changes.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from . import stamp_artifact
from .taxa import TaxonomyDB, extract_lexicon_mentions, extract_taxon_mentions

logger = logging.getLogger(__name__)


def _extract_taxa_and_lexicons(
    chunks_file: Path,
    hash_dir: Path,
    taxonomy_db: Optional[TaxonomyDB],
    lexicons: Optional[Dict[str, Dict[str, Dict]]] = None,
    *,
    taxonomy_fingerprint: Optional[Dict[str, Any]] = None,
    lexicon_fingerprints: Optional[Dict[str, Dict[str, Any]]] = None,
) -> List[Path]:
    """Run taxon + multi-category lexicon extraction on the per-PDF chunks.

    Always-emitted artifacts:
        ``taxa.json`` — when a taxonomy DB is configured
        ``<category>.json`` — one per category in ``lexicons``

    ``lexicons`` is a two-level dict: ``{category: {term: {synonyms,
    translations, description}}}`` (the shape returned by
    :func:`taxa.load_lexicon`). ``lexicon_fingerprints`` carries the
    matching ``{category: {sha256, size, path}}`` so corpus_status can
    detect stale annotations after a section changes.

    Designed to degrade gracefully: missing inputs simply skip their
    output artifact rather than raising.
    """
    out: List[Path] = []
    if not chunks_file.exists():
        return out
    try:
        chunks_data = json.load(chunks_file.open(encoding="utf-8"))
    except Exception as e:
        logger.warning("Skipping taxa/lexicon pass: couldn't read %s: %s", chunks_file, e)
        return out
    chunks = chunks_data.get("chunks") or []

    if taxonomy_db is not None:
        taxa_res = extract_taxon_mentions(chunks, taxonomy_db)
        if taxonomy_fingerprint is not None:
            taxa_res["input_fingerprint"] = taxonomy_fingerprint
        taxa_file = hash_dir / "taxa.json"
        with taxa_file.open("w", encoding="utf-8") as f:
            json.dump(stamp_artifact(taxa_res), f, indent=2, ensure_ascii=False)
        out.append(taxa_file)
        logger.info(
            "Taxon mentions: %d (unique taxa: %d)",
            taxa_res["total_mentions"], taxa_res["unique_taxa"],
        )
    else:
        logger.info("Skipping taxon extraction (no taxonomy DB configured)")

    lexicons = lexicons or {}
    if not lexicons:
        logger.info("Skipping lexicon extraction (no lexicon configured)")
    fingerprints = lexicon_fingerprints or {}
    for category, lex in lexicons.items():
        if not lex:
            continue
        res = extract_lexicon_mentions(chunks, lex, category=category)
        fp = fingerprints.get(category)
        if fp is not None:
            res["input_fingerprint"] = fp
        out_file = hash_dir / f"{category}.json"
        with out_file.open("w", encoding="utf-8") as f:
            json.dump(stamp_artifact(res), f, indent=2, ensure_ascii=False)
        out.append(out_file)
        logger.info(
            "%s mentions: %d (unique terms: %d)",
            category.capitalize(), res["total_mentions"], res["unique_terms"],
        )
    return out
