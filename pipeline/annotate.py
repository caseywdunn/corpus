"""Per-paper taxon + lexicon annotation stage.

Walks the chunks of one paper and emits ``taxa.json`` (when a
taxonomy snapshot is configured) plus one ``<category>.json`` per
configured lexicon (anatomy by default, plus any --lexicon
CATEGORY:PATH overrides from #24).

Each artifact carries an ``input_fingerprint`` (#29) so corpus_status
can detect stale annotations after the input file changes.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from taxa import TaxonomyDB, extract_lexicon_mentions, extract_taxon_mentions

logger = logging.getLogger(__name__)


def _extract_taxa_and_anatomy(
    chunks_file: Path,
    hash_dir: Path,
    taxonomy_db: Optional[TaxonomyDB],
    anatomy_lexicon: Optional[Dict[str, Dict]],
    *,
    taxonomy_fingerprint: Optional[Dict[str, Any]] = None,
    anatomy_fingerprint: Optional[Dict[str, Any]] = None,
    extra_lexicons: Optional[Dict[str, Dict[str, Dict]]] = None,
    extra_fingerprints: Optional[Dict[str, Dict[str, Any]]] = None,
) -> List[Path]:
    """Run taxon + multi-category lexicon extraction on the per-PDF chunks.

    Always-emitted artifacts:
        ``taxa.json`` — when a taxonomy DB is configured
        ``anatomy.json`` — when ``anatomy_lexicon`` is supplied
        ``<category>.json`` — one per ``extra_lexicons`` category (#24)

    Designed to degrade gracefully: missing inputs simply skip their
    output artifact rather than raising.

    Fingerprints (#29) get stamped into each output so corpus_status
    can detect stale annotations after an input changes.
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
            json.dump(taxa_res, f, indent=2, ensure_ascii=False)
        out.append(taxa_file)
        logger.info(
            "Taxon mentions: %d (unique taxa: %d)",
            taxa_res["total_mentions"], taxa_res["unique_taxa"],
        )
    else:
        logger.info("Skipping taxon extraction (no taxonomy DB configured)")

    # Build the unified category → lexicon map. Anatomy is the v0.1
    # default category; extra_lexicons (#24) layers in user-defined
    # ones (biogeography, life_history, …).
    lexicons_to_run: Dict[str, Dict[str, Dict]] = {}
    fingerprints_to_run: Dict[str, Optional[Dict[str, Any]]] = {}
    if anatomy_lexicon:
        lexicons_to_run["anatomy"] = anatomy_lexicon
        fingerprints_to_run["anatomy"] = anatomy_fingerprint
    if extra_lexicons:
        for cat, lex in extra_lexicons.items():
            if cat == "anatomy":
                # Already in the map under the canonical name; don't
                # double-process even if --lexicon anatomy:... is
                # supplied alongside --anatomy-lexicon.
                continue
            if lex:
                lexicons_to_run[cat] = lex
                fingerprints_to_run[cat] = (extra_fingerprints or {}).get(cat)

    if not lexicons_to_run:
        logger.info("Skipping lexicon extraction (no lexicon configured)")
    for category, lex in lexicons_to_run.items():
        res = extract_lexicon_mentions(chunks, lex, category=category)
        fp = fingerprints_to_run.get(category)
        if fp is not None:
            res["input_fingerprint"] = fp
        out_file = hash_dir / f"{category}.json"
        with out_file.open("w", encoding="utf-8") as f:
            json.dump(res, f, indent=2, ensure_ascii=False)
        out.append(out_file)
        logger.info(
            "%s mentions: %d (unique terms: %d)",
            category.capitalize(), res["total_mentions"], res["unique_terms"],
        )
    return out
