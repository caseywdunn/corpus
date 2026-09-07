"""Per-paper taxon + lexicon annotation stage.

Walks the chunks of one paper and emits ``taxa.json`` (when a taxonomy
snapshot is configured) plus one ``<category>.json`` per category in
the configured ``lexicons`` mapping.

Each artifact carries an ``input_fingerprint`` (#29) so corpus_status
can detect stale annotations after the input file changes.
"""
from __future__ import annotations

import json
import hashlib
import logging
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional

from . import stamp_artifact
from .taxa import (RESERVED_ARTIFACT_STEMS, TaxonomyDB, extract_lexicon_mentions,
                   extract_taxon_mentions, validate_category)

logger = logging.getLogger(__name__)


def annotation_outputs_problem(hash_dir):
    """Verify the declared output set, including deliberate absence."""
    try:
        receipt = json.loads((hash_dir / "annotation_outputs.json").read_text())
        if receipt.get("version") != 1 or not isinstance(receipt.get("outputs"), dict):
            return "unverified annotation output receipt"
        for name, digest in receipt["outputs"].items():
            if Path(name).name != name or not name.endswith(".json"):
                return "invalid annotation output path"
            path = hash_dir / name
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
                return f"annotation output changed or missing: {name}"
    except (OSError, ValueError, TypeError, AttributeError):
        return "missing or invalid annotation output receipt"
    return None


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

    Missing configured inputs and unreadable chunks are failures, not empty
    extraction. Removed outputs are archived before publishing a new receipt.
    """
    out: List[Path] = []
    chunks_data = json.loads(chunks_file.read_text(encoding="utf-8"))
    chunks = chunks_data.get("chunks") or []
    lexicons = lexicons or {}
    for category in lexicons:
        validate_category(category)

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
    active = {p.name for p in out}
    for old in sorted(hash_dir.glob("*.json")):
        if old.name in active or old.name == "annotation_outputs.json":
            continue
        if old.stem in RESERVED_ARTIFACT_STEMS and old.name != "taxa.json":
            continue  # Never load large extraction dumps to discover annotations.
        is_annotation = old.name == "taxa.json"
        if not is_annotation:
            try:
                data = json.loads(old.read_text())
                is_annotation = isinstance(data, dict) and data.get("category") == old.stem and "mentions" in data
            except (OSError, ValueError):
                continue
        if is_annotation:
            history = hash_dir / "annotation_history"
            history.mkdir(exist_ok=True)
            old.rename(history / f"{uuid.uuid4().hex}-{old.name}")
    receipt = {"version": 1, "outputs": {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out}}
    temporary = hash_dir / f".annotation-outputs-{uuid.uuid4().hex}.tmp"
    temporary.write_text(json.dumps(receipt, sort_keys=True))
    temporary.replace(hash_dir / "annotation_outputs.json")
    return out
