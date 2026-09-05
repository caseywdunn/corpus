"""Resolved Stage 1 settings and their consumers (#174).

Keep values (not opaque hashes) in receipts so a rerun can name the setting
that changed. Scheduler placement, output paths and unrelated server settings
are not extraction inputs. See OVERVIEW.md, "Corpuscle update contract".
"""

from copy import deepcopy
import json
from pathlib import Path

from .config import _DEFAULT_CONFIG, _deep_merge, load_config


def config_fingerprints(config, *, panel_mode, vision_model=None):
    """Return direct and inherited configuration inputs by stage.

    ``panel_mode`` is the applied CLI/config mode, not a backend availability
    guess. An explicitly requested backend must not stamp success if absent.
    The CPU floor and separately scheduled vision phase have distinct owners.
    """
    cfg = _deep_merge(_DEFAULT_CONFIG, config)

    def select(section, keys):
        return {f"{section}.{key}": deepcopy(cfg.get(section, {}).get(key)) for key in keys}

    scan = select("ocr", (
        "ocr_languages_default", "gibberish_threshold", "visual_script_gibberish_min",
        "reocr_scanned_text_layers", "scan_page_fraction_min", "probe_language_by_ocr",
        "probe_language_min_confidence", "probe_max_languages", "probe_max_gibberish",
        "probe_dpi", "probe_sample_pages",
    ))
    prep = {**scan, **select("ocr", ("optimize_level", "tesseract_page_timeout", "jobs")),
            **select("stage_timeouts", ("ocr", "ocr_per_page"))}
    extract = {**prep, **select("figures", ("resolution_mode", "images_scale", "vector_dpi", "max_dpi")),
               "compute.accelerator": cfg.get("compute", {}).get("accelerator", "auto")}
    chunks = {**extract, **select("chunking", ("max_tokens",))}
    figures = {**extract, "figures.panel_detection": panel_mode}
    if panel_mode.startswith("vision-"):
        figures["figures.model"] = vision_model
    metadata = {**prep,
                "grobid.disable": bool(cfg.get("grobid", {}).get("disable", False)),
                "grobid.producer_id": cfg.get("grobid", {}).get("producer_id"),
                "stage_timeouts.grobid": cfg["stage_timeouts"]["grobid"],
                "grobid.consolidate_header": int(cfg.get("grobid", {}).get("consolidate_header", 1)),
                "grobid.consolidate_citations": int(cfg.get("grobid", {}).get("consolidate_citations", 0))}
    return {
        "scan_detection": scan,
        "pdf_preparation": prep,
        "docling_extraction": extract,
        "metadata_extraction": metadata,
        "text_chunking": chunks,
        "taxa_and_lexicon_extraction": chunks,
        "figure_materialization": figures,
        "figure_crossref": {**chunks, **figures},
        "huge_document_check": select("huge_document", ("max_pages",)),
        "quality_gates": select("quality_gates", _DEFAULT_CONFIG["quality_gates"]),
    }


def configuration_drift(output_dir: Path, config_path: Path):
    """Compare configured settings with receipts, without models or writes.

    This is deliberately not a source-content or full completion audit. CLI
    overrides (including a separately scheduled CPU floor) may legitimately
    differ from config.yaml. Report those differences, don't silently declare
    that the original run used the current file's settings.
    """
    from .config_schema import validate_config
    from .stages import (_load_pipeline_state, _stage_input_changes,
                         _stage_recorded_complete)

    if not config_path.is_file():
        raise FileNotFoundError(config_path)
    config = load_config(config_path)
    validate_config(config)
    figures = config.get("figures", {})
    expected = config_fingerprints(config, panel_mode=figures.get("panel_detection", "ocr"),
                                   vision_model=figures.get("model"))
    affected = {}
    checked = 0
    for hd in sorted((output_dir / "documents").iterdir()):
        if not hd.is_dir():
            continue
        checked += 1
        records = _load_pipeline_state(hd)["stages"]
        changes = {}
        for stage, settings in expected.items():
            # Preserve non-config inputs: checking BibTeX/source inventory and
            # annotation evidence is a separate audit, not assumed current.
            record = records.get(stage) or {}
            fp = {**(record.get("input_fingerprint") or {}), "config": settings}
            if not _stage_recorded_complete(hd, stage, expected_fingerprint=fp):
                changes[stage] = _stage_input_changes(hd, stage, fp)
        if changes:
            affected[hd.name] = changes
    return {"config_path": str(config_path.resolve()), "documents_checked": checked,
            "documents_with_differences": len(affected), "differences": affected,
            "scope": "Stage 1 configuration only; CLI/CPU-floor overrides may differ. "
                     "Does not audit source files, BibTeX, annotation inputs, or external model versions."}


def source_input_drift(output_dir: Path, config_path: Path):
    """Read current PDFs, bibliography and annotations without starting models.

    This complements configuration_drift. Service/model identities and CLI
    overrides cannot be reconstructed from config alone; they are not probed.
    """
    from bib import BibIndex, keeppages_for_pdf, ocrlang_for_pdf, ocrmode_for_pdf
    from .io import find_all_pdfs, get_relative_paths, short_hash
    from .stages import (_expected_fingerprints_for_run, _file_sha256,
                         _load_pipeline_state, _metadata_fingerprint_for_pdf,
                         _stage_input_changes, _stage_recorded_complete)
    from .taxa import lexicon_fingerprints, load_lexicon

    config = load_config(config_path)
    def resolved(value):
        return (config_path.parent / value).resolve() if value else None
    input_dir = resolved(config.get("input_pdfs"))
    if input_dir is None:
        return {"available": False, "scope": "No input_pdfs configured; source inventory not checked."}
    bib_path = resolved(config.get("bib"))
    bib_index = BibIndex.from_path(bib_path) if bib_path else None
    lexicon_path = resolved(config.get("lexicon"))
    if lexicon_path:
        load_lexicon(lexicon_path)  # Validate before treating it as current input.
    lexicons = lexicon_fingerprints(lexicon_path) if lexicon_path else {}
    taxonomy_path = output_dir / "taxonomy.sqlite"
    taxonomy = ({"path": str(taxonomy_path), "sha256": _file_sha256(taxonomy_path),
                 "size": taxonomy_path.stat().st_size} if taxonomy_path.exists() else None)
    taxonomy_source = {"configured": False}
    tx = config.get("taxonomy") or {}
    if tx.get("source"):
        from .taxonomy_ingest import snapshot_matches, source_fingerprint
        fp = source_fingerprint(tx["source"], tx.get("root_id"), resolved(tx.get("path")))
        taxonomy_source = {"configured": True, "current": snapshot_matches(taxonomy_path, fp),
                           "expected_fingerprint": fp,
                           "scope": "WoRMS snapshots refresh only on explicit rebuild; no API probe."}
    inventory = find_all_pdfs(input_dir, exclude_under=output_dir, strict=True)
    current = {short_hash(full): (full, paths) for full, paths in inventory.items()}
    if len(current) != len(inventory):
        raise ValueError("Source inventory has a PDF hash-prefix collision")
    docs = {p.name: p for p in (output_dir / "documents").iterdir() if p.is_dir()}
    differences = {}
    consumed = {"bib_entry_sha256", "filename", "ocrlang", "ocrmode", "keeppages", "taxonomy", "lexicons"}
    for sha in sorted(current.keys() & docs.keys()):
        full, paths = current[sha]
        hd = docs[sha]
        name = paths[0].name
        expected = _expected_fingerprints_for_run(
            config_fingerprints=None,
            metadata_fingerprint=_metadata_fingerprint_for_pdf(bib_index, name),
            taxonomy_fingerprint=taxonomy, lexicon_fingerprints=lexicons,
            ocrlang=ocrlang_for_pdf(bib_index, name),
            ocrmode=ocrmode_for_pdf(bib_index, name),
            keeppages=keeppages_for_pdf(bib_index, name))
        records = _load_pipeline_state(hd)["stages"]
        changes = {}
        for stage in expected.keys() | {"taxa_and_lexicon_extraction", "figure_materialization", "figure_crossref"}:
            old = (records.get(stage) or {}).get("input_fingerprint") or {}
            fp = {k: v for k, v in old.items() if k not in consumed}
            fp.update(expected.get(stage, {}))
            if stage in ("figure_materialization", "figure_crossref"):
                fp.update({k: v for k, v in expected["docling_extraction"].items()
                           if k in {"ocrlang", "ocrmode", "keeppages"}})
            if not _stage_recorded_complete(hd, stage, expected_fingerprint=fp):
                changes[stage] = _stage_input_changes(hd, stage, fp)
        summary_path = hd / "summary.json"
        if not summary_path.is_file():
            changes["source_inventory"] = ["missing summary.json"]
        else:
            summary = json.loads(summary_path.read_text())
            source_changes = []
            if summary.get("pdf_hash_full") not in (None, full):
                source_changes.append("PDF hash-prefix collision")
            if summary.get("relative_paths") != get_relative_paths(paths, input_dir):
                source_changes.append("source paths")
            if summary.get("total_copies_found") != len(paths):
                source_changes.append("source copy count")
            if source_changes:
                changes["source_inventory"] = source_changes
        if changes:
            differences[sha] = changes
    return {"available": True, "input_dir": str(input_dir),
            "current_documents": len(current), "documents_checked": len(current.keys() & docs.keys()),
            "added": sorted(current.keys() - docs.keys()), "removed": sorted(docs.keys() - current.keys()),
            "taxonomy_source": taxonomy_source,
            "documents_with_differences": len(differences), "differences": differences,
            "scope": "Current PDF inventory, resolved BibTeX, OCR/page directives, lexicon and built taxonomy snapshot. "
                     "Configured taxonomy source receipts are checked separately. No external service/model probes; CLI overrides may differ."}
