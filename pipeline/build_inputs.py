"""Resolved Stage 1 settings and their consumers (#174).

Keep values (not opaque hashes) in receipts so a rerun can name the setting
that changed. Scheduler placement, output paths and unrelated server settings
are not extraction inputs. See OVERVIEW.md, "Corpuscle update contract".
"""

from copy import deepcopy
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
            if stage == "taxa_and_lexicon_extraction" and stage not in records:
                continue  # Optional annotations aren't part of this audit.
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
