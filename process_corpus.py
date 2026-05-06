#!/usr/bin/env python3
"""CLI shim for the corpus pipeline (#42 / #45 — split into pipeline/).

The implementation lives in the ``pipeline/`` package now. This file
stays at the repo root so existing ``python process_corpus.py …``
invocations from MCP clients, SLURM batch scripts, and the
``update_corpus.py`` orchestrator keep working unchanged.

Backwards-compat re-exports: every public + private name that used
to live at module level is re-exported here so ``from process_corpus
import _stage`` / ``import _extract_taxa_and_anatomy`` /
``import calculate_pdf_hash`` (etc.) keeps working in test code that
predates the split. The deprecation window for these re-exports
extends until the test suite has migrated to ``from pipeline import …``.
"""
from __future__ import annotations

import sys

# Re-exports — public surface used by tests + downstream tools
from pipeline import config as _pipeline_config  # noqa: F401
from pipeline.annotate import _extract_taxa_and_anatomy  # noqa: F401
from pipeline.chunking import chunk_text, ingest_to_vector_db  # noqa: F401
from pipeline.config import (  # noqa: F401
    _DEFAULT_CONFIG,
    _deep_merge,
    CONFIG,
    classify_section,
    load_config,
)
from pipeline.extract import extract_docling_content  # noqa: F401
from pipeline.figure_passes import (  # noqa: F401
    _crossref_chunks_and_figures,
    _pass25_annotate_figures,
    _pass3a_annotate_rois,
    _pass3b_annotate_rois,
)
from pipeline.io import (  # noqa: F401
    HASH_PREFIX_LEN,
    _verify_or_raise_collision,
    audit_orphans,
    calculate_pdf_hash,
    create_output_structure,
    create_summary_json,
    find_all_pdfs,
    get_relative_paths,
    short_hash,
)
from pipeline.log import per_pdf_file_log, setup_root_logging  # noqa: F401
from pipeline.main import main
from pipeline.metadata import (  # noqa: F401
    _write_placeholder_metadata,
    _year_from_filename,
    extract_metadata,
)
from pipeline.runner import run_pdf_processing_pipeline  # noqa: F401
from pipeline.scan import (  # noqa: F401
    _annotate_pack_availability,
    _available_tesseract_langs,
    _compose_ocr_langs,
    _detect_language,
    _gibberish_score,
    _resolve_tesseract_packs,
    _text_layer_scripts,
    _visual_page_script,
    create_cell_visualizations,
    detect_scan_type,
    prepare_pdf,
)
from pipeline.stages import (  # noqa: F401
    _all_stage_artifacts_complete,
    _classify_exception,
    _file_sha256,
    _HugeDocumentError,
    _load_pipeline_state,
    _pdf_page_count,
    _REASON_CODES,
    _record_stage_completion,
    _run_quality_gates,
    _safe_load_json,
    _save_pipeline_state,
    _should_run_stage,
    _stage,
    _stage_recorded_complete,
    _utcnow_iso,
)


if __name__ == "__main__":
    sys.exit(main())
