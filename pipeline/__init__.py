"""Pipeline package — extracts foundational concerns from process_corpus.py (#42).

This is a partial split: cross-cutting infrastructure (config, logging,
stage runner / failures / quality gates / resume helpers) and per-paper
I/O (hashing, orphan audit, summary writing) live here. The actual stage
implementations (scan detection, OCR, docling, metadata, chunking,
figure passes, taxa/anatomy) remain in ``process_corpus.py`` as a
follow-up.

Submodules:

  config.py   — _DEFAULT_CONFIG, load_config, CONFIG, classify_section
  log.py      — setup_root_logging, per_pdf_file_log
  stages.py   — _stage, _classify_exception, fingerprint helpers,
                _should_run_stage, _all_stage_artifacts_complete,
                _run_quality_gates, _HugeDocumentError, _pdf_page_count
  io.py       — calculate_pdf_hash / short_hash / find_all_pdfs,
                audit_orphans, create_output_structure,
                create_summary_json, _verify_or_raise_collision

The root ``process_corpus.py`` re-exports these names so existing
``from process_corpus import …`` callers (tests, downstream tools)
keep working during the deprecation window.
"""
