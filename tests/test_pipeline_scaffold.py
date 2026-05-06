"""Smoke tests for the pipeline/ package scaffold (#42 partial).

These tests pin the public surface of the new modules so that future
work which migrates process_corpus.py to import from pipeline doesn't
break the entry points by accident.

Behavior is already covered by the existing test_per_stage_resume,
test_quality_gates, and test_stage_recorder suites — those exercise
the same functions through the process_corpus.* re-exports. This file
just confirms the new package paths exist.
"""
from __future__ import annotations

import pytest


def test_pipeline_config_surface():
    from pipeline.config import (
        _DEFAULT_CONFIG,
        _deep_merge,
        CONFIG,
        classify_section,
        load_config,
    )
    assert isinstance(_DEFAULT_CONFIG, dict)
    assert isinstance(CONFIG, dict)
    assert callable(load_config)
    assert callable(_deep_merge)
    # Spot-check classify_section on a known pattern
    assert classify_section(["Materials and Methods"]) == "methods"
    assert classify_section(["Some random heading"]) is None
    assert classify_section([]) is None


def test_pipeline_log_surface():
    from pipeline.log import per_pdf_file_log, setup_root_logging
    assert callable(setup_root_logging)
    assert callable(per_pdf_file_log)


def test_pipeline_stages_surface():
    from pipeline.stages import (
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
    # Closed reason-code set
    assert "timeout" in _REASON_CODES
    assert "too_large" in _REASON_CODES
    assert "quality_gate" in _REASON_CODES

    # Hierarchy: _HugeDocumentError must classify to too_large
    code, _ = _classify_exception(_HugeDocumentError("test"))
    assert code == "too_large"


def test_pipeline_io_surface():
    from pipeline.io import (
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
    assert HASH_PREFIX_LEN == 12
    assert callable(short_hash)
    # 12-char prefix
    assert short_hash("a" * 64) == "a" * 12


def test_process_corpus_still_works():
    """The legacy module-level names that tests already depend on
    should keep working — process_corpus.py is unchanged in this PR.
    """
    import process_corpus  # noqa
    assert callable(process_corpus.load_config)
    assert callable(process_corpus.classify_section)
    assert callable(process_corpus._stage)
    assert callable(process_corpus.calculate_pdf_hash)
