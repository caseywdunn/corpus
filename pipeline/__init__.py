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
                _should_run_stage / _all_stage_artifacts_complete /
                _stage_recorded_complete (pipeline_state.json-backed),
                _run_quality_gates, _HugeDocumentError, _pdf_page_count
  io.py       — calculate_pdf_hash / short_hash / find_all_pdfs,
                audit_orphans, create_output_structure,
                create_summary_json, _verify_or_raise_collision

The root ``process_corpus.py`` re-exports these names so existing
``from process_corpus import …`` callers (tests, downstream tools)
keep working during the deprecation window.
"""

# PIPELINE_VERSION is the value recorded in pipeline_state.json on every
# successful stage. A mismatch on --resume forces the affected stage to
# re-run, so artifacts are never reused across releases. Sourced from
# pipeline/version.py (the single source of truth that also stamps bundle
# manifests, MCP bundle_info, etc.) — the release ritual bumps that
# constant, and per-paper artifacts are invalidated on the same cadence.
from .version import __version__ as PIPELINE_VERSION  # noqa: E402,F401


# Per-artifact schema version. Bumped only when an artifact's on-disk
# JSON shape changes in a backwards-incompatible way. Independent of
# PIPELINE_VERSION (which moves on every release whether or not the
# schema changes). Tools that read artifacts can refuse — or migrate —
# unknown values rather than silently mis-parsing old files.
ARTIFACT_SCHEMA_VERSION = 1


def stamp_artifact(payload: dict, *, schema_version: int = ARTIFACT_SCHEMA_VERSION) -> dict:
    """Prepend ``schema_version`` + ``pipeline_version`` to a per-paper
    artifact payload before writing. The version fields appear first in
    the resulting JSON, so a reader can dispatch on them before
    consuming the rest.

    Idempotent: a payload that already carries ``schema_version`` /
    ``pipeline_version`` from a prior write has those overwritten with
    the current values, so a stage that reads-modifies-writes an
    artifact bumps its stamp to whatever code is producing the new
    version. The caller's payload is not mutated.
    """
    out = {
        "schema_version": schema_version,
        "pipeline_version": PIPELINE_VERSION,
    }
    for k, v in payload.items():
        if k in ("schema_version", "pipeline_version"):
            continue
        out[k] = v
    return out
