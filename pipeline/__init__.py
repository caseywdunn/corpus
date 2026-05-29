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


# Honor a CORPUS_LOG_LEVEL env var if set. logging.basicConfig is a
# no-op after the first call, so configuring at package-import time —
# before any submodule's own basicConfig fires — lets `corpus -q run`
# quiet the entire subprocess tree (orchestrator → pipeline.main /
# embed / status / ...) without each module having to plumb -q through
# its argparse.  When the env var is absent, this block is a no-op and
# each module's own basicConfig wins as before.
import logging as _logging  # noqa: E402
import os as _os  # noqa: E402
import sys as _sys  # noqa: E402


# Platform-portability env-var defaults (#72, KMP duplicate-libomp).
# These must be set *before* any submodule (or anything they pull in)
# imports torch — hence the package __init__. setdefault leaves user
# overrides intact.
#
# TORCH_COMPILE_DISABLE=1 (#72): docling's StandardPdfPipeline triggers
# torch._inductor on certain layout / table code paths, which JIT-shells
# out to `g++` at runtime. On HPC compute nodes without a GCC module
# loaded, the JIT call fails with ``OSError: [Errno 14] Bad address:
# 'g++'`` and the affected papers are silently lost. docling doesn't
# need inductor for correctness, so disabling it is a clean fix that
# also avoids the rare "g++ not on PATH" failure on stripped-down
# laptops + deploy hosts. The (small) compile-time perf gain only
# matters for hot tensor ops we don't hit in the pipeline.
_os.environ.setdefault("TORCH_COMPILE_DISABLE", "1")

# KMP_DUPLICATE_LIB_OK=TRUE on macOS: pip-installed torch ships its own
# libomp.dylib, scikit-learn ships another, and conda-forge ships a third
# via llvm-openmp. Without this, the first ``import torch`` after numpy
# aborts with the duplicate-libomp message. The override is a documented
# OpenMP debug knob — safe in practice for our workloads, and matches the
# behavior tools/run_mcp_server.sh has set for the serve path since v0.2.
if _sys.platform == "darwin":
    _os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

# HF_HUB_DISABLE_IMPLICIT_TOKEN=1 (#97): the three model-load sites
# (SentenceTransformer in embeddings.py, AutoProcessor +
# Qwen2_5_VL in vision.py) pull in huggingface_hub, which logs a noisy
# warning when it falls back to an implicitly-stored token to fetch
# public models we never gate. All our models are public, so disable the
# implicit-token path. Must be set before huggingface_hub is imported —
# the lazy `from sentence_transformers import …` / `from transformers
# import …` inside those functions happens well after this package
# __init__, so setting it here covers every code path (pipeline build +
# the MCP serve path, which imports pipeline.embeddings). setdefault
# leaves an explicit user override intact.
_os.environ.setdefault("HF_HUB_DISABLE_IMPLICIT_TOKEN", "1")

_log_level = _os.environ.get("CORPUS_LOG_LEVEL", "").upper()
if _log_level in {"WARNING", "INFO", "DEBUG"}:
    _logging.basicConfig(
        level=getattr(_logging, _log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    # Stamp the marker that pipeline.log.setup_root_logging looks for
    # so it doesn't add a second StreamHandler on top of ours.
    for _h in _logging.getLogger().handlers:
        if isinstance(_h, _logging.StreamHandler) and not hasattr(_h, "baseFilename"):
            _h._corpus_stream = True


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
