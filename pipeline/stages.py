"""Cross-cutting pipeline infrastructure: stage runner, failures,
quality gates, fingerprints, resume helpers.

* :func:`_stage` — the context manager every pipeline stage wraps.
  Records per-stage timing into ``processing_summary['stage_timings']``
  and on exception appends a structured failure into
  ``processing_summary['stage_failures']`` with a reason code from
  :func:`_classify_exception`. On successful exit, persists a
  completion record (pipeline_version + input_fingerprint) to
  ``pipeline_state.json`` so resume decisions are based on recorded
  state, not on artifact file existence. Re-raises so the
  pipeline-level try/except still catches.
* :func:`_should_run_stage` / :func:`_all_stage_artifacts_complete`
  / :func:`_stage_recorded_complete` — per-stage resume backed by
  ``pipeline_state.json``.
* :func:`_run_quality_gates` — silent-failure detectors (#36).
* :func:`_file_sha256` / :func:`_safe_load_json` — small utilities
  used by other stages and by :func:`_run_quality_gates`.
* :class:`_HugeDocumentError` + :func:`_pdf_page_count` — huge-document
  gate (#35).
"""
from __future__ import annotations

import hashlib
import json
import logging
import subprocess
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from .config import CONFIG

logger = logging.getLogger(__name__)


# Reason codes recorded in stage_failures[]. The set is closed: any
# unmapped exception falls through to "crash". When extending, also
# update tools that group/aggregate failures (corpus_status.py, #40).
_REASON_CODES = (
    "timeout",                # wallclock cap exceeded
    "crash",                  # unhandled exception, subprocess died
    "external_unavailable",   # Grobid / BHL / CrossRef / OpenAlex non-200
    "unsupported_format",     # encrypted, password-protected, etc.
    "corrupted",              # PDF parse error
    "too_large",              # exceeds huge_document.max_pages (#35)
    "quality_gate",           # failed a downstream sanity check (#36)
)


class _HugeDocumentError(Exception):
    """Raised by the huge-document gate (#35) when a PDF's page count
    exceeds ``CONFIG['huge_document']['max_pages']``.

    Caught by the pipeline's outer try/except like any other stage
    error; ``_classify_exception`` maps it to ``reason_code=too_large``
    so corpus_status.py can group these separately from real crashes.
    """


def _pdf_page_count(pdf_path: Path) -> Optional[int]:
    """Return the page count of ``pdf_path`` via PyMuPDF metadata.

    Cheap — opens the document but does not render pages. Returns
    None if the file can't be read (a corrupted PDF will fail later
    with the appropriate reason_code).

    PyMuPDF is a hard dependency (requirements.txt) — its absence
    raises ImportError so the huge-document gate fails as a stage
    error rather than silently bypassing.
    """
    import fitz  # type: ignore
    try:
        with fitz.open(str(pdf_path)) as doc:
            return int(doc.page_count)
    except Exception as e:
        logger.warning("Could not read page count from %s: %s", pdf_path, e)
        return None


def _utcnow_iso() -> str:
    """UTC ISO-8601 timestamp with second precision (no microseconds)."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _classify_exception(e: BaseException) -> Tuple[str, str]:
    """Map an exception to ``(reason_code, short_detail)`` for stage_failures[].

    Conservative — anything not specifically recognized falls through
    to ``crash`` rather than being silently mis-categorized.
    """
    name = type(e).__name__
    if isinstance(e, _HugeDocumentError):
        return "too_large", str(e)
    if isinstance(e, subprocess.TimeoutExpired):
        return "timeout", f"{name}: command timed out after {e.timeout}s"
    # Grobid availability errors are imported lazily — avoid circular import
    # when grobid_client doesn't exist (e.g., in unit tests).
    try:
        from grobid_client import GrobidUnavailableError  # type: ignore
        if isinstance(e, GrobidUnavailableError):
            return "external_unavailable", f"{name}: {e}"
    except ImportError:
        pass
    try:
        import requests
        if isinstance(e, (requests.ConnectionError, requests.Timeout, requests.HTTPError)):
            return "external_unavailable", f"{name}: {e}"
    except ImportError:
        pass
    msg = str(e).lower()
    if "encrypted" in msg or "password" in msg:
        return "unsupported_format", f"{name}: {e}"
    if "corrupt" in msg or "invalid pdf" in msg or "syntax error" in msg:
        return "corrupted", f"{name}: {e}"
    return "crash", f"{name}: {e}"


def _file_sha256(path: Path) -> str:
    """Streaming SHA-256 of ``path``. Used by #29 to fingerprint inputs
    (taxonomy.sqlite, lexicon.yaml) so per-paper annotation
    artifacts can record the exact input version they were built from.

    Streamed in 64 KB chunks — taxonomy.sqlite at corpus scale is tens
    of MB, well within scope for in-process hashing during startup but
    not something to load into memory whole.
    """
    h = hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda: f.read(64 * 1024), b""):
            h.update(buf)
    return h.hexdigest()


def _safe_load_json(path: Path) -> Any:
    """Load JSON; return {} when missing or malformed.

    Used by quality gates so a missing/broken artifact short-circuits
    individual checks rather than crashing the gate runner.
    """
    try:
        if path.exists() and path.stat().st_size > 0:
            with path.open(encoding="utf-8") as f:
                return json.load(f)
    except Exception:
        return {}
    return {}


# ---------------------------------------------------------------------------
# Per-stage resume (#28) — pipeline_state.json
# ---------------------------------------------------------------------------
#
# A stage is "complete" iff pipeline_state.json records a completion
# entry whose pipeline_version matches the current PIPELINE_VERSION
# and whose input_fingerprint matches what the caller now expects.
# File existence on disk is no longer the correctness signal: a paper
# that ran under an older pipeline version, or whose lexicons changed
# since it last ran, is forced to re-execute the affected stages.
#
# The state file is written incrementally by :func:`_stage` on every
# successful stage exit (atomic tmp+rename), so a crash between stages
# loses no completion records.

PIPELINE_STATE_FILE = "pipeline_state.json"


def _load_pipeline_state(hash_dir: Path) -> Dict[str, Any]:
    """Return the parsed pipeline_state.json or an empty skeleton.

    Returns a dict with ``stages: {<stage_name>: {...}}`` even when the
    file is missing or malformed, so callers don't have to guard.
    """
    path = hash_dir / PIPELINE_STATE_FILE
    state = _safe_load_json(path)
    if not isinstance(state, dict):
        state = {}
    state.setdefault("stages", {})
    return state


def _save_pipeline_state(hash_dir: Path, state: Dict[str, Any]) -> None:
    """Atomically persist pipeline_state.json (tmp + rename)."""
    path = hash_dir / PIPELINE_STATE_FILE
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def _record_stage_completion(
    hash_dir: Path,
    stage_name: str,
    *,
    input_fingerprint: Optional[Dict[str, Any]] = None,
) -> None:
    """Record a successful stage completion in pipeline_state.json.

    Captures pipeline_version + completed_at + input_fingerprint so
    resume can detect both code-version and input-version drift.
    """
    from . import PIPELINE_VERSION

    state = _load_pipeline_state(hash_dir)
    state["pipeline_version_latest"] = PIPELINE_VERSION
    state["stages"][stage_name] = {
        "completed_at": _utcnow_iso(),
        "pipeline_version": PIPELINE_VERSION,
        "input_fingerprint": input_fingerprint or {},
    }
    _save_pipeline_state(hash_dir, state)


def _stage_recorded_complete(
    hash_dir: Path,
    stage_name: str,
    *,
    expected_fingerprint: Optional[Dict[str, Any]] = None,
) -> bool:
    """Return True iff pipeline_state.json has a completion record for
    ``stage_name`` whose pipeline_version matches PIPELINE_VERSION and
    whose input_fingerprint matches ``expected_fingerprint`` (when given).
    """
    from . import PIPELINE_VERSION

    rec = _load_pipeline_state(hash_dir)["stages"].get(stage_name)
    if not isinstance(rec, dict):
        return False
    if rec.get("pipeline_version") != PIPELINE_VERSION:
        return False
    if expected_fingerprint is not None:
        if rec.get("input_fingerprint") != expected_fingerprint:
            return False
    return True


def _should_run_stage(
    stage_name: str,
    *,
    hash_dir: Path,
    resume: bool,
    processing_summary: Dict[str, Any],
    expected_fingerprint: Optional[Dict[str, Any]] = None,
) -> bool:
    """Return True if the stage should run.

    With ``resume=True``, skip only when pipeline_state.json records
    this stage as complete under the current PIPELINE_VERSION (and
    matching input_fingerprint, when supplied). On skip, record on
    ``processing_summary['skipped_stages']`` and log.
    """
    if not resume:
        return True
    if not _stage_recorded_complete(hash_dir, stage_name,
                                    expected_fingerprint=expected_fingerprint):
        return True
    logger.info("%s: skipping (recorded complete)", stage_name)
    processing_summary.setdefault("skipped_stages", []).append(stage_name)
    return False


# Stages that every paper produces, regardless of optional inputs.
_CORE_STAGES: Tuple[str, ...] = (
    "scan_detection",
    "pdf_preparation",
    "docling_extraction",
    "metadata_extraction",
    "text_chunking",
)


def _all_stage_artifacts_complete(
    hash_dir: Path,
    *,
    expected_stages: Optional[Iterable[str]] = None,
) -> bool:
    """Return True iff every required stage is recorded as complete in
    pipeline_state.json under the current PIPELINE_VERSION.

    Used by the outer --resume short-circuit. ``expected_stages``
    defaults to ``_CORE_STAGES``; callers extend it with optional
    stages whose presence depends on configured inputs (e.g.
    ``"taxa_anatomy_extraction"`` when a taxonomy DB or lexicon is
    configured).
    """
    stages = tuple(expected_stages) if expected_stages is not None else _CORE_STAGES
    return all(_stage_recorded_complete(hash_dir, s) for s in stages)


def _run_quality_gates(hash_dir: Path) -> List[Dict[str, Any]]:
    """Cheap silent-failure detectors against the produced artifacts (#36).

    Returns a list of ``{gate, severity, detail, metric}`` records — one
    per failed gate. Empty when everything looks clean. Read-only.

    Each gate is informational. Operators decide what to do with flagged
    papers via ``corpus_status.py`` (#40); nothing is rejected here.
    """
    # _gibberish_score lives in process_corpus.py until #45 extracts
    # the scan-detection stage into pipeline.scan. Lazy import so
    # this module's load order doesn't couple to scan.
    try:
        from .scan import _gibberish_score  # type: ignore
    except ImportError:
        from process_corpus import _gibberish_score  # type: ignore

    cfg = CONFIG.get("quality_gates", {})
    flags: List[Dict[str, Any]] = []

    text = _safe_load_json(hash_dir / "text.json")
    chunks = _safe_load_json(hash_dir / "chunks.json")
    figures = _safe_load_json(hash_dir / "figures.json")
    refs = _safe_load_json(hash_dir / "references.json")
    scan = _safe_load_json(hash_dir / "scan_detection.json")

    body = (text.get("text") or "") if isinstance(text, dict) else ""
    pages = int(text.get("pages") or 0) if isinstance(text, dict) else 0
    chunk_list = chunks.get("chunks") or [] if isinstance(chunks, dict) else []
    fig_list = figures.get("figures") or [] if isinstance(figures, dict) else []
    ref_count = int(refs.get("total_references") or 0) if isinstance(refs, dict) else 0
    needs_ocr = bool(scan.get("needs_ocr")) if isinstance(scan, dict) else False

    # empty_text — extracted text is implausibly short
    min_chars = int(cfg.get("empty_text_min_chars", 500))
    if len(body) < min_chars:
        flags.append({
            "gate": "empty_text",
            "severity": "error",
            "detail": f"text.json body has {len(body)} chars (min: {min_chars})",
            "metric": len(body),
        })

    # low_text_density — text/page ratio below threshold
    min_chars_per_page = int(cfg.get("min_chars_per_page", 200))
    if pages > 0:
        density = len(body) / pages
        if density < min_chars_per_page:
            flags.append({
                "gate": "low_text_density",
                "severity": "warn",
                "detail": f"{density:.0f} chars/page over {pages} pages (min: {min_chars_per_page})",
                "metric": round(density, 1),
            })

    # gibberish_after_ocr — OCR'd papers whose final extracted text
    # is still gibberish (silent-failure mode)
    max_gibberish = float(cfg.get("max_gibberish_score", 0.5))
    if needs_ocr and body:
        score = _gibberish_score(body[:50000])  # cap sample for speed
        if score > max_gibberish:
            flags.append({
                "gate": "gibberish_after_ocr",
                "severity": "error",
                "detail": f"gibberish_score={score:.2f} on extracted text (max: {max_gibberish})",
                "metric": round(score, 3),
            })

    # zero_references_unexpected — multi-page paper with no references
    min_pages_for_refs = int(cfg.get("zero_refs_min_pages", 5))
    if pages >= min_pages_for_refs and ref_count == 0:
        flags.append({
            "gate": "zero_references_unexpected",
            "severity": "warn",
            "detail": f"{pages} pages but references.json is empty",
            "metric": pages,
        })

    # single_token_chunks — extraction collapsed; median chunk too short
    min_median_chars = int(cfg.get("min_median_chunk_chars", 50))
    if chunk_list:
        lengths = sorted(len(c.get("text") or "") for c in chunk_list if isinstance(c, dict))
        if lengths:
            median = lengths[len(lengths) // 2]
            if median < min_median_chars:
                flags.append({
                    "gate": "single_token_chunks",
                    "severity": "warn",
                    "detail": f"median chunk length {median} chars over {len(chunk_list)} chunks (min: {min_median_chars})",
                    "metric": median,
                })

    # all_black_figures — extraction artifact where most figures are empty
    if fig_list:
        try:
            from PIL import Image, ImageStat
        except ImportError:
            Image = None  # type: ignore
        if Image is not None:
            mean_threshold = float(cfg.get("min_figure_mean_intensity", 10))
            max_sampled = int(cfg.get("max_figures_sampled", 50))
            n_black = 0
            n_checked = 0
            for fig in fig_list[:max_sampled]:
                if not isinstance(fig, dict):
                    continue
                fname = fig.get("filename")
                if not fname:
                    continue
                fpath = hash_dir / "figures" / fname
                if not fpath.exists():
                    continue
                try:
                    with Image.open(fpath) as im:
                        gray = im.convert("L")
                        mean = ImageStat.Stat(gray).mean[0]
                except Exception:
                    continue
                n_checked += 1
                if mean < mean_threshold:
                    n_black += 1
            if n_checked > 0 and n_black / n_checked > 0.5:
                flags.append({
                    "gate": "all_black_figures",
                    "severity": "warn",
                    "detail": f"{n_black}/{n_checked} figures have mean intensity < {mean_threshold}",
                    "metric": n_black,
                })

    return flags


@contextmanager
def _stage(
    processing_summary: Dict[str, Any],
    name: str,
    *,
    hash_dir: Optional[Path] = None,
    input_fingerprint: Optional[Dict[str, Any]] = None,
):
    """Record per-stage timing into ``stage_timings[]`` and, on exception,
    append a structured failure into ``stage_failures[]``. Re-raises
    so the pipeline-level try/except still catches.

    On successful exit, when ``hash_dir`` is provided, the stage's
    completion is persisted to ``pipeline_state.json`` (with
    ``PIPELINE_VERSION`` + ``input_fingerprint``) so the next ``--resume``
    can decide skip-vs-rerun from a structured signal rather than from
    file existence.
    """
    started_at = _utcnow_iso()
    t0 = time.monotonic()
    err: Optional[BaseException] = None
    try:
        yield
    except Exception as e:
        err = e
        raise
    finally:
        ended_at = _utcnow_iso()
        duration_s = round(time.monotonic() - t0, 3)
        processing_summary.setdefault("stage_timings", []).append({
            "stage": name,
            "started_at": started_at,
            "ended_at": ended_at,
            "duration_s": duration_s,
            "ok": err is None,
        })
        if err is not None:
            reason_code, detail = _classify_exception(err)
            processing_summary.setdefault("stage_failures", []).append({
                "stage": name,
                "reason_code": reason_code,
                "reason_detail": detail[:500],
                "started_at": started_at,
                "ended_at": ended_at,
                "duration_s": duration_s,
            })
        elif hash_dir is not None:
            _record_stage_completion(
                hash_dir, name, input_fingerprint=input_fingerprint,
            )
