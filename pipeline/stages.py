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
import os
import re
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
        from .grobid_client import GrobidUnavailableError  # type: ignore
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
    """Atomically persist pipeline_state.json (per-writer tmp + rename).

    Per-writer tmp filename (``.tmp.<pid>.<ns>``) is belt-and-braces
    against the shared-tmp corruption pattern from #55: if two array
    tasks ever again end up writing the same hash directory in
    parallel, their tmp files don't collide and the atomic rename
    leaves a clean payload from one writer or the other (last write
    wins; never an interleaved-bytes corruption or a missing-tmp
    rename failure).
    """
    from . import stamp_artifact

    path = hash_dir / PIPELINE_STATE_FILE
    tmp = path.with_suffix(f"{path.suffix}.tmp.{os.getpid()}.{time.monotonic_ns()}")
    tmp.write_text(
        json.dumps(stamp_artifact(state), indent=2, sort_keys=True),
        encoding="utf-8",
    )
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
        changes = _stage_input_changes(hash_dir, stage_name, expected_fingerprint)
        logger.info("%s: rerunning (%s)", stage_name, ", ".join(changes))
        processing_summary.setdefault("rerun_reasons", {})[stage_name] = changes
        return True
    logger.info("%s: skipping (recorded complete)", stage_name)
    processing_summary.setdefault("skipped_stages", []).append(stage_name)
    return False


def _stage_input_changes(hash_dir, stage_name, expected_fingerprint):
    """Explain receipt drift without loading a model or changing artifacts."""
    from . import PIPELINE_VERSION
    record = _load_pipeline_state(hash_dir)["stages"].get(stage_name)
    if not isinstance(record, dict):
        return ["missing completion record"]
    if record.get("pipeline_version") != PIPELINE_VERSION:
        return ["pipeline_version"]

    def flatten(value, prefix=""):
        out = {}
        for key, item in value.items():
            name = f"{prefix}.{key}" if prefix else key
            if isinstance(item, dict) and item:
                out.update(flatten(item, name))
            else:
                out[name] = item
        return out

    old = flatten(record.get("input_fingerprint") or {})
    new = flatten(expected_fingerprint or {})
    return sorted(key for key in old.keys() | new.keys()
                  if key not in old or key not in new or old[key] != new[key]) or ["completion record"]


# Every resumable stage that descends from `processed.pdf`, and so is
# invalidated by an OCR-language change (#176). That is all of them: the
# core five plus taxa/lexicon extraction, which reads the chunks. Composite
# figure receipts inherit these directives in the configuration loop below.
_OCR_DEPENDENT_STAGES: Tuple[str, ...] = (
    "scan_detection",
    "pdf_preparation",
    "docling_extraction",
    "metadata_extraction",
    "text_chunking",
    "taxa_and_lexicon_extraction",
)


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
    expected_fingerprints: Optional[Dict[str, Dict[str, Any]]] = None,
) -> bool:
    """Return True iff every required stage is recorded as complete in
    pipeline_state.json under the current PIPELINE_VERSION *and* (when
    applicable) under the matching ``input_fingerprint``.

    Used by the outer --resume short-circuit. ``expected_stages``
    defaults to ``_CORE_STAGES``; callers extend it with optional
    stages whose presence depends on configured inputs (e.g.
    ``"taxa_and_lexicon_extraction"`` when a taxonomy DB or lexicon is
    configured).

    ``expected_fingerprints`` (#56) is ``{stage_name: fingerprint}`` —
    forwarded to :func:`_stage_recorded_complete` per stage so the
    outer gate matches the inner per-stage gate's staleness logic.
    Without this, editing a lexicon never invalidated already-recorded
    stages: the outer gate dropped the doc before the inner gate could
    spot the fingerprint mismatch. Build the mapping with
    :func:`_expected_fingerprints_for_run` so all call sites stay in
    lockstep with what the runner records.
    """
    stages = tuple(expected_stages) if expected_stages is not None else _CORE_STAGES
    fps = expected_fingerprints or {}
    return all(
        _stage_recorded_complete(hash_dir, s, expected_fingerprint=fps.get(s))
        for s in stages
    )


def _metadata_fingerprint_for_pdf(bib_index, filename: str, *,
                                  grobid_context=None, hash_dir=None) -> Dict[str, Any]:
    """Fingerprint the resolved entry, not the whole library bibliography.

    Absence is a value: adding/removing an entry must invalidate metadata.
    Filename is also consumed directly for provenance and fallback title/year,
    even when a rename resolves to the same entry (or to no entry).
    """
    entry = bib_index.lookup(filename) if bib_index is not None else None
    canonical = json.dumps(entry, sort_keys=True, ensure_ascii=False,
                           separators=(",", ":"))
    result = {"bib_entry_sha256": hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
              "filename": filename}
    if grobid_context is not None:
        from .grobid_state import REFERENCE_EVIDENCE_VERSION, grobid_input
        result["grobid"] = grobid_input(grobid_context, hash_dir)
        result["reference_evidence_version"] = REFERENCE_EVIDENCE_VERSION
    return result


def _expected_fingerprints_for_run(
    *,
    config_fingerprints: Optional[Dict[str, Dict[str, Any]]],
    metadata_fingerprint: Optional[Dict[str, Any]],
    ocrlang: Optional[str],
    ocrmode: Optional[str],
    keeppages: Optional[str],
    taxonomy_fingerprint: Optional[Dict[str, Any]] = None,
    lexicon_fingerprints: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Dict[str, Dict[str, Any]]:
    """Build the ``{stage_name: input_fingerprint}`` map both resume
    gates need (#56).

    Metadata consumes a per-paper resolved BibTeX entry and filename;
    annotation consumes taxonomy and lexicon fingerprints. OCR directives
    invalidate the PDF and every descendant. Inputs are required at the call
    boundary so a newly added input cannot silently miss the outer fast path.
    ``metadata_fingerprint=None`` is only for callers without a document
    context (e.g. tests isolating another stage); production resolves it with
    :func:`_metadata_fingerprint_for_pdf`, including absent bib entries.

    Both ``main.py``'s outer fast-path and ``runner.py``'s per-stage
    gate must call this so they agree on what "stale" means; otherwise
    a doc whose stage was recorded under a now-edited lexicon will be
    silently skipped by whichever side forgets the fingerprint.

    ``ocrlang`` (#176) and ``ocrmode`` (#186) are per-*document*, unlike the
    taxonomy and lexicon fingerprints which are per-run — so callers must
    resolve both inside their per-document loop rather than hoisting them.

    It fingerprints *every* resumable stage, not just the two that read
    it. Changing the OCR language or execution mode rewrites
    ``processed.pdf``, and every later stage descends from those bytes —
    docling reads the PDF, Grobid reads it, chunks come from docling's text,
    taxa come from the chunks.
    Fingerprinting only ``scan_detection`` and ``pdf_preparation`` re-OCRs
    the paper and then skips everything that consumes the result, so
    ``text.json`` still holds the old OCR while the log cheerfully reports
    the new ``-l``. That is worse than not re-running at all: it looks
    like it worked.

    Neither OCR directive has a default *on purpose*. The first cut of #176
    gave ``ocrlang`` one, and
    main.py's outer fast path — the gate that actually skips work, and
    which runs before the per-stage gate — silently kept the default:
    a paper with no tag last run recorded ``{}``, the gate expected ``{}``,
    they matched, and the paper was skipped whole before anything could
    see the tag the operator had just added. Requiring the argument turns
    that omission into a TypeError at the call site. Pass ``None``
    explicitly where a per-document value genuinely doesn't apply.

    ``keeppages`` (#188) follows the same contract and for a stronger
    reason: it does not merely change how the PDF is read, it changes which
    pages the PDF *is*. A page-range edit invalidates strictly more than an
    ocrlang edit does, so nothing derived from the old selection may survive
    a resume. It is also required rather than defaulted, for the same
    fast-path reason — a default here reproduces the shipped #176 bug where
    an added directive was skipped before anything could see it.

    Note the empty dict, not None, when no override is set. ``{}`` is what
    :func:`_record_stage_completion` writes for a stage with no
    fingerprint, so an untagged document still compares equal to its
    existing record and is *not* re-OCR'd — while a tag that was added,
    changed, or removed all compare unequal and re-run. Returning None
    here instead would make removal undetectable, since
    :func:`_stage_recorded_complete` skips the comparison entirely on
    None.
    """
    fps: Dict[str, Dict[str, Any]] = {}
    if metadata_fingerprint is not None:
        fps["metadata_extraction"] = dict(metadata_fingerprint)
    taxa_lex_fp: Dict[str, Any] = {}
    if taxonomy_fingerprint is not None:
        taxa_lex_fp["taxonomy"] = taxonomy_fingerprint
    if lexicon_fingerprints:
        taxa_lex_fp["lexicons"] = lexicon_fingerprints
    if taxa_lex_fp:
        fps["taxa_and_lexicon_extraction"] = taxa_lex_fp
    for stage in _OCR_DEPENDENT_STAGES:
        fp = dict(fps.get(stage, {}))
        if ocrlang:
            fp["ocrlang"] = ocrlang
        if ocrmode:
            fp["ocrmode"] = ocrmode
        if keeppages:
            fp["keeppages"] = keeppages
        fps[stage] = fp
    if config_fingerprints is not None:
        for stage, config in config_fingerprints.items():
            fp = fps.setdefault(stage, {})
            fp["config"] = config
            if stage in ("figure_materialization", "figure_crossref"):
                for key, value in (("ocrlang", ocrlang), ("ocrmode", ocrmode), ("keeppages", keeppages)):
                    if value:
                        fp[key] = value
    return fps


_DOCLING_IMAGE_PLACEHOLDER_RE = re.compile(
    r"<!--\s*image\s*-->", re.IGNORECASE,
)


def _meaningful_extracted_text(raw: str) -> str:
    """Text that can count as recovered language for quality gates (#267).

    Docling writes one HTML image placeholder per graphic. Those markers are
    layout metadata, not transcription; counting them made an entirely empty
    figure-heavy document look progressively healthier as it acquired more
    images.
    """
    return _DOCLING_IMAGE_PLACEHOLDER_RE.sub("", raw or "").strip()


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

    raw_body = (text.get("text") or "") if isinstance(text, dict) else ""
    body = _meaningful_extracted_text(raw_body)
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
        # _gibberish_score itself now excludes CJK tokens, whose short
        # runs are ordinary words rather than corruption — see its
        # docstring for the Yamamori 2014 case.
        score = _gibberish_score(body[:50000])  # cap sample for speed
        if score > max_gibberish:
            flags.append({
                "gate": "gibberish_after_ocr",
                "severity": "error",
                "detail": f"gibberish_score={score:.2f} on extracted text (max: {max_gibberish})",
                "metric": round(score, 3),
            })

    # ocr_no_text_recovered — OCR can exit zero without adding text (#268).
    # Read both the explicit current field and the derivable older-artifact
    # shape, so re-running gates against an existing build surfaces it too.
    is_scan = isinstance(scan, dict)
    textless_pages = (scan.get("pages_without_text") or []) if is_scan else []
    ocr_page_count = int((scan.get("page_count") or 0) if is_scan else 0)
    all_ocr_pages_textless = bool(
        ocr_page_count
        and len({int(p) for p in textless_pages if str(p).isdigit()})
        >= ocr_page_count
    )
    if needs_ocr and (
        bool(scan.get("ocr_no_text_recovered")) or all_ocr_pages_textless
    ):
        total = ocr_page_count or pages or "?"
        flags.append({
            "gate": "ocr_no_text_recovered",
            "severity": "error",
            "detail": (
                f"OCR completed but all {total} page(s) have empty text layers"
            ),
            "metric": ocr_page_count or len(textless_pages),
        })

    # ocr_pages_blanked — pages the per-page OCR timeout gave up on and
    # left with no text (#254). Distinct from empty_text: a document can
    # read as "mostly fine" on every other gate while a third of it is
    # gone, which is exactly how ~9.5% of a full cluster build passed
    # unnoticed through two rebuilds. Not a heuristic — this is the
    # intersection of the pages ocrmypdf *named* on stderr with the pages
    # that ended up with no text, so plates and blank versos are excluded
    # by construction rather than by threshold. See pipeline/scan.py
    # `_report_ocr_page_loss`.
    blanked = (scan.get("pages_blanked") or []) if is_scan else []
    if blanked:
        # text.json's page count is the one the operator sees elsewhere;
        # scan_detection's is the fallback when extraction produced nothing.
        total = pages or (scan.get("page_count") if is_scan else None) or "?"
        shown = ", ".join(str(p) for p in blanked[:12])
        flags.append({
            "gate": "ocr_pages_blanked",
            "severity": "error",
            "detail": (
                f"{len(blanked)}/{total} page(s) hit the per-page OCR timeout "
                f"and were left blank: {shown}"
                f"{'…' if len(blanked) > 12 else ''}"
            ),
            "metric": len(blanked),
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


_STAGE_DEPENDENTS = {
    "scan_detection": ("pdf_preparation",),
    "pdf_preparation": ("docling_extraction", "metadata_extraction"),
    "docling_extraction": ("text_chunking", "figure_materialization"),
    "metadata_extraction": ("quality_gates",),
    "text_chunking": ("taxa_and_lexicon_extraction", "figure_crossref", "quality_gates"),
    # These producers share the figure tree. Invalidate both while either is
    # writing; a successful standalone vision overlay restores its CPU-floor
    # receipt, but an interrupted one must not leave that receipt reusable.
    "figure_materialization": ("vision_refresh", "figure_crossref", "quality_gates"),
    "vision_refresh": ("figure_materialization", "figure_crossref", "quality_gates"),
}


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
    if hash_dir is not None:
        # Clear dependent receipts before touching producer artifacts. This
        # also handles interruption between producer commit and consumer work,
        # even when a forced rerun used the same settings.
        state = _load_pipeline_state(hash_dir)
        pending = [name]
        invalidated = set()
        while pending:
            stage = pending.pop()
            if stage not in invalidated:
                invalidated.add(stage)
                pending.extend(_STAGE_DEPENDENTS.get(stage, ()))
        if invalidated.intersection(state["stages"]):
            for stage in invalidated:
                state["stages"].pop(stage, None)
            _save_pipeline_state(hash_dir, state)
    try:
        yield
    except BaseException as e:
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
