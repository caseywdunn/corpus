"""Per-paper pipeline orchestrator.

:func:`run_pdf_processing_pipeline` is the work loop for one PDF —
copies it into a temp dir, runs the huge-document gate, then walks
through scan_detection, pdf_preparation, docling_extraction,
metadata_extraction, text_chunking, the Pass 2.5/3a/3b/3c figure
passes, figure_crossref, taxa_and_lexicon_extraction, figures_report,
and quality_gates. Each stage is wrapped in :func:`_stage` so timing
and structured failures land on processing_summary.

Per-stage resume (#28) skips stages whose artifact is already on
disk; the caller passes ``resume=True`` and per-stage guards do the
checks. summary.json is the nuclear option (delete to force full
reprocessing); deleting any single artifact selectively re-runs only
that stage.
"""
from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path
from typing import Any, Dict, Optional

from bib import BibIndex

from . import stamp_artifact
from .annotate import _extract_taxa_and_lexicons
from .chunking import chunk_text
from .config import CONFIG
from .extract import extract_docling_content
from .figure_passes import (
    _crossref_chunks_and_figures,
    _pass25_annotate_figures,
    _pass3a_annotate_rois,
    _pass3b_annotate_rois,
)
from .figures import generate_figures_report, resolve_compound_figures
from .grobid_client import GrobidClient
from .metadata import extract_metadata
from .scan import detect_scan_type, prepare_pdf
from .taxa import TaxonomyDB
from .stages import (
    _HugeDocumentError,
    _expected_fingerprints_for_run,
    _pdf_page_count,
    _run_quality_gates,
    _should_run_stage,
    _stage,
    _utcnow_iso,
)

logger = logging.getLogger(__name__)


class _PaperLogAdapter(logging.LoggerAdapter):
    """Prepend ``[<pdf_stem>]`` to every record so per-stage progress
    lines stay attributable inside the per-paper log AND in the noisy
    multi-paper console stream (where docling/transformers/httpx output
    drowns out the unprefixed "Chunking text..." line)."""

    def process(self, msg, kwargs):
        return f"[{self.extra['pdf']}] {msg}", kwargs


def run_pdf_processing_pipeline(
    pdf_path: Path,
    hash_dir: Path,
    temp_dir: Path,
    grobid_client: Optional[GrobidClient] = None,
    taxonomy_db: Optional[TaxonomyDB] = None,
    lexicons: Optional[Dict[str, Dict[str, Dict]]] = None,
    content_aware_figures: bool = False,
    vision_backend=None,
    bib_index: Optional[BibIndex] = None,
    resume: bool = False,
    taxonomy_fingerprint: Optional[Dict[str, Any]] = None,
    lexicon_fingerprints: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Dict:
    """Run the per-PDF processing pipeline and return a summary dict.

    Called inside a :func:`per_pdf_file_log` block by the main loop, so all
    ``logging`` calls here are captured to ``<hash_dir>/pipeline.log``.
    """
    pdf_name = pdf_path.stem
    plog = _PaperLogAdapter(logger, {"pdf": pdf_name})
    temp_pdf = temp_dir / f"{pdf_name}.pdf"

    # Copy PDF to temp directory for processing
    shutil.copy2(pdf_path, temp_pdf)

    # Create subdirectories in hash directory
    figures_dir = hash_dir / "figures"
    figures_dir.mkdir(exist_ok=True)
    visualizations_dir = hash_dir / "visualizations"
    visualizations_dir.mkdir(exist_ok=True)

    processing_summary = {
        "original_pdf": str(pdf_path),
        "started_at": _utcnow_iso(),
        "processing_steps": [],   # legacy: stage names, success only
        "files_created": [],
        "errors": [],             # legacy free-text; kept for backwards compat
        "stage_timings": [],      # #34: per-stage timing, success or fail
        "stage_failures": [],     # #34: structured failure records
        "skipped_stages": [],     # #28: stages whose artifact was already present
    }

    try:
        # Huge-document gate (#35) — skip-and-flag PDFs above max_pages
        # before any expensive stage runs.
        with _stage(processing_summary, "huge_document_check", hash_dir=hash_dir):
            max_pages = int(CONFIG.get("huge_document", {}).get("max_pages", 5000))
            n_pages = _pdf_page_count(temp_pdf)
            if n_pages is not None:
                processing_summary["page_count"] = n_pages
                if n_pages > max_pages:
                    raise _HugeDocumentError(
                        f"PDF has {n_pages} pages (max: {max_pages}); skipping. "
                        "Pre-split this paper or raise huge_document.max_pages "
                        "to process it."
                    )
            processing_summary["processing_steps"].append("huge_document_check")

        # ── scan_detection ──────────────────────────────────────────
        detection_file = hash_dir / "scan_detection.json"
        if _should_run_stage("scan_detection", hash_dir=hash_dir,
                             resume=resume, processing_summary=processing_summary):
            with _stage(processing_summary, "scan_detection", hash_dir=hash_dir):
                plog.info("Detecting scan type...")
                detection_result = detect_scan_type(temp_pdf)
                with open(detection_file, "w") as f:
                    json.dump(stamp_artifact(detection_result), f, indent=2)
                processing_summary["files_created"].append(str(detection_file))
                processing_summary["processing_steps"].append("scan_detection")
        else:
            with open(detection_file) as f:
                detection_result = json.load(f)

        # ── pdf_preparation ─────────────────────────────────────────
        processed_pdf = hash_dir / "processed.pdf"
        if _should_run_stage("pdf_preparation", hash_dir=hash_dir,
                             resume=resume, processing_summary=processing_summary):
            with _stage(processing_summary, "pdf_preparation", hash_dir=hash_dir):
                plog.info("Preparing PDF...")
                prepare_pdf(temp_pdf, detection_result, processed_pdf)
                processing_summary["files_created"].append(str(processed_pdf))
                processing_summary["processing_steps"].append("pdf_preparation")

        # ── docling_extraction ─────────────────────────────────────
        text_file = hash_dir / "text.json"
        figures_file = hash_dir / "figures.json"
        docling_doc_file = hash_dir / "docling_doc.json"
        if _should_run_stage("docling_extraction", hash_dir=hash_dir,
                             resume=resume, processing_summary=processing_summary):
            with _stage(processing_summary, "docling_extraction", hash_dir=hash_dir):
                plog.info("Extracting text and figures...")
                extract_docling_content(
                    processed_pdf,
                    text_file,
                    figures_file,
                    figures_dir,
                    visualizations_dir,
                    docling_doc_output=docling_doc_file,
                    scan_file_type=detection_result.get("file_type"),
                )
                processing_summary["files_created"].extend([str(text_file), str(figures_file)])
                if docling_doc_file.exists():
                    processing_summary["files_created"].append(str(docling_doc_file))
                processing_summary["processing_steps"].append("docling_extraction")

        # ── metadata_extraction ─────────────────────────────────────
        metadata_file = hash_dir / "metadata.json"
        references_file = hash_dir / "references.json"
        tei_file = hash_dir / "grobid.tei.xml"
        if _should_run_stage("metadata_extraction", hash_dir=hash_dir,
                             resume=resume, processing_summary=processing_summary):
            with _stage(processing_summary, "metadata_extraction", hash_dir=hash_dir):
                plog.info("Extracting metadata (Grobid)...")
                bib_entry = bib_index.lookup(pdf_path.name) if bib_index is not None else None
                extract_metadata(
                    processed_pdf,
                    metadata_file,
                    references_output=references_file,
                    tei_output=tei_file,
                    grobid_client=grobid_client,
                    original_filename=pdf_path.name,
                    bib_entry=bib_entry,
                )
                processing_summary["files_created"].extend(
                    [str(metadata_file), str(references_file)]
                )
                if tei_file.exists():
                    processing_summary["files_created"].append(str(tei_file))
                processing_summary["processing_steps"].append("metadata_extraction")

        # ── text_chunking ──────────────────────────────────────────
        chunks_file = hash_dir / "chunks.json"
        if _should_run_stage("text_chunking", hash_dir=hash_dir,
                             resume=resume, processing_summary=processing_summary):
            with _stage(processing_summary, "text_chunking", hash_dir=hash_dir):
                plog.info("Chunking text...")
                chunk_text(text_file, chunks_output=chunks_file)
                processing_summary["files_created"].append(str(chunks_file))
                processing_summary["processing_steps"].append("text_chunking")

        with _stage(processing_summary, "figure_pass25_annotation", hash_dir=hash_dir):
            plog.info("Pass 2.5: annotating figures from captions + text...")
            _pass25_annotate_figures(text_file, figures_file)
            processing_summary["processing_steps"].append("figure_pass25_annotation")

        if vision_backend is not None:
            with _stage(processing_summary, "figure_pass3b_rois", hash_dir=hash_dir):
                plog.info("Pass 3b: vision-model-driven panel + compound detection...")
                _pass3b_annotate_rois(figures_file, vision_backend)
                processing_summary["processing_steps"].append("figure_pass3b_rois")
        elif content_aware_figures:
            with _stage(processing_summary, "figure_pass3a_rois", hash_dir=hash_dir):
                plog.info("Pass 3a: OCR-driven panel ROI detection...")
                _pass3a_annotate_rois(figures_file)
                processing_summary["processing_steps"].append("figure_pass3a_rois")

        # Pass 3c — resolve *_compound figures: split their ROIs, match
        # to missing_figures, rename the PNG to range notation. Cheap;
        # worth running whenever 3a or 3b has been run.
        if vision_backend is not None or content_aware_figures:
            with _stage(processing_summary, "figure_pass3c_resolve", hash_dir=hash_dir):
                plog.info("Pass 3c: compound figure resolution + file rename...")
                summary_3c = resolve_compound_figures(figures_file)
                plog.info(
                    "Pass 3c: %d resolved, %d renamed, %d unchanged, %d new records",
                    summary_3c.get("resolved", 0),
                    summary_3c.get("renamed", 0),
                    summary_3c.get("unchanged", 0),
                    summary_3c.get("new_records", 0),
                )
                processing_summary["processing_steps"].append("figure_pass3c_resolve")

        with _stage(processing_summary, "figure_crossref", hash_dir=hash_dir):
            plog.info("Linking chunks to figures...")
            _crossref_chunks_and_figures(figures_file, chunks_file)
            processing_summary["processing_steps"].append("figure_crossref")

        # ── taxa_and_lexicon_extraction ─────────────────────────────────
        # Stage runs only when a taxonomy DB or at least one lexicon
        # category is configured. The input_fingerprint captures the
        # taxonomy + per-category content hashes, so editing one
        # lexicon section forces this stage to re-run on --resume.
        run_taxa_anat = taxonomy_db is not None or bool(lexicons)
        if run_taxa_anat:
            # Build via the shared helper (#56) so the outer per-doc
            # gate in main.py and this inner per-stage gate stay in
            # lockstep — same dict shape, same staleness semantics.
            taxa_anat_fingerprint = _expected_fingerprints_for_run(
                taxonomy_fingerprint=taxonomy_fingerprint if taxonomy_db is not None else None,
                lexicon_fingerprints=lexicon_fingerprints,
            ).get("taxa_and_lexicon_extraction", {})
            if _should_run_stage(
                "taxa_and_lexicon_extraction",
                hash_dir=hash_dir,
                resume=resume,
                processing_summary=processing_summary,
                expected_fingerprint=taxa_anat_fingerprint,
            ):
                with _stage(processing_summary, "taxa_and_lexicon_extraction",
                            hash_dir=hash_dir, input_fingerprint=taxa_anat_fingerprint):
                    plog.info("Extracting taxa + lexicon mentions...")
                    taxa_anat_files = _extract_taxa_and_lexicons(
                        chunks_file,
                        hash_dir,
                        taxonomy_db,
                        lexicons,
                        taxonomy_fingerprint=taxonomy_fingerprint,
                        lexicon_fingerprints=lexicon_fingerprints,
                    )
                    processing_summary["files_created"].extend(str(p) for p in taxa_anat_files)
                    if taxa_anat_files:
                        processing_summary["processing_steps"].append("taxa_and_lexicon_extraction")

        with _stage(processing_summary, "figures_report", hash_dir=hash_dir):
            plog.info("Generating figures report...")
            report_path = generate_figures_report(hash_dir)
            if report_path:
                processing_summary["files_created"].append(str(report_path))
                processing_summary["processing_steps"].append("figures_report")

        # Quality gates (#36) — informational silent-failure detectors.
        # Run after success so artifacts are populated. A failed gate
        # records a quality_flag in summary.json but does not fail the
        # paper; corpus_status.py (#40) rolls these up for review.
        with _stage(processing_summary, "quality_gates", hash_dir=hash_dir):
            qgs = _run_quality_gates(hash_dir)
            processing_summary["quality_flags"] = qgs
            if qgs:
                plog.warning(
                    "Quality gates flagged %d issue(s): %s",
                    len(qgs), ", ".join(g["gate"] for g in qgs),
                )
            processing_summary["processing_steps"].append("quality_gates")

        processing_summary["status"] = "success"

    except Exception as e:
        processing_summary["status"] = "error"
        processing_summary["errors"].append(str(e))
        plog.exception("Error processing PDF: %s", e)

    processing_summary["ended_at"] = _utcnow_iso()
    return processing_summary

