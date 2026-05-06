#!/usr/bin/env python3
"""
Generalized corpus processing script that takes input_dir and output_dir arguments.
Recursively finds PDFs in input_dir, processes them with hash-based organization.
"""

import argparse
import hashlib
import json
import logging
import re
import shutil
import subprocess
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple
import multiprocessing
import tempfile
import os

from grobid_client import (
    GrobidClient,
    GrobidUnavailableError,
    parse_tei_header,
    parse_tei_intext_citations,
    parse_tei_references,
)
from bib import BibIndex, bib_entry_to_metadata
from figures import (
    extract_caption_info,
    parse_figure_number,
    parse_panels_from_caption,
    detect_missing_figures,
    detect_figure_rois,
    detect_figure_rois_via_vision,
    resolve_compound_figures,
    link_chunks_to_figures,
    generate_figures_report,
)
from taxa import (
    TaxonomyDB,
    extract_taxon_mentions,
    extract_lexicon_mentions,
    load_lexicon,
)


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pipeline foundations (#42 / #44 partial split)
# ---------------------------------------------------------------------------
#
# Cross-cutting infrastructure (config, logging, stage runner, failures,
# quality gates, resume helpers, hashing/orphan/summary I/O) lives in
# the ``pipeline/`` package. Imported here so existing references in
# the rest of this file keep working unchanged. Re-exported at the
# module level so callers that do ``from process_corpus import …``
# (tests, downstream tools) keep working during the deprecation
# window. Per-stage extraction is open as #45.
from pipeline.config import (
    _DEFAULT_CONFIG,
    _deep_merge,
    classify_section,
    load_config,
)
from pipeline import config as _pipeline_config
from pipeline.config import CONFIG  # mutable module-level singleton
from pipeline.log import per_pdf_file_log, setup_root_logging
from pipeline.stages import (
    _all_stage_artifacts_complete,
    _classify_exception,
    _file_sha256,
    _HugeDocumentError,
    _pdf_page_count,
    _REASON_CODES,
    _run_quality_gates,
    _safe_load_json,
    _should_run_stage,
    _stage,
    _stage_artifacts_present,
    _utcnow_iso,
)
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
from pipeline.annotate import _extract_taxa_and_anatomy
from pipeline.chunking import chunk_text, ingest_to_vector_db
from pipeline.metadata import _write_placeholder_metadata, _year_from_filename, extract_metadata
from pipeline.figure_passes import (
    _crossref_chunks_and_figures,
    _pass25_annotate_figures,
    _pass3a_annotate_rois,
    _pass3b_annotate_rois,
)
from pipeline.scan import (
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
from pipeline.extract import extract_docling_content
from pipeline.runner import run_pdf_processing_pipeline




def main():
    parser = argparse.ArgumentParser(
        description="Process a corpus of PDFs with hash-based organization"
    )
    parser.add_argument("input_dir", type=Path, help="Input directory containing PDFs")
    parser.add_argument("output_dir", type=Path, help="Output directory for processed files")
    parser.add_argument("--resume", action="store_true", help="Skip documents whose summary.json already exists")
    parser.add_argument("--config", type=Path, default=None, help="Path to config.yaml (defaults to ./config.yaml)")
    parser.add_argument(
        "--grobid-url",
        default=os.environ.get("GROBID_URL", "http://localhost:8070"),
        help="Grobid service URL (default: $GROBID_URL or http://localhost:8070); "
        "pass --grobid-url='' or --no-grobid to skip metadata extraction",
    )
    parser.add_argument(
        "--no-grobid",
        action="store_true",
        help="Skip Grobid even if reachable (useful for dev iteration on non-metadata stages)",
    )
    parser.add_argument(
        "--bib",
        type=Path,
        default=None,
        help="Optional BibTeX file. Records whose 'file = {Foo.pdf}' field "
             "matches an input PDF supply the metadata header (title, authors, "
             "year, journal, DOI) instead of Grobid. References are still "
             "parsed from each PDF by Grobid.",
    )
    parser.add_argument(
        "--taxonomy-db",
        type=Path,
        default=None,
        help="Path to Darwin Core taxonomy SQLite "
             "(default: <output_dir>/taxonomy.sqlite). "
             "Build with: python ingest_taxonomy.py --source <dwc|dwca|worms> ...",
    )
    parser.add_argument(
        "--anatomy-lexicon",
        type=Path,
        default=None,
        help="Path to a domain-specific anatomy lexicon YAML. Optional and "
             "user-supplied — see demo/anatomy_lexicon.yaml for the format. "
             "Without this flag, anatomy extraction is skipped.",
    )
    parser.add_argument(
        "--lexicon",
        action="append",
        default=[],
        metavar="CATEGORY:PATH",
        help="Additional category-tagged lexicon YAML, e.g. "
             "--lexicon biogeography:/path/to/biogeo.yaml. Repeatable. "
             "Same YAML format as --anatomy-lexicon (see "
             "demo/anatomy_lexicon.yaml). Output lands in "
             "<hash>/<CATEGORY>.json; per-stage resume + corpus_status "
             "treat each category independently. (#24)",
    )
    parser.add_argument(
        "--no-taxa",
        action="store_true",
        help="Skip taxon + anatomy extraction",
    )
    parser.add_argument(
        "--content-aware-figures",
        action="store_true",
        help="Run Pass 3a — OCR-driven panel/figure ROI detection on multi-panel "
             "figures (adds ~1-2 s per figure with caption panels; opt-in because "
             "OCR reliability on line-art figures is mixed; see PLAN.md §9)",
    )
    parser.add_argument(
        "--vision-backend",
        choices=["claude", "local"],
        default=None,
        help="Run Pass 3b — vision-LLM-driven panel + compound-figure detection. "
             "'claude' uses the Anthropic API (needs ANTHROPIC_API_KEY); "
             "'local' uses an open-weights VLM on CUDA/MPS (Bouchet production). "
             "Supersedes --content-aware-figures when both are set.",
    )
    parser.add_argument(
        "--vision-model",
        default=None,
        help="Override the per-backend default vision model (e.g. "
             "claude-sonnet-4-6-20251001 for higher quality than Haiku).",
    )
    parser.add_argument(
        "--refresh-vision",
        action="store_true",
        help="With --resume and --vision-backend: instead of skipping hashes whose "
             "summary.json already exists, re-run ONLY Pass 3b on each hash's "
             "existing figures.json. No OCR/Docling/Grobid/chunking is re-done.",
    )
    parser.add_argument(
        "--audit-orphans",
        action="store_true",
        help="Read-only audit. List documents/<HASH>/ directories whose source "
             "PDF is no longer in input_dir, and LanceDB rows whose hash has no "
             "documents/ directory. Re-hashes input PDFs so it's path-independent. "
             "Does not delete anything.",
    )
    parser.add_argument(
        "--strict-network",
        action="store_true",
        help="Fail fast on the first transient external-service failure "
             "(Grobid 5xx, connect error, timeout) instead of retrying. "
             "Use for release-build runs where silent partial data is "
             "worse than aborting.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Discover input PDFs, apply --resume + batch filters, and report "
             "the work plan without processing anything. No external services "
             "are contacted; no files are written.",
    )
    parser.add_argument(
        "--batch-index",
        type=int,
        default=None,
        help="0-based index of the batch to process (use with --batch-size). "
             "Typically set to $SLURM_ARRAY_TASK_ID.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=256,
        help="Number of unique PDFs (hashes) per batch (default: 256).",
    )

    args = parser.parse_args()

    if args.refresh_vision and not args.vision_backend:
        parser.error("--refresh-vision requires --vision-backend to be set")

    setup_root_logging()

    if args.strict_network:
        from external import set_strict_network
        set_strict_network(True)

    # Mutate the pipeline.config singleton in place so all readers
    # (per-stage modules, quality_gates, etc.) see the loaded values.
    # Reassignment via ``global CONFIG = ...`` would only rebind the
    # local re-export — the canonical dict in pipeline.config would
    # stay at defaults.
    loaded = load_config(args.config)
    _pipeline_config.CONFIG.clear()
    _pipeline_config.CONFIG.update(loaded)

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.exists():
        logger.error("Input directory %s does not exist", input_dir)
        sys.exit(1)

    if args.audit_orphans:
        audit_orphans(input_dir, output_dir)
        return

    logger.info("Processing PDFs from: %s", input_dir)
    logger.info("Output directory: %s", output_dir)

    # One-time Grobid health check. If the service is unreachable at
    # startup we log and carry on with placeholder metadata for every
    # document rather than retrying (and logging) per PDF.
    grobid_client: Optional[GrobidClient] = None
    if args.no_grobid or not args.grobid_url:
        logger.info("Grobid skipped (--no-grobid or empty --grobid-url)")
    else:
        probe = GrobidClient(base_url=args.grobid_url)
        if probe.is_alive():
            logger.info("Grobid reachable at %s", args.grobid_url)
            grobid_client = probe
        else:
            logger.warning(
                "Grobid not reachable at %s — metadata will be placeholder. "
                "Start it with: docker compose up -d grobid",
                args.grobid_url,
            )

    # Optional BibTeX-driven metadata override. Loaded once and shared
    # across workers — entries are looked up by PDF basename.
    bib_index: Optional[BibIndex] = None
    if args.bib is not None:
        if not args.bib.exists():
            logger.error("Bib file %s does not exist", args.bib)
            sys.exit(1)
        try:
            bib_index = BibIndex.from_path(args.bib)
        except Exception as e:
            logger.error("Could not parse %s: %s", args.bib, e)
            sys.exit(1)

    # Open taxonomy snapshot and (if supplied) load anatomy lexicon. Both
    # are optional; missing inputs are logged and their output artifacts
    # are skipped.  The lexicon is opt-in via --anatomy-lexicon — there
    # is no default lookup because it's a domain-specific user input.
    taxonomy_db: Optional[TaxonomyDB] = None
    anatomy_lexicon: Optional[Dict[str, Dict]] = None
    taxonomy_fingerprint: Optional[Dict[str, Any]] = None
    anatomy_fingerprint: Optional[Dict[str, Any]] = None
    if not args.no_taxa:
        taxonomy_path = args.taxonomy_db or (args.output_dir / "taxonomy.sqlite")
        if taxonomy_path.exists():
            try:
                taxonomy_db = TaxonomyDB(taxonomy_path)
                # Stamp #29: hash once at startup so per-paper writes are
                # cheap. SHA-256 of taxonomy.sqlite is a stable identifier
                # that survives copy/move and changes any time the DB is
                # rebuilt.
                taxonomy_fingerprint = {
                    "path": str(taxonomy_path),
                    "sha256": _file_sha256(taxonomy_path),
                    "size": taxonomy_path.stat().st_size,
                }
                logger.info(
                    "Taxonomy snapshot loaded from %s (%d names, sha256=%s…)",
                    taxonomy_path, len(taxonomy_db.name_set()),
                    taxonomy_fingerprint["sha256"][:12],
                )
            except Exception as e:
                logger.warning(
                    "Could not open taxonomy snapshot %s: %s", taxonomy_path, e,
                )
        else:
            logger.warning(
                "Taxonomy snapshot %s not found — taxon extraction skipped. "
                "Build it with: python ingest_taxonomy.py %s --source <dwc|dwca|worms> ...",
                taxonomy_path, args.output_dir,
            )

        if args.anatomy_lexicon is not None:
            if args.anatomy_lexicon.exists():
                try:
                    anatomy_lexicon = load_lexicon(args.anatomy_lexicon)
                    anatomy_fingerprint = {
                        "path": str(args.anatomy_lexicon),
                        "sha256": _file_sha256(args.anatomy_lexicon),
                        "size": args.anatomy_lexicon.stat().st_size,
                    }
                    logger.info(
                        "Anatomy lexicon loaded from %s (%d terms, sha256=%s…)",
                        args.anatomy_lexicon, len(anatomy_lexicon),
                        anatomy_fingerprint["sha256"][:12],
                    )
                except Exception as e:
                    logger.warning(
                        "Could not load anatomy lexicon %s: %s",
                        args.anatomy_lexicon, e,
                    )
            else:
                logger.warning(
                    "Anatomy lexicon %s not found — anatomy extraction skipped",
                    args.anatomy_lexicon,
                )

    # Multi-category lexicons (#24). --lexicon CATEGORY:PATH, repeatable.
    extra_lexicons: Dict[str, Dict[str, Dict]] = {}
    extra_fingerprints: Dict[str, Dict[str, Any]] = {}
    for spec in (args.lexicon or []):
        if ":" not in spec:
            logger.error("--lexicon expects CATEGORY:PATH, got %r", spec)
            sys.exit(2)
        category, _, raw_path = spec.partition(":")
        category = category.strip().lower()
        path = Path(raw_path).expanduser()
        if not category or not category.replace("_", "").isalnum():
            logger.error(
                "--lexicon CATEGORY must be alphanumeric/underscore (got %r)",
                category,
            )
            sys.exit(2)
        if not path.exists():
            logger.warning(
                "--lexicon %s:%s not found — %s extraction skipped",
                category, path, category,
            )
            continue
        try:
            lex = load_lexicon(path)
            fp = {
                "path": str(path),
                "sha256": _file_sha256(path),
                "size": path.stat().st_size,
            }
            extra_lexicons[category] = lex
            extra_fingerprints[category] = fp
            logger.info(
                "Lexicon[%s] loaded from %s (%d terms, sha256=%s…)",
                category, path, len(lex), fp["sha256"][:12],
            )
        except Exception as e:
            logger.warning(
                "Could not load --lexicon %s:%s: %s", category, path, e,
            )

    # Vision backend for Pass 3b. Constructed once and reused so the
    # backend can keep long-lived state (API client, loaded model, etc.).
    vision_backend = None
    if args.vision_backend:
        try:
            from vision import get_vision_backend
            kwargs = {}
            if args.vision_model:
                kwargs["model"] = args.vision_model
            vision_backend = get_vision_backend(args.vision_backend, **kwargs)
            logger.info("Vision backend loaded: %s", vision_backend.name)
        except Exception as e:
            logger.error(
                "Could not load vision backend %r: %s — Pass 3b will be skipped",
                args.vision_backend, e,
            )

    # Create output directory structure
    documents_dir, vector_db_dir = create_output_structure(output_dir)

    # Find all PDFs and group by hash
    logger.info("Discovering PDFs...")
    pdf_map = find_all_pdfs(input_dir)

    logger.info("Found %d PDF file(s)", sum(len(paths) for paths in pdf_map.values()))
    logger.info("Unique PDFs (by hash): %d", len(pdf_map))

    # ── Pre-filter completed documents before batch slicing ─────────
    # When --resume is active and we're batching, remove fully-completed
    # hashes *before* dividing into batches. A hash is fully complete iff
    # summary.json exists AND every per-stage artifact is on disk (#28);
    # otherwise it has missing stages and the inner per-stage guards will
    # only re-run those.
    #
    # Exception: --refresh-vision (#27) explicitly targets the
    # already-completed population (it re-runs Pass 3b on existing
    # figures.json files). Skipping the pre-filter for that case
    # preserves the work-set so array tasks aren't silently emptied.
    if args.resume and args.batch_index is not None and not args.refresh_vision:
        expect_taxa = taxonomy_db is not None
        expect_anatomy = bool(anatomy_lexicon)
        before = len(pdf_map)
        kept = {}
        for h, paths in pdf_map.items():
            hd = documents_dir / short_hash(h)
            if (hd / "summary.json").exists() and _all_stage_artifacts_complete(
                hd, expect_taxa=expect_taxa, expect_anatomy=expect_anatomy,
            ):
                continue
            kept[h] = paths
        pdf_map = kept
        skipped = before - len(pdf_map)
        if skipped:
            logger.info(
                "Resume: filtered out %d fully-completed documents "
                "(%d remaining to process)", skipped, len(pdf_map),
            )

    # ── Batch slicing (for SLURM job arrays) ────────────────────────
    if args.batch_index is not None:
        all_hashes = sorted(pdf_map.keys())
        total = len(all_hashes)
        total_batches = -(-total // args.batch_size)  # ceiling division
        start = args.batch_index * args.batch_size
        end = start + args.batch_size
        batch_hashes = all_hashes[start:end]
        logger.info(
            "Batch %d/%d: processing hashes %d–%d (%d of %d total)",
            args.batch_index, total_batches - 1,
            start, min(end, total) - 1,
            len(batch_hashes), total,
        )
        if not batch_hashes:
            logger.warning(
                "Batch %d is empty (only %d hashes exist) — nothing to do",
                args.batch_index, total,
            )
            sys.exit(0)
        pdf_map = {h: pdf_map[h] for h in batch_hashes}

    if args.dry_run:
        n_would_full = n_would_partial = n_would_skip = 0
        expect_taxa = taxonomy_db is not None
        expect_anatomy = bool(anatomy_lexicon)
        for h in pdf_map:
            hd = documents_dir / short_hash(h)
            if not args.resume:
                n_would_full += 1
                continue
            if (hd / "summary.json").exists() and _all_stage_artifacts_complete(
                hd, expect_taxa=expect_taxa, expect_anatomy=expect_anatomy,
            ):
                n_would_skip += 1
            elif (hd / "summary.json").exists():
                # Partial — per-stage guards will run only the missing stages
                n_would_partial += 1
            else:
                n_would_full += 1
        logger.info(
            "Dry-run: %d unique PDF(s) in scope; would full-process %d, "
            "partial-process %d, skip %d (--resume = %s). Vision backend: %s. "
            "Grobid: %s. No files written.",
            len(pdf_map), n_would_full, n_would_partial, n_would_skip,
            "on" if args.resume else "off",
            args.vision_backend or "off",
            "off (--no-grobid or empty URL)" if (args.no_grobid or not args.grobid_url) else args.grobid_url,
        )
        return

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)

        for pdf_hash_full, pdf_paths in pdf_map.items():
            pdf_hash = short_hash(pdf_hash_full)
            hash_dir = documents_dir / pdf_hash

            # Detect prefix collision against any prior run in this dir.
            try:
                prior_matches = _verify_or_raise_collision(hash_dir, pdf_hash_full)
            except RuntimeError as e:
                logger.error(str(e))
                sys.exit(2)

            if args.resume and prior_matches:
                if args.refresh_vision and vision_backend is not None:
                    figures_file = hash_dir / "figures.json"
                    if not figures_file.exists():
                        logger.info(
                            "Skipping %s for vision refresh (no figures.json)", pdf_hash
                        )
                        continue
                    logger.info("Refreshing Pass 3b on %s", pdf_hash)
                    with per_pdf_file_log(hash_dir) as log_path:
                        logger.info("pipeline.log: %s (refresh-vision)", log_path)
                        try:
                            _pass3b_annotate_rois(figures_file, vision_backend)
                        except Exception as e:
                            logger.exception(
                                "Pass 3b refresh failed on %s: %s", pdf_hash, e
                            )
                    continue
                # Per-stage resume (#28): if every artifact is on disk,
                # skip the whole paper (fast path). Otherwise fall through
                # — run_pdf_processing_pipeline runs only the missing
                # stages.
                if _all_stage_artifacts_complete(
                    hash_dir,
                    expect_taxa=taxonomy_db is not None,
                    expect_anatomy=bool(anatomy_lexicon),
                ):
                    logger.info("Skipping %s (all stages complete)", pdf_hash)
                    continue
                logger.info(
                    "Resuming %s (re-running missing stages only)", pdf_hash
                )

            hash_dir.mkdir(exist_ok=True)

            # Use the first copy for processing (they're all identical by hash)
            primary_pdf = pdf_paths[0]

            logger.info(
                "Processing PDF hash %s (%d copies); primary file: %s",
                pdf_hash, len(pdf_paths), primary_pdf.relative_to(input_dir),
            )
            if len(pdf_paths) > 1:
                for path in pdf_paths[1:]:
                    logger.info("  additional copy: %s", path.relative_to(input_dir))

            # Run each document in a child process so that C-level crashes
            # (segfaults in docling/PyMuPDF) don't kill the whole batch.
            # We use the fork start method so the child inherits already-built
            # objects (grobid_client, taxonomy_db, vision_backend, _worker
            # closure) without pickling.
            #
            # macOS exception: PyTorch + Apple's Objective-C frameworks
            # (Metal/MPSGraph) are not fork-safe. The forked child crashes
            # the moment docling loads torch. There is no env-var workaround
            # that holds — OBJC_DISABLE_INITIALIZE_FORK_SAFETY just turns the
            # abort into a SIGSEGV. So on macOS we run inline and rely on
            # --resume to recover from the rare segfault.
            def _worker():
                with per_pdf_file_log(hash_dir) as log_path:
                    logger.info("pipeline.log: %s", log_path)
                    logger.info("pdf_hash_full: %s", pdf_hash_full)

                    processing_summary = run_pdf_processing_pipeline(
                        primary_pdf, hash_dir, temp_dir,
                        grobid_client=grobid_client,
                        taxonomy_db=taxonomy_db,
                        anatomy_lexicon=anatomy_lexicon,
                        content_aware_figures=args.content_aware_figures,
                        vision_backend=vision_backend,
                        bib_index=bib_index,
                        resume=args.resume,
                        taxonomy_fingerprint=taxonomy_fingerprint,
                        anatomy_fingerprint=anatomy_fingerprint,
                        extra_lexicons=extra_lexicons,
                        extra_fingerprints=extra_fingerprints,
                    )

                    summary_file = create_summary_json(
                        pdf_hash_full, pdf_paths, input_dir, hash_dir, processing_summary
                    )
                    logger.info("Created summary: %s", summary_file)

                    if processing_summary.get("status") == "success":
                        chunks_file = hash_dir / "chunks.json"
                        if chunks_file.exists():
                            logger.info("Writing vector-db ingestion marker...")
                            ingest_to_vector_db(chunks_file, vector_db_dir, pdf_hash)

            if sys.platform == "darwin":
                _worker()
            else:
                mp_ctx = multiprocessing.get_context("fork")
                proc = mp_ctx.Process(target=_worker)
                proc.start()
                proc.join()
                if proc.exitcode != 0:
                    if proc.exitcode < 0:
                        logger.error(
                            "Document %s killed by signal %d (segfault?) — skipping",
                            pdf_hash, -proc.exitcode,
                        )
                    else:
                        logger.error(
                            "Document %s worker exited with code %d — skipping",
                            pdf_hash, proc.exitcode,
                        )

    logger.info("Processing complete. Results saved to: %s", output_dir)
    logger.info("  Documents: %s", documents_dir)
    logger.info("  Vector DB: %s", vector_db_dir)


if __name__ == "__main__":
    main()
