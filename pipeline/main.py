"""argparse + main loop for the corpus pipeline.

Parses CLI flags, loads config + taxonomy + lexicons, hashes input
PDFs, slices into batches for SLURM job arrays, and dispatches each
paper through :func:`run_pdf_processing_pipeline`.

Read-only audit (--audit-orphans) and dry-run reporting paths short-
circuit before any processing.
"""
from __future__ import annotations

import argparse
import json
import logging
import multiprocessing
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

from bib import BibIndex

from . import config as _pipeline_config
from .annotate import _extract_taxa_and_lexicons
from .chunking import ingest_to_vector_db
from .config import CONFIG, load_config
from .figure_passes import _pass3b_annotate_rois
from .grobid_client import GrobidClient
from .io import (
    HASH_PREFIX_LEN,
    _verify_or_raise_collision,
    audit_orphans,
    create_output_structure,
    create_summary_json,
    find_all_pdfs,
    short_hash,
)
from .log import per_pdf_file_log, setup_root_logging
from .runner import run_pdf_processing_pipeline
from .taxa import TaxonomyDB, lexicon_fingerprints, load_lexicon
from .stages import (
    _all_stage_artifacts_complete,
    _expected_fingerprints_for_run,
    _file_sha256,
)

logger = logging.getLogger(__name__)


def _slice_hashes_for_batch(
    pdf_map: Dict[str, List[Path]],
    batch_index: int,
    batch_size: int,
) -> tuple[List[str], int, int]:
    """Deterministic batch slice over the sorted hash list (#55).

    The slice is computed on the *unfiltered* hash list. ``find_all_pdfs``
    output depends only on the input directory contents, so every SLURM
    array task — regardless of when it happens to start — produces the
    same sorted hash order and therefore disjoint slices.

    Resume-skip happens per-doc inside the main loop, after slicing.
    Pre-filtering completed hashes before slicing (the prior behavior)
    made the slice depend on disk state at task-start time, which
    differs across array tasks: tasks starting later see fewer remaining
    hashes, so their slice indices land on different members of the
    list. That produced overlapping batches that raced on the same
    summary.json / pipeline_state.json files.

    Returns ``(batch_hashes, total_hashes, total_batches)``.
    """
    all_hashes = sorted(pdf_map.keys())
    total = len(all_hashes)
    total_batches = -(-total // batch_size) if total else 0
    start = batch_index * batch_size
    end = start + batch_size
    return all_hashes[start:end], total, total_batches


def _expected_stages_for_run(
    *,
    taxonomy_db: Any,
    lexicons: Optional[Dict[str, Any]],
) -> List[str]:
    """Return the list of stage names this run is configured to produce.

    Always includes the core stages (scan_detection, pdf_preparation,
    docling_extraction, metadata_extraction, text_chunking).
    ``taxa_and_lexicon_extraction`` is added when a taxonomy DB or any
    lexicon category is configured.
    """
    stages = [
        "scan_detection",
        "pdf_preparation",
        "docling_extraction",
        "metadata_extraction",
        "text_chunking",
    ]
    if taxonomy_db is not None or lexicons:
        stages.append("taxa_and_lexicon_extraction")
    return stages


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
             "Build with: python -m pipeline.taxonomy_ingest --source <dwc|dwca|worms> ...",
    )
    parser.add_argument(
        "--lexicon",
        type=Path,
        default=None,
        help="Path to a multi-category lexicon YAML. Top-level keys are "
             "categories (anatomy, biogeography, …); see demo/lexicon.yaml "
             "for the format. Each category emits its own "
             "<hash>/<category>.json artifact and is fingerprinted "
             "separately, so editing one section only invalidates that "
             "category's annotations on --resume.",
    )
    parser.add_argument(
        "--no-taxa",
        action="store_true",
        help="Skip the taxa_and_lexicon_extraction stage entirely "
             "(taxon mentions and every --lexicon category).",
    )
    parser.add_argument(
        "--content-aware-figures",
        action="store_true",
        help="Run Pass 3a — OCR-driven panel/figure ROI detection on multi-panel "
             "figures (adds ~1-2 s per figure with caption panels; opt-in because "
             "OCR reliability on line-art figures is mixed; see dev_docs/PLAN.md §9)",
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
        from .external import set_strict_network
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

    # Open taxonomy snapshot and (if supplied) the multi-category
    # lexicon. Both are optional; missing inputs are logged and their
    # output artifacts are skipped. The lexicon is opt-in via
    # --lexicon — there is no default lookup because it's a
    # domain-specific user input.
    taxonomy_db: Optional[TaxonomyDB] = None
    lexicons: Dict[str, Dict[str, Dict]] = {}
    taxonomy_fingerprint: Optional[Dict[str, Any]] = None
    lex_fingerprints: Dict[str, Dict[str, Any]] = {}
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
                "Build it with: python -m pipeline.taxonomy_ingest %s --source <dwc|dwca|worms> ...",
                taxonomy_path, args.output_dir,
            )

        if args.lexicon is not None:
            if args.lexicon.exists():
                try:
                    lexicons = load_lexicon(args.lexicon)
                    lex_fingerprints = lexicon_fingerprints(args.lexicon)
                    for category, section in lexicons.items():
                        fp = lex_fingerprints.get(category, {})
                        sha = fp.get("sha256", "?" * 12)
                        logger.info(
                            "Lexicon[%s] loaded from %s (%d terms, sha256=%s…)",
                            category, args.lexicon, len(section), sha[:12],
                        )
                except Exception as e:
                    logger.warning(
                        "Could not load lexicon %s: %s", args.lexicon, e,
                    )
            else:
                logger.warning(
                    "Lexicon %s not found — lexicon extraction skipped",
                    args.lexicon,
                )

    # Vision backend for Pass 3b. Constructed once and reused so the
    # backend can keep long-lived state (API client, loaded model, etc.).
    vision_backend = None
    if args.vision_backend:
        try:
            from .vision import get_vision_backend
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

    # ── Batch slicing (for SLURM job arrays) ────────────────────────
    # Slice BEFORE applying the resume filter (#55). The slice has to be
    # computed on the unfiltered hash list so every array task produces
    # the same partition regardless of when it starts; otherwise tasks
    # starting later see a shorter list and their slice indices shift,
    # producing overlapping batches that race on per-doc writes.
    # Resume skipping happens per-doc inside the main loop.
    if args.batch_index is not None:
        batch_hashes, total, total_batches = _slice_hashes_for_batch(
            pdf_map, args.batch_index, args.batch_size,
        )
        logger.info(
            "Batch %d/%d: processing hashes %d–%d (%d of %d total)",
            args.batch_index, max(total_batches - 1, 0),
            args.batch_index * args.batch_size,
            min(args.batch_index * args.batch_size + args.batch_size, total) - 1,
            len(batch_hashes), total,
        )
        if not batch_hashes:
            logger.warning(
                "Batch %d is empty (only %d hashes exist) — nothing to do",
                args.batch_index, total,
            )
            sys.exit(0)
        pdf_map = {h: pdf_map[h] for h in batch_hashes}

    # Log resume coverage for the (now sliced) work set so ops can see
    # how much of this batch will short-circuit on the per-doc guards.
    # Skipped when --refresh-vision targets the already-completed
    # population on purpose.
    if args.resume and not args.refresh_vision:
        expected_stages = _expected_stages_for_run(
            taxonomy_db=taxonomy_db,
            lexicons=lexicons,
        )
        expected_fingerprints = _expected_fingerprints_for_run(
            taxonomy_fingerprint=taxonomy_fingerprint,
            lexicon_fingerprints=lex_fingerprints,
        )
        completed_in_scope = sum(
            1
            for h in pdf_map
            if (documents_dir / short_hash(h) / "summary.json").exists()
            and _all_stage_artifacts_complete(
                documents_dir / short_hash(h),
                expected_stages=expected_stages,
                expected_fingerprints=expected_fingerprints,
            )
        )
        if completed_in_scope:
            logger.info(
                "Resume: %d of %d documents in scope are already complete "
                "(per-doc guards will skip them)",
                completed_in_scope, len(pdf_map),
            )

    if args.dry_run:
        n_would_full = n_would_partial = n_would_skip = 0
        expected_stages = _expected_stages_for_run(
            taxonomy_db=taxonomy_db,
            lexicons=lexicons,
        )
        expected_fingerprints = _expected_fingerprints_for_run(
            taxonomy_fingerprint=taxonomy_fingerprint,
            lexicon_fingerprints=lex_fingerprints,
        )
        for h in pdf_map:
            hd = documents_dir / short_hash(h)
            if not args.resume:
                n_would_full += 1
                continue
            if (hd / "summary.json").exists() and _all_stage_artifacts_complete(
                hd,
                expected_stages=expected_stages,
                expected_fingerprints=expected_fingerprints,
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

        paper_total = len(pdf_map)
        for paper_idx, (pdf_hash_full, pdf_paths) in enumerate(pdf_map.items(), start=1):
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
                            "[%d/%d] %s — skipping vision refresh (no figures.json; hash %s)",
                            paper_idx, paper_total,
                            pdf_paths[0].name, pdf_hash,
                        )
                        continue
                    logger.info(
                        "[%d/%d] %s — refreshing Pass 3b (hash %s)",
                        paper_idx, paper_total,
                        pdf_paths[0].name, pdf_hash,
                    )
                    with per_pdf_file_log(hash_dir) as log_path:
                        logger.info("pipeline.log: %s (refresh-vision)", log_path)
                        try:
                            _pass3b_annotate_rois(figures_file, vision_backend)
                        except Exception as e:
                            logger.exception(
                                "Pass 3b refresh failed on %s: %s", pdf_hash, e
                            )
                    continue
                # Per-stage resume (#28, #56): if every required stage
                # is recorded as complete in pipeline_state.json under
                # the current PIPELINE_VERSION *and* the matching
                # input_fingerprint (taxonomy + lexicons), skip the whole
                # paper (fast path). Otherwise fall through —
                # run_pdf_processing_pipeline runs only the stages that
                # aren't recorded complete or whose fingerprint is stale.
                if _all_stage_artifacts_complete(
                    hash_dir,
                    expected_stages=_expected_stages_for_run(
                        taxonomy_db=taxonomy_db,
                        lexicons=lexicons,
                    ),
                    expected_fingerprints=_expected_fingerprints_for_run(
                        taxonomy_fingerprint=taxonomy_fingerprint,
                        lexicon_fingerprints=lex_fingerprints,
                    ),
                ):
                    logger.info(
                        "[%d/%d] %s — skipping (all stages complete; hash %s)",
                        paper_idx, paper_total,
                        pdf_paths[0].name, pdf_hash,
                    )
                    continue
                logger.info(
                    "[%d/%d] %s — resuming (re-running missing or stale stages; hash %s)",
                    paper_idx, paper_total,
                    pdf_paths[0].name, pdf_hash,
                )

            hash_dir.mkdir(exist_ok=True)

            # Use the first copy for processing (they're all identical by hash)
            primary_pdf = pdf_paths[0]

            sep = "─" * 72
            logger.info(sep)
            logger.info(
                "[%d/%d] %s  (hash %s, %d cop%s)",
                paper_idx, paper_total,
                primary_pdf.relative_to(input_dir),
                pdf_hash,
                len(pdf_paths),
                "y" if len(pdf_paths) == 1 else "ies",
            )
            logger.info(sep)
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
                        lexicons=lexicons,
                        content_aware_figures=args.content_aware_figures,
                        vision_backend=vision_backend,
                        bib_index=bib_index,
                        resume=args.resume,
                        taxonomy_fingerprint=taxonomy_fingerprint,
                        lexicon_fingerprints=lex_fingerprints,
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
