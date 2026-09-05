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
from typing import Any, Dict, List, Optional, Tuple

from bib import BibIndex, keeppages_for_pdf, ocrlang_for_pdf, ocrmode_for_pdf

from . import config as _pipeline_config
from . import external
from .annotate import _extract_taxa_and_lexicons
from .config import CONFIG, load_config
from .build_inputs import config_fingerprints as _config_fingerprints
from .figure_passes import _crossref_chunks_and_figures, _pass25_annotate_figures, _pass3b_annotate_rois
from .figure_materialization import rebuild_figure_base
from .extract import extract_docling_content
from .figures import resolve_compound_figures
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
    _metadata_fingerprint_for_pdf,
    _file_sha256,
    _stage,
    _record_stage_completion,
    _load_pipeline_state,
    _save_pipeline_state,
    _run_quality_gates,
)

logger = logging.getLogger(__name__)


def _refresh_vision_artifacts(
    figures_file: Path,
    chunks_file: Path,
    vision_backend,
    *,
    reset_base=False,
) -> None:
    """Refresh the complete vision-derived figure layer for one document.

    ``--only vision`` is the independently schedulable form of the figure
    vision phase. It must leave the same artifacts as the inline full-run
    path: Pass 3b evidence, Pass 3c compound materialization, and rebuilt
    chunk/figure links. Keeping those follow-on passes here prevents a GPU
    refresh from writing ROIs that the served bundle cannot actually reach.
    """
    if reset_base:
        rebuild_figure_base(figures_file.parent, extract_docling_content)
        _pass25_annotate_figures(figures_file.parent / "text.json", figures_file)
    _pass3b_annotate_rois(figures_file, vision_backend)
    summary_3c = resolve_compound_figures(figures_file)
    logger.info(
        "Pass 3c: %d resolved, %d renamed, %d unchanged, %d new records",
        summary_3c.get("resolved", 0),
        summary_3c.get("renamed", 0),
        summary_3c.get("unchanged", 0),
        summary_3c.get("new_records", 0),
    )
    _crossref_chunks_and_figures(figures_file, chunks_file)
    if reset_base:
        from .pageselect import annotate_source_pages
        from .figures import generate_figures_report
        scan_path = figures_file.parent / "scan_detection.json"
        scan = json.loads(scan_path.read_text()) if scan_path.exists() else {}
        if scan.get("keeppages_selected"):
            annotate_source_pages([figures_file, chunks_file], scan["keeppages_selected"])
        generate_figures_report(figures_file.parent)


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
    Annotation also records the empty configuration so removing its last
    input retires old outputs instead of skipping them forever.
    """
    stages = [
        "scan_detection",
        "pdf_preparation",
        "docling_extraction",
        "metadata_extraction",
        "text_chunking",
        "figure_materialization",
        "figure_crossref",
        "huge_document_check",
        "quality_gates",
        "taxa_and_lexicon_extraction",
    ]
    return stages


def _audit_corpus_chunks(documents_dir: Path) -> Tuple[int, int, List[str]]:
    """Scan a built ``documents/`` tree and report extraction yield.

    Returns ``(attempted, total_chunks, zero_chunk_hashes)`` where
    ``attempted`` counts per-document dirs that have a ``summary.json``
    (i.e. extraction ran), ``total_chunks`` sums ``chunks.json`` chunk
    counts across them, and ``zero_chunk_hashes`` lists the dirs that
    produced none. Missing/malformed ``chunks.json`` counts as zero.
    Used by the silent-failure guard (#99)."""
    attempted = 0
    total_chunks = 0
    zero_chunk_hashes: List[str] = []
    if not documents_dir.exists():
        return attempted, total_chunks, zero_chunk_hashes
    for hash_dir in sorted(documents_dir.iterdir()):
        if not hash_dir.is_dir() or not (hash_dir / "summary.json").exists():
            continue
        attempted += 1
        n_chunks = 0
        chunks_path = hash_dir / "chunks.json"
        if chunks_path.exists():
            try:
                n_chunks = len(json.loads(chunks_path.read_text()).get("chunks", []))
            except (OSError, json.JSONDecodeError):
                n_chunks = 0
        total_chunks += n_chunks
        if n_chunks == 0:
            zero_chunk_hashes.append(hash_dir.name)
    return attempted, total_chunks, zero_chunk_hashes


def _panels_to_legacy(mode: str) -> Tuple[bool, Optional[str]]:
    """Map the ``--figure-panels`` selector (#102) to the legacy
    ``(content_aware_figures, vision_backend_name)`` pair the runner and
    its refresh/pass-3c gates still consume.

    ``ocr`` → Pass 3a (content-aware, no vision backend); ``vision-local``
    / ``vision-claude`` → Pass 3b with the named backend; ``off`` →
    neither pass runs.
    """
    content_aware = mode == "ocr"
    backend = {"vision-local": "local", "vision-claude": "claude"}.get(mode)
    return content_aware, backend


def main():
    parser = argparse.ArgumentParser(
        description="Process a corpus of PDFs with hash-based organization"
    )
    parser.add_argument("input_dir", type=Path, help="Input directory containing PDFs")
    parser.add_argument("output_dir", type=Path, help="Output directory for processed files")
    parser.add_argument("--resume", action="store_true", help="Skip stages with current producer and input fingerprints")
    parser.add_argument("--config", type=Path, default=None, help="Path to config.yaml (defaults to ./config.yaml)")
    parser.add_argument(
        "--grobid-url",
        default=None,
        help="Grobid service URL (default: $GROBID_URL, config, then http://localhost:8070); "
        "pass --grobid-url='' or --no-grobid to disable Grobid-derived data",
    )
    parser.add_argument(
        "--no-grobid",
        action="store_true",
        help="Disable Grobid-derived metadata/references (archive active TEI; retain BibTeX headers)",
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
        "--require-taxonomy",
        action="store_true",
        help="Fail instead of warning when the taxonomy snapshot is "
             "missing (#139). Passed automatically by `corpus run` when "
             "the corpuscle configures taxonomy.source, so a batch job "
             "cannot quietly produce a corpus with empty taxa.json for "
             "every paper.",
    )
    parser.add_argument(
        "--figure-panels",
        choices=["ocr", "vision-local", "vision-claude", "off"],
        default=None,
        help="Panel-ROI detection mode (#102). 'ocr' (default) runs Pass 3a "
             "— OCR-driven panel ROIs, CPU-only, self-gated to multi-panel "
             "figures. 'vision-local' / 'vision-claude' run Pass 3b "
             "(open-weights VLM on CUDA/MPS, or the Anthropic API — needs "
             "ANTHROPIC_API_KEY). 'off' disables panel-ROI detection.",
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
        help="With --resume and --figure-panels vision-*: instead of skipping "
             "hashes whose summary.json already exists, re-run ONLY Pass 3b on "
             "each hash's existing figures.json. No OCR/Docling/Grobid/chunking "
             "is re-done.",
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

    setup_root_logging()

    if args.strict_network:
        from .external import set_strict_network
        set_strict_network(True)

    # Mutate the pipeline.config singleton in place so all readers
    # (per-stage modules, quality_gates, etc.) see the loaded values.
    # Reassignment via ``global CONFIG = ...`` would only rebind the
    # local re-export — the canonical dict in pipeline.config would
    # stay at defaults.
    # An explicitly requested config that cannot be read is an error, not a
    # fallback (#210). `load_config` treats a missing file as "use defaults",
    # which is right for the implicit ./config.yaml but silently discards
    # everything the operator tuned when they named a file and it was not
    # found. That is how a relative --config path used to fail: forwarded
    # verbatim to a subprocess running from a different directory, missed, and
    # replaced by defaults with a single INFO line as the only trace.
    if args.config is not None and not Path(args.config).exists():
        parser.error(
            f"--config {args.config}: no such file. Every setting in it would "
            f"be silently replaced by built-in defaults, so this is refused "
            f"rather than run."
        )
    loaded = load_config(args.config)
    grobid_config = loaded.get("grobid", {})
    if args.grobid_url is None:
        args.grobid_url = os.environ.get("GROBID_URL", grobid_config.get("url", "http://localhost:8070"))
    args.no_grobid = args.no_grobid or grobid_config.get("disable", False) or not args.grobid_url
    loaded["grobid"] = {**grobid_config, "disable": bool(args.no_grobid)}
    _pipeline_config.CONFIG.clear()
    _pipeline_config.CONFIG.update(loaded)
    figures_config = loaded.get("figures", {})
    args.figure_panels = args.figure_panels or figures_config.get("panel_detection", "ocr")
    if args.figure_panels not in ("off", "ocr", "vision-local", "vision-claude"):
        parser.error('figures.panel_detection must be off, ocr, vision-local or vision-claude (quote "off" in YAML)')
    args.vision_model = args.vision_model or figures_config.get("model")
    args.content_aware_figures, args.vision_backend = _panels_to_legacy(args.figure_panels)
    if args.refresh_vision and not args.vision_backend:
        parser.error("--refresh-vision requires --figure-panels vision-local|vision-claude")
    run_config_fingerprints = _config_fingerprints(
        loaded, panel_mode=args.figure_panels, vision_model=args.vision_model)

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
    grobid_context = {"enabled": not args.no_grobid, "available": False,
                      "service_version": None}
    if args.no_grobid or not args.grobid_url or args.dry_run:
        logger.info("Grobid skipped (--no-grobid or empty --grobid-url)")
    else:
        probe = GrobidClient(base_url=args.grobid_url,
                             timeout=loaded["stage_timeouts"]["grobid"])
        if probe.is_alive():
            logger.info("Grobid reachable at %s", args.grobid_url)
            grobid_client = probe
            grobid_context.update(available=True, service_version=probe.get_version())
        else:
            if external.STRICT_NETWORK:
                logger.error("Requested Grobid service is unavailable (--strict-network)")
                return 1
            logger.warning(
                "Grobid not reachable at %s — metadata will be placeholder. "
                "Start it with `docker compose up -d grobid` (laptop/Docker "
                "host), or on a cluster `sbatch slurm/batch_grobid.sh` then "
                "`export GROBID_URL=http://<node>:8070`.",
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
    # lexicon. Both are optional, but a configured unreadable source is
    # a failure, not permission to retire previous outputs. Lexicon is opt-in via
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
                logger.error(
                    "Could not open taxonomy snapshot %s: %s", taxonomy_path, e,
                )
                return 1
        elif args.require_taxonomy and args.dry_run:
            # A dry-run writes nothing, so it cannot produce the empty
            # taxa.json that #139 guards against. On a corpuscle that has
            # not been built yet the snapshot is *always* missing at this
            # point — ingest_taxonomy dry-runs too — so failing here made
            # `corpus run --dry-run` unusable as a first-run sanity check.
            logger.info(
                "Taxonomy snapshot %s not found; the ingest_taxonomy step "
                "builds it on a real run. Continuing the dry-run.",
                taxonomy_path,
            )
        elif args.require_taxonomy or args.taxonomy_db is not None:
            # #139 — the corpuscle configures a taxonomy, so a missing
            # snapshot is a hard error rather than a skip. This is the
            # layer where the damage actually happens: the first full
            # Bouchet production run put 1763 papers through with empty
            # taxa.json, and the only signal was the WARNING below.
            # `corpus run` pre-checks the same condition in
            # orchestrator._check_taxonomy_available; this is the
            # defense-in-depth copy that also covers direct
            # `python -m pipeline.main` invocations and a snapshot that
            # disappears between the pre-check and the work.
            logger.error(
                "Taxonomy snapshot %s not found, but taxonomy was explicitly "
                "configured — refusing to run and silently "
                "produce empty taxa.json for every paper.\n"
                "  Build it first:  corpus taxonomy ingest --source <dwc|dwca|worms> ...\n"
                "  WoRMS needs outbound internet: build the snapshot once "
                "up front, or export a DwC-A snapshot and switch config to "
                "`source: dwca` — which is also faster and version-pinned "
                "(see INSTALL.md and dev_docs/BOUCHET.md).\n"
                "  Genuinely don't want taxon extraction? Pass --no-taxa.",
                taxonomy_path,
            )
            return 1
        else:
            logger.warning(
                "Taxonomy snapshot %s not found — taxon extraction skipped. "
                "Build it with: corpus taxonomy ingest --source <dwc|dwca|worms> ...",
                taxonomy_path,
            )

        if args.lexicon is not None:
            if args.lexicon.exists():
                # Distinguish deliberate removal from an unreadable input.
                import yaml as _yaml
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
                except (OSError, _yaml.YAMLError) as e:
                    logger.error(
                        "Could not load lexicon %s: %s", args.lexicon, e,
                    )
                    return 1
            else:
                logger.error(
                    "Configured lexicon %s not found — refusing to change annotations",
                    args.lexicon,
                )
                return 1

    # Vision backend for Pass 3b. Constructed once and reused so the
    # backend can keep long-lived state (API client, loaded model, etc.).
    vision_backend = None
    if args.vision_backend and not args.dry_run:
        try:
            from .vision import get_vision_backend
            kwargs = {}
            if args.vision_model:
                kwargs["model"] = args.vision_model
            vision_backend = get_vision_backend(args.vision_backend, **kwargs)
            logger.info("Vision backend loaded: %s", vision_backend.name)
        except Exception as e:
            logger.error(
                "Could not load requested vision backend %r: %s — refusing to mark it complete",
                args.vision_backend, e,
            )
            return 1

    # Create output directory structure. A dry-run promises "No files
    # written", so it must not leave a half-scaffolded corpuscle behind
    # either — the paths are still computed, just not created.
    if args.dry_run:
        documents_dir = output_dir / "documents"
        vector_db_dir = output_dir / "vector_db"
    else:
        documents_dir, vector_db_dir = create_output_structure(output_dir)

    # Find all PDFs and group by hash. Exclude anything under output_dir
    # so re-runs don't re-ingest per-paper processed.pdf artifacts (which
    # have different SHA-256s than the originals because OCR adds a text
    # layer) — only matters when input_dir is an ancestor of output_dir
    # (the demo's `input_pdfs: .` case).
    logger.info("Discovering PDFs...")
    pdf_map = find_all_pdfs(input_dir, exclude_under=output_dir)

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
        # #176 — ocrlang is per-document, so the fingerprint map has to be
        # rebuilt inside the loop rather than hoisted out of it.
        completed_in_scope = sum(
            1
            for h in pdf_map
            if (documents_dir / short_hash(h) / "summary.json").exists()
            and _all_stage_artifacts_complete(
                documents_dir / short_hash(h),
                expected_stages=expected_stages,
                expected_fingerprints=_expected_fingerprints_for_run(
                    config_fingerprints=run_config_fingerprints,
                    metadata_fingerprint=_metadata_fingerprint_for_pdf(
                        bib_index, pdf_map[h][0].name, grobid_context=grobid_context,
                        hash_dir=documents_dir / short_hash(h)),
                    taxonomy_fingerprint=taxonomy_fingerprint,
                    lexicon_fingerprints=lex_fingerprints,
                    ocrlang=ocrlang_for_pdf(bib_index, pdf_map[h][0].name),
                    ocrmode=ocrmode_for_pdf(bib_index, pdf_map[h][0].name),
                    keeppages=keeppages_for_pdf(bib_index, pdf_map[h][0].name),
                ),
            )
        )
        if completed_in_scope:
            logger.info(
                "Resume: %d of %d documents in scope are already complete "
                "(per-doc guards will skip them)",
                completed_in_scope, len(pdf_map),
            )

    if args.dry_run:
        would_full: List[Tuple[str, str]] = []     # (filename, short_hash)
        would_partial: List[Tuple[str, str]] = []
        would_skip: List[Tuple[str, str]] = []
        expected_stages = _expected_stages_for_run(
            taxonomy_db=taxonomy_db,
            lexicons=lexicons,
        )
        for h, paths in pdf_map.items():
            sh = short_hash(h)
            hd = documents_dir / sh
            # Use the first PDF path's basename as the human-readable label;
            # in the rare multi-copy case the others are dedup-equivalent.
            label = paths[0].name
            entry = (label, sh)
            if not args.resume:
                would_full.append(entry)
                continue
            if (hd / "summary.json").exists() and _all_stage_artifacts_complete(
                hd,
                expected_stages=expected_stages,
                expected_fingerprints=_expected_fingerprints_for_run(
                    config_fingerprints=run_config_fingerprints,
                    metadata_fingerprint=_metadata_fingerprint_for_pdf(
                        bib_index, label, grobid_context=grobid_context, hash_dir=hd),
                    taxonomy_fingerprint=taxonomy_fingerprint,
                    lexicon_fingerprints=lex_fingerprints,
                    ocrlang=ocrlang_for_pdf(bib_index, label),
                    ocrmode=ocrmode_for_pdf(bib_index, label),
                    keeppages=keeppages_for_pdf(bib_index, label),
                ),
            ):
                would_skip.append(entry)
            elif (hd / "summary.json").exists():
                # Partial — per-stage guards will run only the missing stages
                would_partial.append(entry)
            else:
                would_full.append(entry)
        logger.info(
            "Dry-run: %d unique PDF(s) in scope; would full-process %d, "
            "partial-process %d, skip %d (--resume = %s). Vision backend: %s. "
            "Grobid: %s. No files written.",
            len(pdf_map),
            len(would_full), len(would_partial), len(would_skip),
            "on" if args.resume else "off",
            args.vision_backend or "off",
            "off (--no-grobid or empty URL)" if (args.no_grobid or not args.grobid_url) else args.grobid_url,
        )
        # Print per-paper buckets so the operator can see exactly which
        # files would be touched. Cap each bucket at 20 entries to keep
        # large-corpus output bounded; sort alphabetically by filename
        # for stable, scan-friendly output.
        def _emit_bucket(label: str, entries: List[Tuple[str, str]]) -> None:
            if not entries:
                return
            logger.info("  %s (%d):", label, len(entries))
            cap = 20
            for fn, sh in sorted(entries)[:cap]:
                logger.info("    %s (%s)", fn, sh)
            if len(entries) > cap:
                logger.info("    ... and %d more", len(entries) - cap)
        _emit_bucket("would full-process", would_full)
        _emit_bucket("would partial-process", would_partial)
        _emit_bucket("would skip", would_skip)
        return

    # Track papers whose extraction worker crashed (#issue: stage-1
    # child-process failures used to be log-only and the pipeline still
    # exited 0, so segfaulting papers vanished from the build with no
    # signal to the operator). Each entry: {hash, pdf, exitcode, signal}.
    worker_failures: List[Dict[str, Any]] = []

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
                        logger.error(
                            "[%d/%d] %s (%s) — cannot refresh vision (no figures.json)",
                            paper_idx, paper_total,
                            pdf_paths[0].name, pdf_hash,
                        )
                        worker_failures.append({"pdf_hash": pdf_hash, "exitcode": 1,
                                                "signal": None, "reason": "missing figures.json"})
                        continue
                    logger.info(
                        "[%d/%d] %s (%s) — refreshing Pass 3b",
                        paper_idx, paper_total,
                        pdf_paths[0].name, pdf_hash,
                    )
                    with per_pdf_file_log(hash_dir) as log_path:
                        logger.info("pipeline.log: %s (refresh-vision)", log_path)
                        try:
                            baseline = _load_pipeline_state(hash_dir)["stages"].get("figure_materialization")
                            summary_path = hash_dir / "summary.json"
                            summary = json.loads(summary_path.read_text())
                            processing = summary.setdefault("processing_summary", {})
                            with _stage({}, "vision_refresh", hash_dir=hash_dir,
                                        input_fingerprint={"config": run_config_fingerprints["figure_materialization"]}):
                                _refresh_vision_artifacts(
                                    figures_file, hash_dir / "chunks.json", vision_backend,
                                    reset_base=True,
                                )
                                with _stage(processing, "quality_gates", hash_dir=hash_dir,
                                            input_fingerprint={"config": run_config_fingerprints["quality_gates"]}):
                                    processing["quality_flags"] = _run_quality_gates(hash_dir)
                                    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8",
                                                                     dir=hash_dir, prefix=".summary-",
                                                                     delete=False) as output:
                                        temporary = Path(output.name)
                                        try:
                                            json.dump(summary, output, indent=2)
                                            output.close()
                                            temporary.replace(summary_path)
                                        finally:
                                            temporary.unlink(missing_ok=True)
                                _record_stage_completion(hash_dir, "figure_crossref",
                                                         input_fingerprint={"config": run_config_fingerprints["figure_crossref"]})
                                if baseline is not None:
                                    state = _load_pipeline_state(hash_dir)
                                    state["stages"]["figure_materialization"] = baseline
                                    _save_pipeline_state(hash_dir, state)
                        except Exception as e:
                            logger.exception(
                                "Vision refresh failed on %s: %s", pdf_hash, e
                            )
                            worker_failures.append({"pdf_hash": pdf_hash, "exitcode": 1,
                                                    "signal": None, "reason": "vision refresh failed"})
                    continue
                # Per-stage resume (#28, #56): if every required stage
                # is recorded as complete in pipeline_state.json under
                # the current PIPELINE_VERSION *and* the matching
                # input_fingerprint (taxonomy + lexicons + ocrlang), skip
                # the whole paper (fast path). Otherwise fall through —
                # run_pdf_processing_pipeline runs only the stages that
                # aren't recorded complete or whose fingerprint is stale.
                #
                # ocrlang must be passed here, not just to the per-stage
                # gate (#176). This is the gate that actually skips work,
                # and it runs first: omit the tag and a paper that had no
                # ocrlang last run compares {} == {} and is skipped whole,
                # so the per-stage gate never sees the tag the operator
                # just added.
                if _all_stage_artifacts_complete(
                    hash_dir,
                    expected_stages=_expected_stages_for_run(
                        taxonomy_db=taxonomy_db,
                        lexicons=lexicons,
                    ),
                    expected_fingerprints=_expected_fingerprints_for_run(
                        config_fingerprints=run_config_fingerprints,
                        metadata_fingerprint=_metadata_fingerprint_for_pdf(
                            bib_index, pdf_paths[0].name, grobid_context=grobid_context, hash_dir=hash_dir),
                        taxonomy_fingerprint=taxonomy_fingerprint,
                        lexicon_fingerprints=lex_fingerprints,
                        ocrlang=ocrlang_for_pdf(bib_index, pdf_paths[0].name),
                        ocrmode=ocrmode_for_pdf(bib_index, pdf_paths[0].name),
                        keeppages=keeppages_for_pdf(bib_index, pdf_paths[0].name),
                    ),
                ):
                    from .io import refresh_summary_source_paths
                    if refresh_summary_source_paths(pdf_hash_full, pdf_paths, input_dir, hash_dir):
                        logger.info("%s: refreshed source paths (extraction unchanged)", pdf_hash)
                    logger.info(
                        "[%d/%d] %s (%s) — skipping (all stages complete)",
                        paper_idx, paper_total,
                        pdf_paths[0].name, pdf_hash,
                    )
                    continue
                logger.info(
                    "[%d/%d] %s (%s) — resuming (re-running missing or stale stages)",
                    paper_idx, paper_total,
                    pdf_paths[0].name, pdf_hash,
                )

            hash_dir.mkdir(exist_ok=True)

            # Use the first copy for processing (they're all identical by hash)
            primary_pdf = pdf_paths[0]

            sep = "─" * 72
            logger.info(sep)
            logger.info(
                "[%d/%d] %s (%s)%s",
                paper_idx, paper_total,
                primary_pdf.relative_to(input_dir),
                pdf_hash,
                "" if len(pdf_paths) == 1 else f" — {len(pdf_paths)} copies",
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
                        grobid_context=grobid_context,
                        taxonomy_db=taxonomy_db,
                        lexicons=lexicons,
                        content_aware_figures=args.content_aware_figures,
                        run_config_fingerprints=run_config_fingerprints,
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

            if sys.platform == "darwin":
                _worker()
            else:
                mp_ctx = multiprocessing.get_context("fork")
                proc = mp_ctx.Process(target=_worker)
                proc.start()
                proc.join()
                if proc.exitcode != 0:
                    if proc.exitcode < 0:
                        signal_num = -proc.exitcode
                        logger.error(
                            "Document %s killed by signal %d (segfault?) — skipping",
                            pdf_hash, signal_num,
                        )
                        worker_failures.append({
                            "pdf_hash": pdf_hash,
                            "pdf_hash_full": pdf_hash_full,
                            "pdf_path": str(primary_pdf),
                            "exitcode": proc.exitcode,
                            "signal": signal_num,
                        })
                    else:
                        logger.error(
                            "Document %s worker exited with code %d — skipping",
                            pdf_hash, proc.exitcode,
                        )
                        worker_failures.append({
                            "pdf_hash": pdf_hash,
                            "pdf_hash_full": pdf_hash_full,
                            "pdf_path": str(primary_pdf),
                            "exitcode": proc.exitcode,
                            "signal": None,
                        })

            # The per-paper runner catches stage errors to preserve a useful
            # summary. A normal worker exit therefore does not imply success.
            summary_path = hash_dir / "summary.json"
            if summary_path.exists():
                summary = json.loads(summary_path.read_text())
                if (summary.get("processing_summary", {}).get("status") != "success"
                        and not any(f["pdf_hash"] == pdf_hash for f in worker_failures)):
                    worker_failures.append({"pdf_hash": pdf_hash,
                                            "pdf_path": str(primary_pdf),
                                            "exitcode": 1, "signal": None,
                                            "reason": "per-paper stage failure"})

    logger.info("Processing complete. Results saved to: %s", output_dir)
    logger.info("  Documents: %s", documents_dir)
    logger.info("  Vector DB: %s", vector_db_dir)

    if worker_failures:
        # Persist a structured failure record so downstream tooling
        # (orchestrator, `corpus status --report`) can see what dropped
        # out of this batch without grepping the log stream.
        from .version import __version__ as _pipeline_version
        failure_summary = {
            "pipeline_version": _pipeline_version,
            "n_failed": len(worker_failures),
            "failures": worker_failures,
        }
        failures_path = output_dir / "stage1_failures.json"
        try:
            failures_path.write_text(json.dumps(failure_summary, indent=2))
            logger.error(
                "%d document(s) crashed during extraction — see %s",
                len(worker_failures), failures_path,
            )
        except OSError as e:
            logger.error(
                "%d document(s) crashed during extraction; could not write "
                "%s (%s)",
                len(worker_failures), failures_path, e,
            )
        return 1

    # Silent-failure guard (#99): a run can "complete" with every document
    # producing zero chunks — e.g. when the extraction backend (docling)
    # regresses on a platform and returns empty documents without crashing.
    # On macOS the per-document worker runs inline, so such a document lands
    # a non-success summary.json but no subprocess failure; the run would
    # otherwise exit 0 and ship an empty served bundle. Treat a corpus-wide
    # zero-chunk result as a hard error, and warn on individual empties.
    attempted, total_chunks, zero_chunk_hashes = _audit_corpus_chunks(documents_dir)

    if attempted and total_chunks == 0:
        logger.error(
            "Extraction produced zero chunks across all %d processed "
            "document(s) — refusing to package an empty corpus. The "
            "extraction backend likely failed silently (e.g. docling "
            "returning empty output); inspect pipeline.log in each document "
            "dir under %s.",
            attempted, documents_dir,
        )
        return 1
    if zero_chunk_hashes:
        logger.warning(
            "%d of %d document(s) produced zero chunks: %s",
            len(zero_chunk_hashes), attempted, ", ".join(zero_chunk_hashes),
        )

    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
