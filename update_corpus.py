#!/usr/bin/env python3
"""Run the pipeline + post-pipeline scripts in dependency order (#32).

After adding new papers (or after a taxonomy / lexicon edit), the
corpus needs to walk through six tools in a specific order:

  1. process_corpus.py        — Stage 1: extraction + annotation
  2. embed_chunks.py          — Stage 2: embeddings (GPU)
  3. build_biblio_authority.py — bibliographic authority DB
  4. build_taxon_mentions.py   — taxon mention SQLite
  5. backfill_intext_citations.py — TEI body → intext_citations.json
  6. reconcile_corpus_to_biblio.py — merge ghost cited-references

Forgetting a step produces a silently inconsistent corpus. This
wrapper runs all six in order, with ``--resume`` semantics throughout
so the loop is "add papers / change inputs, run this command, done."

Each step is a subprocess; failure fails the whole run fast (no
attempt to recover by skipping later steps that might depend on the
failed one).

Pairs with #28 (per-stage resume in process_corpus.py) and #30
(idempotency audit of the post-pipeline scripts) — the underlying
guarantees that make a single re-run loop safe.

For Bouchet runs use ``slurm/batch_pipeline.sh`` instead; the GPU
embedding step (2) needs a different partition. ``update_corpus.py``
is the local-dev / small-corpus path.

Usage::

    python update_corpus.py <input_dir> <output_dir>
    python update_corpus.py <input_dir> <output_dir> --bib siphonophores.bib
    python update_corpus.py <input_dir> <output_dir> --skip-pipeline
    python update_corpus.py <input_dir> <output_dir> --from build_biblio
    python update_corpus.py <input_dir> <output_dir> --dry-run
"""
from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger("update_corpus")

REPO_ROOT = Path(__file__).resolve().parent


@dataclass
class Step:
    """One subprocess invocation in the dependency chain."""
    name: str
    script: str
    description: str

    def argv(self, args: argparse.Namespace) -> List[str]:
        """Build the argv for this step from CLI args."""
        cmd = [sys.executable, str(REPO_ROOT / self.script)]
        if self.name == "extract":
            cmd += [str(args.input_dir), str(args.output_dir)]
            if args.resume:
                cmd.append("--resume")
            if args.dry_run:
                cmd.append("--dry-run")
            if args.config:
                cmd += ["--config", str(args.config)]
            if args.bib:
                cmd += ["--bib", str(args.bib)]
            if args.anatomy_lexicon:
                cmd += ["--anatomy-lexicon", str(args.anatomy_lexicon)]
            if args.taxonomy_db:
                cmd += ["--taxonomy-db", str(args.taxonomy_db)]
            for spec in args.lexicon or []:
                cmd += ["--lexicon", spec]
            if args.no_taxa:
                cmd.append("--no-taxa")
            if args.no_grobid:
                cmd.append("--no-grobid")
            if args.grobid_url is not None:
                cmd += ["--grobid-url", args.grobid_url]
            if args.strict_network:
                cmd.append("--strict-network")
            if args.content_aware_figures:
                cmd.append("--content-aware-figures")
            if args.vision_backend:
                cmd += ["--vision-backend", args.vision_backend]
            if args.vision_model:
                cmd += ["--vision-model", args.vision_model]
            if args.refresh_vision:
                cmd.append("--refresh-vision")
        elif self.name == "embed":
            cmd += [str(args.output_dir)]
            if args.resume:
                cmd.append("--resume")
            if args.dry_run:
                cmd.append("--dry-run")
        elif self.name in ("build_biblio", "build_taxa", "backfill_intext"):
            cmd += [str(args.output_dir)]
            if args.dry_run:
                cmd.append("--dry-run")
        elif self.name == "reconcile":
            cmd += [str(args.output_dir)]
            if args.dry_run:
                cmd.append("--dry-run")
        return cmd


# Ordered. Names are stable identifiers for --from.
STEPS: List[Step] = [
    Step("extract",         "process_corpus.py",
         "Stage 1: extraction, OCR, docling, metadata, chunking, annotation"),
    Step("embed",           "embed_chunks.py",
         "Stage 2: BGE-M3 embeddings → LanceDB"),
    Step("build_biblio",    "build_biblio_authority.py",
         "Bibliographic authority DB"),
    Step("build_taxa",      "build_taxon_mentions.py",
         "Taxon mention rollup"),
    Step("backfill_intext", "backfill_intext_citations.py",
         "TEI body → intext_citations.json"),
    Step("reconcile",       "reconcile_corpus_to_biblio.py",
         "Merge ghost cited-references onto corpus papers"),
]

# Steps that depend on Stage 1 outputs only — usable with --skip-pipeline.
POST_PIPELINE_STEPS = {"build_biblio", "build_taxa", "backfill_intext", "reconcile"}


def _delete_stale_artifacts(args: argparse.Namespace) -> int:
    """Delete anatomy.json / taxa.json for papers whose stored
    input_fingerprint disagrees with the current input file (#33).

    The pipeline's per-stage resume (#28) regenerates the deleted
    artifact on the next run. In --dry-run mode, deletions are
    reported but not performed.
    """
    from corpus_status import detect_stale_fingerprints

    documents_dir = args.output_dir / "documents"
    if not documents_dir.is_dir():
        logger.warning(
            "--re-annotate-stale: no documents/ at %s, nothing to scan",
            args.output_dir,
        )
        return 0

    # Resolve current taxonomy / anatomy paths
    taxonomy_path = args.taxonomy_db
    if taxonomy_path is None:
        candidate = args.output_dir / "taxonomy.sqlite"
        if candidate.exists():
            taxonomy_path = candidate

    mode = args.re_annotate_stale
    want_taxonomy = mode in ("taxonomy", "any")
    want_anatomy = mode in ("anatomy", "any")

    if want_taxonomy and taxonomy_path is None:
        logger.error(
            "--re-annotate-stale=%s needs taxonomy.sqlite (not at the "
            "default <output>/taxonomy.sqlite; pass --taxonomy-db).", mode,
        )
        return 1
    if want_anatomy and args.anatomy_lexicon is None:
        logger.error(
            "--re-annotate-stale=%s needs --anatomy-lexicon to compare "
            "against.", mode,
        )
        return 1

    stale = detect_stale_fingerprints(
        documents_dir,
        taxonomy_path=taxonomy_path if want_taxonomy else None,
        anatomy_lexicon_path=args.anatomy_lexicon if want_anatomy else None,
    )

    n_taxa = len(stale.get("taxonomy_stale", []))
    n_anat = len(stale.get("anatomy_stale", []))
    logger.info(
        "Stale fingerprints: %d taxonomy_stale, %d anatomy_stale",
        n_taxa, n_anat,
    )

    if not (n_taxa or n_anat):
        logger.info("No stale annotations — nothing to delete.")
        return 0

    deleted = 0
    for h in stale.get("taxonomy_stale", []):
        path = documents_dir / h / "taxa.json"
        if args.dry_run:
            logger.info("would delete %s", path)
        else:
            try:
                path.unlink()
                deleted += 1
            except FileNotFoundError:
                pass
    for h in stale.get("anatomy_stale", []):
        path = documents_dir / h / "anatomy.json"
        if args.dry_run:
            logger.info("would delete %s", path)
        else:
            try:
                path.unlink()
                deleted += 1
            except FileNotFoundError:
                pass

    if args.dry_run:
        logger.info(
            "Dry-run: would delete %d artifact(s).", n_taxa + n_anat,
        )
    else:
        logger.info(
            "Deleted %d annotation artifact(s); per-stage resume will "
            "regenerate them.", deleted,
        )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input_dir", type=Path, help="Input directory of PDFs")
    parser.add_argument("output_dir", type=Path, help="Corpus output directory")
    parser.add_argument(
        "--bib", type=Path, default=None,
        help="Optional BibTeX for metadata override (passed to process_corpus.py).",
    )
    parser.add_argument(
        "--anatomy-lexicon", type=Path, default=None,
        help="Anatomy lexicon YAML (passed to process_corpus.py).",
    )
    parser.add_argument(
        "--taxonomy-db", type=Path, default=None,
        help="Taxonomy SQLite path (passed to process_corpus.py).",
    )
    parser.add_argument(
        "--lexicon", action="append", default=[], metavar="CATEGORY:PATH",
        help="Additional category-tagged lexicon YAML; repeatable "
             "(passed to process_corpus.py).",
    )
    parser.add_argument(
        "--no-taxa", action="store_true",
        help="Skip taxon + anatomy extraction (passed to process_corpus.py).",
    )
    parser.add_argument(
        "--no-grobid", action="store_true",
        help="Skip Grobid even if reachable (passed to process_corpus.py).",
    )
    parser.add_argument(
        "--grobid-url", default=None,
        help="Grobid service URL (passed to process_corpus.py).",
    )
    parser.add_argument(
        "--strict-network", action="store_true",
        help="Fail fast on the first transient external-service failure "
             "(passed to process_corpus.py). Recommended for release runs.",
    )
    parser.add_argument(
        "--content-aware-figures", action="store_true",
        help="Run Pass 3a OCR-driven panel ROI detection "
             "(passed to process_corpus.py).",
    )
    parser.add_argument(
        "--vision-backend", choices=["claude", "local"], default=None,
        help="Run Pass 3b vision-LLM-driven panel detection "
             "(passed to process_corpus.py).",
    )
    parser.add_argument(
        "--vision-model", default=None,
        help="Override the per-backend default vision model "
             "(passed to process_corpus.py).",
    )
    parser.add_argument(
        "--refresh-vision", action="store_true",
        help="With --resume + --vision-backend, re-run only Pass 3b on "
             "existing figures.json (passed to process_corpus.py).",
    )
    parser.add_argument(
        "--config", type=Path, default=None,
        help="Path to config.yaml (passed to process_corpus.py).",
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Pass --resume to the pipeline scripts (idempotent re-run).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Pass --dry-run to every step. No artifacts are written.",
    )
    parser.add_argument(
        "--skip-pipeline", action="store_true",
        help="Skip Stage 1 + Stage 2 (extract + embed). Only run the "
             "post-pipeline scripts. Use when no per-paper artifacts "
             "need to be (re)generated — e.g. you only want to rebuild "
             "the biblio/taxon DBs from existing taxa.json/anatomy.json. "
             "Incompatible with --re-annotate-stale, which requires the "
             "extract stage to regenerate the artifacts it deletes.",
    )
    parser.add_argument(
        "--from", dest="from_step", default=None,
        choices=[s.name for s in STEPS],
        help="Start from this step (skipping all earlier steps).",
    )
    parser.add_argument(
        "--re-annotate-stale", choices=["anatomy", "taxonomy", "any"],
        default=None,
        help="Before running the pipeline, scan for papers whose stored "
             "input_fingerprint disagrees with the current "
             "anatomy_lexicon / taxonomy.sqlite, then delete the relevant "
             "annotation artifact so per-stage resume regenerates it. "
             "Requires --resume and the matching input flag(s) "
             "(--anatomy-lexicon / --taxonomy-db).",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if not args.input_dir.exists():
        logger.error("Input directory %s does not exist", args.input_dir)
        return 1

    # Re-annotate-stale (#33): delete annotation artifacts whose stored
    # input_fingerprint disagrees with the current input. Pipeline's
    # per-stage resume (#28) will then regenerate them.
    if args.re_annotate_stale:
        if not args.resume:
            parser.error("--re-annotate-stale requires --resume")
        if args.skip_pipeline:
            parser.error(
                "--re-annotate-stale and --skip-pipeline are incompatible: "
                "the deleted artifacts are only regenerated by the extract "
                "stage that --skip-pipeline drops."
            )
        rc = _delete_stale_artifacts(args)
        if rc != 0:
            return rc

    # Determine the step subset
    if args.skip_pipeline and args.from_step:
        parser.error("--skip-pipeline and --from are mutually exclusive")

    selected = list(STEPS)
    if args.skip_pipeline:
        selected = [s for s in selected if s.name in POST_PIPELINE_STEPS]
    elif args.from_step is not None:
        idx = next(i for i, s in enumerate(STEPS) if s.name == args.from_step)
        selected = STEPS[idx:]

    logger.info("Running %d step(s):", len(selected))
    for s in selected:
        logger.info("  %-20s — %s", s.name, s.description)
    logger.info("Resume: %s. Dry-run: %s.",
                "on" if args.resume else "off",
                "on" if args.dry_run else "off")

    # Run each step in order. Real runs fail fast on any non-zero exit;
    # dry-run failures are downgraded to warnings because later steps'
    # dry-runs can't observe earlier steps' (un-written) outputs and
    # legitimately fail with "DB not found" or similar.
    overall_t0 = time.monotonic()
    n_warnings = 0
    for s in selected:
        cmd = s.argv(args)
        logger.info("═══ %s — %s ═══", s.name, " ".join(cmd[1:]))
        t0 = time.monotonic()
        try:
            subprocess.run(cmd, check=True, cwd=REPO_ROOT)
        except subprocess.CalledProcessError as e:
            if args.dry_run:
                logger.warning(
                    "%s failed in dry-run (exit %d) — likely a downstream "
                    "step that needs prior steps' real outputs. Continuing.",
                    s.name, e.returncode,
                )
                n_warnings += 1
                continue
            logger.error("%s failed (exit %d)", s.name, e.returncode)
            return e.returncode
        elapsed = time.monotonic() - t0
        logger.info("✓ %s completed in %.1fs", s.name, elapsed)

    elapsed_all = time.monotonic() - overall_t0
    if n_warnings:
        logger.info(
            "═══ All %d step(s) ran in %.1fs (%d dry-run warning(s)) ═══",
            len(selected), elapsed_all, n_warnings,
        )
    else:
        logger.info("═══ All %d step(s) succeeded in %.1fs ═══", len(selected), elapsed_all)
    return 0


if __name__ == "__main__":
    sys.exit(main())
