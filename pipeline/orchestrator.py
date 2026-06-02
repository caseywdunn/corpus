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
    """One subprocess invocation in the dependency chain.

    ``module`` is a dotted module path invoked via ``python -m``
    (the v0.3 packaging — root-level scripts are gone per #60).
    """
    name: str
    module: str
    description: str

    def argv(self, args: argparse.Namespace) -> List[str]:
        """Build the argv for this step from CLI args."""
        cmd = [sys.executable, "-m", self.module]
        if self.name == "ingest_taxonomy":
            cmd += [str(args.output_dir), "--source", args.taxonomy_source]
            if args.taxonomy_root_id is not None:
                cmd += ["--root-id", str(args.taxonomy_root_id)]
            if args.taxonomy_path is not None:
                cmd += ["--input", str(args.taxonomy_path)]
            if args.dry_run:
                cmd.append("--dry-run")
        elif self.name == "extract":
            cmd += [str(args.input_dir), str(args.output_dir)]
            if args.resume:
                cmd.append("--resume")
            if args.dry_run:
                cmd.append("--dry-run")
            if args.config:
                cmd += ["--config", str(args.config)]
            if args.bib:
                cmd += ["--bib", str(args.bib)]
            if args.taxonomy_db:
                cmd += ["--taxonomy-db", str(args.taxonomy_db)]
            if args.lexicon:
                cmd += ["--lexicon", str(args.lexicon)]
            if args.no_taxa:
                cmd.append("--no-taxa")
            if args.no_grobid:
                cmd.append("--no-grobid")
            if args.grobid_url is not None:
                cmd += ["--grobid-url", args.grobid_url]
            if args.strict_network:
                cmd.append("--strict-network")
            if args.figure_panels:
                cmd += ["--figure-panels", args.figure_panels]
            if args.vision_model:
                cmd += ["--vision-model", args.vision_model]
            if args.refresh_vision:
                cmd.append("--refresh-vision")
            cmd += _batch_flags(args)
        elif self.name == "vision":
            # HPC `--only vision` phase: re-run only Pass 3b (vision ROI
            # detection) on existing figures.json, on the GPU partition.
            # The extract phase may have run the CPU ocr floor, so force a
            # vision backend here. No Grobid / taxa work — figures only.
            cmd += [str(args.input_dir), str(args.output_dir), "--resume",
                    "--refresh-vision", "--no-grobid", "--no-taxa"]
            if args.dry_run:
                cmd.append("--dry-run")
            if args.config:
                cmd += ["--config", str(args.config)]
            fp = (args.figure_panels
                  if args.figure_panels in ("vision-local", "vision-claude")
                  else "vision-local")
            cmd += ["--figure-panels", fp]
            if args.vision_model:
                cmd += ["--vision-model", args.vision_model]
            cmd += _batch_flags(args)
        elif self.name == "embed":
            cmd += [str(args.output_dir)]
            if args.resume:
                cmd.append("--resume")
            if args.dry_run:
                cmd.append("--dry-run")
        elif self.name == "build_biblio":
            cmd += [str(args.output_dir)]
            if args.dry_run:
                cmd.append("--dry-run")
            if args.force_rebuild_biblio or args.force_rebuild:
                cmd.append("--rebuild")
            if args.enrich_bhl:
                cmd.append("--enrich-bhl")
            if args.config:
                # #51 — bib.authority reads `licensing.pd_cutoff_years`
                # from the per-corpuscle config for publishable derivation.
                cmd += ["--config", str(args.config)]
        elif self.name == "build_taxa":
            cmd += [str(args.output_dir)]
            if args.dry_run:
                cmd.append("--dry-run")
            if args.force_rebuild_taxon_mentions or args.force_rebuild:
                cmd.append("--rebuild")
        elif self.name == "backfill_intext":
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
    Step("ingest_taxonomy", "pipeline.taxonomy_ingest",
         "Auto-build taxonomy SQLite (#64; conditional on --taxonomy-source)"),
    Step("extract",         "pipeline.main",
         "Stage 1: extraction, OCR, docling, metadata, chunking, annotation"),
    Step("embed",           "pipeline.embed",
         "Stage 2: BGE-M3 embeddings → LanceDB"),
    Step("build_biblio",    "bib.authority",
         "Bibliographic authority DB"),
    Step("build_taxa",      "pipeline.taxon_mentions",
         "Taxon mention rollup"),
    Step("backfill_intext", "pipeline.intext_citations",
         "TEI body → intext_citations.json"),
    Step("reconcile",       "bib.reconcile",
         "Merge ghost cited-references onto corpus papers"),
]

# Steps that depend on Stage 1 outputs only — usable with --skip-pipeline.
POST_PIPELINE_STEPS = {"build_biblio", "build_taxa", "backfill_intext", "reconcile"}

# The GPU Pass-3b refresh, runnable as its own phase (`--only vision`).
# Not part of the default chain — when running the whole pipeline on one
# node, Pass 3b happens inline inside `extract` (figures.panel_detection
# = vision-*). It's a separate phase only so HPC can put it on a GPU
# partition after a CPU extract pass.
VISION_STEP = Step("vision", "pipeline.main",
                   "Pass 3b: vision ROI refresh (GPU; --only vision)")

# `--only <phase>` → the step name(s) that phase runs. ``post`` is the
# four cross-paper builds; ``bundle`` is handled in the CLI (served-bundle
# distill), not here.
ONLY_PHASES = {
    "extract": ["extract"],
    "vision": ["vision"],
    "embed": ["embed"],
    "post": ["build_biblio", "build_taxa", "backfill_intext", "reconcile"],
}


def select_steps(args: argparse.Namespace) -> List[Step]:
    """Resolve the ordered step list to run from the CLI args.

    Precedence: ``--only <phase>`` (single HPC phase, no chaining/
    taxonomy-autobuild) > ``--skip-pipeline`` (post-pipeline only) >
    ``--from <step>`` (suffix) > the full chain. The ingest_taxonomy
    autobuild step is dropped unless it's actually needed (#64).
    """
    if args.only is not None:
        names = ONLY_PHASES[args.only]
        catalog = {s.name: s for s in [*STEPS, VISION_STEP]}
        return [catalog[n] for n in names]

    if args.skip_pipeline:
        selected = [s for s in STEPS if s.name in POST_PIPELINE_STEPS]
    elif args.from_step is not None:
        idx = next(i for i, s in enumerate(STEPS) if s.name == args.from_step)
        selected = list(STEPS[idx:])
    else:
        selected = list(STEPS)

    # #64: drop ingest_taxonomy unless we actually need to (re)build the
    # snapshot. The other DBs handle the same via their own --rebuild.
    if _should_skip_taxonomy_ingest(args):
        selected = [s for s in selected if s.name != "ingest_taxonomy"]
    return selected


def _batch_flags(args: argparse.Namespace) -> List[str]:
    """`--batch-index/--batch-size` passthrough for HPC array slicing
    (#corpus-run-hpc). Both must be set; otherwise no slicing."""
    idx = getattr(args, "batch_index", None)
    size = getattr(args, "batch_size", None)
    if idx is not None and size is not None:
        return ["--batch-index", str(idx), "--batch-size", str(size)]
    return []


def _check_taxonomy_available(
    args: argparse.Namespace, selected: List[Step]
) -> Optional[str]:
    """Return an error string if taxonomy is configured but unavailable.

    Two failure modes:

    * ``ingest_taxonomy`` is *not* in the selected steps (e.g. ``--only
      extract`` on an HPC compute node) but ``taxonomy.sqlite`` is absent
      — the extract step will silently skip taxon annotation. Fail before
      any work starts so the operator knows to pre-build the taxonomy on a
      network-connected node first.

    * Source is ``dwca`` or ``dwc`` and the local archive/directory path
      does not exist — ``ingest_taxonomy`` would fail immediately; surface
      this sooner with a more descriptive message.
    """
    if not args.taxonomy_source:
        return None

    step_names = {s.name for s in selected}

    if "ingest_taxonomy" in step_names:
        # ingest_taxonomy will run — only check that local path exists for
        # file-based sources (worms reaches out to the network, which the
        # step itself will fail on if unavailable).
        if args.taxonomy_source in ("dwca", "dwc"):
            tx_path = getattr(args, "taxonomy_path", None)
            if tx_path is None or not Path(tx_path).exists():
                path_str = str(tx_path) if tx_path is not None else "(not set)"
                return (
                    f"taxonomy.source={args.taxonomy_source!r} requires a local "
                    f"archive but taxonomy.path={path_str!r} does not exist. "
                    "Provide a valid DwC-A archive or DwC directory."
                )
        return None  # ingest_taxonomy will handle the rest

    # ingest_taxonomy is NOT in the selected steps — taxonomy.sqlite must
    # already exist for taxon annotation to work.
    db_path = getattr(args, "taxonomy_db", None) or (args.output_dir / "taxonomy.sqlite")
    if Path(db_path).exists():
        return None

    if args.taxonomy_source == "worms":
        hint = (
            "WoRMS requires internet access. HPC compute nodes are "
            "network-restricted, so the taxonomy must be pre-built on a "
            "login node before submitting the array job. Options:\n"
            "  1. Run `corpus run --only ingest_taxonomy` on the login node.\n"
            "  2. Export a WoRMS subtree as a DwC-A snapshot, copy it to "
            "the project dir, and switch config to source: dwca."
        )
    else:  # dwca / dwc
        tx_path = getattr(args, "taxonomy_path", None)
        path_str = str(tx_path) if tx_path is not None else "(path not set)"
        hint = (
            f"Run `corpus run --only ingest_taxonomy` to build it from "
            f"{path_str}."
        )
    return (
        f"taxonomy.source={args.taxonomy_source!r} is configured but "
        f"{db_path} does not exist. {hint}"
    )


def _should_skip_taxonomy_ingest(args: argparse.Namespace) -> bool:
    """#64: skip ingest_taxonomy when (a) no source set, or (b) the
    taxonomy SQLite already exists and the operator hasn't asked for
    a rebuild via --force-rebuild or --force-rebuild-taxonomy.
    """
    if not args.taxonomy_source:
        return True
    db = args.output_dir / "taxonomy.sqlite"
    if db.exists() and not (args.force_rebuild or args.force_rebuild_taxonomy):
        return True
    return False


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
        "--taxonomy-db", type=Path, default=None,
        help="Taxonomy SQLite path (passed to process_corpus.py).",
    )
    parser.add_argument(
        "--lexicon", type=Path, default=None,
        help="Multi-category lexicon YAML (passed to process_corpus.py). "
             "Top-level keys are categories; see demo/lexicon.yaml.",
    )
    parser.add_argument(
        "--no-taxa", action="store_true",
        help="Skip the taxa_and_lexicon_extraction stage entirely "
             "(taxon mentions and every --lexicon category; passed to "
             "process_corpus.py).",
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
        "--figure-panels",
        choices=["ocr", "vision-local", "vision-claude", "off"],
        default="ocr",
        help="Panel-ROI detection mode (#102): ocr = Pass 3a OCR floor "
             "(default), vision-local / vision-claude = Pass 3b vision, "
             "off = none (passed to process_corpus.py).",
    )
    parser.add_argument(
        "--vision-model", default=None,
        help="Override the per-backend default vision model "
             "(passed to process_corpus.py).",
    )
    parser.add_argument(
        "--refresh-vision", action="store_true",
        help="With --resume + --figure-panels vision-*, re-run only Pass 3b "
             "on existing figures.json (passed to process_corpus.py).",
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
             "the biblio/taxon DBs from existing per-paper artifacts.",
    )
    parser.add_argument(
        "--from", dest="from_step", default=None,
        choices=[s.name for s in STEPS],
        help="Start from this step (skipping all earlier steps).",
    )
    parser.add_argument(
        "--only", choices=list(ONLY_PHASES), default=None,
        help="Run only one phase, instead of the full chain — for HPC "
             "where each phase maps to a separate SLURM job/partition: "
             "extract (CPU), vision (GPU Pass 3b refresh), embed (GPU), "
             "post (cross-paper DBs). Omit to run the whole pipeline.",
    )
    parser.add_argument(
        "--batch-index", type=int, default=None,
        help="0-based array-task index; with --batch-size, processes a "
             "deterministic slice of the sorted hash list (SLURM job "
             "arrays). Applies to the extract + vision phases.",
    )
    parser.add_argument(
        "--batch-size", type=int, default=None,
        help="PDFs per array task (pairs with --batch-index).",
    )
    # #64 — auto-build cross-paper databases.
    parser.add_argument(
        "--taxonomy-source", choices=["worms", "dwc", "dwca"], default=None,
        help="Source for the taxonomy SQLite (auto-built into "
             "<output_dir>/taxonomy.sqlite if missing). #64",
    )
    parser.add_argument(
        "--taxonomy-root-id", type=int, default=None,
        help="WoRMS AphiaID root for --taxonomy-source=worms.",
    )
    parser.add_argument(
        "--taxonomy-path", type=Path, default=None,
        help="Local DwC archive path for --taxonomy-source in "
             "{dwc, dwca}.",
    )
    parser.add_argument(
        "--enrich-bhl", action="store_true",
        help="Enrich pre-DOI references against the Biodiversity "
             "Heritage Library (slow, rate-limited). #64",
    )
    parser.add_argument(
        "--force-rebuild", action="store_true",
        help="Rebuild every cross-paper DB (taxonomy, biblio_authority, "
             "taxon_mentions) regardless of input-hash state. #64",
    )
    parser.add_argument("--force-rebuild-taxonomy", action="store_true")
    parser.add_argument("--force-rebuild-biblio", action="store_true")
    parser.add_argument("--force-rebuild-taxon-mentions", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if not args.input_dir.exists():
        logger.error("Input directory %s does not exist", args.input_dir)
        return 1

    # Lexicon / taxonomy edits no longer need a separate flag: the
    # per-paper pipeline_state.json records the input_fingerprint of
    # every stage, and --resume re-runs whichever stages' fingerprint
    # disagrees with the current input. Edit your inputs, run with
    # --resume, done.

    # Determine the step subset
    if sum(bool(x) for x in (args.skip_pipeline, args.from_step, args.only)) > 1:
        parser.error("--only, --skip-pipeline, and --from are mutually exclusive")
    if args.batch_index is not None and args.only not in ("extract", "vision"):
        parser.error("--batch-index/--batch-size apply only to "
                     "--only extract or --only vision")

    selected = select_steps(args)

    # Fail early if taxonomy is configured but unavailable for the selected
    # steps — avoids a silent "no taxon annotations" corpus on HPC where
    # compute nodes lack internet and --only skips ingest_taxonomy.
    taxonomy_err = _check_taxonomy_available(args, selected)
    if taxonomy_err:
        logger.error(taxonomy_err)
        return 1

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
