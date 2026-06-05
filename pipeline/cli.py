"""`corpus` CLI entry point.

Subcommand router for v0.3 (#60 + #61). Replaces the 13 root-level
scripts that v0.2 shipped. Subcommands:

  corpus init                  scaffold config.yaml in cwd (#59)
  corpus run                   full pipeline + post-pipeline + bundle
  corpus status                build-state rollup + report (#57)
  corpus serve                 MCP server
  corpus bib export|import     BibTeX round-trip
  corpus taxonomy ingest       build taxonomy.sqlite from WoRMS / DwC-A / DwC
  corpus taxonomy export       dump taxonomy.sqlite to a Darwin Core Archive
  corpus check                 environment + config pre-flight (stub, #62)
  corpus completion <shell>    shell completion script (#61)

Global options (#61):
  --config, -c PATH            override the per-corpuscle config.yaml
                               (env: CORPUS_CONFIG)
  --version, -V                package version + git SHA when sourced
                               from a development clone
  --cite [bibtex|plain]        print the citation (default: plain text;
                               source of truth: CITATION.cff)
  -v, -vv                      verbosity: INFO → DEBUG → all chatty
  -q, --quiet                  WARNING and above only

Defined exit codes (#61):
  0 success
  1 generic failure
  2 config error (missing/invalid config.yaml, schema mismatch,
                  unresolvable input path)
  3 precondition not met (`corpus check` fail, missing bundle,
                  port in use)

Each non-trivial verb dispatches via ``python -m <module>`` to the
moved private modules (pipeline.orchestrator, pipeline.status,
mcpsrv.main, etc.). The repo root has zero Python files in v0.3 —
the unified ``corpus`` binary is the only operator entry point.
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import textwrap
from importlib import resources
from pathlib import Path
from typing import List, Optional, Tuple

import yaml

from .config_schema import CorpuscleConfig, ValidationError, validate_config
from .console import console, print_status
from .version import __version__


_TEMPLATE_PACKAGE = "pipeline"
_TEMPLATE_NAME = "config.template.yaml"

# Defined exit codes (#61). Used by both the CLI and `corpus check`.
EXIT_OK = 0
EXIT_GENERIC = 1
EXIT_CONFIG_ERROR = 2
EXIT_PRECONDITION = 3


# ---------------------------------------------------------------------------
# Config-resolution helpers
# ---------------------------------------------------------------------------


def _resolve_config_path(flag_value: Optional[Path]) -> Optional[Path]:
    """Resolve which config.yaml to load.

    Precedence (#61): --config flag > $CORPUS_CONFIG env var > ./config.yaml.
    Returns None if no config file is found anywhere; the caller decides
    whether that's an error (corpus run errors; corpus init doesn't).
    """
    if flag_value is not None:
        return flag_value
    env = os.environ.get("CORPUS_CONFIG")
    if env:
        return Path(env)
    default = Path.cwd() / "config.yaml"
    if default.exists():
        return default
    return None


def _load_validated(config_path: Path) -> CorpuscleConfig:
    """Load config.yaml + validate against :class:`CorpuscleConfig`.

    Field-level errors print one line per failing key + bail with exit
    code 2 (config error, per #61 exit-code convention).
    """
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    try:
        return validate_config(raw)
    except ValidationError as e:
        print_status(f"config error in {config_path}", status="fail")
        for err in e.errors():
            loc = ".".join(str(x) for x in err["loc"])
            print(f"  {loc}: {err['msg']}", file=sys.stderr)
        sys.exit(2)


def _resolve_against(config_path: Path, value: Optional[Path]) -> Optional[Path]:
    """Resolve a path field in config.yaml against the config file's parent (#61).

    So `cd demo && corpus run` and `corpus --config demo/config.yaml run`
    from anywhere give identical results.
    """
    if value is None:
        return None
    if value.is_absolute():
        return value
    return (config_path.parent / value).resolve()


# ---------------------------------------------------------------------------
# `corpus init`  (#59)
# ---------------------------------------------------------------------------


def _cmd_init(args: argparse.Namespace) -> int:
    target = Path(args.path) if args.path else Path.cwd() / "config.yaml"
    if target.exists() and not args.force:
        print(
            f"corpus init: {target} already exists; pass --force to overwrite",
            file=sys.stderr,
        )
        return 2
    src = resources.files(_TEMPLATE_PACKAGE).joinpath(_TEMPLATE_NAME)
    target.parent.mkdir(parents=True, exist_ok=True)
    with resources.as_file(src) as src_path:
        shutil.copy(src_path, target)
    print(f"wrote {target}")
    print("Edit the input paths (input_pdfs, taxonomy, optional bib + lexicon)")
    print("then run `corpus run`.")
    return 0


# ---------------------------------------------------------------------------
# `corpus run`  (#60)
# ---------------------------------------------------------------------------


def _run_preconditions(
    cfg: CorpuscleConfig,
    config_path: Path,
    args: argparse.Namespace,
) -> int:
    """Fail-fast check before invoking the orchestrator.

    Two cliffs we surface here rather than letting the run silently
    degrade:
      - Grobid unreachable (when not disabled): every paper would get
        placeholder metadata, biblio_authority would be empty, hours of
        OCR/embed wasted on a useless bibliography pass.
      - input_pdfs path missing or empty: corpus run would happily
        complete in zero seconds with no work done.

    Returns ``EXIT_PRECONDITION`` on failure, ``EXIT_OK`` otherwise.
    Operators who know what they're doing pass ``--skip-checks``.
    """
    failures: List[str] = []

    # 1. input_pdfs reachable + non-empty
    if cfg.input_pdfs is None:
        failures.append(
            "config.yaml has no `input_pdfs:` set; nothing to process. "
            "Edit the config or pass --skip-checks to override."
        )
    else:
        input_dir = _resolve_against(config_path, cfg.input_pdfs)
        if not input_dir.exists():
            failures.append(
                f"input_pdfs path does not exist: {input_dir}. "
                f"Check the path in {config_path} or pass --skip-checks."
            )
        elif input_dir.is_dir():
            n_pdfs = sum(1 for _ in input_dir.rglob("*.pdf"))
            if n_pdfs == 0:
                failures.append(
                    f"input_pdfs path has zero PDFs under {input_dir}. "
                    "Drop PDFs in or pass --skip-checks to run anyway."
                )

    # 2. Grobid reachable (unless disabled). Skipping Grobid yields
    # working OCR + embedding + chunking but a bibliography graph that
    # can't link cited references to corpus papers. The deliberate-skip
    # path is `grobid.disable: true` in config; that opts out of the
    # check entirely.
    if not cfg.grobid.disable:
        ok, detail = _ping_grobid(cfg.grobid.url)
        if not ok:
            failures.append(
                f"Grobid not reachable at {cfg.grobid.url} ({detail}). "
                "Start it with `docker compose up -d grobid` (it persists "
                "across runs). To deliberately run without metadata "
                "extraction, set `grobid: {disable: true}` in config.yaml. "
                "To bypass this check just for this run, pass --skip-checks."
            )

    if failures:
        for f in failures:
            print_status(f, status="fail")
        print()
        print_status(
            f"{len(failures)} precondition(s) failed before the orchestrator "
            "started. No compute was spent on this run. See "
            "`corpus check` for the full pre-flight surface.",
            status="fail",
        )
        return EXIT_PRECONDITION
    return EXIT_OK


def _invalidate_flagged(
    cfg: CorpuscleConfig,
    config_path: Path,
    gate: str,
) -> None:
    """Delete pipeline_state.json for every paper carrying ``gate`` (#54).

    Per-stage resume re-runs whichever stages no longer have a state
    record, so removing the file forces a full re-extraction. We don't
    touch any other artifact — figures.json, summary.json, etc. stay so
    the operator can compare before/after.
    """
    output_dir = _resolve_against(config_path, cfg.output_dir)
    documents_dir = output_dir / "documents"
    if not documents_dir.is_dir():
        return
    from .status import aggregate
    rollup = aggregate(documents_dir)
    hashes = rollup["papers_by_gate"].get(gate, set())
    if not hashes:
        print_status(
            f"--re-process-flagged: no papers carry gate {gate!r}; "
            f"nothing to invalidate.",
            status="warn",
        )
        return
    n = 0
    for h in hashes:
        state_path = documents_dir / h / "pipeline_state.json"
        if state_path.exists():
            state_path.unlink()
            n += 1
    print_status(
        f"--re-process-flagged {gate!r}: invalidated pipeline_state.json "
        f"for {n} paper(s); they will re-extract on this run.",
        status="info",
    )


def _prune_orphans(
    cfg: CorpuscleConfig,
    config_path: Path,
    args: argparse.Namespace,
) -> int:
    """Run #66's orphan cleanup before dispatching to the orchestrator."""
    if cfg.input_pdfs is None:
        return EXIT_OK  # nothing to compare against; let extract handle it
    input_dir = _resolve_against(config_path, cfg.input_pdfs)
    output_dir = _resolve_against(config_path, cfg.output_dir)
    if not (output_dir / "documents").is_dir():
        return EXIT_OK  # fresh corpuscle — nothing to prune

    from .io import audit_orphans, prune_orphans

    if args.no_prune:
        # Read-only audit (matches #31 behavior).
        n = audit_orphans(input_dir, output_dir)
        print_status(f"--no-prune: audit reported {n} orphan(s)", status="info")
        return EXIT_OK

    try:
        result = prune_orphans(
            input_dir,
            output_dir,
            dry_run=args.dry_run,
            force=args.force_prune,
        )
    except RuntimeError as e:
        print_status(str(e), status="fail")
        return EXIT_PRECONDITION

    if args.dry_run:
        if result["doc_pruned"]:
            print_status(
                f"dry-run: would prune {result['doc_pruned']} of "
                f"{result['doc_total']} hash dir(s) + LanceDB rows",
                status="warn",
            )
        else:
            print_status("dry-run: no orphans to prune", status="ok")
    else:
        if result["doc_pruned"]:
            print_status(
                f"pruned {result['doc_pruned']} of {result['doc_total']} "
                f"orphan hash dir(s) + {result['vec_pruned']} LanceDB row(s)",
                status="warn",
            )
    return EXIT_OK


def _build_orchestrator_argv(
    cfg: CorpuscleConfig,
    config_path: Path,
    args: argparse.Namespace,
) -> List[str]:
    """Translate validated config + CLI flags → orchestrator argv."""
    if cfg.input_pdfs is None:
        print_status(
            "config error: `input_pdfs` is required in config.yaml for `corpus run`",
            status="fail",
        )
        sys.exit(2)
    input_dir = _resolve_against(config_path, cfg.input_pdfs)
    output_dir = _resolve_against(config_path, cfg.output_dir)

    sub_argv: List[str] = [str(input_dir), str(output_dir)]

    # Forward the per-corpuscle config so pipeline.main + downstream
    # steps read system tuning (ocr, chunking, stage_timeouts, ...)
    # from the same file the corpus router validated.
    sub_argv += ["--config", str(config_path)]

    # Implicit resume per #60. Operators opt out with --force-rebuild.
    if not args.force_rebuild:
        sub_argv.append("--resume")
    if args.dry_run:
        sub_argv.append("--dry-run")

    bib = _resolve_against(config_path, cfg.bib)
    if bib is not None:
        sub_argv += ["--bib", str(bib)]
    lexicon = _resolve_against(config_path, cfg.lexicon)
    if lexicon is not None:
        sub_argv += ["--lexicon", str(lexicon)]

    if cfg.grobid.disable:
        sub_argv.append("--no-grobid")
    if cfg.grobid.url:
        sub_argv += ["--grobid-url", cfg.grobid.url]

    # Figure panel-ROI detection (#102). Effective mode precedence:
    # explicit `--figure-panels` > `--no-vision` (deprecated alias for
    # `ocr`) > the corpuscle's `figures.panel_detection` config key.
    if args.figure_panels is not None:
        panel_mode = args.figure_panels
    elif args.no_vision:
        print_status(
            "--no-vision is a deprecated alias for `--figure-panels ocr` "
            "(CPU OCR-panel floor)",
            status="info",
        )
        panel_mode = "ocr"
    else:
        panel_mode = cfg.figures.panel_detection
    # Capability detection (#65) downgrades an unusable vision backend to
    # the OCR floor with a one-line nudge rather than hard-failing Pass 3b
    # deep into the run.
    if panel_mode in ("vision-local", "vision-claude"):
        skip_reason = _vision_skip_reason(panel_mode)
        if skip_reason is not None:
            print_status(
                f"vision panel pass downgraded to the OCR floor: {skip_reason}",
                status="warn",
            )
            panel_mode = "ocr"
    sub_argv += ["--figure-panels", panel_mode]
    if panel_mode in ("vision-local", "vision-claude") and cfg.figures.model:
        sub_argv += ["--vision-model", cfg.figures.model]

    # #64: auto-build cross-paper databases.
    if cfg.taxonomy.source is not None:
        sub_argv += ["--taxonomy-source", cfg.taxonomy.source]
        if cfg.taxonomy.root_id is not None:
            sub_argv += ["--taxonomy-root-id", str(cfg.taxonomy.root_id)]
        if cfg.taxonomy.path is not None:
            tx_path = _resolve_against(config_path, cfg.taxonomy.path)
            sub_argv += ["--taxonomy-path", str(tx_path)]
    if args.enrich_bhl or cfg.bibliography.enrich_bhl:
        sub_argv.append("--enrich-bhl")
    if args.force_rebuild:
        sub_argv.append("--force-rebuild")
    if args.force_rebuild_taxonomy:
        sub_argv.append("--force-rebuild-taxonomy")
    if args.force_rebuild_biblio:
        sub_argv.append("--force-rebuild-biblio")
    if args.force_rebuild_taxon_mentions:
        sub_argv.append("--force-rebuild-taxon-mentions")
    # HPC single-phase + array slicing (#corpus-run-hpc). `bundle` is a
    # CLI-level phase (handled before the orchestrator is invoked), so
    # it's never forwarded here.
    only = getattr(args, "only", None)
    if only is not None and only != "bundle":
        sub_argv += ["--only", only]
    if getattr(args, "batch_index", None) is not None:
        sub_argv += ["--batch-index", str(args.batch_index)]
    if getattr(args, "batch_size", None) is not None:
        sub_argv += ["--batch-size", str(args.batch_size)]
    return sub_argv


def _cmd_run(args: argparse.Namespace) -> int:
    config_path = _resolve_config_path(args.config)
    if config_path is None:
        print(
            "corpus run: no config.yaml in cwd; run `corpus init` to scaffold one, "
            "or pass --config PATH",
            file=sys.stderr,
        )
        return EXIT_CONFIG_ERROR
    if not config_path.exists():
        print_status(f"config not found: {config_path}", status="fail")
        return EXIT_CONFIG_ERROR
    cfg = _load_validated(config_path)
    # $GROBID_URL overrides the configured Grobid address (#138). On HPC,
    # Grobid runs on a dynamically-allocated compute node whose hostname
    # isn't known until submit time, so the static config.yaml value can't
    # name it; slurm/batch_pipeline.sh discovers the node and exports
    # GROBID_URL. Matches the old `pipeline.main` default-from-$GROBID_URL
    # behavior. Applied before the precondition ping and the argv build so
    # both see the override.
    grobid_url_env = os.environ.get("GROBID_URL")
    if grobid_url_env:
        cfg.grobid.url = grobid_url_env
    output_dir = _resolve_against(config_path, cfg.output_dir)
    # Single-phase HPC selector (#corpus-run-hpc): None = full run on one
    # node; extract/vision/embed/post run one orchestrator phase; bundle
    # distills the served bundle only. Each maps to its own SLURM job.
    only = getattr(args, "only", None)
    if (args.batch_index is None) != (args.batch_size is None):
        print_status("--batch-index and --batch-size must be given together",
                     status="fail")
        return EXIT_CONFIG_ERROR
    if args.batch_index is not None and only not in ("extract", "vision"):
        print_status(
            "--batch-index/--batch-size apply only to --only extract|vision",
            status="fail")
        return EXIT_CONFIG_ERROR

    # --only bundle: (re)distill the served bundle from existing
    # output_dir artifacts. No preconditions / prune / pipeline.
    if only == "bundle":
        if args.dry_run:
            print_status(
                f"dry-run: would distill served bundle into {output_dir / '_serve'}",
                status="info")
            return EXIT_OK
        if args.no_bundle:
            print_status("--only bundle with --no-bundle is a no-op", status="warn")
            return EXIT_OK
        return _distill_bundle(output_dir)

    # Preconditions + orphan prune + flag-invalidation are extract-phase
    # concerns — skip them for the GPU (vision/embed) and post phases,
    # which scan no input and need no Grobid.
    if only in (None, "extract"):
        # Fail-fast preconditions: catch the two ergonomic cliffs that
        # would otherwise burn hours (Grobid unreachable → placeholder
        # metadata; input_pdfs missing → silent zero-paper run). The
        # vision/ANTHROPIC_API_KEY paths are warn-and-skip (#65).
        if not args.skip_checks:
            rc = _run_preconditions(cfg, config_path, args)
            if rc != EXIT_OK:
                return rc
        # #66: orphan cleanup before extract so later stages don't consume
        # orphan hashes. Default prune; --no-prune = audit; --force-prune
        # bypasses the safety rail.
        rc = _prune_orphans(cfg, config_path, args)
        if rc != EXIT_OK:
            return rc
        # #54: --re-process-flagged invalidates pipeline_state.json for
        # papers carrying the named quality_flag so resume re-extracts them.
        if args.re_process_flagged:
            _invalidate_flagged(cfg, config_path, args.re_process_flagged)

    sub_argv = _build_orchestrator_argv(cfg, config_path, args)
    cmd = [sys.executable, "-m", "pipeline.orchestrator", *sub_argv]
    print_status(f"corpus run → {' '.join(cmd)}", status="info")
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        return rc

    # Shepherd <corpuscle>/instructions.md into <output_dir>/ before
    # bundle distill. The README documents instructions.md as a
    # corpuscle-root file (next to config.yaml — the demo follows
    # that pattern), but every downstream reader (mcpsrv.main's
    # InitializeResult.instructions, mcpsrv.bundle's served-bundle
    # whitelist) looks for it under <output_dir>/. Without this
    # forwarding step the user's instructions never reach the served
    # bundle — exactly the warning the smoke test surfaced:
    # `Expected top-level file instructions.md missing from <output_dir>`.
    # We copy on mtime change, so re-runs are idempotent and edits
    # propagate.
    if not args.dry_run and config_path is not None:
        src_instructions = config_path.parent / "instructions.md"
        if src_instructions.exists():
            dst_instructions = output_dir / "instructions.md"
            output_dir.mkdir(parents=True, exist_ok=True)
            needs_copy = (
                not dst_instructions.exists()
                or src_instructions.stat().st_mtime > dst_instructions.stat().st_mtime
            )
            if needs_copy:
                shutil.copy2(src_instructions, dst_instructions)

    # #60 — bundle distillation in line, on a full run only. The served
    # bundle lands at `<output_dir>/_serve/`. Sub-phases (extract / vision
    # / embed / post) skip it — the HPC flow distills once via a final
    # `--only bundle` job; `--only bundle` is handled at the top.
    if only is None:
        if args.no_bundle:
            print_status(
                "skipping served-bundle distillation (--no-bundle)",
                status="info",
            )
        elif args.dry_run:
            print_status(
                f"dry-run: would distill served bundle into {output_dir / '_serve'}",
                status="info",
            )
        else:
            rc = _distill_bundle(output_dir)
            if rc != 0:
                print_status(
                    f"served-bundle distillation failed (exit {rc}); the build "
                    "bundle is intact but no served bundle was produced. Re-run "
                    "`python -m mcpsrv.bundle <output_dir> <serve_dir> "
                    "--version vX.Y.Z` to retry, or pass --no-bundle to skip.",
                    status="fail",
                )
                # Propagate: _serve/ is the deployable artifact (DEPLOY.md);
                # a failed distill must not stamp the run as successful.
                return rc

    if args.dry_run:
        print_status(
            "dry-run complete. No artifacts written — re-run without "
            "`--dry-run` to execute the plan above.",
            status="ok",
        )
    elif only not in (None, "post"):
        # A sub-phase (extract / vision / embed) — its job is done; the
        # run-log + report come from the full run or the post phase.
        print_status(f"phase '{only}' complete.", status="ok")
    else:
        run_log = _write_run_log(output_dir, cfg=cfg, argv=sys.argv)
        if run_log is None:
            print_status(
                "run complete. Try `corpus status --report` and `corpus serve` next.",
                status="ok",
            )
        else:
            print_status(
                f"run complete. Report → {run_log} (same as `corpus status "
                f"--report`). `corpus serve` to serve the bundle.",
                status="ok",
            )
    return EXIT_OK


def _cmd_debug_pdf(args: argparse.Namespace) -> int:
    """Run the per-paper pipeline on a single PDF with verbose stage
    tracing and no bundle / cross-paper steps (#92).

    A focused troubleshooting runner for hard documents. It drives
    exactly the per-PDF stages ``run_pdf_processing_pipeline`` runs
    (huge-doc gate, scan detection, PDF prep, docling extraction,
    chunking, figures — plus the bibliography stage when ``--grobid-url``
    is given), writes the intermediate artifacts under a debug directory
    you can inspect, and prints a per-stage ok/fail + timing table so a
    failure is immediately attributable to a stage. It deliberately runs
    inline (no fork) and with ``resume=False`` so every stage executes
    every time. The taxa/lexicon stages — which need corpus-level
    ``taxonomy.sqlite`` / lexicon inputs — are out of scope; use a full
    `corpus run` for those.
    """
    import tempfile

    from .grobid_client import GrobidClient
    from .io import calculate_pdf_hash, create_summary_json, short_hash
    from .log import per_pdf_file_log, setup_root_logging
    from .runner import run_pdf_processing_pipeline

    pdf_path: Path = args.pdf
    if not pdf_path.is_file():
        print_status(f"not a PDF file: {pdf_path}", status="fail")
        return EXIT_CONFIG_ERROR

    # This is a debug command — default to verbose so the stage trace is
    # visible without requiring -v. -q/-vv still take precedence.
    setup_root_logging(logging.DEBUG if args.verbose >= 2 else logging.INFO)
    log = logging.getLogger("pipeline.debug_pdf")

    debug_root = (args.output_dir or (Path.cwd() / "debug-pdf")).resolve()
    pdf_hash_full = calculate_pdf_hash(pdf_path)
    hash_dir = debug_root / short_hash(pdf_hash_full)
    hash_dir.mkdir(parents=True, exist_ok=True)
    print_status(
        f"debug-pdf: {pdf_path.name} ({short_hash(pdf_hash_full)}) → {hash_dir}",
        status="info",
    )

    grobid_client: Optional[GrobidClient] = None
    if args.grobid_url:
        grobid_client = GrobidClient(base_url=args.grobid_url)
        print_status(f"bibliography stage enabled via Grobid at {args.grobid_url}",
                     status="info")
    else:
        print_status(
            "bibliography stage skipped (pass --grobid-url to enable); "
            "taxa/lexicon stages are out of scope for debug-pdf",
            status="info",
        )

    with tempfile.TemporaryDirectory() as temp_dir:
        with per_pdf_file_log(hash_dir) as log_path:
            log.info("pipeline.log: %s", log_path)
            log.info("pdf_hash_full: %s", pdf_hash_full)
            summary = run_pdf_processing_pipeline(
                pdf_path, hash_dir, Path(temp_dir),
                grobid_client=grobid_client,
                taxonomy_db=None,
                lexicons=None,
                bib_index=None,
                resume=False,
            )
        create_summary_json(
            pdf_hash_full, [pdf_path], pdf_path.parent, hash_dir, summary,
        )

    _print_debug_pdf_trace(summary, hash_dir)
    return EXIT_OK if summary.get("status") == "success" else EXIT_GENERIC


def _print_debug_pdf_trace(summary: dict, hash_dir: Path) -> None:
    """Print a per-stage ok/fail + timing table for `corpus debug-pdf`."""
    timings = summary.get("stage_timings") or []
    failures = {f.get("stage"): f for f in (summary.get("stage_failures") or [])}
    skipped = set(summary.get("skipped_stages") or [])

    print()
    print(f"Stage trace ({len(timings)} stages run):")
    for t in timings:
        stage = t.get("stage", "?")
        ok = t.get("ok")
        dur = t.get("duration_s")
        mark = "✓" if ok else "✗"
        dur_s = f"{dur:6.2f}s" if isinstance(dur, (int, float)) else "     ?"
        line = f"  {mark} {stage:<28s} {dur_s}"
        if not ok and stage in failures:
            line += f"   {failures[stage].get('error', '')}"
        print(line)
    if skipped:
        print(f"  (skipped, artifact already present: {', '.join(sorted(skipped))})")

    status = summary.get("status", "unknown")
    n_files = len(summary.get("files_created") or [])
    print()
    print_status(
        f"status={status}; {n_files} artifact(s) written under {hash_dir}",
        status="ok" if status == "success" else "fail",
    )
    print(f"Full log: {hash_dir / 'pipeline.log'}")


def _write_run_log(
    output_dir: Path,
    cfg: Optional[CorpuscleConfig] = None,
    argv: Optional[List[str]] = None,
) -> Optional[Path]:
    """Write the end-of-run report (#57, #90).

    Captures a self-contained, reproducible record of the run: the
    ``corpus status --report`` pass-rate + artifact rollup, plus (#90)
    the invocation's argv, the *resolved* config, the dependency-stack
    versions that produced the output, and aggregate stage success /
    failure counts.

    Two copies are written:

    * ``<output_dir>/runs/<timestamp>/run.log`` — the per-invocation
      archive, so a sequence of runs leaves an auditable trail rather
      than overwriting one another.
    * ``<output_dir>/run.log`` — the "latest" copy, preserved at the
      top level for #57 back-compat (``corpus status``, the served
      bundle, and existing tooling read it there).

    Returns the top-level path, or ``None`` if there's nothing to
    report (no ``documents/`` tree yet — e.g. the orchestrator bailed
    before extract ran).
    """
    from datetime import datetime
    from .status import aggregate, render_artifacts, render_text

    documents_dir = output_dir / "documents"
    if not documents_dir.is_dir():
        return None
    rollup = aggregate(documents_dir)

    now = datetime.now()
    header_lines = [
        "# corpus run report",
        f"# generated at {now.isoformat(timespec='seconds')}",
        f"# corpus version {__version__}",
    ]
    sha = _git_sha_or_none()
    if sha:
        header_lines.append(f"# git {sha}")
    header_lines.append(f"# output_dir {output_dir.resolve()}")
    if argv is not None:
        header_lines.append(f"# argv {' '.join(argv)}")
    header = "\n".join(header_lines) + "\n\n"

    sections = [render_text(rollup), _render_run_summary(rollup)]
    if cfg is not None:
        sections.append(_render_resolved_config(cfg))
    sections.append(_render_dependency_versions())
    sections.append(render_artifacts(output_dir))
    content = header + "\n\n".join(s.rstrip() for s in sections) + "\n"

    # Per-invocation archive (#90) — second-resolution timestamp; full
    # pipeline runs take minutes, so collisions don't occur in practice.
    runs_dir = output_dir / "runs" / now.strftime("%Y%m%dT%H%M%S")
    runs_dir.mkdir(parents=True, exist_ok=True)
    (runs_dir / "run.log").write_text(content, encoding="utf-8")

    # Top-level "latest" copy (#57 back-compat).
    log_path = output_dir / "run.log"
    log_path.write_text(content, encoding="utf-8")
    return log_path


def _render_run_summary(rollup: dict) -> str:
    """Aggregate stage success/failure counts across all stages (#90)."""
    stages = rollup.get("stages", {}) or {}
    total_ok = sum(r.get("ok", 0) for r in stages.values())
    total_fail = sum(r.get("fail", 0) for r in stages.values())
    return "\n".join([
        "Run summary:",
        f"  documents:          {rollup.get('total_documents', 0)}",
        f"  stage runs ok:      {total_ok}",
        f"  stage runs failed:  {total_fail}",
    ])


def _render_resolved_config(cfg: CorpuscleConfig) -> str:
    """Serialize the validated, resolved config (#90). ``mode='json'``
    coerces Path / enum fields to strings so the dump is plain YAML."""
    try:
        dumped = yaml.safe_dump(
            cfg.model_dump(mode="json"), sort_keys=False, default_flow_style=False,
        )
    except Exception as e:  # never let a serialization quirk fail the run
        dumped = f"(could not serialize config: {e})"
    return "Resolved config:\n" + textwrap.indent(dumped.rstrip(), "  ")


def _render_dependency_versions() -> str:
    """Versions of the stack that produced this output (#90). Centered
    on the ML pins (#98) whose silent drift broke a v0.5 build, plus the
    serving deps, so a run.log pins exactly what to reproduce against."""
    from importlib.metadata import PackageNotFoundError, version

    pkgs = [
        "docling", "torch", "transformers", "sentence-transformers",
        "lancedb", "pyarrow", "mcp", "anthropic", "pydantic",
    ]
    lines = ["Dependency versions:", f"  {'python':<22s} {sys.version.split()[0]}"]
    for p in pkgs:
        try:
            v = version(p)
        except PackageNotFoundError:
            v = "(not installed)"
        lines.append(f"  {p:<22s} {v}")
    return "\n".join(lines)


def _distill_bundle(output_dir: Path) -> int:
    """Invoke mcpsrv.bundle to produce <output_dir>/_serve/ (#60)."""
    serve_dir = output_dir / "_serve"
    cmd = [
        sys.executable, "-m", "mcpsrv.bundle",
        str(output_dir), str(serve_dir),
        "--version", __version__,
    ]
    print_status(
        f"distilling served bundle → {serve_dir} (v{__version__})",
        status="info",
    )
    return subprocess.run(cmd).returncode


# ---------------------------------------------------------------------------
# `corpus status`  /  `corpus serve`  /  `corpus bib …`  (#60 passthroughs)
# ---------------------------------------------------------------------------


def _passthrough(module: str, extra_argv: List[str]) -> int:
    """Spawn ``python -m <module> <extra_argv...>`` and return its exit code."""
    cmd = [sys.executable, "-m", module, *extra_argv]
    return subprocess.run(cmd).returncode


def _cmd_status(args: argparse.Namespace) -> int:
    output_dir = _resolve_output_dir(args)
    return _passthrough("pipeline.status", [str(output_dir), *args.passthrough])


def _resolve_output_dir(args: argparse.Namespace) -> Path:
    """Pick output_dir for a verb that needs it: --output-dir flag wins; else
    pull from the resolved config.yaml's `output_dir`.
    """
    if args.output_dir is not None:
        return args.output_dir
    config_path = _resolve_config_path(args.config)
    if config_path is None or not config_path.exists():
        print_status(
            "no --output-dir and no config.yaml found; pass --output-dir PATH "
            "or `cd <corpuscle>`",
            status="fail",
        )
        sys.exit(2)
    cfg = _load_validated(config_path)
    return _resolve_against(config_path, cfg.output_dir)


def _cmd_serve(args: argparse.Namespace) -> int:
    output_dir = _resolve_output_dir(args)
    return _passthrough("mcpsrv.main", [str(output_dir), *args.passthrough])


def _cmd_bib_export(args: argparse.Namespace) -> int:
    output_dir = _resolve_output_dir(args)
    return _passthrough("bib.export", [str(output_dir), *args.passthrough])


def _cmd_bib_import(args: argparse.Namespace) -> int:
    output_dir = _resolve_output_dir(args)
    return _passthrough("bib.importer", [str(output_dir), *args.passthrough])


def _cmd_taxonomy_export(args: argparse.Namespace) -> int:
    output_dir = _resolve_output_dir(args)
    return _passthrough("pipeline.taxonomy_export",
                        [str(output_dir), *args.passthrough])


def _cmd_taxonomy_ingest(args: argparse.Namespace) -> int:
    output_dir = _resolve_output_dir(args)
    return _passthrough("pipeline.taxonomy_ingest",
                        [str(output_dir), *args.passthrough])


# ---------------------------------------------------------------------------
# `corpus check`  (#62 stub)  /  `corpus completion`  (#61 stub)
# ---------------------------------------------------------------------------


def _cmd_check(args: argparse.Namespace) -> int:
    """Environment + config pre-flight (#62).

    Distinct from ``corpus run --dry-run`` (which prints the *plan*)
    and ``corpus status`` (postmortem against a built corpuscle):
    answers "can the next run actually succeed on this host?" without
    consulting the output tree.

    Exit codes (per #61): 0 on green, 2 on config-error, 3 on
    environment-precondition failure.
    """
    # `corpus check` exists to print its findings. The default verbosity
    # is WARNING (see _setup_logging), which would silently drop every
    # `ok` line via print_status's level gate — leaving the operator
    # with just the warnings and a green exit code, unable to tell
    # what was actually checked. Force every status line in this
    # function to print regardless of root log level. `corpus -q check`
    # still works (the gate is bypassed, not the log threshold), since
    # check's whole purpose is the report.
    from functools import partial
    pstatus = partial(print_status, force=True)

    config_path = _resolve_config_path(args.config)
    failures: List[str] = []  # precondition (exit 3)
    config_failures: List[str] = []  # config (exit 2)

    # 1. config.yaml resolution + schema validation
    if config_path is None:
        pstatus(
            "no config.yaml in cwd; pass --config PATH or run `corpus init`",
            status="fail",
        )
        return EXIT_CONFIG_ERROR
    if not config_path.exists():
        pstatus(f"config not found: {config_path}", status="fail")
        return EXIT_CONFIG_ERROR
    try:
        raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        cfg = validate_config(raw)
        pstatus(f"config.yaml valid ({config_path})", status="ok")
    except ValidationError as e:
        pstatus(f"config error in {config_path}", status="fail")
        for err in e.errors():
            loc = ".".join(str(x) for x in err["loc"])
            print(f"  {loc}: {err['msg']}", file=sys.stderr)
        return EXIT_CONFIG_ERROR

    # 2. GPU availability
    gpu = _detect_accelerator()
    if gpu == "cuda":
        pstatus("GPU: CUDA available", status="ok")
    elif gpu == "mps":
        pstatus("GPU: MPS available (Apple Silicon)", status="ok")
    else:
        pstatus("GPU: none (CPU-only; embedding pass will be slow)", status="warn")

    # 3. Figure panel-ROI detection (#102). Under #65 an unusable vision
    # backend is a warn (the run downgrades to the OCR floor with the
    # same message), not a hard precondition failure.
    panel_mode = cfg.figures.panel_detection
    if panel_mode == "off":
        pstatus("Figure panels: disabled in config (panel_detection: off)",
                status="warn")
    elif panel_mode == "ocr":
        pstatus("Figure panels: OCR floor (panel_detection: ocr; CPU-only)",
                status="ok")
    else:
        vision_skip = _vision_skip_reason(panel_mode)
        if vision_skip is None:
            pstatus(f"Figure panels: `{panel_mode}` ready", status="ok")
        else:
            pstatus(f"Figure panels: would downgrade to OCR floor — {vision_skip}",
                    status="warn")

    # 4. Grobid reachability
    if cfg.grobid.disable:
        pstatus(f"Grobid: disabled in config (header metadata will use --bib only)", status="warn")
    else:
        ok, detail = _ping_grobid(cfg.grobid.url)
        if ok:
            pstatus(f"Grobid: reachable at {cfg.grobid.url}", status="ok")
        else:
            failures.append(
                f"Grobid not reachable at {cfg.grobid.url} ({detail}). "
                "Start it with: docker compose up -d grobid"
            )
            pstatus(f"Grobid: unreachable at {cfg.grobid.url} ({detail})", status="fail")

    # 5. Output disk space (warn if < 5 GB free)
    output_dir = _resolve_against(config_path, cfg.output_dir)
    free_gb = _free_disk_gb(output_dir)
    if free_gb is None:
        pstatus(f"output_dir disk: cannot stat {output_dir}", status="warn")
    elif free_gb < 5:
        pstatus(f"output_dir disk: {free_gb:.1f} GB free (low — recommend ≥ 20 GB for a real corpus)", status="warn")
    else:
        pstatus(f"output_dir disk: {free_gb:.1f} GB free", status="ok")

    # 6. input_pdfs reachable
    if cfg.input_pdfs is None:
        pstatus("input_pdfs: not set (required for `corpus run`)", status="warn")
    else:
        input_dir = _resolve_against(config_path, cfg.input_pdfs)
        if input_dir.exists():
            # When input_pdfs is a parent of output_dir (the demo's
            # `input_pdfs: .` is the canonical case), naive rglob picks
            # up every output/documents/<HASH>/processed.pdf and reports
            # twice the real count. Exclude paths under output_dir; the
            # pipeline already dedupes by SHA-256 so processing is
            # correct either way, but the operator-facing count must
            # reflect actual source PDFs.
            output_resolved = output_dir.resolve() if input_dir.is_dir() else None
            def _is_under_output(p: Path) -> bool:
                if output_resolved is None:
                    return False
                try:
                    p.resolve().relative_to(output_resolved)
                    return True
                except ValueError:
                    return False
            n_pdfs = (
                sum(1 for p in input_dir.rglob("*.pdf") if not _is_under_output(p))
                if input_dir.is_dir() else 0
            )
            pstatus(f"input_pdfs: {input_dir} ({n_pdfs} PDFs)", status="ok")
        else:
            failures.append(f"input_pdfs path does not exist: {input_dir}")
            pstatus(f"input_pdfs: {input_dir} not found", status="fail")

    # 7. Taxonomy source availability
    if cfg.taxonomy.source is not None:
        db_path = output_dir / "taxonomy.sqlite"
        if cfg.taxonomy.source == "worms":
            # WoRMS reaches out to the network. On an internet-connected node
            # a first `corpus run` will build taxonomy.sqlite automatically.
            # On HPC compute nodes (no outbound internet), the sqlite must be
            # pre-built on a login node first — `corpus run --only extract`
            # will fail loudly if it is absent.
            if db_path.exists():
                pstatus(
                    f"Taxonomy: {db_path.name} present (WoRMS pre-built)",
                    status="ok",
                )
            else:
                pstatus(
                    "Taxonomy: taxonomy.sqlite not yet built — WoRMS will be "
                    "fetched on first run (requires internet; pre-build on a "
                    "login node before submitting HPC array jobs)",
                    status="warn",
                )
        else:  # dwca / dwc
            tx_path = (
                _resolve_against(config_path, cfg.taxonomy.path)
                if cfg.taxonomy.path is not None else None
            )
            if tx_path is None:
                failures.append(
                    f"taxonomy.source={cfg.taxonomy.source!r} but "
                    "taxonomy.path is not set in config."
                )
                pstatus(
                    f"Taxonomy: source={cfg.taxonomy.source!r}, path not set",
                    status="fail",
                )
            elif not tx_path.exists():
                failures.append(
                    f"taxonomy.source={cfg.taxonomy.source!r} but "
                    f"taxonomy.path={tx_path} does not exist."
                )
                pstatus(f"Taxonomy: {tx_path} not found", status="fail")
            elif db_path.exists():
                pstatus(
                    f"Taxonomy: {db_path.name} present "
                    f"(source: {cfg.taxonomy.source})",
                    status="ok",
                )
            else:
                # Archive/dir is present but sqlite not yet built — the next
                # full `corpus run` will build it via ingest_taxonomy.
                pstatus(
                    f"Taxonomy: {tx_path.name} found, "
                    f"{db_path.name} not yet built (will be created on first run)",
                    status="warn",
                )

    # 8. macOS Python arch — Rosetta'd Python on Apple Silicon traps the
    # env in the unsupported macOS x86_64 matrix (Apple dropped Intel-mac
    # torch wheels after 2.2, breaking docling + transformers ≥ 5). Hard
    # fail loud rather than letting `corpus run` discover it deep in a
    # model load. Linux is not checked: there is no Rosetta equivalent.
    arch, arch_failure = _check_python_arch()
    if arch:
        if arch_failure is None:
            pstatus(f"Python arch: {arch} (native)", status="ok")
        else:
            failures.append(arch_failure)
            pstatus(f"Python arch: {arch} (Rosetta — unsupported)", status="fail")

    if failures:
        print()
        pstatus(f"{len(failures)} precondition(s) failed:", status="fail")
        for f in failures:
            print(f"  - {f}")
        return EXIT_PRECONDITION
    print()
    pstatus("ready: `corpus run` should succeed on this host.", status="ok")
    return EXIT_OK


def _check_python_arch() -> Tuple[str, Optional[str]]:
    """Check the active Python's architecture on macOS.

    Returns ``(arch_label, failure_message)``:
    - On non-darwin: ``("", None)`` — caller skips the check entirely.
    - On darwin/arm64: ``("arm64", None)``.
    - On darwin/<other>: ``(arch, "...explanatory failure...")``.

    Separated from ``_cmd_check`` so the contract is unit-testable
    via mocked ``sys.platform`` + ``platform.machine`` without
    standing up a full config + grobid + disk environment.
    """
    if sys.platform != "darwin":
        return "", None
    import platform as _platform
    arch = _platform.machine()
    if arch == "arm64":
        return arch, None
    return arch, (
        f"Python arch is {arch} (Rosetta on Apple Silicon). corpus "
        f"requires an arm64-native conda env; see README §Supported "
        f"platforms for the fix."
    )


def _detect_accelerator() -> Optional[str]:
    """Return 'cuda', 'mps', or None. Uses torch if importable; otherwise None."""
    try:
        import torch
        if torch.cuda.is_available():
            return "cuda"
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return "mps"
    except Exception:
        pass
    return None


def _vision_skip_reason(mode: str) -> Optional[str]:
    """Return a human-readable reason the ``vision-*`` panel mode isn't
    usable on this host, or None if it is (#65, #102).

    Gracefully avoids loading a multi-GB open-weights model on a CPU-only
    host or hitting the Anthropic API without credentials, both of which
    hard-fail Pass 3b deep into the pipeline. The caller downgrades to
    the OCR panel floor (``ocr``) when this returns a reason.
    """
    if mode == "vision-local":
        if _detect_accelerator() is None:
            return (
                "figures.panel_detection=vision-local but no CUDA/MPS "
                "detected on this host. Use `figures.panel_detection: "
                "vision-claude` (and export ANTHROPIC_API_KEY), or `ocr`, "
                "or run on a GPU host."
            )
    elif mode == "vision-claude":
        if not os.environ.get("ANTHROPIC_API_KEY"):
            return (
                "figures.panel_detection=vision-claude but ANTHROPIC_API_KEY "
                "is unset. Either export ANTHROPIC_API_KEY, or use "
                "`figures.panel_detection: vision-local` / `ocr`."
            )
    return None


def _ping_grobid(url: str) -> tuple[bool, str]:
    """Probe Grobid's /api/isalive endpoint with a 2s timeout."""
    try:
        import requests
        r = requests.get(url.rstrip("/") + "/api/isalive", timeout=2)
        if r.status_code == 200 and r.text.strip().lower() == "true":
            return True, "alive"
        return False, f"HTTP {r.status_code}: {r.text[:40]}"
    except Exception as e:
        return False, type(e).__name__


def _free_disk_gb(path: Path) -> Optional[float]:
    """Free disk-space for the partition holding ``path`` (creates parents if missing)."""
    try:
        # Walk up to the first existing parent if the path itself is missing.
        p = path
        while not p.exists() and p != p.parent:
            p = p.parent
        return shutil.disk_usage(p).free / (1024 ** 3)
    except Exception:
        return None


_BASH_COMPLETION = """\
# corpus(1) bash completion — generated by `corpus completion bash`
_corpus_complete() {
    local cur prev cmd
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Top-level verbs.
    local verbs="init run debug-pdf status serve bib check completion"
    if [[ ${COMP_CWORD} -eq 1 ]] || [[ "$prev" == "-c" || "$prev" == "--config" ]]; then
        if [[ "$prev" == "-c" || "$prev" == "--config" ]]; then
            COMPREPLY=( $(compgen -f -- "$cur") )
        else
            COMPREPLY=( $(compgen -W "$verbs --help --version --cite --config -c -V -h -v -vv -q --quiet" -- "$cur") )
        fi
        return 0
    fi

    # Per-verb flag completions.
    cmd="${COMP_WORDS[1]}"
    case "$cmd" in
        run)        COMPREPLY=( $(compgen -W "--dry-run --force-rebuild --figure-panels --no-vision --no-bundle --enrich-bhl -h --help" -- "$cur") ) ;;
        debug-pdf)  COMPREPLY=( $(compgen -W "--output-dir --grobid-url -h --help" -- "$cur") ) ;;
        status)     COMPREPLY=( $(compgen -W "--output-dir --report --json --list-hashes --filter-stage --filter-reason --filter-gate -v -h --help" -- "$cur") ) ;;
        serve)      COMPREPLY=( $(compgen -W "--output-dir --transport --host --port --auth-token-file -h --help" -- "$cur") ) ;;
        bib)        COMPREPLY=( $(compgen -W "export import -h --help" -- "$cur") ) ;;
        init)       COMPREPLY=( $(compgen -W "--path --force -h --help" -- "$cur") ) ;;
        check)      COMPREPLY=( $(compgen -W "-h --help" -- "$cur") ) ;;
        completion) COMPREPLY=( $(compgen -W "bash zsh fish -h --help" -- "$cur") ) ;;
    esac
    return 0
}
complete -o default -F _corpus_complete corpus
"""

_ZSH_COMPLETION = """\
#compdef corpus
# corpus(1) zsh completion — generated by `corpus completion zsh`
_corpus() {
    local -a verbs
    verbs=(
        'init:scaffold config.yaml in cwd'
        'run:full pipeline + post-pipeline + bundle'
        'debug-pdf:run the per-paper pipeline on one PDF with stage tracing'
        'status:build-state rollup + report'
        'serve:MCP server'
        'bib:BibTeX round-trip'
        'check:environment + config pre-flight'
        'completion:generate a shell completion script'
    )
    if (( CURRENT == 2 )); then
        _describe 'command' verbs
    elif (( CURRENT >= 3 )); then
        case "$words[2]" in
            run)        _arguments '--dry-run' '--force-rebuild' '--figure-panels:mode:(ocr vision-local vision-claude off)' '--no-vision' '--no-bundle' '--enrich-bhl' ;;
            debug-pdf)  _arguments '--output-dir:DIR:_files -/' '--grobid-url:URL:' '1:PDF:_files' ;;
            status)     _arguments '--output-dir:DIR:_files -/' '--report' '--json' '--list-hashes' ;;
            serve)      _arguments '--output-dir:DIR:_files -/' '--transport' '--host' '--port' '--auth-token-file:FILE:_files' ;;
            bib)        _values 'subcommand' export import ;;
            completion) _values 'shell' bash zsh fish ;;
        esac
    fi
}
_corpus "$@"
"""

_FISH_COMPLETION = """\
# corpus(1) fish completion — generated by `corpus completion fish`
complete -c corpus -n "__fish_use_subcommand" -a init        -d 'scaffold config.yaml in cwd'
complete -c corpus -n "__fish_use_subcommand" -a run         -d 'full pipeline + post-pipeline + bundle'
complete -c corpus -n "__fish_use_subcommand" -a debug-pdf   -d 'run the per-paper pipeline on one PDF with stage tracing'
complete -c corpus -n "__fish_use_subcommand" -a status      -d 'build-state rollup + report'
complete -c corpus -n "__fish_use_subcommand" -a serve       -d 'MCP server'
complete -c corpus -n "__fish_use_subcommand" -a bib         -d 'BibTeX round-trip'
complete -c corpus -n "__fish_use_subcommand" -a check       -d 'environment + config pre-flight'
complete -c corpus -n "__fish_use_subcommand" -a completion  -d 'generate a shell completion script'
complete -c corpus -s c -l config -r -d 'Path to config.yaml'
complete -c corpus -s V -l version
complete -c corpus -l cite -d 'Print the citation (plain | bibtex)'
complete -c corpus -s v -d 'verbose (-vv = debug all)'
complete -c corpus -s q -l quiet
complete -c corpus -n "__fish_seen_subcommand_from run" -l dry-run -l force-rebuild -l no-vision -l no-bundle -l enrich-bhl
complete -c corpus -n "__fish_seen_subcommand_from run" -l figure-panels -x -a 'ocr vision-local vision-claude off'
complete -c corpus -n "__fish_seen_subcommand_from debug-pdf" -l output-dir -l grobid-url
complete -c corpus -n "__fish_seen_subcommand_from completion" -a 'bash zsh fish'
"""


def _cmd_completion(args: argparse.Namespace) -> int:
    snippets = {"bash": _BASH_COMPLETION, "zsh": _ZSH_COMPLETION, "fish": _FISH_COMPLETION}
    sys.stdout.write(snippets[args.shell])
    return EXIT_OK


# ---------------------------------------------------------------------------
# Parser construction
# ---------------------------------------------------------------------------


def _git_sha_or_none() -> Optional[str]:
    """Short git SHA when the install was sourced from a development clone."""
    repo_root = Path(__file__).resolve().parent.parent
    if not (repo_root / ".git").exists():
        return None
    try:
        out = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True, timeout=2,
        )
        return out.stdout.strip() or None
    except (subprocess.SubprocessError, FileNotFoundError):
        return None


def _version_string() -> str:
    sha = _git_sha_or_none()
    if sha:
        return f"corpus {__version__} ({sha})"
    return f"corpus {__version__}"


def _read_citation_cff() -> dict:
    """Read CITATION.cff. The canonical copy at the repo root is shipped
    via package_data into ``pipeline/CITATION.cff`` so wheel installs
    (`pip install git+...@vX.Y.Z`) can find it via importlib.resources.
    Falls back to the repo-root copy for editable installs that
    haven't been re-installed since CITATION.cff was added.
    """
    # Packaged copy first (works for both editable and wheel installs).
    try:
        src = resources.files(_TEMPLATE_PACKAGE).joinpath("CITATION.cff")
        with resources.as_file(src) as cff_path:
            if cff_path.exists():
                return yaml.safe_load(cff_path.read_text(encoding="utf-8"))
    except (ModuleNotFoundError, FileNotFoundError):
        pass
    # Fallback: repo root (editable install pre-dating the package_data move).
    repo_root = Path(__file__).resolve().parent.parent
    cff_path = repo_root / "CITATION.cff"
    if cff_path.exists():
        return yaml.safe_load(cff_path.read_text(encoding="utf-8"))
    print_status(
        f"CITATION.cff not found in package data or at {cff_path}. "
        "Re-install with `pip install -e .` (or `pip install git+...`) "
        "to pick up the packaged copy.",
        status="fail",
    )
    sys.exit(EXIT_GENERIC)


def _format_citation_plain(cff: dict) -> str:
    pref = cff.get("preferred-citation") or cff
    authors = pref.get("authors") or []
    parts = []
    for a in authors:
        family = a.get("family-names", "")
        given = a.get("given-names", "")
        # "C. W." style initials
        initials = ". ".join(p[0] for p in given.split()) + "." if given else ""
        parts.append(f"{family}, {initials}".strip(", "))
    if len(parts) > 1:
        author_str = ", ".join(parts[:-1]) + ", & " + parts[-1]
    else:
        author_str = parts[0] if parts else ""
    title = pref.get("title", "")
    year = pref.get("year", "")
    doi = pref.get("doi", "")
    url = pref.get("url") or (f"https://doi.org/{doi}" if doi else "")
    return f"{author_str} ({year}). {title}. {url}".rstrip(" .")


def _format_citation_bibtex(cff: dict) -> str:
    pref = cff.get("preferred-citation") or cff
    authors = pref.get("authors") or []
    bib_authors = " and ".join(
        f"{a.get('family-names', '')}, {a.get('given-names', '')}".strip(", ")
        for a in authors
    )
    title = pref.get("title", "").replace("{", "").replace("}", "")
    year = pref.get("year", "")
    doi = pref.get("doi", "")
    url = pref.get("url") or (f"https://doi.org/{doi}" if doi else "")
    return (
        "@misc{corpus,\n"
        f"  author = {{{bib_authors}}},\n"
        f"  title  = {{{title}}},\n"
        f"  year   = {{{year}}},\n"
        f"  doi    = {{{doi}}},\n"
        f"  url    = {{{url}}},\n"
        "}\n"
    )


def _print_citation(fmt: str) -> int:
    cff = _read_citation_cff()
    if fmt == "bibtex":
        sys.stdout.write(_format_citation_bibtex(cff))
    else:
        sys.stdout.write(_format_citation_plain(cff) + "\n")
    return EXIT_OK


def _setup_logging(args: argparse.Namespace) -> None:
    """Wire root logging from -v/-vv/-q.

    Also exports CORPUS_LOG_LEVEL into the environment so subprocesses
    (orchestrator + downstream modules, mcpsrv, bib) honor the same
    threshold via each package's __init__.py — without this, `-q`
    quiets only this CLI process and not the subprocess tree it spawns.
    """
    if args.quiet:
        level = logging.WARNING
        env_level = "WARNING"
    elif args.verbose >= 2:
        level = logging.DEBUG
        env_level = "DEBUG"
    elif args.verbose == 1:
        level = logging.INFO
        env_level = "INFO"
    else:
        level = logging.WARNING
        env_level = None  # leave subprocesses at their existing default
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    if env_level is not None:
        os.environ["CORPUS_LOG_LEVEL"] = env_level


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="corpus",
        description=(
            "Top-level orchestrator for a corpus instance — a 'corpuscle' "
            "assembled from a folder of scientific PDFs. Scaffolds the "
            "corpuscle's config, runs the full extraction + embedding + "
            "cross-paper pipeline, and serves the result as a Model Context "
            "Protocol (MCP) server that LLM clients (Claude Desktop, Claude "
            "Code, claude.ai, Cursor) can query. "
            "Run `corpus <command> --help` for details on any subcommand."
        ),
    )
    # --version (#58 + #61): includes git SHA on dev clones.
    p.add_argument(
        "-V", "--version", action="version", version=_version_string(),
    )
    # --cite (#61): plain by default, BibTeX with --cite=bibtex.
    p.add_argument(
        "--cite",
        nargs="?",
        const="plain",
        choices=["plain", "bibtex"],
        default=None,
        help="Print the citation. Source of truth: CITATION.cff. "
        "`--cite` for plain text; `--cite=bibtex` for BibTeX.",
    )
    # Global --config (pre-verb, git/cargo/kubectl style; #60 + #61).
    p.add_argument(
        "-c", "--config",
        type=Path,
        default=None,
        help="Path to config.yaml. Default: $CORPUS_CONFIG, else "
             "./config.yaml in cwd.",
    )
    # Verbosity (#61). -v INFO, -vv DEBUG, -q WARNING.
    p.add_argument(
        "-v", "--verbose",
        action="count",
        default=0,
        help="Raise log threshold (-v=INFO, -vv=DEBUG).",
    )
    p.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Drop log threshold to WARNING.",
    )

    sub = p.add_subparsers(dest="command", metavar="<command>")

    # --- init ---
    init_p = sub.add_parser(
        "init",
        help="Scaffold config.yaml in cwd from the bundled template",
        description=(
            "Scaffold config.yaml in the current directory from the "
            "bundled template. Run once per new corpuscle, then edit the "
            "file to point at your PDFs, optional BibTeX, optional "
            "lexicon, and a taxonomy source (WoRMS / DwC-A). After that, "
            "`corpus check` validates the config and `corpus run` executes "
            "the pipeline."
        ),
    )
    init_p.add_argument("--path", type=Path, default=None)
    init_p.add_argument("--force", action="store_true")
    init_p.set_defaults(func=_cmd_init)

    # --- run ---
    run_p = sub.add_parser(
        "run",
        help="Full pipeline: extract → embed → cross-paper builds → bundle",
        description=(
            "Run the full corpus build pipeline end-to-end against the "
            "PDFs in `input_pdfs`. Stages: per-paper extraction (OCR if "
            "needed, docling layout, Grobid metadata, chunking, taxa + "
            "lexicon tagging) → BGE-M3 embeddings into LanceDB → cross-"
            "paper bibliography reconciliation + taxon-mention rollup → "
            "distillation of a served bundle at <output_dir>/_serve/. "
            "Idempotent: a re-run only re-processes papers whose inputs "
            "have changed."
        ),
    )
    run_p.add_argument("--dry-run", action="store_true",
                       help="Plan only — no artifacts written")
    run_p.add_argument("--force-rebuild", action="store_true",
                       help="Rebuild every cross-paper DB (#64)")
    run_p.add_argument("--force-rebuild-taxonomy", action="store_true",
                       help="Rebuild only taxonomy.sqlite (#64)")
    run_p.add_argument("--force-rebuild-biblio", action="store_true",
                       help="Rebuild only biblio_authority.sqlite (#64)")
    run_p.add_argument("--force-rebuild-taxon-mentions", action="store_true",
                       help="Rebuild only taxon_mentions.sqlite (#64)")
    # Figure panel-ROI detection (#102). The per-run override of the
    # corpuscle's `figures.panel_detection` config key. `--no-vision` is
    # the deprecated alias for `--figure-panels ocr`; the two are
    # mutually exclusive.
    panels_grp = run_p.add_mutually_exclusive_group()
    panels_grp.add_argument(
        "--figure-panels",
        choices=["ocr", "vision-local", "vision-claude", "off"],
        default=None,
        help="Override config `figures.panel_detection` for this run: "
             "ocr (CPU OCR-panel floor), vision-local / vision-claude "
             "(Pass 3b vision), or off. Omit to use config.yaml.",
    )
    panels_grp.add_argument(
        "--no-vision", action="store_true",
        help="Deprecated alias for `--figure-panels ocr`: skip the vision "
             "pass (Pass 3b) and use the CPU OCR-panel floor instead. Does "
             "NOT disable panels — pass `--figure-panels off` for that.",
    )
    run_p.add_argument("--no-bundle", action="store_true",
                       help="Skip the served-bundle distillation step. "
                       "Build artifacts in `<output_dir>/` are unaffected; "
                       "only the `<output_dir>/_serve/` distillation is "
                       "skipped. (#60)")
    run_p.add_argument("--no-prune", action="store_true",
                       help="Skip orphan cleanup; run a read-only audit "
                       "instead (matches #31 behavior). #66")
    run_p.add_argument("--force-prune", action="store_true",
                       help="Bypass the safety rail (#66) that refuses to "
                       "prune when more than 25%% of hash dirs would be "
                       "removed (likely a config mistake or unmounted volume).")
    run_p.add_argument("--re-process-flagged", default=None, metavar="GATE",
                       help="Clear pipeline_state.json for every paper "
                       "carrying this quality_flag (e.g. gibberish_after_ocr, "
                       "empty_text), then run the pipeline. After the "
                       "operator fixes the root cause (installs a language "
                       "pack, etc.), this re-extracts just the affected "
                       "papers without touching the rest. (#54)")
    run_p.add_argument("--skip-checks", action="store_true",
                       help="Bypass the fail-fast preconditions (Grobid "
                       "reachability + input_pdfs presence). Use when you "
                       "know what you're doing — for example, running a "
                       "vision-only refresh on an existing build.")
    run_p.add_argument("--enrich-bhl", action="store_true",
                       help="Enrich pre-DOI references against the Biodiversity "
                       "Heritage Library (slow, rate-limited; #64)")
    # HPC / job-array support: run one phase per SLURM job, on the right
    # partition, and slice the per-paper stages into array tasks. Omit
    # --only to run the whole pipeline on one node (the default).
    run_p.add_argument(
        "--only", choices=["extract", "vision", "embed", "post", "bundle"],
        default=None,
        help="Run only one phase (for HPC, where each maps to a separate "
             "SLURM job/partition): extract (CPU), vision (GPU Pass 3b "
             "refresh), embed (GPU), post (cross-paper DBs), bundle (served-"
             "bundle distill). Omit to run the full pipeline on one node.")
    run_p.add_argument(
        "--batch-index", type=int, default=None,
        help="0-based array-task index; with --batch-size, processes a "
             "deterministic slice of the sorted hash list (SLURM job "
             "arrays). Only with --only extract|vision.")
    run_p.add_argument(
        "--batch-size", type=int, default=None,
        help="PDFs per array task (pairs with --batch-index).")
    run_p.set_defaults(func=_cmd_run)

    # --- debug-pdf ---
    debug_pdf_p = sub.add_parser(
        "debug-pdf",
        help="Run the per-paper pipeline on a single PDF with stage tracing.",
        description=(
            "Troubleshooting runner for hard documents (#92). Drives one "
            "PDF through the per-paper extraction pipeline "
            "(huge-doc gate → scan detection → PDF prep → docling "
            "extraction → chunking → figures, plus bibliography when "
            "--grobid-url is given), writes the intermediate artifacts "
            "under a debug directory for inspection, and prints a "
            "per-stage ok/fail + timing table so a failure is immediately "
            "attributable to a stage. Runs inline with resume disabled — "
            "every stage executes every time. No bundle, no cross-paper "
            "steps. The taxa/lexicon stages need corpus-level inputs and "
            "are out of scope; use a full `corpus run` for those."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    debug_pdf_p.add_argument(
        "pdf", type=Path, help="Path to the single PDF to debug.",
    )
    debug_pdf_p.add_argument(
        "--output-dir", type=Path, default=None,
        help="Where to write the debug artifacts (default: ./debug-pdf/). "
             "A <hash>/ subdir is created per PDF, mirroring documents/.",
    )
    debug_pdf_p.add_argument(
        "--grobid-url", default=None, metavar="URL",
        help="Enable the bibliography stage against a Grobid server "
             "(e.g. http://localhost:8070). Skipped when omitted.",
    )
    debug_pdf_p.set_defaults(func=_cmd_debug_pdf)

    # --- status ---
    status_p = sub.add_parser(
        "status",
        help="Postmortem report on the build state of a corpuscle.",
        description=(
            "Read-only postmortem of a built corpuscle: stage success "
            "rates, per-paper failures grouped by reason code, quality-"
            "flag rollup, cross-paper artifact presence (taxonomy / "
            "biblio / taxon-mentions SQLite), and triage modes. Reads "
            "<output_dir>/documents/*/summary.json and cross-paper "
            "SQLite; writes nothing."
        ),
        epilog=(
            "Forwarded flags (see `python -m pipeline.status --help`):\n"
            "  --report               full text report (default rollup is short)\n"
            "  --json                 emit the rollup as JSON\n"
            "  --sort-by <metric>     worst-N papers; pair with --tail\n"
            "  --propose-skips        BibTeX-ready `serve = {false}` candidates\n"
            "  --skipped              currently-excluded papers, grouped by reason\n"
            "  --list-hashes          one short hash per line (xargs-friendly)\n"
            "  --filter-{stage,reason,gate} <name>   narrow the set"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    status_p.add_argument(
        "--output-dir", type=Path, default=None,
        help="Override the corpuscle output dir from config.yaml.",
    )
    status_p.set_defaults(func=_cmd_status)

    # --- serve ---
    serve_p = sub.add_parser(
        "serve",
        help="Serve the built corpus as a Model Context Protocol server.",
        description=(
            "Serve the built corpuscle as a Model Context Protocol (MCP) "
            "server. Default transport is stdio — meant for an MCP "
            "client (Claude Desktop, Claude Code, Cursor) that launches "
            "this process itself. For remote deploys or manual curl "
            "testing, pass `--transport sse --host <ip> --port <port> "
            "--auth-token-file <file>`. The server reads the build "
            "bundle at `output_dir` (from config.yaml) by default; point "
            "at `<output_dir>/_serve/` to serve a distilled bundle "
            "directly."
        ),
        epilog=(
            "Forwarded flags (see `python -m mcpsrv.main --help`):\n"
            "  --transport {stdio,sse}        default: stdio\n"
            "  --host <ip>                    SSE bind address; default 127.0.0.1\n"
            "  --port <n>                     SSE bind port; default 8080\n"
            "  --auth-token-file <path>       bearer token (required for SSE)\n"
            "  --instructions <path>          override the per-corpuscle nudges\n"
            "  --taxonomy-db / --biblio-sqlite / --taxon-mention-sqlite <path>\n"
            "                                 per-component DB overrides\n"
            "  --check                        non-binding pre-flight + exit\n"
            "  --name <id>                    server name reported to MCP clients\n"
            "  --default-profile <name>       fallback output profile for calls\n"
            "                                 without profile= (report|manuscript|\n"
            "                                 presentation; default report)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    serve_p.add_argument("--output-dir", type=Path, default=None)
    serve_p.set_defaults(func=_cmd_serve)

    # --- bib (export | import) ---
    bib_p = sub.add_parser(
        "bib",
        help="Export / import the corpus bibliography for hand curation.",
        description=(
            "Round-trip the corpus bibliography through BibTeX so you "
            "can curate it in your editor of choice. `corpus bib export` "
            "writes the current `biblio_authority.sqlite` to a `.bib` "
            "file with stable `corpus_hash` keys; edit titles, authors, "
            "years, or the `serve = {false}` skip flag as you would any "
            "BibTeX file; `corpus bib import` reads the edits back into "
            "the authority DB. Use this to iteratively refine "
            "bibliographic metadata when Grobid mis-parses an unusual "
            "layout."
        ),
    )
    bib_sub = bib_p.add_subparsers(dest="bib_command", metavar="{export,import}")

    bib_export_p = bib_sub.add_parser(
        "export",
        help="Dump biblio_authority.sqlite to a .bib file.",
        description=(
            "Dump the current `biblio_authority.sqlite` to a BibTeX file. "
            "Forwards to `python -m bib.export`; pass `-o <path>` to "
            "choose the output file (default: stdout)."
        ),
    )
    bib_export_p.add_argument("--output-dir", type=Path, default=None)
    bib_export_p.set_defaults(func=_cmd_bib_export)

    bib_import_p = bib_sub.add_parser(
        "import",
        help="Re-apply hand-edited BibTeX entries back into the authority DB.",
        description=(
            "Read a hand-curated BibTeX file and apply the edits back "
            "into `biblio_authority.sqlite`. Matches entries by their "
            "`corpus_hash` field (DOI as a fallback) so renamed entries "
            "and re-builds don't break the round-trip. Pass the .bib "
            "file as the first positional argument."
        ),
    )
    bib_import_p.add_argument("--output-dir", type=Path, default=None)
    bib_import_p.set_defaults(func=_cmd_bib_import)

    bib_p.set_defaults(func=lambda a: (bib_p.print_help() or 2))

    # --- taxonomy (export | ingest) ---
    tax_p = sub.add_parser(
        "taxonomy",
        help="Build or share the corpus taxonomy via Darwin Core.",
        description=(
            "Round-trip the corpus taxonomy through a Darwin Core "
            "Archive. `corpus taxonomy ingest` builds "
            "`taxonomy.sqlite` from a WoRMS subtree, a DwC-A .zip, or "
            "a directory of Darwin Core .tsv files. `corpus taxonomy "
            "export` dumps the current `taxonomy.sqlite` back out as "
            "a DwC-A — useful for sharing a snapshot without forcing "
            "the recipient to re-walk the WoRMS API, and for "
            "committing small fixtures into a repo so CI exercises "
            "the dwca ingest path without external network calls."
        ),
    )
    tax_sub = tax_p.add_subparsers(dest="taxonomy_command",
                                   metavar="{export,ingest}")

    tax_export_p = tax_sub.add_parser(
        "export",
        help="Dump taxonomy.sqlite to a Darwin Core Archive (.zip).",
        description=(
            "Dump the current `taxonomy.sqlite` to a DwC-A ZIP "
            "containing meta.xml + taxon.tsv (and vernacularname.tsv "
            "when vernacular names are present). Forwards to "
            "`python -m pipeline.taxonomy_export`; pass `-o <path>` "
            "to choose the output ZIP (required)."
        ),
    )
    tax_export_p.add_argument("--output-dir", type=Path, default=None)
    tax_export_p.set_defaults(func=_cmd_taxonomy_export)

    tax_ingest_p = tax_sub.add_parser(
        "ingest",
        help="Build taxonomy.sqlite from WoRMS / DwC-A / DwC files.",
        description=(
            "Build `taxonomy.sqlite` from one of three sources. "
            "`--source worms --root-id <AphiaID>` walks the WoRMS REST "
            "API. `--source dwca --input <path.zip>` ingests a Darwin "
            "Core Archive (use this with a fixture produced by "
            "`corpus taxonomy export`). `--source dwc --input "
            "<Taxon.tsv>` ingests a single DwC Taxon file. Forwards "
            "to `python -m pipeline.taxonomy_ingest`."
        ),
    )
    tax_ingest_p.add_argument("--output-dir", type=Path, default=None)
    tax_ingest_p.set_defaults(func=_cmd_taxonomy_ingest)

    tax_p.set_defaults(func=lambda a: (tax_p.print_help() or 2))

    # --- check (#62) ---
    check_p = sub.add_parser(
        "check",
        help="Pre-flight: can this host run `corpus run` successfully?",
        description=(
            "Pre-flight check: does this host have what `corpus run` "
            "will need? Validates config.yaml against the schema, "
            "probes for Grobid reachability, GPU/MPS availability, "
            "ANTHROPIC_API_KEY presence, and disk space. Exits 0 on "
            "green, 2 on config error, 3 on environment precondition "
            "failure. Distinct from `corpus serve --check` (which is a "
            "serve-time bundle pre-flight) and `corpus run --dry-run` "
            "(which prints the work plan)."
        ),
    )
    check_p.set_defaults(func=_cmd_check)

    # --- completion (#61) ---
    comp_p = sub.add_parser(
        "completion",
        help="Emit a shell-completion script for bash/zsh/fish.",
        description=(
            "Emit a shell-completion script to stdout for the given "
            "shell. Pipe it into your dotfiles for tab-completion of "
            "corpus subcommands and their options. Example:\n\n"
            "  corpus completion bash > ~/.local/share/corpus.bash\n"
            "  echo 'source ~/.local/share/corpus.bash' >> ~/.bashrc"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    comp_p.add_argument("shell", choices=["bash", "zsh", "fish"])
    comp_p.set_defaults(func=_cmd_completion)

    return p


_PASSTHROUGH_VERBS = {"status", "serve", "bib", "taxonomy"}


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    # Use parse_known_args so passthrough subcommands (status, serve,
    # bib) can forward flags like --transport, --report, --json to the
    # downstream module without argparse rejecting them as unrecognized.
    args, extras = parser.parse_known_args(argv)
    if args.cite is not None:
        return _print_citation(args.cite)
    if not args.command:
        parser.print_help()
        return EXIT_OK
    if extras and args.command not in _PASSTHROUGH_VERBS:
        parser.error(f"unrecognized arguments: {' '.join(extras)}")
    # Strip a leading `--` separator if present. Operators following the
    # documented `corpus serve --output-dir <bundle> -- --transport sse ...`
    # form would otherwise forward the literal `--` to the downstream
    # argparse, which interprets it as end-of-options and demotes
    # subsequent flags to positional args (which then fail).
    if extras and extras[0] == "--":
        extras = extras[1:]
    args.passthrough = extras
    _setup_logging(args)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
