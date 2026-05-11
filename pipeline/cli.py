"""`corpus` CLI entry point.

Subcommand router for v0.3 (#60 + #61). Replaces the 13 root-level
scripts that v0.2 shipped. Subcommands:

  corpus init                  scaffold config.yaml in cwd (#59)
  corpus run                   full pipeline + post-pipeline + bundle
  corpus status                build-state rollup + report (#57)
  corpus serve                 MCP server
  corpus bib export|import     BibTeX round-trip
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
from importlib import resources
from pathlib import Path
from typing import List, Optional

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

    # Vision pass (#65): opt-out via --no-vision flag, vision.backend=none,
    # OR capability detection — skip gracefully (with a one-line nudge)
    # when the configured backend isn't usable on this host.
    if not args.no_vision and cfg.vision.backend != "none":
        skip_reason = _vision_skip_reason(cfg.vision.backend)
        if skip_reason is not None:
            print_status(
                f"vision pass skipped: {skip_reason}",
                status="warn",
            )
        else:
            sub_argv += ["--vision-backend", cfg.vision.backend]
            if cfg.vision.model:
                sub_argv += ["--vision-model", cfg.vision.model]

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

    # Fail-fast preconditions: catch the two ergonomic cliffs that would
    # otherwise burn hours of compute (Grobid unreachable → placeholder
    # metadata across the whole bibliography pass; input_pdfs missing →
    # silent zero-paper run). The vision/ANTHROPIC_API_KEY paths are
    # warn-and-skip (#65) — graceful degradation, not a fail-fast case.
    if not args.skip_checks:
        rc = _run_preconditions(cfg, config_path, args)
        if rc != EXIT_OK:
            return rc

    # #66: orphan cleanup runs *before* the extract step so subsequent
    # stages don't consume the orphan hashes. Default = prune; --no-prune
    # = read-only audit only; --force-prune bypasses the safety rail.
    rc = _prune_orphans(cfg, config_path, args)
    if rc != EXIT_OK:
        return rc

    # #54: --re-process-flagged invalidates pipeline_state.json for every
    # paper carrying the named quality_flag, so the per-stage resume
    # logic re-extracts them on the upcoming run.
    if args.re_process_flagged:
        _invalidate_flagged(cfg, config_path, args.re_process_flagged)

    sub_argv = _build_orchestrator_argv(cfg, config_path, args)
    cmd = [sys.executable, "-m", "pipeline.orchestrator", *sub_argv]
    print_status(f"corpus run → {' '.join(cmd)}", status="info")
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        return rc

    # #60 — bundle distillation in line. Produce a ready-to-ship served
    # bundle alongside the build bundle unless --no-bundle. The served
    # bundle lands at `<output_dir>/_serve/` and is the directory remote
    # deploys (DEPLOY.md) ship to S3 / EC2.
    output_dir = _resolve_against(config_path, cfg.output_dir)
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
                status="warn",
            )
            # Don't propagate — the build bundle is valid; bundle is a
            # ship-to-host convenience.

    if args.dry_run:
        print_status(
            "dry-run complete. No artifacts written — re-run without "
            "`--dry-run` to execute the plan above.",
            status="ok",
        )
    else:
        print_status(
            "run complete. Try `corpus status --report` and `corpus serve` next.",
            status="ok",
        )
    return EXIT_OK


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
    config_path = _resolve_config_path(args.config)
    failures: List[str] = []  # precondition (exit 3)
    config_failures: List[str] = []  # config (exit 2)

    # 1. config.yaml resolution + schema validation
    if config_path is None:
        print_status(
            "no config.yaml in cwd; pass --config PATH or run `corpus init`",
            status="fail",
        )
        return EXIT_CONFIG_ERROR
    if not config_path.exists():
        print_status(f"config not found: {config_path}", status="fail")
        return EXIT_CONFIG_ERROR
    try:
        raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        cfg = validate_config(raw)
        print_status(f"config.yaml valid ({config_path})", status="ok")
    except ValidationError as e:
        print_status(f"config error in {config_path}", status="fail")
        for err in e.errors():
            loc = ".".join(str(x) for x in err["loc"])
            print(f"  {loc}: {err['msg']}", file=sys.stderr)
        return EXIT_CONFIG_ERROR

    # 2. GPU availability
    gpu = _detect_accelerator()
    if gpu == "cuda":
        print_status("GPU: CUDA available", status="ok")
    elif gpu == "mps":
        print_status("GPU: MPS available (Apple Silicon)", status="ok")
    else:
        print_status("GPU: none (CPU-only; embedding pass will be slow)", status="warn")

    # 3. Vision backend usability — under #65 a missing capability is
    # a warn (the pass auto-skips at run time with the same message),
    # not a hard precondition failure.
    vision_skip = _vision_skip_reason(cfg.vision.backend)
    if cfg.vision.backend == "none":
        print_status("Vision pass: disabled in config", status="warn")
    elif vision_skip is None:
        print_status(f"Vision pass: backend `{cfg.vision.backend}` ready", status="ok")
    else:
        print_status(f"Vision pass: would skip — {vision_skip}", status="warn")

    # 4. Grobid reachability
    if cfg.grobid.disable:
        print_status(f"Grobid: disabled in config (header metadata will use --bib only)", status="warn")
    else:
        ok, detail = _ping_grobid(cfg.grobid.url)
        if ok:
            print_status(f"Grobid: reachable at {cfg.grobid.url}", status="ok")
        else:
            failures.append(
                f"Grobid not reachable at {cfg.grobid.url} ({detail}). "
                "Start it with: docker compose up -d grobid"
            )
            print_status(f"Grobid: unreachable at {cfg.grobid.url} ({detail})", status="fail")

    # 5. Output disk space (warn if < 5 GB free)
    output_dir = _resolve_against(config_path, cfg.output_dir)
    free_gb = _free_disk_gb(output_dir)
    if free_gb is None:
        print_status(f"output_dir disk: cannot stat {output_dir}", status="warn")
    elif free_gb < 5:
        print_status(f"output_dir disk: {free_gb:.1f} GB free (low — recommend ≥ 20 GB for a real corpus)", status="warn")
    else:
        print_status(f"output_dir disk: {free_gb:.1f} GB free", status="ok")

    # 6. input_pdfs reachable
    if cfg.input_pdfs is None:
        print_status("input_pdfs: not set (required for `corpus run`)", status="warn")
    else:
        input_dir = _resolve_against(config_path, cfg.input_pdfs)
        if input_dir.exists():
            n_pdfs = len(list(input_dir.rglob("*.pdf"))) if input_dir.is_dir() else 0
            print_status(f"input_pdfs: {input_dir} ({n_pdfs} PDFs)", status="ok")
        else:
            failures.append(f"input_pdfs path does not exist: {input_dir}")
            print_status(f"input_pdfs: {input_dir} not found", status="fail")

    if failures:
        print()
        print_status(f"{len(failures)} precondition(s) failed:", status="fail")
        for f in failures:
            print(f"  - {f}")
        return EXIT_PRECONDITION
    print()
    print_status("ready: `corpus run` should succeed on this host.", status="ok")
    return EXIT_OK


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


def _vision_skip_reason(backend: str) -> Optional[str]:
    """Return a human-readable reason to skip the vision pass on this host
    given the configured ``backend``, or None if the backend is usable.

    #65 — capability-aware skip. Gracefully avoids loading a multi-GB
    open-weights model on a CPU-only host or hitting the Anthropic API
    without credentials, both of which hard-fail Pass 3b deep into the
    pipeline today.
    """
    if backend == "local":
        if _detect_accelerator() is None:
            return (
                "vision.backend=local but no CUDA/MPS detected on this host. "
                "Either set `vision.backend: claude` (and export "
                "ANTHROPIC_API_KEY) or `vision.backend: none`, or run on a "
                "GPU host."
            )
    elif backend == "claude":
        if not os.environ.get("ANTHROPIC_API_KEY"):
            return (
                "vision.backend=claude but ANTHROPIC_API_KEY is unset. "
                "Either export ANTHROPIC_API_KEY, or set "
                "`vision.backend: local` / `vision.backend: none`."
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
    local verbs="init run status serve bib check completion"
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
        run)        COMPREPLY=( $(compgen -W "--dry-run --force-rebuild --no-vision --no-bundle --enrich-bhl -h --help" -- "$cur") ) ;;
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
            run)        _arguments '--dry-run' '--force-rebuild' '--no-vision' '--no-bundle' '--enrich-bhl' ;;
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
    """Wire root logging from -v/-vv/-q (#61)."""
    if args.quiet:
        level = logging.WARNING
    elif args.verbose >= 2:
        level = logging.DEBUG
    elif args.verbose == 1:
        level = logging.INFO
    else:
        level = logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


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
    run_p.add_argument("--no-vision", action="store_true",
                       help="Skip the vision pass (Pass 3b)")
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
    run_p.set_defaults(func=_cmd_run)

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
            "  --allow-unpublishable          bypass figure-licensing gate"
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


_PASSTHROUGH_VERBS = {"status", "serve", "bib"}


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
