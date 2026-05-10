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

    # Vision pass: opt-out via --no-vision flag or vision.backend=none.
    # #65 adds capability detection to skip gracefully when the backend
    # isn't usable on this host.
    if not args.no_vision and cfg.vision.backend != "none":
        sub_argv += ["--vision-backend", cfg.vision.backend]
        if cfg.vision.model:
            sub_argv += ["--vision-model", cfg.vision.model]
    return sub_argv


def _cmd_run(args: argparse.Namespace) -> int:
    config_path = _resolve_config_path(args.config)
    if config_path is None:
        print(
            "corpus run: no config.yaml in cwd; run `corpus init` to scaffold one, "
            "or pass --config PATH",
            file=sys.stderr,
        )
        return 2
    if not config_path.exists():
        print_status(f"config not found: {config_path}", status="fail")
        return 2
    cfg = _load_validated(config_path)

    sub_argv = _build_orchestrator_argv(cfg, config_path, args)
    cmd = [sys.executable, "-m", "pipeline.orchestrator", *sub_argv]
    print_status(f"corpus run → {' '.join(cmd)}", status="info")
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        return rc

    # Success summary + bundle distillation are deferred to follow-up
    # work on #60 (run.log writer + invoking mcpsrv.bundle in line).
    # Status report and bundle remain explicit subcommands until then.
    print_status(
        "run complete. Try `corpus status` and `corpus serve` next.",
        status="ok",
    )
    return 0


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
    print_status(
        "corpus check is not yet implemented (#62). For now: validate "
        "config.yaml manually with `python -c \"import yaml; "
        "from pipeline.config_schema import validate_config; "
        "validate_config(yaml.safe_load(open('config.yaml')))\"`.",
        status="warn",
    )
    return 0


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
    """Read CITATION.cff from the repo root (editable / git+ install)."""
    repo_root = Path(__file__).resolve().parent.parent
    cff_path = repo_root / "CITATION.cff"
    if not cff_path.exists():
        print_status(
            f"CITATION.cff not found at {cff_path}. Available only in "
            "editable / git+ installs of corpus until PyPI publishing "
            "lands (PLAN §1).",
            status="fail",
        )
        sys.exit(EXIT_GENERIC)
    return yaml.safe_load(cff_path.read_text(encoding="utf-8"))


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
        description="MCP server for taxonomic literature corpora.",
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
    )
    init_p.add_argument("--path", type=Path, default=None)
    init_p.add_argument("--force", action="store_true")
    init_p.set_defaults(func=_cmd_init)

    # --- run ---
    run_p = sub.add_parser(
        "run",
        help="Full pipeline: extract → embed → cross-paper builds → bundle",
    )
    run_p.add_argument("--dry-run", action="store_true",
                       help="Plan only — no artifacts written")
    run_p.add_argument("--force-rebuild", action="store_true",
                       help="Disable implicit resume; re-run every stage")
    run_p.add_argument("--no-vision", action="store_true",
                       help="Skip the vision pass (Pass 3b)")
    run_p.add_argument("--no-bundle", action="store_true",
                       help="Skip the served-bundle distillation step "
                       "(deferred follow-up on #60)")
    run_p.add_argument("--enrich-bhl", action="store_true",
                       help="Enrich pre-DOI references against the Biodiversity "
                       "Heritage Library (slow, rate-limited; #64)")
    run_p.set_defaults(func=_cmd_run)

    # --- status ---
    status_p = sub.add_parser(
        "status",
        help="Build-state rollup + report (#57). Extra flags forwarded "
        "to pipeline.status (--report, --json, --sort-by, etc.).",
    )
    status_p.add_argument(
        "--output-dir", type=Path, default=None,
        help="Override the corpuscle output dir from config.yaml.",
    )
    status_p.set_defaults(func=_cmd_status)

    # --- serve ---
    serve_p = sub.add_parser(
        "serve",
        help="MCP server. Extra flags forwarded to mcpsrv.main "
        "(--transport, --host, --port, --auth-token-file, etc.).",
    )
    serve_p.add_argument("--output-dir", type=Path, default=None)
    serve_p.set_defaults(func=_cmd_serve)

    # --- bib (export | import) ---
    bib_p = sub.add_parser("bib", help="BibTeX round-trip")
    bib_sub = bib_p.add_subparsers(dest="bib_command", metavar="{export,import}")

    bib_export_p = bib_sub.add_parser("export")
    bib_export_p.add_argument("--output-dir", type=Path, default=None)
    bib_export_p.set_defaults(func=_cmd_bib_export)

    bib_import_p = bib_sub.add_parser("import")
    bib_import_p.add_argument("--output-dir", type=Path, default=None)
    bib_import_p.set_defaults(func=_cmd_bib_import)

    bib_p.set_defaults(func=lambda a: (bib_p.print_help() or 2))

    # --- check (stub) ---
    check_p = sub.add_parser(
        "check",
        help="Environment + config pre-flight (stub, #62)",
    )
    check_p.set_defaults(func=_cmd_check)

    # --- completion (#61) ---
    comp_p = sub.add_parser(
        "completion",
        help="Generate a shell completion script — pipe into your dotfiles "
        "(e.g. `corpus completion bash > ~/.local/share/corpus.bash; "
        "echo 'source ~/.local/share/corpus.bash' >> ~/.bashrc`)",
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
    args.passthrough = extras
    _setup_logging(args)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
