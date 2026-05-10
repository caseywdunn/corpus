"""`corpus` CLI entry point (#58 packaging, #59 init, #60 router).

Currently dispatches `corpus init` (#59) and prints help / version.
The remaining verbs (`run`, `check`, `status`, `serve`, `bib`,
`completion`) land across #60–#62. This file grows as those issues
close; the dispatch table is the only thing that changes.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from importlib import resources
from pathlib import Path

from .version import __version__


_TEMPLATE_PACKAGE = "pipeline"
_TEMPLATE_NAME = "config.template.yaml"


def _cmd_init(args: argparse.Namespace) -> int:
    """Scaffold a fresh ``config.yaml`` in cwd from the bundled template (#59)."""
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
    print("then run `corpus run` (lands in #60).")
    return 0


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="corpus",
        description="MCP server for taxonomic literature corpora.",
    )
    p.add_argument("-V", "--version", action="version", version=f"corpus {__version__}")
    sub = p.add_subparsers(dest="command", metavar="<command>")

    init_p = sub.add_parser(
        "init",
        help="Scaffold config.yaml in cwd from the bundled template",
    )
    init_p.add_argument(
        "--path",
        type=Path,
        default=None,
        help="Where to write the config (default: ./config.yaml in cwd)",
    )
    init_p.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing config.yaml",
    )
    init_p.set_defaults(func=_cmd_init)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
