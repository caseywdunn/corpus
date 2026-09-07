"""Hard F821 gate: undefined names cannot be suppressed in the source tree.

A failure of this test is by definition a NameError waiting to fire at
runtime — Ruff's Pyflakes-compatible rules check syntax without executing code.
Two real instances surfaced in the platform-portability smoke iteration
that this test would have caught instantly:

  bib/authority.py:1471   undefined name 'db_path'
                          ('Could not apply staged bib overrides:
                          name "db_path" is not defined' on every run)

  mcpsrv/tools/chunks.py  undefined name 'json'   (translate_chunk
                          undefined name 'EmbeddingError'  embedding
                          error handler in get_chunks_for_topic)

The configured lint gate checks the complete ``F`` family. This deliberately
separate assertion runs isolated with ``--ignore-noqa`` so an accidental
``# noqa: F821`` cannot turn a guaranteed NameError into a green build (#259).

Run locally with::

    python -m pytest tests/test_no_undefined_names.py
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
# `tools/` is in the list even though it ships no importable package: its
# scripts are run by operators at release time, where a NameError costs a whole
# manual run rather than a fast test failure (#75, #193).
SOURCE_DIRS = ("pipeline", "mcpsrv", "bib", "tools")


def _which_ruff() -> str:
    """Ruff is a declared dev dependency; absence is a broken test env."""
    executable = shutil.which("ruff")
    assert executable is not None, (
        "ruff not on PATH — update the development environment from "
        "environment.yaml or install the project's dev extra"
    )
    return executable


def test_no_undefined_names():
    """Ruff reports zero F821 findings in every tracked source directory.

    An undefined name is a guaranteed runtime NameError — the only
    reason it doesn't already crash is that the affected code path
    hasn't been exercised by tests. Treat every finding as a hard
    failure; tests cannot retroactively cover every code path.
    """
    ruff = _which_ruff()
    targets = [str(REPO_ROOT / d) for d in SOURCE_DIRS]
    result = subprocess.run(
        [
            ruff, "check", "--isolated", "--select", "F821",
            "--ignore-noqa", *targets,
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        "Ruff reported undefined names — these will NameError at runtime. "
        "Fix or add the missing import; F821 suppressions are intentionally "
        "ignored by this gate.\n\n" + result.stdout + result.stderr
    )
