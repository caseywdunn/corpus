"""Static-analysis gate: zero ``undefined name`` findings in the source tree.

A failure of this test is by definition a NameError waiting to fire at
runtime — pyflakes is a syntactic check that doesn't execute code.
Two real instances surfaced in the platform-portability smoke iteration
that this test would have caught instantly:

  bib/authority.py:1471   undefined name 'db_path'
                          ('Could not apply staged bib overrides:
                          name "db_path" is not defined' on every run)

  mcpsrv/tools/chunks.py  undefined name 'json'   (translate_chunk
                          undefined name 'EmbeddingError'  embedding
                          error handler in get_chunks_for_topic)

The same scan also turned up a handful of bare ``Dict``/``List``
references in modules that already use ``from __future__ import
annotations`` — harmless because the annotations are strings, but
noisy. This test deliberately fails *only* on outright undefined
names; cleanup of the unused-import + bare-annotation noise is a
separate quality-of-life item.

Run locally with::

    python -m pytest tests/test_no_undefined_names.py
"""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SOURCE_DIRS = ("pipeline", "mcpsrv", "bib")


# pyflakes' undefined-name line looks like:
#   path/to/file.py:LINE:COL: undefined name 'SYMBOL'
_UNDEF_RE = re.compile(r": undefined name '([^']+)'")


def _which_pyflakes() -> str:
    """Locate the pyflakes binary; skip the test if pyflakes isn't installed."""
    py = shutil.which("pyflakes")
    if py is None:
        pytest.skip(
            "pyflakes not on PATH — install it via the conda env "
            "(`conda env update -f environment.yaml`) or "
            "`pip install pyflakes`."
        )
    return py


def test_no_undefined_names():
    """``pyflakes`` reports zero ``undefined name`` lines in any tracked source dir.

    An undefined name is a guaranteed runtime NameError — the only
    reason it doesn't already crash is that the affected code path
    hasn't been exercised by tests. Treat every finding as a hard
    failure; tests cannot retroactively cover every code path.
    """
    pyflakes = _which_pyflakes()
    targets = [str(REPO_ROOT / d) for d in SOURCE_DIRS]
    result = subprocess.run(
        [pyflakes, *targets],
        capture_output=True,
        text=True,
    )
    # pyflakes prints every finding (one per line) to stdout and exits
    # nonzero when there's anything to report. We only care about
    # undefined names; other categories (unused imports, etc.) are
    # noisy cleanup, not bug classes.
    offending = [
        line for line in result.stdout.splitlines() if _UNDEF_RE.search(line)
    ]
    assert not offending, (
        "pyflakes reported undefined names — these will NameError at "
        "runtime the moment the affected code path executes. Fix or "
        "add the missing import.\n\n"
        + "\n".join(offending)
    )
