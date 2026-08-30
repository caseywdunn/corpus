"""Nothing junk at the top level of the repo (#257).

A 9-byte file named `6` — the `%PDF-1.4` header a test stub writes — sat
next to `README.md`, tracked, through four green CI runs and two release
PRs. It was caught by eye at the v1.2.0 tag boundary, one merge away from
being in a Zenodo archive permanently, because the GitHub release archives
whatever is in the tree under a DOI that `CITATION.cff`'s concept DOI
resolves to.

Nothing looked. T0 lints `pipeline/`, `mcpsrv/`, `bib/` and `tools/` for
undefined names; no check had an opinion about the top level. Same shape as
the `tools/` pyflakes gap (#75, #193): the check that would have caught it
did not exist because nobody had been bitten yet.

The name varies by machine, which is what makes an allowlist the right
shape rather than a pattern match — it is the resolved OCR job count, so it
was `6` here and `3` on the 4-CPU allocation that reported #257.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent

# Every file that belongs at the top level. Adding one here should be a
# deliberate act with a reviewer, which is the entire point.
ALLOWED = {
    ".gitignore",
    ".mcp.json",
    "AGENTS.md",
    "CHANGELOG.md",
    "CITATION.cff",
    "CLAUDE.md",
    "CONTRIBUTING.md",
    "DEPLOY.md",
    "INSTALL.md",
    "LICENSE",
    "README.md",
    "docker-compose.yml",
    "environment.yaml",
    "pyproject.toml",
    "pytest.ini",
    "requirements.txt",
}


def _tracked_root_files() -> set[str]:
    out = subprocess.run(
        ["git", "ls-files", "-z"], cwd=REPO,
        capture_output=True, text=True, check=True,
    ).stdout
    return {p for p in out.split("\0") if p and "/" not in p}


@pytest.mark.skipif(
    not (REPO / ".git").exists(), reason="not a git checkout (sdist/wheel)"
)
def test_no_stray_files_are_tracked_at_the_repo_root():
    unexpected = _tracked_root_files() - ALLOWED
    assert not unexpected, (
        f"unexpected tracked file(s) at the repo root: {sorted(unexpected)}. "
        "If one of these belongs, add it to ALLOWED with a reason; otherwise "
        "`git rm` it before it reaches a release archive (#257)."
    )


@pytest.mark.skipif(
    not (REPO / ".git").exists(), reason="not a git checkout (sdist/wheel)"
)
def test_the_allowlist_has_not_rotted():
    """A name removed from the repo but left here would silently stop
    guarding anything."""
    stale = ALLOWED - _tracked_root_files()
    assert not stale, f"ALLOWED names files that no longer exist: {sorted(stale)}"
