"""A configured-but-missing taxonomy must fail, not skip silently (#139).

The reported failure: the first full production run on Bouchet put 1763
papers through with empty `taxa.json` and an empty
`taxon_mentions.sqlite`, because `taxonomy.source: worms` needs outbound
internet that compute nodes don't have, and the per-paper stage merely
logged a WARNING and carried on.

Two layers now guard it, and both matter:

* `orchestrator._check_taxonomy_available` pre-flights before any work
  starts (already present).
* `pipeline.main --require-taxonomy` fails at the layer that would
  otherwise write the empty artifacts — which also covers direct
  `python -m pipeline.main` invocations and a snapshot that disappears
  after the pre-check.

The negative case is equally important: a corpuscle with no `taxonomy:`
block at all is a legitimate, documented configuration and must keep
working. `--no-taxa` and "no taxonomy configured" must not be turned into
errors by this change.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def _run_main(tmp_path, *extra_args):
    in_dir = tmp_path / "in"
    out_dir = tmp_path / "out"
    in_dir.mkdir(exist_ok=True)
    out_dir.mkdir(exist_ok=True)
    return subprocess.run(
        [sys.executable, "-m", "pipeline.main", str(in_dir), str(out_dir),
         *extra_args],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=300,
    )


def test_missing_taxonomy_with_require_flag_fails(tmp_path):
    proc = _run_main(
        tmp_path, "--require-taxonomy",
        "--taxonomy-db", str(tmp_path / "absent.sqlite"),
    )
    out = proc.stdout + proc.stderr
    assert proc.returncode == 1, out
    assert "refusing to run" in out, out
    # The message must be actionable — this is the whole point of #139.
    assert "--no-taxa" in out, out
    assert "dwca" in out, out


def test_missing_taxonomy_without_require_flag_only_warns(tmp_path):
    """A corpuscle with no taxonomy configured is a supported setup."""
    proc = _run_main(
        tmp_path, "--taxonomy-db", str(tmp_path / "absent.sqlite"),
    )
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0, out
    assert "taxon extraction skipped" in out, out
    assert "refusing to run" not in out


def test_no_taxa_wins_over_require_taxonomy(tmp_path):
    """An explicit --no-taxa must not be overridden into a failure — the
    orchestrator only passes --require-taxonomy when --no-taxa is absent,
    but the per-paper stage should be robust to receiving both."""
    proc = _run_main(
        tmp_path, "--no-taxa", "--require-taxonomy",
        "--taxonomy-db", str(tmp_path / "absent.sqlite"),
    )
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0, out
    assert "refusing to run" not in out


def test_orchestrator_passes_require_taxonomy_when_configured():
    """The flag is only useful if `corpus run` actually sends it. Asserted
    against the dry-run's echoed sub-command for the demo corpuscle, which
    configures taxonomy.source: dwca."""
    proc = subprocess.run(
        # --skip-checks: this asserts on argv construction, not on
        # pre-flight. Without it `corpus run --dry-run` probes Grobid and
        # exits 3 in the T0 tier, which has no Grobid service.
        [sys.executable, "-m", "pipeline.cli", "run", "--dry-run",
         "--skip-checks", "--only", "extract"],
        cwd=REPO_ROOT / "demo", capture_output=True, text=True, timeout=300,
    )
    out = proc.stdout + proc.stderr
    assert "-m pipeline.main" in out, out
    assert "--require-taxonomy" in out, out


def test_orchestrator_omits_require_taxonomy_with_no_taxa():
    proc = subprocess.run(
        [sys.executable, "-m", "pipeline.cli", "run", "--dry-run",
         "--skip-checks", "--only", "extract", "--no-taxa"],
        cwd=REPO_ROOT / "demo", capture_output=True, text=True, timeout=300,
    )
    out = proc.stdout + proc.stderr
    assert "--no-taxa" in out, out
    assert "--require-taxonomy" not in out, out
