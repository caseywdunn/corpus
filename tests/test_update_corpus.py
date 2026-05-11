"""Tests for the orchestrator (#32; moved from ``update_corpus.py`` to
``pipeline.orchestrator`` per #60)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


def _run(*argv: str, **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "pipeline.orchestrator", *argv],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        **kwargs,
    )


def test_help_lists_all_steps():
    r = _run("--help")
    assert r.returncode == 0
    for step in ("extract", "embed", "build_biblio", "build_taxa",
                 "backfill_intext", "reconcile"):
        assert step in r.stdout


def test_skip_pipeline_and_from_are_mutually_exclusive(tmp_path):
    inp = tmp_path / "in"
    inp.mkdir()
    out = tmp_path / "out"
    r = _run(str(inp), str(out), "--skip-pipeline", "--from", "build_biblio")
    assert r.returncode != 0
    assert "mutually exclusive" in r.stderr


def test_skip_pipeline_only_runs_post_pipeline(tmp_path):
    inp = tmp_path / "in"
    inp.mkdir()
    out = tmp_path / "out"
    (out / "documents").mkdir(parents=True)
    r = _run(str(inp), str(out), "--skip-pipeline", "--dry-run")
    # Should NOT see the pipeline steps in the run plan
    assert "extract" not in r.stderr.split("Running")[-1].split("═══")[0]
    # Should see the post-pipeline steps
    for step in ("build_biblio", "build_taxa", "backfill_intext", "reconcile"):
        assert step in r.stderr


def test_from_build_biblio_skips_extract_and_embed(tmp_path):
    inp = tmp_path / "in"
    inp.mkdir()
    out = tmp_path / "out"
    (out / "documents").mkdir(parents=True)
    r = _run(str(inp), str(out), "--from", "build_biblio", "--dry-run")
    plan_section = r.stderr.split("Running")[-1].split("═══")[0]
    assert "extract" not in plan_section
    assert "embed" not in plan_section
    assert "build_biblio" in plan_section


def test_dry_run_completes_with_warnings_on_empty_corpus(tmp_path):
    """An empty output dir + dry-run should not abort — later steps
    legitimately fail (e.g. reconcile needs the SQLite that build_biblio
    would have written), but dry-run downgrades those to warnings.
    """
    inp = tmp_path / "in"
    inp.mkdir()
    # Use a fake PDF so process_corpus has something to discover but
    # nothing to actually open downstream
    out = tmp_path / "out"
    r = _run(str(inp), str(out), "--dry-run", "--resume")
    assert r.returncode == 0
    assert "warning(s)" in r.stderr or "succeeded" in r.stderr


def test_missing_input_dir_errors(tmp_path):
    r = _run(str(tmp_path / "nonexistent"), str(tmp_path / "out"), "--dry-run")
    assert r.returncode == 1
    assert "does not exist" in r.stderr


# --re-annotate-stale was removed in v0.2: pipeline_state.json
# (per-stage completion records stamped with PIPELINE_VERSION +
# input_fingerprint) lets a plain --resume re-run only the stages
# whose fingerprint disagrees with the current input. The flag and
# its tests are no longer needed.
