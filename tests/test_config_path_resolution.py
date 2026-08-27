"""`--config` with a relative path must survive the trip to a subprocess (#210).

Reported by @ejedwards against 1.1.1. `corpus --config <relative-path> run`
silently discarded every system-tuning block — `ocr`, `chunking`,
`stage_timeouts`, `huge_document`, `quality_gates` — while the run looked
entirely normal.

Two independent causes, both fixed here:

1. `_resolve_config_path` returned the flag value verbatim. It is forwarded to
   each stage subprocess as `--config <path>`, and the orchestrator runs those
   from `REPO_ROOT`, not the operator's directory — so a relative path missed.
2. `load_config` treats a missing file as "use built-in defaults", which is
   right for the implicit `./config.yaml` and wrong for a file the operator
   named. The subprocess therefore fell back silently, leaving one INFO line
   as the only trace.

`input_pdfs`, `bib`, `lexicon` and `taxonomy` survived either way, because the
CLI resolves those itself and passes them as absolute arguments — which is
what made the failure so quiet, and what made it mask #209: an `ocr.jobs`
setting appeared to have no effect because the file carrying it was never read.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from pipeline.cli import _resolve_config_path

REPO_ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture
def cfg_dir(tmp_path, monkeypatch):
    (tmp_path / "tuned.yaml").write_text(
        "output_dir: ./output\nocr:\n  tesseract_page_timeout: 111\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    return tmp_path


# --- resolution ---------------------------------------------------------------


def test_a_relative_flag_becomes_absolute(cfg_dir):
    """The whole bug in one assertion: the returned path is handed to a
    subprocess with a different working directory."""
    got = _resolve_config_path(Path("tuned.yaml"))
    assert got.is_absolute()
    assert got == cfg_dir / "tuned.yaml"


def test_the_resolved_path_is_findable_from_the_subprocess_cwd(cfg_dir):
    """`orchestrator` runs each stage with `cwd=REPO_ROOT`. A path that only
    exists relative to the operator's directory is the failure."""
    got = _resolve_config_path(Path("tuned.yaml"))
    assert (REPO_ROOT / got).exists() or got.exists()
    # And explicitly: the *unresolved* form would not have been findable.
    assert not (REPO_ROOT / "tuned.yaml").exists()


def test_the_env_var_is_resolved_too(cfg_dir, monkeypatch):
    """`CORPUS_CONFIG` had the same hole and is reached by the same code path."""
    monkeypatch.setenv("CORPUS_CONFIG", "tuned.yaml")
    got = _resolve_config_path(None)
    assert got.is_absolute() and got == cfg_dir / "tuned.yaml"


def test_a_tilde_path_is_expanded(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    (tmp_path / "c.yaml").write_text("output_dir: ./o\n", encoding="utf-8")
    got = _resolve_config_path(Path("~/c.yaml"))
    assert got.is_absolute() and got.exists()


def test_the_implicit_default_is_absolute(cfg_dir):
    (cfg_dir / "config.yaml").write_text("output_dir: ./output\n", encoding="utf-8")
    got = _resolve_config_path(None)
    assert got.is_absolute() and got.name == "config.yaml"


def test_no_config_anywhere_is_still_none(tmp_path, monkeypatch):
    """Absent config is not an error here — `corpus init` relies on it."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("CORPUS_CONFIG", raising=False)
    assert _resolve_config_path(None) is None


# --- the silent fallback ------------------------------------------------------


def test_a_named_config_that_does_not_exist_is_refused(tmp_path):
    """The second cause. Even with paths resolved, a typo would otherwise
    replace every tuned setting with defaults and run anyway."""
    proc = subprocess.run(
        [sys.executable, "-m", "pipeline.main", str(tmp_path), str(tmp_path),
         "--config", str(tmp_path / "nope.yaml")],
        capture_output=True, text=True, cwd=REPO_ROOT,
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT)},
    )
    assert proc.returncode != 0
    combined = proc.stdout + proc.stderr
    assert "no such file" in combined
    # The message has to say *why* it refuses, or it reads as pedantry.
    assert "silently replaced" in combined


def test_an_absent_implicit_config_still_uses_defaults():
    """`load_config(None)` must keep falling back — only a *named* file that
    is missing is an error."""
    from pipeline.config import load_config
    cfg = load_config(Path("/nonexistent/config.yaml"))
    assert isinstance(cfg, dict) and cfg, "defaults must still be returned"
