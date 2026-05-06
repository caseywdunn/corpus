"""Tests for the update_corpus.py orchestrator (#32)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


def _run(*argv: str, **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "update_corpus.py", *argv],
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


# ---------------------------------------------------------------------------
# --re-annotate-stale (#33)
# ---------------------------------------------------------------------------


import hashlib
import json


def _make_corpus(tmp_path: Path, anatomy_sha: str) -> Path:
    """Synthetic 3-paper corpus with anatomy.json files stamped with
    the given fingerprint sha. Two papers stale, one fresh.
    """
    inp = tmp_path / "in"
    inp.mkdir()
    out = tmp_path / "out"
    docs = out / "documents"
    docs.mkdir(parents=True)

    fresh_sha = anatomy_sha  # matches current
    stale_sha = hashlib.sha256(b"old version").hexdigest()

    for h, sha in [("aaa", fresh_sha), ("bbb", stale_sha), ("ccc", stale_sha)]:
        hd = docs / h
        hd.mkdir()
        (hd / "anatomy.json").write_text(json.dumps({
            "total_mentions": 0, "unique_terms": 0, "mentions": [],
            "input_fingerprint": {"path": "...", "sha256": sha, "size": 1},
        }))
    return out


def test_re_annotate_stale_requires_resume(tmp_path):
    inp = tmp_path / "in"
    inp.mkdir()
    lexicon = tmp_path / "lex.yaml"
    lexicon.write_text("terms: []")
    r = _run(
        str(inp), str(tmp_path / "out"),
        "--re-annotate-stale", "anatomy",
        "--anatomy-lexicon", str(lexicon),
    )
    assert r.returncode != 0
    assert "requires --resume" in r.stderr


def test_re_annotate_stale_dry_run_lists_but_does_not_delete(tmp_path):
    lexicon = tmp_path / "lex.yaml"
    lexicon.write_text("terms: ['nectophore']")
    sha = hashlib.sha256(lexicon.read_bytes()).hexdigest()
    out = _make_corpus(tmp_path, sha)
    inp = tmp_path / "in"

    r = _run(
        str(inp), str(out),
        "--resume", "--dry-run",
        "--re-annotate-stale", "anatomy",
        "--anatomy-lexicon", str(lexicon),
    )
    assert r.returncode == 0, r.stderr
    # Two stale papers — should report them
    assert "2 anatomy_stale" in r.stderr
    assert "would delete" in r.stderr
    # Files still present (dry-run)
    for h in ("aaa", "bbb", "ccc"):
        assert (out / "documents" / h / "anatomy.json").exists()


def _call_delete_stale(*, output_dir: Path, anatomy_lexicon: Path,
                       dry_run: bool = False, mode: str = "anatomy",
                       taxonomy_db: Path = None):
    """Invoke ``update_corpus._delete_stale_artifacts`` directly.

    The two tests that exercise the deletion logic need to verify
    on-disk effects without running the rest of the wrapper, since
    --re-annotate-stale + --skip-pipeline is now mutually exclusive
    and the rest of the post-pipeline chain pulls in optional deps
    (lancedb etc.) that aren't required just to test deletion.
    """
    import argparse
    import sys
    sys.path.insert(0, str(REPO_ROOT))
    try:
        import update_corpus as uc
    finally:
        sys.path.pop(0)
    args = argparse.Namespace(
        output_dir=output_dir,
        anatomy_lexicon=anatomy_lexicon,
        taxonomy_db=taxonomy_db,
        re_annotate_stale=mode,
        dry_run=dry_run,
    )
    return uc._delete_stale_artifacts(args)


def test_re_annotate_stale_deletes_only_stale_artifacts(tmp_path):
    lexicon = tmp_path / "lex.yaml"
    lexicon.write_text("terms: ['nectophore']")
    sha = hashlib.sha256(lexicon.read_bytes()).hexdigest()
    out = _make_corpus(tmp_path, sha)

    rc = _call_delete_stale(output_dir=out, anatomy_lexicon=lexicon)
    assert rc == 0
    # aaa was fresh — anatomy.json must remain
    assert (out / "documents" / "aaa" / "anatomy.json").exists()
    # bbb + ccc were stale — anatomy.json must be gone
    assert not (out / "documents" / "bbb" / "anatomy.json").exists()
    assert not (out / "documents" / "ccc" / "anatomy.json").exists()


def test_re_annotate_stale_no_change_when_nothing_stale(tmp_path):
    lexicon = tmp_path / "lex.yaml"
    lexicon.write_text("terms: ['x']")
    sha = hashlib.sha256(lexicon.read_bytes()).hexdigest()

    out = tmp_path / "out"
    docs = out / "documents"
    docs.mkdir(parents=True)
    hd = docs / "fresh"
    hd.mkdir()
    (hd / "anatomy.json").write_text(json.dumps({
        "total_mentions": 0, "unique_terms": 0, "mentions": [],
        "input_fingerprint": {"sha256": sha, "size": 1, "path": "..."},
    }))

    rc = _call_delete_stale(output_dir=out, anatomy_lexicon=lexicon)
    assert rc == 0
    assert (hd / "anatomy.json").exists()
