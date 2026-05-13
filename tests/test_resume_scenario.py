"""T3 — 4 + 1 implicit-resume scenario (issue #75).

Pytest port of ``dev_docs/clean_install_walkthrough.sh`` §8. Builds the
demo corpuscle from scratch on its 4 bundled PDFs, then copies
``tests/fixtures/round2_paper/Siebert_etal2011.pdf`` into ``demo/`` and
re-runs ``corpus run``. The second run must skip everything from
round 1 and embed only the new paper — the canonical regression check
for implicit-resume (#71 broke this; pyflakes / unit tests can't catch
it because it needs a real LanceDB on disk to surface).

Marker: ``@pytest.mark.resume_scenario``. Deselected by default so the
local ``pytest`` run stays fast; T3 in ``.github/workflows/integration.yml``
opts in with ``pytest -m resume_scenario tests/test_resume_scenario.py``.

Wall time: ~2 min on Linux CI runners (round 1 ~80s, round 2 ~30s).
Requires Grobid reachable at ``http://localhost:8070`` (the demo's
config.yaml points at the default port).
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.resume_scenario

REPO_ROOT = Path(__file__).resolve().parent.parent
DEMO_DIR = REPO_ROOT / "demo"
ROUND2_PDF = REPO_ROOT / "tests" / "fixtures" / "round2_paper" / "Siebert_etal2011.pdf"

# pipeline/embed.py prints exactly this line at the end of the embed stage:
#   "Done. embedded=%d, skipped=%d, failed=%d"
_EMBED_TALLY_RE = re.compile(
    r"Done\.\s+embedded=(\d+),\s+skipped=(\d+),\s+failed=(\d+)"
)


def _grobid_alive(url: str = "http://localhost:8070/api/isalive") -> bool:
    import urllib.error
    import urllib.request

    try:
        with urllib.request.urlopen(url, timeout=2) as resp:
            return resp.read().strip() == b"true"
    except (urllib.error.URLError, OSError):
        return False


def _run_corpus(cwd: Path) -> subprocess.CompletedProcess:
    """Invoke ``corpus run --no-vision`` in ``cwd`` and capture combined output."""
    return subprocess.run(
        ["corpus", "-v", "run", "--no-vision"],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def _parse_embed_tally(output: str) -> tuple[int, int, int]:
    """Return (embedded, skipped, failed) from the last embed-stage tally line."""
    matches = _EMBED_TALLY_RE.findall(output)
    assert matches, (
        "embed-stage tally line not found in run output. "
        "Expected a line like 'Done. embedded=N, skipped=N, failed=N'."
    )
    last = matches[-1]
    return int(last[0]), int(last[1]), int(last[2])


def _load_manifest(output_dir: Path) -> dict:
    manifest_path = output_dir / "_serve" / "bundle_manifest.json"
    assert manifest_path.exists(), f"bundle_manifest.json missing at {manifest_path}"
    return json.loads(manifest_path.read_text())


@pytest.fixture
def clean_demo_workspace(tmp_path):
    """Copy demo/ into a tmp_path and yield it; tear down on exit.

    A temp copy keeps the source tree clean (the round-2 PDF would
    otherwise have to be removed from demo/ after every run) and lets
    parallel CI matrix legs share the same working directory model.
    """
    workspace = tmp_path / "demo"
    shutil.copytree(DEMO_DIR, workspace, ignore=shutil.ignore_patterns("output"))
    yield workspace


@pytest.mark.skipif(
    shutil.which("corpus") is None,
    reason="corpus CLI not on PATH; run `pip install -e .` first",
)
@pytest.mark.skipif(
    not _grobid_alive(),
    reason="Grobid not reachable at http://localhost:8070 — start it with "
           "`docker compose up -d grobid` before running this test",
)
@pytest.mark.skipif(
    not ROUND2_PDF.exists(),
    reason=f"round-2 fixture missing at {ROUND2_PDF}",
)
def test_resume_scenario(clean_demo_workspace):
    """Round-1 build on 4 papers, round-2 build adds the 5th and skips the rest."""
    workspace = clean_demo_workspace
    output_dir = workspace / "output"

    # ── Round 1: 4 papers, cold build. ────────────────────────────
    initial_pdfs = sorted(p.name for p in workspace.glob("*.pdf"))
    assert len(initial_pdfs) == 4, (
        f"demo/ should ship 4 PDFs; found {len(initial_pdfs)}: {initial_pdfs}"
    )
    assert ROUND2_PDF.name not in initial_pdfs, (
        "round-2 fixture leaked into demo/ — check tests/fixtures/round2_paper/"
    )

    r1 = _run_corpus(workspace)
    if r1.returncode != 0:
        pytest.fail(
            f"corpus run (round 1) exited {r1.returncode}\n"
            f"--- stdout (tail) ---\n{r1.stdout[-2000:]}\n"
            f"--- stderr (tail) ---\n{r1.stderr[-2000:]}"
        )
    r1_log = r1.stdout + r1.stderr
    embedded, skipped, failed = _parse_embed_tally(r1_log)
    assert (embedded, skipped, failed) == (4, 0, 0), (
        f"round 1: expected embedded=4 skipped=0 failed=0, "
        f"got embedded={embedded} skipped={skipped} failed={failed}"
    )

    m1 = _load_manifest(output_dir)
    assert m1["paper_count"] == 4, f"round 1 paper_count = {m1['paper_count']}"
    r1_chunks = m1["chunk_count"]
    assert r1_chunks > 0, "round 1 chunk_count should be positive"

    # Sanity: the round-2 PDF is held OUTSIDE input_pdfs, so it must
    # not have been picked up on round 1. This guards against a future
    # refactor that accidentally widens demo/'s input_pdfs scope.
    documents_dir = output_dir / "documents"
    round1_doc_dirs = {d.name for d in documents_dir.iterdir() if d.is_dir()}
    assert len(round1_doc_dirs) == 4, (
        f"round 1: expected 4 document dirs, got {len(round1_doc_dirs)}"
    )

    # ── Round 2: add the held-back paper, re-run. ─────────────────
    shutil.copy(ROUND2_PDF, workspace / ROUND2_PDF.name)
    assert (workspace / ROUND2_PDF.name).exists()

    r2 = _run_corpus(workspace)
    if r2.returncode != 0:
        pytest.fail(
            f"corpus run (round 2) exited {r2.returncode}\n"
            f"--- stdout (tail) ---\n{r2.stdout[-2000:]}\n"
            f"--- stderr (tail) ---\n{r2.stderr[-2000:]}"
        )
    r2_log = r2.stdout + r2.stderr
    embedded, skipped, failed = _parse_embed_tally(r2_log)
    assert (embedded, skipped, failed) == (1, 4, 0), (
        f"round 2: expected embedded=1 skipped=4 failed=0, "
        f"got embedded={embedded} skipped={skipped} failed={failed}\n"
        f"This is the canonical #71 (lancedb 0.30.x list_tables shape) "
        f"regression — a green round 1 followed by round 2 re-embedding "
        f"all 5 papers means implicit-resume against an existing LanceDB "
        f"is broken."
    )

    m2 = _load_manifest(output_dir)
    assert m2["paper_count"] == 5, f"round 2 paper_count = {m2['paper_count']}"
    assert m2["chunk_count"] > r1_chunks, (
        f"round 2 chunk_count ({m2['chunk_count']}) should exceed "
        f"round 1 ({r1_chunks}) — the new paper produced no chunks"
    )

    # #70 audit clean — the path-scrub line is the proof that the
    # bundle-distillation scrubbers covered every JSON written this
    # round, including the round-2 additions.
    assert "audit clean" in r2_log, (
        "round 2 run log missing 'audit clean' from the served-bundle "
        "path-scrub stage (#70 regression — a new JSON path is leaking "
        "absolute paths into the bundle)"
    )
