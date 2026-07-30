"""Model prefetch + the new `corpus check` probes (#159, #160).

The point of both features is that a *missing* prerequisite is reported
at pre-flight rather than discovered mid-run, so the interesting
assertions are the negative ones: no models cached, no tesseract
language packs, no OCR binaries on PATH.

`pipeline.prefetch`'s reporting half must never touch the network —
`corpus check` calls it, and `corpus check` has to work on an air-gapped
HPC compute node (#139, #153). The tests below monkeypatch the cache
scan rather than mocking HTTP, so a regression that introduced a network
call would surface as a hang or an error here rather than silently
passing.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from pipeline import prefetch as pf

REPO_ROOT = Path(__file__).resolve().parent.parent


# --- pure reporting helpers ---------------------------------------------------


def test_required_repos_covers_docling_and_embedding():
    repos = pf.required_repos()
    # The two docling models a run always needs...
    assert "docling-project/docling-layout-heron" in repos
    assert "docling-project/docling-models" in repos
    # ...plus the embedding model.
    assert "BAAI/bge-m3" in repos
    # The ~16 GB local VLM is opt-in and must NOT be implied by default.
    assert pf.VISION_REPO not in repos


def test_vision_repo_included_only_on_request():
    assert pf.VISION_REPO in pf.required_repos(include_vision=True)


def test_embedding_model_override_is_honored():
    repos = pf.required_repos(embedding_model="BAAI/bge-base-en-v1.5")
    assert "BAAI/bge-base-en-v1.5" in repos
    assert "BAAI/bge-m3" not in repos


def test_model_report_all_absent(monkeypatch):
    """The fresh-install case: nothing cached, so nothing is present."""
    monkeypatch.setattr(pf, "cached_repos", lambda: {})
    rows, all_present = pf.model_report()
    assert all_present is False
    assert rows and all(r["cached"] is False for r in rows)
    assert all(r["size_bytes"] == 0 for r in rows)


def test_model_report_all_present(monkeypatch):
    monkeypatch.setattr(
        pf, "cached_repos",
        lambda: {r: 1024 for r in pf.required_repos()},
    )
    rows, all_present = pf.model_report()
    assert all_present is True
    assert all(r["cached"] for r in rows)


def test_model_report_partial(monkeypatch):
    """One model cached, the rest not — the common mid-download state."""
    monkeypatch.setattr(
        pf, "cached_repos", lambda: {"BAAI/bge-m3": 4_600_000_000},
    )
    rows, all_present = pf.model_report()
    assert all_present is False
    cached = [r for r in rows if r["cached"]]
    assert [r["repo"] for r in cached] == ["BAAI/bge-m3"]


def test_cached_repos_survives_an_unreadable_cache(monkeypatch):
    """A fresh install has no HF cache at all; that is not an error."""
    def _boom():
        raise OSError("no such cache")
    monkeypatch.setattr(
        "huggingface_hub.scan_cache_dir", lambda *a, **k: _boom(),
    )
    assert pf.cached_repos() == {}


def test_human_size_units():
    assert pf.human_size(0) == "0 B"
    assert pf.human_size(2048).endswith("KB")
    assert pf.human_size(5 * 1024 ** 3).endswith("GB")


def test_prefetch_retries_then_raises(monkeypatch):
    """A throttled host (HF 429) must retry and then fail loudly, not
    silently proceed to a run that will die mid-extract."""
    calls = {"n": 0}

    def _always_429():
        calls["n"] += 1
        raise RuntimeError("429 Too Many Requests")

    monkeypatch.setattr(pf.time, "sleep", lambda s: None)  # no real backoff
    with pytest.raises(RuntimeError, match="failed after 3 attempts"):
        pf._with_retry("thing", _always_429, attempts=3)
    assert calls["n"] == 3


def test_prefetch_succeeds_on_a_later_attempt(monkeypatch):
    calls = {"n": 0}

    def _flaky():
        calls["n"] += 1
        if calls["n"] < 3:
            raise RuntimeError("429 Too Many Requests")

    monkeypatch.setattr(pf.time, "sleep", lambda s: None)
    pf._with_retry("thing", _flaky, attempts=5)
    assert calls["n"] == 3


def test_docling_prefetch_does_not_warm_the_picture_classifier(monkeypatch):
    """#140 removed do_picture_classification from the real pipeline, so
    prefetch must not re-introduce the download it eliminated."""
    seen = {}

    class _FakeOpts:
        def __init__(self, **kw):
            seen.update(kw)

    class _FakeConv:
        def __init__(self, **kw):
            pass

        def convert(self, path):
            return None

    monkeypatch.setitem(
        sys.modules, "docling.datamodel.pipeline_options",
        type(sys)("docling.datamodel.pipeline_options"),
    )
    sys.modules["docling.datamodel.pipeline_options"].PdfPipelineOptions = _FakeOpts
    mod_dc = type(sys)("docling.document_converter")
    mod_dc.DocumentConverter = _FakeConv
    mod_dc.PdfFormatOption = lambda **kw: None
    monkeypatch.setitem(sys.modules, "docling.document_converter", mod_dc)
    mod_bm = type(sys)("docling.datamodel.base_models")

    class _IF:
        PDF = "pdf"
    mod_bm.InputFormat = _IF
    monkeypatch.setitem(sys.modules, "docling.datamodel.base_models", mod_bm)

    pf.prefetch_docling(Path("whatever.pdf"), attempts=1)
    assert seen.get("do_picture_classification") is False
    assert seen.get("do_table_structure") is True


# --- `corpus check` probes ----------------------------------------------------
#
# Driven as a subprocess against the demo corpuscle, because the probes
# are wired into _cmd_check's output and exit code, which is the contract
# a user (and CI) actually sees.


def _run_check(env_overrides: dict, expect_rc=None) -> subprocess.CompletedProcess:
    import os
    env = dict(os.environ)
    env.update(env_overrides)
    demo = REPO_ROOT / "demo"
    proc = subprocess.run(
        [sys.executable, "-m", "pipeline.cli", "check"],
        cwd=demo, env=env, capture_output=True, text=True, timeout=300,
    )
    if expect_rc is not None:
        assert proc.returncode == expect_rc, (
            f"rc={proc.returncode}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc


@pytest.mark.corpus_required
def test_check_warns_when_no_tesseract_language_packs(tmp_path):
    """#160: skipping `tools/install_tessdata.sh` must be a pre-flight
    warning, not a silently-wrong OCR result later."""
    empty = tmp_path / "empty_tessdata"
    empty.mkdir()
    proc = _run_check({"TESSDATA_PREFIX": str(empty)})
    out = proc.stdout + proc.stderr
    assert "Tesseract language packs" in out
    assert "install_tessdata.sh" in out, out
    # A missing pack degrades quality, it does not prevent a run — so it
    # must be a warning, never a precondition failure. Asserted on the
    # failure list rather than the exit code, because the exit code also
    # reflects Grobid reachability, which this test doesn't control.
    assert "precondition(s) failed" not in out or "tessdata" not in out.split(
        "precondition(s) failed")[1], out


@pytest.mark.corpus_required
def test_check_warns_when_no_models_cached(tmp_path):
    """#159: an empty model cache is reported without any network call."""
    proc = _run_check({"HF_HOME": str(tmp_path / "empty_hf")})
    out = proc.stdout + proc.stderr
    assert "Models:" in out
    assert "not in the local" in out, out
    assert "corpus prefetch" in out, out
    # Absent models are a warning: an online host just downloads mid-run.
    # Same reasoning as above for not asserting the exit code directly.
    assert "precondition(s) failed" not in out or "Models" not in out.split(
        "precondition(s) failed")[1], out


@pytest.mark.corpus_required
def test_check_reports_ocr_toolchain_ok_on_a_working_env():
    proc = _run_check({})
    out = proc.stdout + proc.stderr
    assert "OCR toolchain" in out
    assert "Models:" in out


# --- shell completions vs the real parser -------------------------------------


def test_completions_offer_every_verb():
    """The completion snippets are hand-maintained string constants, so they
    drift silently when a verb is added. `prefetch` (#159) exposed that
    `taxonomy` had been missing from all three for some time."""
    import argparse
    import re

    from pipeline.cli import (
        _BASH_COMPLETION, _FISH_COMPLETION, _ZSH_COMPLETION, _build_parser,
    )

    verbs = set()
    for action in _build_parser()._actions:
        if isinstance(action, argparse._SubParsersAction):
            verbs |= set(action.choices)
    assert "prefetch" in verbs  # guard against the test silently passing

    bash = set(re.search(r'local verbs="([^"]+)"', _BASH_COMPLETION).group(1).split())
    assert bash == verbs, f"bash completion drift: {verbs ^ bash}"

    for name, snippet in (("zsh", _ZSH_COMPLETION), ("fish", _FISH_COMPLETION)):
        missing = [v for v in verbs if v not in snippet]
        assert not missing, f"{name} completion missing verbs: {missing}"


# --- `corpus check` must not be more permissive than `corpus run` ------------


def test_check_fails_on_zero_pdfs_like_run_does(tmp_path):
    """`corpus check`'s promise is "can the next run succeed on this host?",
    and `corpus run` refuses a zero-PDF input dir outright. `check` used to
    report "ready" for exactly that condition — while run's refusal tells
    the user to "See `corpus check` for the full pre-flight surface", which
    was then less complete than the thing it deferred to.

    Found by the T4 operator walkthrough (dev_docs/clean_install_walkthrough.sh)
    on a freshly `corpus init`-ed corpuscle.
    """
    import os
    import textwrap

    corpuscle = tmp_path / "empty_corpuscle"
    (corpuscle / "pdfs").mkdir(parents=True)
    (corpuscle / "config.yaml").write_text(textwrap.dedent("""
        input_pdfs: ./pdfs
        output_dir: ./output
        grobid:
          disable: true
    """).lstrip())

    env = dict(os.environ)
    proc = subprocess.run(
        [sys.executable, "-m", "pipeline.cli", "check"],
        cwd=corpuscle, env=env, capture_output=True, text=True, timeout=300,
    )
    out = proc.stdout + proc.stderr
    assert "0 PDFs" in out, out
    assert proc.returncode == 3, (
        f"expected exit 3 (precondition), got {proc.returncode} — check is "
        f"greener than run:\n{out}"
    )
    # Same remedy wording as the orchestrator's guard, so the two agree.
    assert "--skip-checks" in out, out


# --- the prefetch code path must bind to the REAL backend API ----------------


def test_prefetch_calls_a_method_the_embedding_backend_actually_has():
    """`corpus prefetch` called `emb.encode(...)`, but the EmbeddingBackend
    ABC defines `embed`. It shipped because every unit test here mocks the
    backend, and T3 exercises the cold download through `corpus run` rather
    than through `corpus prefetch` — so nothing bound the call to the real
    API. The failure only appeared on a genuine cold prefetch, *after* the
    4.8 GB download had already succeeded.

    Asserted against the abstract base class, so it stays true for any
    backend, and without constructing one (which would download a model).
    """
    import inspect

    from pipeline import prefetch as pf
    from pipeline.embeddings import EmbeddingBackend

    src = inspect.getsource(pf.prefetch_embedding)
    called = {m for m in ("embed", "encode") if f"emb.{m}(" in src}
    assert called, "prefetch_embedding no longer touches the model at all"
    for method in called:
        assert hasattr(EmbeddingBackend, method), (
            f"prefetch_embedding calls emb.{method}(), which EmbeddingBackend "
            f"does not define. Available: "
            f"{sorted(m for m in vars(EmbeddingBackend) if not m.startswith('_'))}"
        )


def test_programming_errors_are_not_retried(monkeypatch):
    """Retrying an AttributeError costs ~5 min of backoff and then reports
    'failed after 6 attempts', which reads like a flaky network."""
    calls = {"n": 0}

    def _bug():
        calls["n"] += 1
        raise AttributeError("'LocalBackend' object has no attribute 'encode'")

    from pipeline import prefetch as pf
    monkeypatch.setattr(pf.time, "sleep", lambda s: None)
    with pytest.raises(RuntimeError, match="bug in corpus"):
        pf._with_retry("thing", _bug, attempts=6)
    assert calls["n"] == 1, "a programming error must fail on the first attempt"
