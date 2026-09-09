"""`compute.accelerator: require` — fail instead of degrading to CPU (#270).

`auto` falling back to CPU is right on a workstation, where the
alternative is every docling page dying with "no kernel image is
available". Inside a scheduler allocation it is the wrong answer and an
expensive one: a 2026-08-31 embed job resolved to CPU on an allocated
RTX 5000 Ada, logged one WARNING into a stderr stream that is mostly
HuggingFace chatter, ran 181 of 1,775 documents in 75 minutes against a
4-hour wall, and was cancelled by the cluster's 0%-GPU-utilization
policy. Nothing distinguished it from a slow-but-fine run — it logged
steady per-document progress the whole time.

Pinning `cuda` was not a way to say this: a pinned value is honoured
verbatim, which *disables* the capability check rather than enforcing it.
"""
from __future__ import annotations

import argparse

import pytest

from pipeline.accelerator import AcceleratorUnavailable, resolve_device
from pipeline.orchestrator import _GPU_STEPS, STEPS, _report_accelerator


def _unusable(monkeypatch, reason="RTX 5000 Ada has CUDA compute capability "
                                 "8.9, ... Falling back to CPU."):
    monkeypatch.setattr("pipeline.accelerator.unsupported_cuda_reason",
                        lambda: reason)


def _no_gpu(monkeypatch):
    monkeypatch.setattr("pipeline.accelerator.unsupported_cuda_reason",
                        lambda: None)
    monkeypatch.setattr("pipeline.accelerator._detect_device", lambda: "cpu")


# --- resolve_device -----------------------------------------------------


def test_require_raises_on_a_visible_but_unusable_gpu(monkeypatch):
    _unusable(monkeypatch)
    with pytest.raises(AcceleratorUnavailable) as exc:
        resolve_device("require")
    message = str(exc.value)
    # The diagnosis unsupported_cuda_reason already builds is the whole
    # value of failing here rather than letting torch throw per launch.
    assert "compute capability 8.9" in message
    # ...but not its "Falling back to CPU", which is what is *not* happening.
    assert "Falling back to CPU" not in message
    assert "failure rather than a fallback" in message


def test_require_raises_when_there_is_no_gpu_at_all(monkeypatch):
    """A different cause needing a different action: this is a submission
    problem, not a torch/hardware mismatch."""
    _no_gpu(monkeypatch)
    with pytest.raises(AcceleratorUnavailable) as exc:
        resolve_device("require")
    message = str(exc.value)
    assert "no CUDA or MPS device" in message
    assert "did not request one" in message


def test_require_returns_the_device_when_one_is_usable(monkeypatch):
    monkeypatch.setattr("pipeline.accelerator.unsupported_cuda_reason",
                        lambda: None)
    monkeypatch.setattr("pipeline.accelerator._detect_device", lambda: "cuda")
    assert resolve_device("require") == "cuda"


def test_auto_still_falls_back(monkeypatch):
    """The workstation behaviour #198 added must not change."""
    _unusable(monkeypatch)
    assert resolve_device("auto") == "cpu"


def test_require_is_a_valid_config_value():
    from pipeline.config_schema import ComputeConfig
    assert ComputeConfig(accelerator="require").accelerator == "require"
    with pytest.raises(Exception):
        ComputeConfig(accelerator="gpu-please")


# --- the run banner and the abort ---------------------------------------


def _args(require_gpu=False):
    return argparse.Namespace(require_gpu=require_gpu)


def _steps(*names):
    by_name = {s.name: s for s in STEPS}
    return [by_name[n] for n in names]


def test_a_cpu_fallback_is_reported_in_the_first_screen(monkeypatch, caplog):
    """The 2026-08-31 warning existed; it was ~200 lines into stderr. This
    one is next to the step list, before any step runs."""
    _unusable(monkeypatch)
    with caplog.at_level("INFO"):
        assert _report_accelerator(_steps("embed")) is None
    assert "Accelerator: cpu" in caplog.text
    assert "GPU steps will run on CPU" in caplog.text
    assert "#270" in caplog.text


def test_require_gpu_aborts_before_any_step_runs(monkeypatch):
    _unusable(monkeypatch)
    err = _report_accelerator(_steps("embed"), require=True)
    assert err is not None
    assert "No usable GPU" in err


def test_a_run_with_no_gpu_step_is_not_blocked(monkeypatch):
    """`--only post` is SQLite and JSON. Failing it for want of a GPU would
    make the flag unusable in the chained-job layout the SLURM scripts use."""
    _unusable(monkeypatch)
    steps = _steps("build_biblio", "build_taxa")
    assert not any(s.name in _GPU_STEPS for s in steps)
    assert _report_accelerator(steps, require=True) is None


def test_require_gpu_says_so_when_a_pin_makes_it_inert(monkeypatch, caplog):
    """A pinned device is honoured verbatim, so the flag cannot enforce
    anything on top of it. Silently doing nothing is the failure mode this
    whole issue is about."""
    monkeypatch.setitem(
        __import__("pipeline.config", fromlist=["CONFIG"]).CONFIG,
        "compute", {"accelerator": "cuda"},
    )
    with caplog.at_level("WARNING"):
        _report_accelerator(_steps("embed"), require=True)
    assert "--require-gpu has no effect" in caplog.text


# --- wiring -------------------------------------------------------------


def test_the_gpu_steps_are_the_ones_that_use_a_gpu():
    assert _GPU_STEPS == {"extract", "vision", "embed"}


def test_prefetch_downgrades_require_so_a_login_node_can_warm_the_cache():
    """Prefetch runs where there is *network* access, typically a login node
    with no GPU. Failing it would make the flag a hazard."""
    import inspect

    from pipeline import prefetch
    src = inspect.getsource(prefetch)
    assert '"auto" if _configured == "require" else _configured' in src


def test_the_gpu_slurm_scripts_ask_for_the_failure():
    """A GPU allocation is a statement that CPU is not acceptable."""
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent
    for script in ("batch_embed.sh", "batch_pass3b.sh"):
        text = (root / "slurm" / script).read_text()
        assert "--require-gpu" in text, script


def test_corpus_run_forwards_the_flag():
    import inspect

    from pipeline import cli
    src = inspect.getsource(cli)
    assert 'sub_argv.append("--require-gpu")' in src


def test_the_run_args_fixtures_cover_every_flag_the_builder_reads():
    """A guard for how adding this flag broke twelve unrelated tests.

    `_build_orchestrator_argv` reads `args.<flag>` directly, and two test
    modules hand-build the Namespace it receives. Adding a flag to the real
    parser leaves those fixtures a field short, and the failure surfaces as
    an AttributeError in tests about figure panels and HPC batching — a
    long way from the change that caused it. Assert the coupling instead.
    """
    import inspect
    import re

    from pipeline import cli
    from tests import test_corpus_run_hpc, test_figure_panels_flag

    src = inspect.getsource(cli._build_orchestrator_argv)
    # A name reached through `getattr(args, "x", default)` anywhere in the
    # function is already defended — the author made it optional on purpose,
    # and the bare `args.x` that follows such a guard only runs once it has
    # passed. Only unguarded attribute access couples the fixture.
    defended = set(re.findall(r"getattr\(\s*args\s*,\s*[\"']([a-z_]+)", src))
    read = set(re.findall(r"\bargs\.([a-z_]+)", src)) - defended
    assert "require_gpu" in read, "guard is checking the wrong function"
    assert defended, "expected the HPC flags to be read defensively"
    for module in (test_figure_panels_flag, test_corpus_run_hpc):
        provided = set(vars(module._run_args()))
        missing = read - provided
        assert not missing, (
            f"{module.__name__}._run_args() is missing {sorted(missing)}, "
            f"which _build_orchestrator_argv reads"
        )
