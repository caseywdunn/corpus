"""Accelerator capability checking (#198).

`torch.cuda.is_available()` answers "is there a visible NVIDIA GPU", not "can
this torch build run kernels on it". The two differ on exactly the hardware a
lab workstation has, and when they differ docling fails every page rather than
falling back:

    2026-08-06  Accelerator device: 'cpu'      -> clean build
    2026-08-26  Accelerator device: 'cuda:0'   -> CUDA error: no kernel image
                                                  is available for execution

Same machine, same pinned `torch==2.12.0`/`docling==2.94.0`, 20 days apart. A
GTX 1080 (compute capability 6.1) became visible to torch, whose pinned build
ships `sm_75` and up. Nothing in the project changed — which is the point, and
what the #98 pins exist to prevent.
"""
from __future__ import annotations

from pipeline.accelerator import (
    cuda_capability_supported,
    resolve_device,
    unsupported_cuda_reason,
)

# What the pinned torch build actually reports.
PINNED = ["sm_75", "sm_80", "sm_86", "sm_90", "sm_100", "sm_120"]


def test_the_gtx_1080_case_that_prompted_this():
    """Capability 6.1 against a build starting at sm_75 — the observed failure."""
    assert cuda_capability_supported((6, 1), PINNED) is False


def test_a_supported_card_is_not_rejected():
    """The check must not push working hardware onto the CPU."""
    assert cuda_capability_supported((8, 0), PINNED) is True
    assert cuda_capability_supported((12, 0), PINNED) is True


def test_an_exact_binary_kernel_counts():
    assert cuda_capability_supported((6, 1), ["sm_61"]) is True


def test_ptx_jits_forward_but_never_backward():
    """`compute_75` JITs onto anything at or above 7.5, and nothing below it.

    Both directions matter: allowing only exact `sm_` matches would reject
    cards a PTX-carrying build serves fine, and allowing PTX backwards would
    re-admit the very failure this module exists to prevent.
    """
    assert cuda_capability_supported((8, 6), ["compute_75"]) is True
    assert cuda_capability_supported((7, 5), ["compute_75"]) is True
    assert cuda_capability_supported((6, 1), ["compute_75"]) is False


def test_architecture_suffixes_are_tolerated():
    """torch reports `sm_90a` for architecture-specific builds."""
    assert cuda_capability_supported((9, 0), ["sm_90a"]) is True


def test_a_cpu_only_torch_reports_no_architectures():
    """Empty arch list means a CPU-only build; callers fall through to CPU."""
    assert cuda_capability_supported((6, 1), []) is False


def test_unparseable_entries_are_ignored_not_crashed_on():
    assert cuda_capability_supported((8, 0), ["nonsense", "sm_80"]) is True
    assert cuda_capability_supported((8, 0), ["nonsense"]) is False


# --- resolve_device -----------------------------------------------------------


def test_a_pinned_device_is_honoured_verbatim(monkeypatch):
    """Second-guessing an explicit pin would make the knob useless for the
    case it exists for — taking the decision out of play on a host where
    detection is the thing that is wrong."""
    called = []
    monkeypatch.setattr("pipeline.accelerator.unsupported_cuda_reason",
                        lambda: called.append(1))
    assert resolve_device("cpu") == "cpu"
    assert resolve_device("cuda") == "cuda"
    assert not called, "a pinned device must not consult the capability probe"


def test_auto_falls_back_to_cpu_when_the_gpu_is_unusable(monkeypatch):
    monkeypatch.setattr("pipeline.accelerator.unsupported_cuda_reason",
                        lambda: "GTX 1080 has capability 6.1 ...")
    assert resolve_device("auto") == "cpu"


def test_this_host_resolves_to_something_real():
    """Whatever the host is, the resolver must name a device docling accepts."""
    assert resolve_device("auto") in {"cuda", "mps", "cpu"}


def test_reason_is_none_when_there_is_no_cuda_at_all(monkeypatch):
    """"No GPU" is a different question, answered elsewhere; this probe only
    reports a GPU that is present and unusable."""
    import sys, types
    fake = types.ModuleType("torch")
    fake.cuda = types.SimpleNamespace(is_available=lambda: False)
    monkeypatch.setitem(sys.modules, "torch", fake)
    assert unsupported_cuda_reason() is None


def test_reason_names_the_card_and_the_architectures(monkeypatch):
    """The message has to be actionable — an operator needs to know which
    device was rejected and what the build supports."""
    import sys, types
    fake = types.ModuleType("torch")
    fake.cuda = types.SimpleNamespace(
        is_available=lambda: True,
        get_device_capability=lambda i: (6, 1),
        get_arch_list=lambda: PINNED,
        get_device_name=lambda i: "NVIDIA GeForce GTX 1080",
    )
    monkeypatch.setitem(sys.modules, "torch", fake)
    reason = unsupported_cuda_reason()
    assert "GTX 1080" in reason
    assert "6.1" in reason
    assert "sm_75" in reason
    assert "no kernel image is available" in reason


# --- config -------------------------------------------------------------------


def test_config_exposes_the_knob():
    from pipeline.config_schema import CorpuscleConfig, validate_config
    assert CorpuscleConfig().compute.accelerator == "auto"
    assert validate_config({"compute": {"accelerator": "cpu"}}).compute.accelerator == "cpu"


def test_config_rejects_a_device_docling_cannot_take():
    import pytest
    from pipeline.config_schema import validate_config
    with pytest.raises(Exception):
        validate_config({"compute": {"accelerator": "tpu"}})


def test_extract_pins_the_device_rather_than_leaving_docling_on_auto():
    """Guard the wiring: docling's default is `auto`, which is the bug."""
    import inspect
    from pipeline import extract
    src = inspect.getsource(extract)
    assert "accelerator_options=AcceleratorOptions(" in src
    assert "resolve_device(" in src


def test_the_vision_gate_does_not_trust_an_unusable_gpu():
    """`vision-local` is routed by _detect_accelerator; a GPU that cannot run
    kernels must not make that path look available."""
    import inspect
    from pipeline import cli
    assert "unsupported_cuda_reason" in inspect.getsource(cli._detect_accelerator)
