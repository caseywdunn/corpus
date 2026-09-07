"""Which compute device the pipeline may actually use (#198).

``torch.cuda.is_available()`` answers "is there a visible NVIDIA GPU", which
is not the question anything here needs. The question is "can this torch build
run kernels on it", and the two differ on exactly the hardware a lab
workstation has.

The failure this module exists to prevent, observed on the same machine with
an unchanged dependency set 20 days apart::

    2026-08-06  INFO docling.utils.accelerator_utils: Accelerator device: 'cpu'
    2026-08-26  INFO docling.utils.accelerator_utils: Accelerator device: 'cuda:0'
                ERROR docling.pipeline.standard_pdf_pipeline: Stage layout
                  failed: CUDA error: no kernel image is available for
                  execution on the device

The card is a GTX 1080 — compute capability 6.1 — and the pinned torch ships
``['sm_75', 'sm_80', 'sm_86', 'sm_90', 'sm_100', 'sm_120']``. Nothing in the
project changed; a driver became visible, ``is_available()`` flipped to true,
and every page failed. **A driver appearing turned a working install into a
broken one**, which is precisely the reproducibility the #98 version pins
exist to provide.

So capability is checked, not availability, and the resolved device is
recorded rather than left implicit — two corpuscles that differ because one
ran on CPU and one on CUDA were previously indistinguishable after the fact.
"""
from __future__ import annotations

import logging
import re
from typing import Optional, Tuple

logger = logging.getLogger(__name__)

# `sm_90a` and friends: architecture-specific suffixes on the numeric part.
_ARCH_RE = re.compile(r"^(sm|compute)_(\d+)([a-z]*)$")


def _parse_arch(name: str) -> Optional[Tuple[str, int]]:
    """``"sm_86"`` → ``("sm", 86)``. ``None`` for anything unrecognised."""
    m = _ARCH_RE.match(name.strip())
    return (m.group(1), int(m.group(2))) if m else None


def cuda_capability_supported(capability, arch_list) -> bool:
    """Can a torch built for ``arch_list`` run on a device of ``capability``?

    ``capability`` is torch's ``(major, minor)``; ``arch_list`` is
    ``torch.cuda.get_arch_list()``.

    Two ways a device can be served, and both must be allowed or this rejects
    working hardware:

    * a **binary** kernel compiled for exactly this architecture — ``sm_86``
      for a (8, 6) device;
    * **PTX** for an *older* architecture — ``compute_75`` JITs forward onto
      any device at or above 7.5. Forward only: PTX never runs on a device
      older than it was emitted for, which is the GTX 1080 case.

    An empty ``arch_list`` means a CPU-only torch build, where the question
    does not arise; treated as unsupported so callers fall through to CPU.
    """
    if not arch_list:
        return False
    device = capability[0] * 10 + capability[1]
    for name in arch_list:
        parsed = _parse_arch(name)
        if parsed is None:
            continue
        kind, arch = parsed
        if kind == "sm" and arch == device:
            return True
        if kind == "compute" and arch <= device:
            return True
    return False


def unsupported_cuda_reason() -> Optional[str]:
    """Why the visible CUDA device cannot run this torch build, or ``None``.

    ``None`` also when there is no CUDA device at all — the caller's question
    is "should I avoid CUDA", and "there isn't any" is answered elsewhere.
    """
    try:
        import warnings
        with warnings.catch_warnings():
            # Same reasoning as cli.py's detector: a driver/runtime mismatch
            # makes torch write two lines of internals to stderr ahead of any
            # corpus output, on every invocation including `--help`.
            warnings.simplefilter("ignore")
            import torch
            if not torch.cuda.is_available():
                return None
            capability = torch.cuda.get_device_capability(0)
            arch_list = torch.cuda.get_arch_list()
            if cuda_capability_supported(capability, arch_list):
                return None
            name = torch.cuda.get_device_name(0)
    except Exception:
        # A torch that cannot answer the question cannot be trusted to run
        # kernels either, but that is not this function's call to make.
        return None
    return (
        f"{name} has CUDA compute capability "
        f"{capability[0]}.{capability[1]}, and this torch build ships kernels "
        f"for {', '.join(arch_list)}. Every kernel launch would fail with "
        f"'no kernel image is available for execution on the device'. "
        f"Falling back to CPU."
    )


class AcceleratorUnavailable(RuntimeError):
    """``compute.accelerator: require`` was set and no GPU is usable.

    Raised instead of returning ``"cpu"``, because inside a scheduler
    allocation a CPU fallback is not a degraded success — it is a failure
    that costs the whole allocation and looks healthy while doing it.
    """


def _detect_device() -> str:
    """The best device this torch build can actually run kernels on."""
    reason = unsupported_cuda_reason()
    if reason:
        logger.warning("Ignoring the visible GPU: %s", reason)
        return "cpu"
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import torch
            if torch.cuda.is_available():
                return "cuda"
            if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
                return "mps"
    except Exception:
        pass
    return "cpu"


def resolve_device(configured: str = "auto") -> str:
    """Return the device to use: ``"cuda"``, ``"mps"`` or ``"cpu"``.

    ``configured`` comes from ``compute.accelerator`` in ``config.yaml``.

    * ``auto`` detects, and falls back to CPU when the visible GPU is one
      this torch build ships no kernels for. Right on a workstation, where
      the alternative is every docling page failing.
    * ``require`` detects the same way and **raises**
      :class:`AcceleratorUnavailable` rather than falling back. For a
      scheduler allocation, where CPU is not an acceptable outcome (#270).
    * anything else — ``cpu``, ``cuda``, ``mps`` — is honoured verbatim. An
      operator pinning a device is making a deliberate choice, and
      second-guessing it would make the knob useless for the case it exists
      for. Note pinning ``cuda`` *disables* the capability check rather than
      enforcing it, which is why it is not a way to say "require".
    """
    if configured == "require":
        device = _detect_device()
        if device == "cpu":
            raise AcceleratorUnavailable(_require_failure_message())
        return device
    if configured and configured != "auto":
        return configured
    return _detect_device()


def _require_failure_message() -> str:
    """Why ``require`` could not be satisfied, in the operator's terms.

    Two distinct causes, and conflating them sends the reader in the wrong
    direction. A visible-but-unusable GPU is a torch/hardware mismatch and
    :func:`unsupported_cuda_reason` already builds exactly that diagnosis;
    no GPU at all is a submission problem — the wrong partition, or a
    missing ``--gpus``.
    """
    reason = unsupported_cuda_reason()
    if reason:
        # unsupported_cuda_reason ends with "Falling back to CPU.", which is
        # precisely what is *not* happening here.
        reason = reason.replace(
            "Falling back to CPU.",
            "compute.accelerator is 'require', so this is a failure rather "
            "than a fallback.",
        )
        return (
            f"No usable GPU: {reason} Either install a torch build with "
            f"kernels for this card, request a card this build supports "
            f"(the SLURM scripts pin a type for this reason), or set "
            f"compute.accelerator to 'auto' to accept CPU."
        )
    return (
        "No usable GPU: compute.accelerator is 'require' but torch reports "
        "no CUDA or MPS device. Inside a scheduler allocation this usually "
        "means the job did not request one, or requested it on a partition "
        "that has none. Set compute.accelerator to 'auto' to accept CPU."
    )
