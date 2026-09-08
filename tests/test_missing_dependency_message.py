"""A dependency error must name the module that actually failed (#258 follow-up).

The local VLM backend reported "transformers >= 4.45 is required" on a Mac
where `torch` was the missing package and transformers was fine. A confident
wrong instruction is worse than a vague one: following it leaves the install
broken and costs the reader a cycle.
"""
from __future__ import annotations

import pytest

from pipeline.optional_deps import missing_dependency_message


def _import_error(name):
    exc = ImportError(f"No module named {name!r}")
    exc.name = name
    return exc


# ── the case that was wrong ────────────────────────────────────────────

def test_a_broken_dependency_chain_names_the_broken_link():
    msg = missing_dependency_message(
        _import_error("torch"),
        feature="the local VLM backend", module="transformers",
        requirement="transformers >= 4.45",
        install="transformers>=4.45 qwen-vl-utils torch accelerate")
    assert "'torch' is not installed" in msg
    # And it must say plainly that transformers is not the thing to fix,
    # because that is the instruction that wasted the reader's time.
    assert "not transformers" in msg


def test_the_asked_for_package_is_still_named_when_it_is_the_one_missing():
    msg = missing_dependency_message(
        _import_error("transformers"),
        feature="the local VLM backend", module="transformers",
        requirement="transformers >= 4.45",
        install="transformers>=4.45 qwen-vl-utils torch accelerate")
    assert msg.startswith("transformers >= 4.45 is required")
    assert "pip install transformers>=4.45" in msg


def test_a_missing_symbol_reads_as_the_version_bound_it_is():
    """`from transformers import <new class>` on an old transformers sets
    `name` to 'transformers', so the version bound in `requirement` is
    what the reader needs — not a claim about a different package."""
    exc = ImportError("cannot import name 'Qwen2_5_VLForConditionalGeneration'")
    exc.name = "transformers"
    msg = missing_dependency_message(
        exc, feature="the local VLM backend", module="transformers",
        requirement="transformers >= 4.45", install="transformers>=4.45")
    assert "transformers >= 4.45 is required" in msg


# ── shape of the API ───────────────────────────────────────────────────

def test_an_importerror_with_no_name_falls_back_rather_than_saying_none():
    """Not every ImportError carries `name`; `'None' is not installed`
    would be worse than the old message."""
    msg = missing_dependency_message(
        ImportError("something odd"), feature="the local backend",
        module="sentence_transformers", requirement="sentence-transformers",
        install="sentence-transformers")
    assert "None" not in msg
    assert msg.startswith("sentence-transformers is required")


def test_the_import_name_and_the_pip_name_can_differ():
    """`sentence_transformers` is imported with an underscore and
    installed with a hyphen; the reader needs the hyphen."""
    msg = missing_dependency_message(
        _import_error("sentence_transformers"),
        feature="the local embedding backend",
        module="sentence_transformers", requirement="sentence-transformers",
        install="sentence-transformers")
    assert "sentence-transformers is required" in msg
    assert "sentence_transformers is required" not in msg


def test_requirement_defaults_to_the_module_name():
    msg = missing_dependency_message(
        _import_error("anthropic"), feature="the Claude vision backend",
        module="anthropic", install="anthropic")
    assert msg == ("anthropic is required for the Claude vision backend "
                   "(pip install anthropic)")


# ── the real call sites use it ─────────────────────────────────────────

@pytest.mark.parametrize("module,missing,expect", [
    ("transformers", "torch", "'torch' is not installed"),
    ("transformers", "transformers", "transformers >= 4.45 is required"),
])
def test_the_local_vlm_backend_reports_the_real_culprit(
        module, missing, expect, monkeypatch):
    """End to end: the backend's own handler, not just the helper."""
    import builtins

    from pipeline.vision import LocalVLMBackend, VisionBackendError

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == module:
            raise _import_error(missing)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(VisionBackendError) as caught:
        LocalVLMBackend(device="cpu")
    assert expect in str(caught.value)
