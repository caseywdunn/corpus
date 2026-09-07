"""Build memory is boundable from configuration (#182).

corpus exposed no way to bound the memory a build uses — not in
config.yaml, not on any CLI, not via an environment variable. On a
12-core / 32 GB CPU-only host, 7 concurrent `--only extract` workers put
three docling processes in flight at once, one of them on a 314-page
scan, and the build died mid-stage with `oom_kill 2` in /proc/vmstat and
no error in any worker log. Dropping to 4 workers removed the symptom,
but nothing in the product said so.

Two things are guarded here. The docling options that bound extraction,
which is where the memory actually goes; and `embeddings.batch_size`,
which `LocalBackend` has always taken and whose docstring described it as
the memory lever, while nothing could reach it.
"""
from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest
import yaml

from pipeline.config import load_config
from pipeline.config_schema import CorpuscleConfig


def _cfg(tmp_path: Path, **blocks) -> Path:
    p = tmp_path / "config.yaml"
    p.write_text(yaml.safe_dump({
        "input_pdfs": "pdfs", "output_dir": "out", **blocks,
    }))
    return p


# --- the schema ---------------------------------------------------------


def test_the_docling_bounds_are_configurable(tmp_path):
    cfg = CorpuscleConfig(**yaml.safe_load(_cfg(tmp_path, docling={
        "queue_max_size": 8, "layout_batch_size": 2,
        "ocr_batch_size": 2, "table_batch_size": 2, "document_timeout": 900,
    }).read_text()))
    assert cfg.docling.queue_max_size == 8
    assert cfg.docling.document_timeout == 900


def test_they_default_to_leaving_doclings_own_defaults_alone(tmp_path):
    """Picking numbers for everyone would be guessing at hardware we
    cannot see, and would change every existing build's behaviour."""
    cfg = CorpuscleConfig(**yaml.safe_load(_cfg(tmp_path).read_text()))
    assert cfg.docling.model_dump() == {
        "queue_max_size": None, "layout_batch_size": None,
        "ocr_batch_size": None, "table_batch_size": None,
        "document_timeout": None,
    }
    assert cfg.embeddings.batch_size is None
    assert cfg.compute.num_threads is None


@pytest.mark.parametrize("block,value", [
    ("docling", {"queue_max_size": 0}),
    ("docling", {"layout_batch_size": 0}),
    ("docling", {"document_timeout": 0}),
    ("embeddings", {"batch_size": 0}),
    ("compute", {"num_threads": 0}),
])
def test_nonsense_values_are_refused(tmp_path, block, value):
    with pytest.raises(Exception):
        CorpuscleConfig(**yaml.safe_load(_cfg(tmp_path, **{block: value}).read_text()))


@pytest.mark.parametrize("block", ["docling", "embeddings"])
def test_a_typo_in_the_new_blocks_is_refused(tmp_path, block):
    """`extra="forbid"`, so a misspelled knob is an error rather than a
    setting that silently does nothing — which is the shape of the defect
    this issue is about."""
    with pytest.raises(Exception):
        CorpuscleConfig(**yaml.safe_load(
            _cfg(tmp_path, **{block: {"queue_max_sixe": 4}}).read_text()))


def test_the_blocks_survive_load_config(tmp_path):
    """The schema accepting them is not the same as the loader keeping
    them; the pipeline reads the loader's dict."""
    loaded = load_config(_cfg(
        tmp_path,
        docling={"queue_max_size": 8},
        embeddings={"batch_size": 7},
        compute={"num_threads": 3},
    ))
    assert loaded["docling"]["queue_max_size"] == 8
    assert loaded["embeddings"]["batch_size"] == 7
    assert loaded["compute"]["num_threads"] == 3


# --- extraction ---------------------------------------------------------


def test_only_the_options_that_were_set_are_passed_to_docling(monkeypatch):
    """Omitted keys must not be passed as docling's documented default:
    that would silently pin a value the operator did not choose, and drift
    if docling ever changes it."""
    import inspect

    from pipeline import extract
    src = inspect.getsource(extract)
    assert "if value is not None" in src
    for key in ("queue_max_size", "layout_batch_size", "ocr_batch_size",
                "table_batch_size", "document_timeout"):
        assert key in src, key
    assert 'accel_kwargs["num_threads"]' in src


def test_the_option_names_match_the_installed_docling():
    """A silently-ignored option is worse than none, and these are
    keyword arguments to a pydantic model with extra fields forbidden —
    so a rename upstream must fail here, not in a build."""
    from docling.datamodel.pipeline_options import (
        AcceleratorOptions, PdfPipelineOptions,
    )
    for key in ("queue_max_size", "layout_batch_size", "ocr_batch_size",
                "table_batch_size", "document_timeout"):
        assert key in PdfPipelineOptions.model_fields, key
    assert "num_threads" in AcceleratorOptions.model_fields


# --- embeddings ---------------------------------------------------------


def _stub_sentence_transformers(monkeypatch):
    """Let LocalBackend construct without downloading 600 MB."""
    class _Model:
        def __init__(self, *a, **k):
            pass

        def get_sentence_embedding_dimension(self):
            return 8

        def half(self):
            return self

        def to(self, *_a):
            return self

        def eval(self):
            return self

    mod = types.ModuleType("sentence_transformers")
    mod.SentenceTransformer = _Model
    monkeypatch.setitem(sys.modules, "sentence_transformers", mod)


def test_batch_size_reaches_the_backend(monkeypatch):
    """The knob whose docstring documented it and whose plumbing did not
    exist: get_embedder forwards **kwargs, and its only caller built them
    from --device alone."""
    _stub_sentence_transformers(monkeypatch)
    from pipeline.embeddings import get_embedder

    assert get_embedder(None, device="cpu", batch_size=7).batch_size == 7
    # Unset still means the backend's own default.
    assert get_embedder(None, device="cpu").batch_size == 32


def test_embed_has_a_batch_size_flag_and_reads_the_config():
    import inspect

    from pipeline import embed
    src = inspect.getsource(embed)
    assert '"--batch-size"' in src
    assert '"--config"' in src
    # Flag beats config, per the project's usual precedence.
    assert "batch_size = args.batch_size" in src
    assert 'embedder_kwargs["batch_size"]' in src


def test_embed_pushes_the_config_into_the_module_dict():
    """Reading it locally is not enough: `embeddings.py` resolves its
    device from `CONFIG["compute"]["accelerator"]`, and because this
    module took no --config at all, that lookup always saw an empty CONFIG
    and fell back to "auto". A corpuscle pinning `compute.accelerator:
    cpu` was silently ignored by Stage 2 while being honoured by Stage 1 —
    and `embeddings.py` carried a comment saying the config reached the
    encoder here."""
    import inspect

    from pipeline import embed
    src = inspect.getsource(embed)
    assert "_pipeline_config.CONFIG.clear()" in src
    assert "_pipeline_config.CONFIG.update(_config)" in src


def test_the_orchestrator_forwards_config_to_the_embed_step():
    import argparse
    import inspect

    from pipeline.orchestrator import STEPS
    step = next(s for s in STEPS if s.name == "embed")
    args = argparse.Namespace(
        output_dir=Path("out"), resume=True, dry_run=False,
        config=Path("config.yaml"),
    )
    argv = step.argv(args)
    assert "--config" in argv and "config.yaml" in argv
    assert inspect.isfunction(type(step).argv)


def test_a_named_config_that_does_not_exist_is_refused():
    """Same rule as pipeline.main (#210): silently replacing every tuned
    setting with defaults is worse than stopping."""
    import inspect

    from pipeline import embed
    src = inspect.getsource(embed)
    assert "no such file" in src


# --- the external cap ---------------------------------------------------


def test_the_cgroup_cap_is_documented():
    """The outer bound needs no code and carries most of the value, so it
    has to be written down where an operator will find it. MemoryHigh
    throttles by reclaim and only MemoryMax kills, so a capped build is
    squeezed first and any kill lands inside its own cgroup rather than on
    a bystander — which on the affected host was a tmux server hosting
    unrelated work."""
    root = Path(__file__).resolve().parent.parent
    text = (root / "INSTALL.md").read_text()
    assert "systemd-run" in text
    assert "MemoryMax" in text
    assert "MemoryHigh" in text
