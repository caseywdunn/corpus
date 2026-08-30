"""`corpus check` must not promise a taxonomy that `--only` won't build (#251).

`check` reported a missing `taxonomy.sqlite` as "will be created on first
run". True for a full `corpus run`; false for the phase-split
`run --only <phase>` path that `slurm/` uses, where
`_check_taxonomy_available` hard-errors before any work starts. A corpuscle
that passed `check` therefore killed every Stage 1 array task about a minute
in, and `afterok` took Pass 3b, Embed and Finalize down with it. Two
consecutive siphonophore builds lost their chain this way.

The WoRMS branch had said the right thing since #139 — "pre-build it once
before submitting HPC array jobs" — and the dwca/dwc branch forty lines below
was never brought into line. These tests hold both branches to it.
"""
from __future__ import annotations

import argparse
import inspect
from pathlib import Path

import pytest

from pipeline import cli

REPO = Path(__file__).resolve().parent.parent


# --- what `check` tells the operator -----------------------------------------


def test_check_no_longer_promises_a_build_that_only_wont_do():
    src = inspect.getsource(cli)
    assert "will be created on first run" not in src, (
        "the dwca/dwc branch is promising the optimistic full-run case again"
    )


def test_both_taxonomy_branches_tell_the_operator_to_pre_build():
    """Wording will drift; the instruction must not."""
    src = inspect.getsource(cli)
    i = src.index("# 8. Taxonomy source availability")
    section = src[i:i + 4000]
    assert section.count("pre-build") + section.count("Pre-build") >= 2, (
        "each of the worms and dwca/dwc branches should tell the operator to "
        "pre-build when the sqlite is absent"
    )


def test_the_config_template_does_not_repeat_the_promise():
    """`corpus init` copies this verbatim, so the operator would be told twice
    before ever submitting."""
    tmpl = (REPO / "pipeline" / "config.template.yaml").read_text()
    assert "on the first `corpus run`" not in tmpl
    assert "corpus taxonomy ingest" in tmpl


# --- and what the submit path does about it ----------------------------------


def test_the_slurm_entry_point_pre_builds_before_submitting():
    """The wording fix only helps an operator who reads `check`. Both builds
    that were lost failed with `check` output on screen."""
    sh = (REPO / "slurm" / "batch_pipeline.sh").read_text()
    assert "corpus taxonomy ingest" in sh
    assert sh.index("corpus taxonomy ingest") < sh.index("Submitting Stage 1"), (
        "the pre-build has to happen before the Stage 1 sbatch, not after"
    )


def test_the_fatal_precondition_reaches_stdout():
    """SLURM sends the logger's stderr to a sibling .err file, so a chain
    dying on this looked like slow first documents from every vantage point
    an operator actually uses."""
    from pipeline import orchestrator
    src = inspect.getsource(orchestrator)
    i = src.index("taxonomy_err = _check_taxonomy_available")
    assert "print(" in src[i:i + 900]


# --- driving the ingest from config ------------------------------------------


@pytest.fixture
def cfg(tmp_path):
    (tmp_path / "tax.zip").write_bytes(b"PK\x03\x04")
    (tmp_path / "config.yaml").write_text(
        "input_pdfs: .\n"
        f"output_dir: {tmp_path / 'out'}\n"
        "taxonomy:\n"
        "  source: dwca\n"
        "  path: ./tax.zip\n"
    )
    return tmp_path / "config.yaml"


def _argv(monkeypatch, cfg, passthrough):
    """Capture the argv `corpus taxonomy ingest` would forward."""
    seen = {}
    monkeypatch.setattr(cli, "_passthrough",
                        lambda mod, argv: seen.setdefault("argv", argv) or 0)
    args = argparse.Namespace(output_dir=None, config=cfg, passthrough=passthrough)
    cli._cmd_taxonomy_ingest(args)
    return seen["argv"]


def test_source_and_input_come_from_the_config(monkeypatch, cfg):
    """Without this the SLURM guard would have to parse YAML in bash."""
    argv = _argv(monkeypatch, cfg, [])
    assert argv[argv.index("--source") + 1] == "dwca"
    assert argv[argv.index("--input") + 1].endswith("tax.zip")


def test_the_input_path_is_resolved_against_the_config_not_the_cwd(monkeypatch, cfg):
    """`path: ./tax.zip` is relative to config.yaml. A submit script runs from
    the repo, not the corpuscle."""
    argv = _argv(monkeypatch, cfg, [])
    assert Path(argv[argv.index("--input") + 1]).is_absolute()
    assert Path(argv[argv.index("--input") + 1]).exists()


def test_an_explicit_source_still_wins(monkeypatch, cfg):
    argv = _argv(monkeypatch, cfg, ["--source", "worms", "--root-id", "1371"])
    assert argv.count("--source") == 1
    assert argv[argv.index("--source") + 1] == "worms"
