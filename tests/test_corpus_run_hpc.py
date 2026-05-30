"""`corpus run` HPC phases — single-phase + array slicing.

Makes `corpus run` usable on HPC (one phase per SLURM job/partition,
job-array slicing) so the production path is the same CLI users run,
instead of a divergent raw-`pipeline.main` reimplementation.

Covers:
- `corpus run --only {extract,vision,embed,post,bundle}` + `--batch-index`
  / `--batch-size` parse, and their passthrough to the orchestrator.
- The orchestrator's `select_steps` phase resolution.
- The orchestrator extract/vision step argv (batch flags, vision forces
  a vision backend).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pytest

import pipeline.cli as cli
import pipeline.orchestrator as orch
from pipeline.config_schema import validate_config


# --- CLI: parse + _build_orchestrator_argv passthrough -----------------------


def _cfg(tmp_path: Path):
    (tmp_path / "pdfs").mkdir(exist_ok=True)
    return validate_config({
        "input_pdfs": "./pdfs", "output_dir": "./output",
        "grobid": {"disable": True, "url": ""},
    })


def _run_args(**overrides):
    base = dict(
        force_rebuild=False, dry_run=False, no_vision=False, figure_panels=None,
        enrich_bhl=False, force_rebuild_taxonomy=False, force_rebuild_biblio=False,
        force_rebuild_taxon_mentions=False, only=None, batch_index=None,
        batch_size=None,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _argv(tmp_path, **arg_overrides):
    return cli._build_orchestrator_argv(
        _cfg(tmp_path), tmp_path / "config.yaml", _run_args(**arg_overrides))


def test_run_parses_only_and_batch_flags():
    p = cli._build_parser()
    a = p.parse_args(["run", "--only", "extract", "--batch-index", "3",
                      "--batch-size", "256"])
    assert a.only == "extract" and a.batch_index == 3 and a.batch_size == 256
    assert p.parse_args(["run"]).only is None  # default = full run


def test_build_argv_forwards_only_and_batch(tmp_path):
    argv = _argv(tmp_path, only="extract", batch_index=2, batch_size=128)
    assert "--only" in argv and argv[argv.index("--only") + 1] == "extract"
    assert "--batch-index" in argv and argv[argv.index("--batch-index") + 1] == "2"
    assert "--batch-size" in argv and argv[argv.index("--batch-size") + 1] == "128"


def test_build_argv_omits_only_for_full_run(tmp_path):
    assert "--only" not in _argv(tmp_path)  # no --only → full chain


def test_build_argv_never_forwards_bundle(tmp_path):
    # bundle is a CLI-level phase, handled before the orchestrator runs.
    assert "--only" not in _argv(tmp_path, only="bundle")


# --- orchestrator: select_steps phase resolution -----------------------------


def _orch_args(**overrides):
    base = dict(only=None, skip_pipeline=False, from_step=None,
                taxonomy_source=None, force_rebuild=False,
                force_rebuild_taxonomy=False, output_dir=Path("/nope"),
                batch_index=None, batch_size=None)
    base.update(overrides)
    return argparse.Namespace(**base)


@pytest.mark.parametrize("phase,expected", [
    ("extract", ["extract"]),
    ("vision", ["vision"]),
    ("embed", ["embed"]),
    ("post", ["build_biblio", "build_taxa", "backfill_intext", "reconcile"]),
])
def test_select_steps_only_phase(phase, expected):
    names = [s.name for s in orch.select_steps(_orch_args(only=phase))]
    assert names == expected


def test_select_steps_full_run_has_no_vision_step():
    # vision runs inline in extract on a full run; it's a separate step
    # only via --only vision.
    names = [s.name for s in orch.select_steps(_orch_args())]
    assert "vision" not in names
    assert "extract" in names and "embed" in names


# --- orchestrator: step argv (batch slicing + vision backend) ----------------


def _argv_args(**overrides):
    base = dict(
        input_dir=Path("/in"), output_dir=Path("/out"), resume=True,
        dry_run=False, config=None, bib=None, taxonomy_db=None, lexicon=None,
        no_taxa=False, no_grobid=False, grobid_url=None, strict_network=False,
        figure_panels="ocr", vision_model=None, refresh_vision=False,
        batch_index=None, batch_size=None,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def test_extract_step_carries_batch_flags():
    step = next(s for s in orch.STEPS if s.name == "extract")
    argv = step.argv(_argv_args(batch_index=4, batch_size=64))
    assert "--batch-index" in argv and "4" in argv
    assert "--batch-size" in argv and "64" in argv


def test_vision_step_forces_vision_backend_and_refresh():
    argv = orch.VISION_STEP.argv(_argv_args(figure_panels="ocr",
                                            batch_index=1, batch_size=8))
    assert "--refresh-vision" in argv
    # extract phase used the ocr floor; vision phase must force a backend
    assert argv[argv.index("--figure-panels") + 1] == "vision-local"
    assert "--no-grobid" in argv and "--no-taxa" in argv
    assert "--batch-index" in argv


def test_vision_step_honors_configured_vision_backend():
    argv = orch.VISION_STEP.argv(_argv_args(figure_panels="vision-claude"))
    assert argv[argv.index("--figure-panels") + 1] == "vision-claude"


def test_batch_flags_need_both():
    assert orch._batch_flags(_argv_args(batch_index=1)) == []  # size missing
    assert orch._batch_flags(_argv_args(batch_index=1, batch_size=2)) == \
        ["--batch-index", "1", "--batch-size", "2"]
