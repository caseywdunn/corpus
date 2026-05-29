"""`--figure-panels` selector wiring (#102).

The v0.5 `--vision-backend` + `--content-aware-figures` pair (config key
`vision.backend`) was replaced by a single `--figure-panels
{ocr,vision-local,vision-claude,off}` selector (config key
`figures.panel_detection`, default `ocr`). These tests pin the two ends:

* `_build_orchestrator_argv` translates `figures.panel_detection` into a
  single `--figure-panels <mode>`, with `--no-vision` and capability
  detection downgrading a vision mode to the OCR floor.
* `pipeline.main`'s arg parser derives the legacy
  (content_aware_figures, vision_backend) pair the runner still threads.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pytest

import pipeline.cli as cli
from pipeline.config_schema import validate_config


def _cfg(tmp_path: Path, panel_detection="ocr"):
    (tmp_path / "pdfs").mkdir(exist_ok=True)
    return validate_config({
        "input_pdfs": "./pdfs",
        "output_dir": "./output",
        "figures": {"panel_detection": panel_detection},
        # disable grobid so the argv stays minimal/deterministic
        "grobid": {"disable": True, "url": ""},
    })


def _run_args(**overrides):
    base = dict(
        force_rebuild=False, dry_run=False, no_vision=False, enrich_bhl=False,
        force_rebuild_taxonomy=False, force_rebuild_biblio=False,
        force_rebuild_taxon_mentions=False,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _argv(tmp_path, panel_detection="ocr", **arg_overrides):
    cfg = _cfg(tmp_path, panel_detection)
    return cli._build_orchestrator_argv(cfg, tmp_path / "config.yaml", _run_args(**arg_overrides))


def test_default_emits_figure_panels_ocr(tmp_path):
    argv = _argv(tmp_path, "ocr")
    assert "--figure-panels" in argv
    assert argv[argv.index("--figure-panels") + 1] == "ocr"
    # legacy flags are gone
    assert "--content-aware-figures" not in argv
    assert "--vision-backend" not in argv


def test_off_mode_passes_through(tmp_path):
    argv = _argv(tmp_path, "off")
    assert argv[argv.index("--figure-panels") + 1] == "off"


def test_no_vision_downgrades_vision_mode_to_ocr(tmp_path, monkeypatch):
    # vision-claude would otherwise need ANTHROPIC_API_KEY; --no-vision
    # downgrades to the OCR floor before any capability check.
    argv = _argv(tmp_path, "vision-claude", no_vision=True)
    assert argv[argv.index("--figure-panels") + 1] == "ocr"


def test_unusable_vision_backend_downgrades_to_ocr(tmp_path, monkeypatch):
    # No ANTHROPIC_API_KEY → vision-claude isn't usable → OCR floor.
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    argv = _argv(tmp_path, "vision-claude")
    assert argv[argv.index("--figure-panels") + 1] == "ocr"


def test_usable_vision_mode_is_forwarded(tmp_path, monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-test")
    argv = _argv(tmp_path, "vision-claude")
    assert argv[argv.index("--figure-panels") + 1] == "vision-claude"


# --- pipeline.main derivation -------------------------------------------------


@pytest.mark.parametrize("mode,expect_caf,expect_backend", [
    ("ocr", True, None),
    ("vision-local", False, "local"),
    ("vision-claude", False, "claude"),
    ("off", False, None),
])
def test_main_derives_legacy_pair(mode, expect_caf, expect_backend):
    """main._panels_to_legacy maps --figure-panels → the
    (content_aware_figures, vision_backend) pair the runner consumes."""
    from pipeline.main import _panels_to_legacy

    content_aware, backend = _panels_to_legacy(mode)
    assert content_aware is expect_caf
    assert backend == expect_backend
