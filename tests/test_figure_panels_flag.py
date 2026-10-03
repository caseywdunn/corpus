"""`--figure-panels` selector wiring (#102).

The v0.5 `--vision-backend` + `--content-aware-figures` pair (config key
`vision.backend`) was replaced by a single `--figure-panels
{ocr,vision-local,vision-claude,off}` selector (config key
`figures.panel_detection`, default `ocr`). These tests pin the two ends:

* `_build_orchestrator_argv` translates `figures.panel_detection` into a
  single `--figure-panels <mode>`, with `--no-vision` and capability
  detection deferring unavailable vision during extract-only runs, with an
  OCR fallback for standalone runs.
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
        force_rebuild=False, dry_run=False, no_vision=False, figure_panels=None,
        enrich_bhl=False, force_rebuild_taxonomy=False, force_rebuild_biblio=False,
        force_rebuild_taxon_mentions=False, require_gpu=False,
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


# --- corpus run --figure-panels override (#102 follow-up) ---------------------


def test_run_figure_panels_overrides_config(tmp_path, monkeypatch):
    """An explicit `corpus run --figure-panels` wins over the config key."""
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-test")
    # config says ocr; the CLI override asks for vision-claude.
    argv = _argv(tmp_path, "ocr", figure_panels="vision-claude")
    assert argv[argv.index("--figure-panels") + 1] == "vision-claude"


def test_run_figure_panels_off_overrides_config(tmp_path):
    argv = _argv(tmp_path, "vision-claude", figure_panels="off")
    assert argv[argv.index("--figure-panels") + 1] == "off"


def test_no_vision_is_alias_for_ocr_even_from_off(tmp_path):
    """`--no-vision` ≡ `--figure-panels ocr`, so it forces the OCR floor
    even when the config selected `off`."""
    argv = _argv(tmp_path, "off", no_vision=True)
    assert argv[argv.index("--figure-panels") + 1] == "ocr"


def test_run_rejects_figure_panels_with_no_vision(tmp_path):
    """The two are mutually exclusive in the `corpus run` parser."""
    parser = cli._build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["run", "--figure-panels", "ocr", "--no-vision"])


def test_run_accepts_figure_panels_choice(tmp_path):
    parser = cli._build_parser()
    args = parser.parse_args(["run", "--figure-panels", "vision-local"])
    assert args.figure_panels == "vision-local"
    assert args.no_vision is False


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


# --- #263: the capability check only runs on phases that use it ---------
#
# In the standard HPC chain the finalize job runs `corpus run --only post`
# on a CPU node, so the "vision panel pass downgraded to the OCR floor"
# warning fired on every build — while Pass 3b had already completed on a
# GPU an hour earlier and its ROIs were intact (3,465 ROIs across 288
# documents, 100% `source: vision:qwen2.5-vl-7b-instruct`, zero OCR-floor).
#
# These phases must not change a configured mode or probe unused hardware.


@pytest.mark.parametrize("phase", ["post", "embed", "bundle"])
def test_no_vision_warning_on_a_phase_that_never_runs_vision(
    tmp_path, monkeypatch, capsys, phase,
):
    monkeypatch.setattr(cli, "_vision_skip_reason",
                        lambda _mode: "no CUDA/MPS detected on this host")
    argv = _argv(tmp_path, "vision-local", only=phase)
    out = capsys.readouterr().out
    assert "downgraded to the OCR floor" not in out
    # The mode is forwarded untouched, so nothing downstream sees `ocr`
    # either — the phase simply ignores it.
    assert "vision-local" in argv


@pytest.mark.parametrize("phase", ["vision", None])
def test_the_warning_survives_where_it_matters(
    tmp_path, monkeypatch, capsys, phase,
):
    """#65's behaviour, preserved: these phases do run the vision pass."""
    monkeypatch.setattr(cli, "_vision_skip_reason",
                        lambda _mode: "no CUDA/MPS detected on this host")
    argv = _argv(tmp_path, "vision-local", only=phase)
    assert "downgraded to the OCR floor" in capsys.readouterr().out
    assert "ocr" in argv and "vision-local" not in argv


@pytest.mark.parametrize('override', [None, 'vision-local'])
def test_extract_defers_unavailable_local_vision_and_later_vision_runs(tmp_path, monkeypatch, capsys, override):
    from pipeline.main import _panels_to_legacy
    monkeypatch.setattr(cli, '_detect_accelerator', lambda: None)
    argv = _argv(tmp_path, 'vision-local', only='extract', figure_panels=override)
    mode = argv[argv.index('--figure-panels')+1]
    assert mode == 'off' and _panels_to_legacy(mode) == (False, None)
    warning = capsys.readouterr().out
    assert 'deferred for extract-only run' in warning and 'downgraded to the OCR floor' not in warning
    assert 'corpus run --only vision' in warning
    # The same unchanged config on the later accelerator phase selects the
    # actual vision backend; extraction did not persist an "off" override.
    monkeypatch.setattr(cli, '_detect_accelerator', lambda: 'cuda')
    later = _argv(tmp_path, 'vision-local', only='vision', figure_panels=override)
    mode = later[later.index('--figure-panels')+1]
    assert mode == 'vision-local' and _panels_to_legacy(mode) == (False, 'local')


@pytest.mark.parametrize('arguments', [{'figure_panels': 'ocr'}, {'no_vision': True}, {}])
def test_explicit_or_configured_ocr_still_runs_during_extract(tmp_path, monkeypatch, arguments):
    from pipeline.main import _panels_to_legacy
    def unused(_mode):
        raise AssertionError('OCR requested; no vision capability probe is needed')
    monkeypatch.setattr(cli, '_vision_skip_reason', unused)
    config_mode = 'vision-local' if arguments else 'ocr'
    argv = _argv(tmp_path, config_mode, only='extract', **arguments)
    mode = argv[argv.index('--figure-panels')+1]
    assert mode == 'ocr' and _panels_to_legacy(mode) == (True, None)


@pytest.mark.parametrize('accelerator', ['cuda', 'mps'])
def test_extract_keeps_available_local_vision(tmp_path, monkeypatch, accelerator):
    monkeypatch.setattr(cli, '_detect_accelerator', lambda: accelerator)
    argv = _argv(tmp_path, 'vision-local', only='extract')
    assert argv[argv.index('--figure-panels')+1] == 'vision-local'


@pytest.mark.parametrize('has_key,expected', [(True, 'vision-claude'), (False, 'off')])
def test_extract_cloud_vision_uses_credentials_not_local_accelerator(tmp_path, monkeypatch, has_key, expected):
    if has_key:
        monkeypatch.setenv('ANTHROPIC_API_KEY', 'test-only-credential')
    else:
        monkeypatch.delenv('ANTHROPIC_API_KEY', raising=False)
    def unused():
        raise AssertionError('Cloud vision must not probe local accelerators')
    monkeypatch.setattr(cli, '_detect_accelerator', unused)
    argv = _argv(tmp_path, 'vision-claude', only='extract')
    assert argv[argv.index('--figure-panels')+1] == expected


@pytest.mark.parametrize("phase", ["post", "embed", "bundle"])
def test_a_usable_backend_is_not_probed_on_those_phases(
    tmp_path, monkeypatch, phase,
):
    """Not just the warning — the capability probe itself is skipped, so a
    CPU finalize node does not pay for a CUDA/MPS detection it cannot use."""
    def should_not_run(_mode):
        raise AssertionError("capability probe ran on a non-vision phase")

    monkeypatch.setattr(cli, "_vision_skip_reason", should_not_run)
    _argv(tmp_path, "vision-local", only=phase)
