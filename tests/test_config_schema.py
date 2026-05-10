"""Tests for the per-corpuscle config.yaml pydantic schema (#59).

Locks in the field-level validation contract: missing required keys,
wrong types, and unknown enum values surface with errors that point
at the exact key + value.
"""
from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from pipeline.config_schema import (
    CorpuscleConfig,
    ValidationError,
    validate_config,
)

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_empty_config_uses_defaults():
    cfg = CorpuscleConfig()
    assert cfg.vision.backend == "none"
    assert cfg.licensing.pd_cutoff_years == 95
    assert cfg.bibliography.enrich_bhl is False
    assert cfg.grobid.url.startswith("http")


def test_bundled_template_validates():
    template = REPO_ROOT / "pipeline" / "config.template.yaml"
    raw = yaml.safe_load(template.read_text(encoding="utf-8"))
    cfg = validate_config(raw)
    # Template ships with WoRMS taxonomy and demo-friendly defaults.
    assert cfg.taxonomy.source == "worms"
    assert cfg.taxonomy.root_id == 1267


def test_demo_config_validates():
    demo_cfg = REPO_ROOT / "demo" / "config.yaml"
    raw = yaml.safe_load(demo_cfg.read_text(encoding="utf-8"))
    cfg = validate_config(raw)
    assert cfg.taxonomy.source == "worms"
    assert cfg.bib == Path("siphonophores.bib")
    assert cfg.lexicon == Path("lexicon.yaml")


def test_bad_enum_points_at_key():
    with pytest.raises(ValidationError) as exc:
        validate_config({"vision": {"backend": "claud"}})
    locs = [".".join(str(x) for x in err["loc"]) for err in exc.value.errors()]
    assert "vision.backend" in locs


def test_unknown_top_level_key_rejected():
    with pytest.raises(ValidationError) as exc:
        validate_config({"visions": {"backend": "claude"}})
    msgs = [err["msg"] for err in exc.value.errors()]
    assert any("Extra inputs" in m for m in msgs)


def test_taxonomy_worms_requires_root_id():
    with pytest.raises(ValidationError) as exc:
        validate_config({"taxonomy": {"source": "worms"}})
    msgs = [err["msg"] for err in exc.value.errors()]
    assert any("root_id" in m for m in msgs)


def test_taxonomy_dwca_requires_path():
    with pytest.raises(ValidationError) as exc:
        validate_config({"taxonomy": {"source": "dwca"}})
    msgs = [err["msg"] for err in exc.value.errors()]
    assert any("path" in m for m in msgs)


def test_pd_cutoff_must_be_non_negative():
    with pytest.raises(ValidationError):
        validate_config({"licensing": {"pd_cutoff_years": -1}})
