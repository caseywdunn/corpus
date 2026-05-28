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
    # Template ships with the taxonomy block commented out so a fresh
    # `corpus init` doesn't lock new users into a Siphonophorae walk.
    # Operators uncomment and pick a source (worms / dwca / dwc).
    assert cfg.taxonomy.source is None
    assert cfg.taxonomy.root_id is None


def test_demo_config_validates():
    demo_cfg = REPO_ROOT / "demo" / "config.yaml"
    raw = yaml.safe_load(demo_cfg.read_text(encoding="utf-8"))
    cfg = validate_config(raw)
    # demo ships the Siphonophorae taxonomy as a pre-built DwC-A
    # (commit 447437f) so `corpus run` doesn't walk the WoRMS REST API
    # on every fresh build. Schema must accept the dwca path + a
    # relative `taxonomy.path`.
    assert cfg.taxonomy.source == "dwca"
    assert cfg.taxonomy.path == Path("./taxonomy.zip")
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


def test_taxonomy_root_id_accepts_integer():
    """Regression guard for the integer-AphiaID path that worms uses."""
    cfg = validate_config({"taxonomy": {"source": "worms", "root_id": 1371}})
    assert cfg.taxonomy.root_id == 1371


def test_taxonomy_root_id_accepts_lsid_string(tmp_path):
    """#78: DwC-A snapshots from WoRMS use LSIDs (e.g.
    ``urn:lsid:marinespecies.org:taxname:558``) as the ``taxonID``
    field. The CLI's ``--root-id`` flag is ``type=str``, and
    ``prune_to_subgraph`` matches strings against the ``taxonID``
    column, but the schema previously rejected anything non-integer
    — leaving LSID pruning unreachable via config.
    """
    dwca = tmp_path / "archive.zip"
    dwca.write_bytes(b"placeholder")
    cfg = validate_config({
        "taxonomy": {
            "source": "dwca",
            "path": str(dwca),
            "root_id": "urn:lsid:marinespecies.org:taxname:558",
        },
    })
    assert cfg.taxonomy.root_id == "urn:lsid:marinespecies.org:taxname:558"


def test_taxonomy_dwca_requires_path():
    with pytest.raises(ValidationError) as exc:
        validate_config({"taxonomy": {"source": "dwca"}})
    msgs = [err["msg"] for err in exc.value.errors()]
    assert any("path" in m for m in msgs)


def test_pd_cutoff_must_be_non_negative():
    with pytest.raises(ValidationError):
        validate_config({"licensing": {"pd_cutoff_years": -1}})
