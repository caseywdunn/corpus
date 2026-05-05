"""Tests for multi-category lexicons (#24).

The matcher and extractor functions were generalized from anatomy-only
to category-agnostic. These tests confirm:
  1. Backwards-compat aliases (load_anatomy_lexicon, extract_anatomy_mentions)
     keep working.
  2. extract_lexicon_mentions stamps the category onto its output.
  3. _extract_taxa_and_anatomy writes <category>.json per configured lexicon.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from taxa import (
    extract_anatomy_mentions,
    extract_lexicon_mentions,
    load_anatomy_lexicon,
    load_lexicon,
)


# ---------------------------------------------------------------------------
# Aliases
# ---------------------------------------------------------------------------


def test_load_anatomy_lexicon_alias_is_load_lexicon():
    assert load_anatomy_lexicon is load_lexicon


def test_extract_anatomy_mentions_alias_is_extract_lexicon_mentions():
    assert extract_anatomy_mentions is extract_lexicon_mentions


# ---------------------------------------------------------------------------
# Category labelling
# ---------------------------------------------------------------------------


@pytest.fixture
def biogeo_lexicon() -> dict:
    return {
        "pelagic": {
            "synonyms": ["open_water", "open water"],
            "translations": {},
            "description": "Free-swimming, off the bottom.",
        },
        "benthic": {
            "synonyms": ["seafloor", "bottom-dwelling"],
            "translations": {},
            "description": "Bottom-associated.",
        },
    }


def test_extract_lexicon_mentions_stamps_category(biogeo_lexicon):
    chunks = [{"chunk_id": "c0", "text": "found in pelagic waters"}]
    res = extract_lexicon_mentions(chunks, biogeo_lexicon, category="biogeography")
    assert res["category"] == "biogeography"
    assert res["total_mentions"] == 1
    assert res["mentions"][0]["canonical"] == "pelagic"


def test_extract_lexicon_mentions_no_category_when_unset(biogeo_lexicon):
    chunks = [{"chunk_id": "c0", "text": "found in pelagic waters"}]
    res = extract_lexicon_mentions(chunks, biogeo_lexicon)
    assert "category" not in res


def test_empty_lexicon_returns_empty_with_category():
    res = extract_lexicon_mentions([{"chunk_id": "c0", "text": "x"}], {}, category="ecology")
    assert res["category"] == "ecology"
    assert res["total_mentions"] == 0


# ---------------------------------------------------------------------------
# _extract_taxa_and_anatomy with extra_lexicons
# ---------------------------------------------------------------------------


def test_extract_writes_extra_category_json(tmp_path, biogeo_lexicon):
    from process_corpus import _extract_taxa_and_anatomy

    hd = tmp_path / "abc"
    hd.mkdir()
    (hd / "chunks.json").write_text(json.dumps({
        "chunks": [
            {"chunk_id": "c0", "text": "Specimens collected in pelagic waters."},
            {"chunk_id": "c1", "text": "Benthic specimens were rare."},
        ],
    }))

    out = _extract_taxa_and_anatomy(
        chunks_file=hd / "chunks.json",
        hash_dir=hd,
        taxonomy_db=None,
        anatomy_lexicon=None,
        extra_lexicons={"biogeography": biogeo_lexicon},
    )
    assert (hd / "biogeography.json") in out
    data = json.loads((hd / "biogeography.json").read_text())
    assert data["category"] == "biogeography"
    assert data["total_mentions"] == 2
    canonicals = sorted(m["canonical"] for m in data["mentions"])
    assert canonicals == ["benthic", "pelagic"]


def test_extract_anatomy_lexicon_under_anatomy_takes_precedence_over_extra(tmp_path):
    """When --anatomy-lexicon is set AND --lexicon anatomy:... is also
    set (unlikely but possible), the explicit --anatomy-lexicon wins
    and the duplicate is silently ignored.
    """
    from process_corpus import _extract_taxa_and_anatomy

    hd = tmp_path / "abc"
    hd.mkdir()
    (hd / "chunks.json").write_text(json.dumps({
        "chunks": [{"chunk_id": "c0", "text": "nectophore observed"}],
    }))

    real = {"nectophore": {"synonyms": [], "translations": {}, "description": "real"}}
    duplicate = {"foo": {"synonyms": [], "translations": {}, "description": "should not run"}}

    _extract_taxa_and_anatomy(
        chunks_file=hd / "chunks.json",
        hash_dir=hd,
        taxonomy_db=None,
        anatomy_lexicon=real,
        extra_lexicons={"anatomy": duplicate},
    )
    data = json.loads((hd / "anatomy.json").read_text())
    assert data["mentions"], "anatomy.json should reflect --anatomy-lexicon, not the duplicate"
    assert data["mentions"][0]["canonical"] == "nectophore"


def test_extract_with_multiple_categories(tmp_path, biogeo_lexicon):
    """All configured categories are emitted as their own JSON files."""
    from process_corpus import _extract_taxa_and_anatomy

    hd = tmp_path / "abc"
    hd.mkdir()
    (hd / "chunks.json").write_text(json.dumps({
        "chunks": [
            {"chunk_id": "c0", "text": "A pelagic colony bears one nectophore."},
        ],
    }))
    anatomy = {"nectophore": {"synonyms": [], "translations": {}, "description": ""}}

    _extract_taxa_and_anatomy(
        chunks_file=hd / "chunks.json",
        hash_dir=hd,
        taxonomy_db=None,
        anatomy_lexicon=anatomy,
        extra_lexicons={"biogeography": biogeo_lexicon},
    )

    anat = json.loads((hd / "anatomy.json").read_text())
    biogeo = json.loads((hd / "biogeography.json").read_text())
    assert anat["category"] == "anatomy"
    assert biogeo["category"] == "biogeography"
    assert anat["total_mentions"] == 1
    assert biogeo["total_mentions"] == 1


def test_fingerprint_stamped_on_extra_category(tmp_path, biogeo_lexicon):
    from process_corpus import _extract_taxa_and_anatomy

    hd = tmp_path / "abc"
    hd.mkdir()
    (hd / "chunks.json").write_text(json.dumps({"chunks": []}))

    fp = {"path": "/x/biogeo.yaml", "sha256": "deadbeef" * 8, "size": 100}
    _extract_taxa_and_anatomy(
        chunks_file=hd / "chunks.json",
        hash_dir=hd,
        taxonomy_db=None,
        anatomy_lexicon=None,
        extra_lexicons={"biogeography": biogeo_lexicon},
        extra_fingerprints={"biogeography": fp},
    )
    data = json.loads((hd / "biogeography.json").read_text())
    assert data["input_fingerprint"] == fp
