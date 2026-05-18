"""Tests for the multi-category lexicon format.

The lexicon YAML is two-level: top-level keys are categories
(``anatomy``, ``biogeography``, …); each value is the canonical-term
map that drives extraction. ``load_lexicon`` returns that two-level
dict, and ``lexicon_fingerprints`` returns per-category content hashes
so resume can detect which categories changed.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from pipeline.taxa import (
    extract_lexicon_mentions,
    lexicon_fingerprints,
    load_lexicon,
)


# ---------------------------------------------------------------------------
# YAML format
# ---------------------------------------------------------------------------


@pytest.fixture
def two_category_yaml(tmp_path: Path) -> Path:
    p = tmp_path / "lex.yaml"
    p.write_text(
        yaml.safe_dump({
            "anatomy": {
                "nectophore": {
                    "synonyms": ["nectophores"],
                    "translations": {"de": ["Schwimmglocke"]},
                    "description": "Medusoid swimming bell.",
                },
            },
            "biogeography": {
                "pelagic": {
                    "synonyms": ["open water"],
                    "translations": {},
                    "description": "Open-ocean.",
                },
                "benthic": {
                    "synonyms": ["seafloor"],
                    "translations": {},
                    "description": "Bottom-associated.",
                },
            },
        }),
        encoding="utf-8",
    )
    return p


def test_load_lexicon_returns_two_level_dict(two_category_yaml):
    lex = load_lexicon(two_category_yaml)
    assert set(lex.keys()) == {"anatomy", "biogeography"}
    assert "nectophore" in lex["anatomy"]
    assert "pelagic" in lex["biogeography"]
    assert lex["biogeography"]["pelagic"]["synonyms"] == ["open water"]


def test_load_lexicon_normalizes_category_to_lowercase(tmp_path):
    p = tmp_path / "lex.yaml"
    p.write_text("Anatomy:\n  foo:\n    synonyms: []\n", encoding="utf-8")
    lex = load_lexicon(p)
    assert "anatomy" in lex
    assert "Anatomy" not in lex


def test_load_lexicon_rejects_non_mapping_root(tmp_path):
    p = tmp_path / "bad.yaml"
    p.write_text("- this is a list\n", encoding="utf-8")
    with pytest.raises(ValueError, match="root must be a mapping"):
        load_lexicon(p)


def test_load_lexicon_rejects_non_mapping_category(tmp_path):
    p = tmp_path / "bad.yaml"
    p.write_text("anatomy:\n  - foo\n", encoding="utf-8")
    with pytest.raises(ValueError, match="must be a mapping"):
        load_lexicon(p)


def test_load_lexicon_rejects_list_term_entry(tmp_path):
    """The natural shortcut `term: [synonym, synonym]` used to crash
    silently with `AttributeError: 'list' object has no attribute
    'get'`, which main.py's `except Exception` swallowed into a
    warning — leaving the whole annotation pass as a no-op. Now
    raises a ValueError with an actionable message showing the user
    the right shape.
    """
    p = tmp_path / "bad.yaml"
    p.write_text(
        "anatomy:\n"
        "  pneumatophore:\n"
        "    - pneumatophores\n"
        "    - float\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError) as exc:
        load_lexicon(p)
    msg = str(exc.value)
    assert "category 'anatomy'" in msg
    assert "term 'pneumatophore'" in msg
    assert "list" in msg          # surfaces the offending type
    assert "synonyms" in msg      # shows the corrected shape


def test_demo_lexicon_loads(tmp_path):
    """The shipped demo lexicon round-trips through load_lexicon."""
    demo = Path(__file__).resolve().parent.parent / "demo" / "lexicon.yaml"
    if not demo.exists():
        pytest.skip("demo/lexicon.yaml not present")
    lex = load_lexicon(demo)
    assert "anatomy" in lex
    assert len(lex["anatomy"]) > 10


# ---------------------------------------------------------------------------
# Fingerprints
# ---------------------------------------------------------------------------


def test_lexicon_fingerprints_per_category(two_category_yaml):
    fps = lexicon_fingerprints(two_category_yaml)
    assert set(fps.keys()) == {"anatomy", "biogeography"}
    assert fps["anatomy"]["sha256"] != fps["biogeography"]["sha256"]
    assert fps["anatomy"]["path"] == str(two_category_yaml)


def test_lexicon_fingerprints_unchanged_when_unrelated_category_edits(tmp_path):
    """A change to one section's content must not perturb another's hash."""
    p1 = tmp_path / "v1.yaml"
    p1.write_text(yaml.safe_dump({
        "anatomy": {"a": {"synonyms": []}},
        "biogeo": {"x": {"synonyms": []}},
    }), encoding="utf-8")
    p2 = tmp_path / "v2.yaml"
    p2.write_text(yaml.safe_dump({
        "anatomy": {"a": {"synonyms": []}},
        "biogeo": {"x": {"synonyms": ["NEW"]}},  # only this section changed
    }), encoding="utf-8")
    fp1, fp2 = lexicon_fingerprints(p1), lexicon_fingerprints(p2)
    assert fp1["anatomy"]["sha256"] == fp2["anatomy"]["sha256"]
    assert fp1["biogeo"]["sha256"] != fp2["biogeo"]["sha256"]


# ---------------------------------------------------------------------------
# Category labelling on extraction output
# ---------------------------------------------------------------------------


def test_extract_lexicon_mentions_stamps_category():
    section = {"pelagic": {"synonyms": [], "translations": {}, "description": ""}}
    chunks = [{"chunk_id": "c0", "text": "found in pelagic waters"}]
    res = extract_lexicon_mentions(chunks, section, category="biogeography")
    assert res["category"] == "biogeography"
    assert res["total_mentions"] == 1
    assert res["mentions"][0]["canonical"] == "pelagic"


def test_extract_lexicon_mentions_no_category_when_unset():
    section = {"pelagic": {"synonyms": [], "translations": {}, "description": ""}}
    chunks = [{"chunk_id": "c0", "text": "pelagic"}]
    res = extract_lexicon_mentions(chunks, section)
    assert "category" not in res


# ---------------------------------------------------------------------------
# _extract_taxa_and_lexicons with the new lexicons param
# ---------------------------------------------------------------------------


def test_extract_writes_one_json_per_category(tmp_path):
    from pipeline.annotate import _extract_taxa_and_lexicons

    hd = tmp_path / "abc"
    hd.mkdir()
    (hd / "chunks.json").write_text(json.dumps({
        "chunks": [
            {"chunk_id": "c0", "text": "Specimens collected in pelagic waters."},
            {"chunk_id": "c1", "text": "A nectophore was observed."},
        ],
    }))

    lexicons = {
        "anatomy": {
            "nectophore": {"synonyms": [], "translations": {}, "description": ""},
        },
        "biogeography": {
            "pelagic": {"synonyms": [], "translations": {}, "description": ""},
        },
    }
    out = _extract_taxa_and_lexicons(
        chunks_file=hd / "chunks.json",
        hash_dir=hd,
        taxonomy_db=None,
        lexicons=lexicons,
    )
    assert (hd / "anatomy.json") in out
    assert (hd / "biogeography.json") in out

    anat = json.loads((hd / "anatomy.json").read_text())
    biogeo = json.loads((hd / "biogeography.json").read_text())
    assert anat["category"] == "anatomy"
    assert biogeo["category"] == "biogeography"
    assert anat["total_mentions"] == 1
    assert biogeo["total_mentions"] == 1


def test_extract_skips_empty_category(tmp_path):
    """An empty section is silently skipped — no <category>.json written."""
    from pipeline.annotate import _extract_taxa_and_lexicons

    hd = tmp_path / "abc"
    hd.mkdir()
    (hd / "chunks.json").write_text(json.dumps({"chunks": []}))

    out = _extract_taxa_and_lexicons(
        chunks_file=hd / "chunks.json",
        hash_dir=hd,
        taxonomy_db=None,
        lexicons={"empty_cat": {}, "anatomy": {
            "x": {"synonyms": [], "translations": {}, "description": ""},
        }},
    )
    assert (hd / "anatomy.json") in out
    assert not (hd / "empty_cat.json").exists()


def test_fingerprint_stamped_per_category(tmp_path):
    from pipeline.annotate import _extract_taxa_and_lexicons

    hd = tmp_path / "abc"
    hd.mkdir()
    (hd / "chunks.json").write_text(json.dumps({"chunks": []}))

    fp = {"path": "/x/lex.yaml", "sha256": "deadbeef" * 8, "size": 100}
    _extract_taxa_and_lexicons(
        chunks_file=hd / "chunks.json",
        hash_dir=hd,
        taxonomy_db=None,
        lexicons={"biogeography": {
            "pelagic": {"synonyms": [], "translations": {}, "description": ""},
        }},
        lexicon_fingerprints={"biogeography": fp},
    )
    data = json.loads((hd / "biogeography.json").read_text())
    assert data["input_fingerprint"] == fp
