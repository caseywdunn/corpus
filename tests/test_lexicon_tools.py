"""Tests for lexicon_matrix and get_lexicon_term_dossier (#76, priority 4/6).

The two tools share scoring contract: column / row ordering is
deterministic so prompt caches stay warm across runs, sparse rows
(papers with no <category>.json) encode as zero-filled vectors,
and unknown-category / unknown-term return structured errors with
the available list (#76 principle 7 — category-agnostic).
"""
from __future__ import annotations

import json
import types
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools.lexicon import get_lexicon_term_dossier, lexicon_matrix


def _write_paper(
    docs: Path,
    hash_id: str,
    *,
    title: str,
    year: int,
    chunks: List[Dict],
    lexicon: Dict[str, Dict[str, int]] | None = None,
    mentions: Dict[str, List[Dict]] | None = None,
    descriptions: Dict[str, Dict[str, str]] | None = None,
) -> Path:
    """Write per-paper artifacts under ``docs/<hash_id>/`` with a
    minimal shape sufficient for the two lexicon tools.

    ``lexicon`` maps category → {term: mention_count}; one
    <category>.json is written per category. ``mentions`` is the
    per-category list of {chunk_id, canonical} (defaults to empty if
    unset — the matrix tool only needs term-rollups, not mentions).
    """
    hash_dir = docs / hash_id
    hash_dir.mkdir(parents=True)
    (hash_dir / "metadata.json").write_text(json.dumps(
        {"title": title, "year": year},
    ))
    (hash_dir / "chunks.json").write_text(json.dumps(
        {"chunks": chunks, "total_chunks": len(chunks)},
    ))
    for category, term_counts in (lexicon or {}).items():
        cat_mentions = (mentions or {}).get(category, [])
        cat_descs = (descriptions or {}).get(category, {})
        (hash_dir / f"{category}.json").write_text(json.dumps({
            "category": category,
            "terms": [
                {
                    "canonical": term,
                    "mention_count": count,
                    "description": cat_descs.get(term, ""),
                }
                for term, count in term_counts.items()
            ],
            "mentions": cat_mentions,
        }))
    return hash_dir


@pytest.fixture
def corpus(tmp_path: Path):
    """Three papers across two years, one shared lexicon category."""
    docs = tmp_path / "output" / "documents"
    docs.mkdir(parents=True)

    _write_paper(
        docs, "aaaaaaaaaaaa",
        title="Paper A on Marrus", year=2005,
        chunks=[{"chunk_id": "c0", "section_class": "description",
                 "headings": ["Description"]}],
        lexicon={"anatomy": {"pneumatophore": 5, "nectophore": 2}},
        mentions={"anatomy": [
            {"chunk_id": "c0", "canonical": "pneumatophore"},
            {"chunk_id": "c0", "canonical": "nectophore"},
        ]},
        descriptions={"anatomy": {
            "pneumatophore": "Apical gas-filled float.",
            "nectophore": "A swimming bell.",
        }},
    )

    _write_paper(
        docs, "bbbbbbbbbbbb",
        title="Paper B on Erenna", year=2010,
        chunks=[{"chunk_id": "c0", "section_class": "abstract",
                 "headings": ["Abstract"]}],
        lexicon={"anatomy": {"pneumatophore": 1, "tentacle": 3}},
        mentions={"anatomy": [
            {"chunk_id": "c0", "canonical": "pneumatophore"},
            {"chunk_id": "c0", "canonical": "tentacle"},
        ]},
    )

    # Paper with no <category>.json at all — must surface as a
    # zero-filled row in the matrix, not crash.
    _write_paper(
        docs, "cccccccccccc",
        title="Paper C, no anatomy section", year=2010,
        chunks=[{"chunk_id": "c0"}],
    )

    papers = {
        "aaaaaaaaaaaa": {"title": "Paper A on Marrus", "year": 2005,
                         "hash_dir": str(docs / "aaaaaaaaaaaa")},
        "bbbbbbbbbbbb": {"title": "Paper B on Erenna", "year": 2010,
                         "hash_dir": str(docs / "bbbbbbbbbbbb")},
        "cccccccccccc": {"title": "Paper C, no anatomy section",
                         "year": 2010,
                         "hash_dir": str(docs / "cccccccccccc")},
    }
    lexicon_to_papers = defaultdict(lambda: defaultdict(list))
    lexicon_to_papers["anatomy"]["pneumatophore"] = ["aaaaaaaaaaaa", "bbbbbbbbbbbb"]
    lexicon_to_papers["anatomy"]["nectophore"] = ["aaaaaaaaaaaa"]
    lexicon_to_papers["anatomy"]["tentacle"] = ["bbbbbbbbbbbb"]
    lexicon_mention_counts = defaultdict(lambda: defaultdict(dict))
    lexicon_mention_counts["anatomy"]["pneumatophore"] = {
        "aaaaaaaaaaaa": 5, "bbbbbbbbbbbb": 1,
    }
    lexicon_mention_counts["anatomy"]["nectophore"] = {"aaaaaaaaaaaa": 2}
    lexicon_mention_counts["anatomy"]["tentacle"] = {"bbbbbbbbbbbb": 3}

    fake = types.SimpleNamespace(
        papers=papers,
        lexicon_to_papers=lexicon_to_papers,
        lexicon_mention_counts=lexicon_mention_counts,
    )
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    yield fake
    mcp_app.set_index(original)


# ---------------------------------------------------------------------------
# lexicon_matrix
# ---------------------------------------------------------------------------


def test_matrix_default_columns_are_top_n_by_mention_count(corpus):
    out = lexicon_matrix("anatomy", top_n=2)
    # pneumatophore: 5+1=6, tentacle: 3, nectophore: 2 → top-2 is
    # pneumatophore, tentacle.
    assert out["terms"] == ["pneumatophore", "tentacle"]


def test_matrix_caller_columns_preserve_order_and_drop_unknown(corpus):
    out = lexicon_matrix("anatomy", terms=["tentacle", "no_such_term", "nectophore"])
    assert out["terms"] == ["tentacle", "nectophore"]


def test_matrix_rows_sorted_by_year_desc_then_hash(corpus):
    out = lexicon_matrix("anatomy", terms=["pneumatophore"])
    row_hashes = [r["paper_hash"] for r in out["rows"]]
    # 2010 papers first (b, c by hash), then 2005 (a).
    assert row_hashes == ["bbbbbbbbbbbb", "cccccccccccc", "aaaaaaaaaaaa"]


def test_matrix_paper_with_no_category_file_is_zero_row(corpus):
    """Paper cccc has no anatomy.json. The contract is that it still
    appears in the matrix as a zero-filled row — papers shouldn't
    silently disappear from coverage tables just because they have
    no hits in the chosen category."""
    out = lexicon_matrix("anatomy", terms=["pneumatophore", "tentacle"])
    row_c = next(r for r in out["rows"] if r["paper_hash"] == "cccccccccccc")
    assert row_c["counts"] == [0, 0]


def test_matrix_counts_match_per_paper_mention_counts(corpus):
    out = lexicon_matrix("anatomy", terms=["pneumatophore", "nectophore"])
    row_a = next(r for r in out["rows"] if r["paper_hash"] == "aaaaaaaaaaaa")
    assert row_a["counts"] == [5, 2]


def test_matrix_year_filter_restricts_rows(corpus):
    out = lexicon_matrix("anatomy", terms=["pneumatophore"],
                         year_from=2010, year_to=2010)
    hashes = {r["paper_hash"] for r in out["rows"]}
    assert hashes == {"bbbbbbbbbbbb", "cccccccccccc"}


def test_matrix_paper_hashes_filter_restricts_rows(corpus):
    out = lexicon_matrix("anatomy", terms=["pneumatophore"],
                         paper_hashes=["aaaaaaaaaaaa", "bogus"])
    hashes = {r["paper_hash"] for r in out["rows"]}
    assert hashes == {"aaaaaaaaaaaa"}


def test_matrix_unknown_category_error(corpus):
    out = lexicon_matrix("not_a_category")
    assert out["error"] == "unknown_category"
    assert "anatomy" in out["available"]


# ---------------------------------------------------------------------------
# get_lexicon_term_dossier
# ---------------------------------------------------------------------------


def test_term_dossier_carries_rollup_counts(corpus):
    out = get_lexicon_term_dossier("anatomy", "pneumatophore")
    assert out["n_papers"] == 2
    assert out["n_mentions_total"] == 6


def test_term_dossier_papers_sorted_by_n_mentions_desc(corpus):
    out = get_lexicon_term_dossier("anatomy", "pneumatophore")
    hashes = [p["hash"] for p in out["papers"]]
    assert hashes == ["aaaaaaaaaaaa", "bbbbbbbbbbbb"]


def test_term_dossier_chunk_examples_carry_section_and_headings(corpus):
    out = get_lexicon_term_dossier("anatomy", "pneumatophore")
    # First example is from paper aaa (more mentions), chunk c0.
    e0 = out["chunk_examples"][0]
    assert e0["paper_hash"] == "aaaaaaaaaaaa"
    assert e0["section_class"] == "description"


def test_term_dossier_description_pulled_from_first_paper(corpus):
    out = get_lexicon_term_dossier("anatomy", "pneumatophore")
    assert out["description"] == "Apical gas-filled float."


def test_term_dossier_unknown_term_error(corpus):
    out = get_lexicon_term_dossier("anatomy", "no_such_term")
    assert out["error"] == "unknown_term"
    assert out["queried_term"] == "no_such_term"
    assert out["category"] == "anatomy"


def test_term_dossier_unknown_category_error(corpus):
    out = get_lexicon_term_dossier("not_a_category", "pneumatophore")
    assert out["error"] == "unknown_category"
    assert "anatomy" in out["available"]


def test_term_dossier_chunk_examples_capped(corpus):
    out = get_lexicon_term_dossier("anatomy", "pneumatophore",
                                   max_chunk_examples=1)
    assert len(out["chunk_examples"]) == 1
