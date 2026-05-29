"""Tests for get_taxon_lexicon_slice (#76, priority 3/6).

The slice tool's value over `corpus_summary` / `get_taxon_dossier`'s
`lexicon_aggregated` is **chunk-level co-occurrence**: terms only
count when they appear in the same chunk as the taxon. These tests
pin that semantic — a term that appears in the *paper* but not in
the chunks where the taxon is mentioned must NOT show up.
"""
from __future__ import annotations

import json
import types
from pathlib import Path
from typing import Dict, List

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools.taxonomy import get_taxon_lexicon_slice


def _write_paper(
    docs_dir: Path,
    hash_id: str,
    *,
    chunks: List[Dict],
    taxa_mentions: List[Dict],
    lexicon_mentions: Dict[str, List[Dict]],
) -> None:
    """Minimal per-paper writer — only the fields slice() reads."""
    hash_dir = docs_dir / hash_id
    hash_dir.mkdir(parents=True)
    (hash_dir / "metadata.json").write_text(json.dumps({"title": hash_id}))
    (hash_dir / "chunks.json").write_text(json.dumps({
        "chunks": chunks, "total_chunks": len(chunks),
    }))
    (hash_dir / "taxa.json").write_text(json.dumps({
        "mentions": taxa_mentions,
        "taxa": [], "unique_taxa": 0,
    }))
    for category, mentions in lexicon_mentions.items():
        (hash_dir / f"{category}.json").write_text(json.dumps({
            "category": category,
            "mentions": mentions,
            "terms": [],
        }))


class _FakeTaxonomyDB:
    def __init__(self, by_name: Dict[str, Dict]):
        self._by_name = {k.lower(): v for k, v in by_name.items()}

    def lookup(self, name: str):
        return self._by_name.get((name or "").strip().lower())


@pytest.fixture
def corpus(tmp_path: Path):
    """One paper with the key signal: taxon Marrus appears in chunks
    c0 + c1. Anatomy term `pneumatophore` appears in c0 (same as
    Marrus) AND c3 (without Marrus). The slice should count c0 only —
    chunk-level co-occurrence is the whole point.
    """
    output_dir = tmp_path / "output"
    docs = output_dir / "documents"
    docs.mkdir(parents=True)

    _write_paper(
        docs, "aaaaaaaaaaaa",
        chunks=[
            {"chunk_id": "c0", "text": ""},
            {"chunk_id": "c1", "text": ""},
            {"chunk_id": "c2", "text": ""},
            {"chunk_id": "c3", "text": ""},
        ],
        taxa_mentions=[
            {"accepted_taxon_id": "t:marrus", "chunk_id": "c0"},
            {"accepted_taxon_id": "t:marrus", "chunk_id": "c1"},
            # Note: Marrus does NOT appear in c2 or c3.
        ],
        lexicon_mentions={
            "anatomy": [
                # Same-chunk co-occurrence with Marrus
                {"chunk_id": "c0", "canonical": "pneumatophore"},
                {"chunk_id": "c1", "canonical": "nectophore"},
                # Token in chunks without Marrus — must NOT count
                {"chunk_id": "c2", "canonical": "pneumatophore"},
                {"chunk_id": "c3", "canonical": "tentacle"},
            ],
        },
    )

    # Second paper for the n_papers aggregation: Marrus in c0,
    # pneumatophore co-occurring.
    _write_paper(
        docs, "bbbbbbbbbbbb",
        chunks=[{"chunk_id": "c0", "text": ""}],
        taxa_mentions=[{"accepted_taxon_id": "t:marrus", "chunk_id": "c0"}],
        lexicon_mentions={
            "anatomy": [
                {"chunk_id": "c0", "canonical": "pneumatophore"},
            ],
        },
    )

    papers = {
        "aaaaaaaaaaaa": {"hash_dir": str(docs / "aaaaaaaaaaaa")},
        "bbbbbbbbbbbb": {"hash_dir": str(docs / "bbbbbbbbbbbb")},
    }
    taxon_to_papers = {"t:marrus": ["aaaaaaaaaaaa", "bbbbbbbbbbbb"]}
    # lexicon_to_papers tells us which categories the bundle declares —
    # the slice tool uses it for the unknown-category check.
    lexicon_to_papers = {"anatomy": {}}
    taxonomy_db = _FakeTaxonomyDB({
        "Marrus": {"accepted_taxon_id": "t:marrus",
                   "accepted_name": "Marrus", "rank": "genus"},
    })

    fake = types.SimpleNamespace(
        papers=papers,
        taxon_to_papers=taxon_to_papers,
        lexicon_to_papers=lexicon_to_papers,
        taxonomy_db=taxonomy_db,
        taxon_mention_counts={},
        taxon_display={},
    )
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    yield fake
    mcp_app.set_index(original)


# ---------------------------------------------------------------------------
# Same-chunk co-occurrence: the key semantic
# ---------------------------------------------------------------------------


def test_only_same_chunk_terms_count(corpus):
    """tentacle appears in c3 (no Marrus). Must not surface."""
    out = get_taxon_lexicon_slice("Marrus", "anatomy")
    terms = {t["term"] for t in out["terms"]}
    assert "tentacle" not in terms


def test_paper_local_co_occurrence_only(corpus):
    """pneumatophore appears in c2 of paper aaa (without Marrus) AND
    in c0 of paper aaa (with Marrus) AND in c0 of paper bbb (with
    Marrus). Only the two with-Marrus chunks should count: 2 chunks,
    2 papers."""
    out = get_taxon_lexicon_slice("Marrus", "anatomy")
    pneu = next(t for t in out["terms"] if t["term"] == "pneumatophore")
    assert pneu["n_chunks"] == 2
    assert pneu["n_papers"] == 2


def test_terms_sorted_by_chunk_count_desc(corpus):
    out = get_taxon_lexicon_slice("Marrus", "anatomy")
    counts = [t["n_chunks"] for t in out["terms"]]
    assert counts == sorted(counts, reverse=True)


def test_paper_examples_capped(corpus):
    out = get_taxon_lexicon_slice("Marrus", "anatomy", max_paper_examples=1)
    pneu = next(t for t in out["terms"] if t["term"] == "pneumatophore")
    assert len(pneu["paper_examples"]) == 1


# ---------------------------------------------------------------------------
# Headline counts
# ---------------------------------------------------------------------------


def test_n_chunks_total_counts_taxon_chunks_across_papers(corpus):
    out = get_taxon_lexicon_slice("Marrus", "anatomy")
    # Paper aaa: c0 + c1 = 2 chunks. Paper bbb: c0 = 1 chunk. Total 3.
    assert out["n_chunks_total"] == 3


def test_n_papers_total(corpus):
    out = get_taxon_lexicon_slice("Marrus", "anatomy")
    assert out["n_papers_total"] == 2


# ---------------------------------------------------------------------------
# Category-agnostic error shapes
# ---------------------------------------------------------------------------


def test_unknown_category_returns_available_list(corpus):
    out = get_taxon_lexicon_slice("Marrus", "made_up")
    assert out["error"] == "unknown_category"
    assert out["queried_category"] == "made_up"
    assert "anatomy" in out["available"]


def test_taxon_not_found(corpus):
    out = get_taxon_lexicon_slice("nonexistent_genus", "anatomy")
    assert out["not_found"] is True


def test_taxon_known_but_no_corpus_papers_returns_empty_terms(corpus):
    """An empty corpus for the taxon shouldn't crash — return the
    same shape with zero rollup counts."""
    corpus.taxon_to_papers["t:orphan"] = []
    corpus.taxonomy_db._by_name["orphan"] = {
        "accepted_taxon_id": "t:orphan", "accepted_name": "Orphan",
    }
    out = get_taxon_lexicon_slice("Orphan", "anatomy")
    assert out["terms"] == []
    assert out["n_chunks_total"] == 0
    assert out["n_papers_total"] == 0
