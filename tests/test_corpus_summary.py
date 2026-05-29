"""Tests for the corpus_summary MCP tool (#76).

Pins the shape + ordering contract the LLM relies on: stable sort
(count desc, name asc), bounded top-N lists, sparse encoding of
nulls. A drift here directly invalidates the prompt-cache hits that
make the tool useful in the first place.

The tool's whole value is "one call instead of 12× paginated
list_papers" — so the test exercises it against a hand-built
CorpusIndex with realistic shape (papers across multiple decades,
multiple lexicon categories, taxa with varying mention counts)
rather than trying to load a real build.
"""
from __future__ import annotations

import types
from collections import defaultdict
from typing import Dict, List

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools.papers import corpus_summary


def _make_index(papers, lexicon_to_papers, lexicon_mention_counts,
                taxon_to_papers, taxon_mention_counts, taxon_display,
                bundle_manifest=None):
    """Hand-roll a CorpusIndex-shaped namespace. We only populate the
    fields corpus_summary touches; everything else stays empty."""
    return types.SimpleNamespace(
        papers=papers,
        lexicon_to_papers=lexicon_to_papers,
        lexicon_mention_counts=lexicon_mention_counts,
        taxon_to_papers=taxon_to_papers,
        taxon_mention_counts=taxon_mention_counts,
        taxon_display=taxon_display,
        bundle_manifest=bundle_manifest,
    )


@pytest.fixture
def index():
    """Three papers, two decades, two lexicon categories, three taxa."""
    papers: Dict[str, Dict] = {
        "aaaaaaaaaaaa": {"year": 2005, "n_figures": 5},
        "bbbbbbbbbbbb": {"year": 2010, "n_figures": 3},
        "cccccccccccc": {"year": 2010, "n_figures": 0},
    }

    # anatomy category: pneumatophore (most-mentioned), nectophore
    lexicon_to_papers = defaultdict(lambda: defaultdict(list))
    lexicon_to_papers["anatomy"]["pneumatophore"] = ["aaaaaaaaaaaa", "bbbbbbbbbbbb"]
    lexicon_to_papers["anatomy"]["nectophore"] = ["aaaaaaaaaaaa"]
    lexicon_to_papers["biogeography"]["pelagic"] = ["bbbbbbbbbbbb", "cccccccccccc"]

    lexicon_mention_counts = defaultdict(lambda: defaultdict(dict))
    lexicon_mention_counts["anatomy"]["pneumatophore"] = {
        "aaaaaaaaaaaa": 7, "bbbbbbbbbbbb": 3,
    }
    lexicon_mention_counts["anatomy"]["nectophore"] = {"aaaaaaaaaaaa": 2}
    lexicon_mention_counts["biogeography"]["pelagic"] = {
        "bbbbbbbbbbbb": 4, "cccccccccccc": 1,
    }

    # Three taxa with varying mention counts to test ordering.
    taxon_to_papers = {
        "t:marrus":   ["aaaaaaaaaaaa", "cccccccccccc"],
        "t:agalma":   ["aaaaaaaaaaaa"],
        "t:hippopodius": ["bbbbbbbbbbbb"],
    }
    taxon_mention_counts = {
        "t:marrus":      {"aaaaaaaaaaaa": 10, "cccccccccccc": 5},
        "t:agalma":      {"aaaaaaaaaaaa": 2},
        "t:hippopodius": {"bbbbbbbbbbbb": 4},
    }
    taxon_display = {
        "t:marrus":      {"accepted_name": "Marrus", "rank": "genus"},
        "t:agalma":      {"accepted_name": "Agalma", "rank": "genus"},
        "t:hippopodius": {"accepted_name": "Hippopodius", "rank": "genus"},
    }

    fake = _make_index(papers, lexicon_to_papers, lexicon_mention_counts,
                       taxon_to_papers, taxon_mention_counts, taxon_display)
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    yield fake
    mcp_app.set_index(original)


# ---------------------------------------------------------------------------
# Headline counts
# ---------------------------------------------------------------------------


def test_n_papers(index):
    assert corpus_summary()["n_papers"] == 3


def test_year_range(index):
    out = corpus_summary()
    assert out["year_range"] == {"min": 2005, "max": 2010}


def test_by_decade_sorted_ascending(index):
    out = corpus_summary()
    assert out["by_decade"] == [
        {"decade": 2000, "n_papers": 1},
        {"decade": 2010, "n_papers": 2},
    ]


def test_n_figures_total_sums_across_papers(index):
    assert corpus_summary()["n_figures_total"] == 8  # 5 + 3 + 0


def test_n_unique_taxa(index):
    assert corpus_summary()["n_unique_taxa"] == 3


# ---------------------------------------------------------------------------
# Lexicon coverage
# ---------------------------------------------------------------------------


def test_lexicon_categories_sorted(index):
    out = corpus_summary()
    assert out["lexicon_categories"] == ["anatomy", "biogeography"]


def test_lexicon_coverage_anatomy(index):
    out = corpus_summary()
    cov = out["lexicon_coverage"]["anatomy"]
    assert cov["n_terms_hit"] == 2
    assert cov["n_papers_with_hits"] == 2  # aaaa + bbbb
    assert cov["n_mentions_total"] == 12   # 7 + 3 + 2


def test_lexicon_top_terms_sorted_by_mention_count_desc(index):
    out = corpus_summary()
    anatomy_terms = [t["term"] for t in out["lexicon_coverage"]["anatomy"]["top_terms"]]
    assert anatomy_terms == ["pneumatophore", "nectophore"]


def test_top_terms_cap_respected(index):
    out = corpus_summary(top_terms_per_category=1)
    anatomy_terms = out["lexicon_coverage"]["anatomy"]["top_terms"]
    assert len(anatomy_terms) == 1
    assert anatomy_terms[0]["term"] == "pneumatophore"


def test_top_terms_zero_omits(index):
    """Caller asking for headline counts only — top_terms is the
    bulky per-category field that should drop out cleanly."""
    out = corpus_summary(top_terms_per_category=0)
    assert out["lexicon_coverage"]["anatomy"]["top_terms"] == []
    # Headline counts still present.
    assert out["lexicon_coverage"]["anatomy"]["n_terms_hit"] == 2


# ---------------------------------------------------------------------------
# Top taxa
# ---------------------------------------------------------------------------


def test_top_taxa_sorted_by_mention_count_desc(index):
    out = corpus_summary()
    names = [t["name"] for t in out["top_taxa"]]
    assert names == ["Marrus", "Hippopodius", "Agalma"]  # 15, 4, 2


def test_top_taxa_cap_respected(index):
    out = corpus_summary(top_taxa=1)
    assert len(out["top_taxa"]) == 1
    assert out["top_taxa"][0]["name"] == "Marrus"


def test_top_taxa_zero_omits(index):
    out = corpus_summary(top_taxa=0)
    assert out["top_taxa"] == []
    # Headline still present.
    assert out["n_unique_taxa"] == 3


def test_top_taxa_carries_paper_and_mention_counts(index):
    out = corpus_summary()
    marrus = next(t for t in out["top_taxa"] if t["name"] == "Marrus")
    assert marrus["n_papers"] == 2
    assert marrus["n_mentions"] == 15
    assert marrus["rank"] == "genus"


def test_top_taxa_sparse_encoding_drops_null_rank(index):
    """If a taxon has no rank, the field is omitted from the entry
    rather than emitted as null. Saves tokens in dense top-N lists."""
    # Re-build index with one rankless taxon.
    index.taxon_display["t:agalma"] = {"accepted_name": "Agalma"}  # no rank
    out = corpus_summary()
    agalma = next(t for t in out["top_taxa"] if t["name"] == "Agalma")
    assert "rank" not in agalma


# ---------------------------------------------------------------------------
# Bundle identity
# ---------------------------------------------------------------------------


def test_bundle_version_null_for_build_output(index):
    out = corpus_summary()
    assert out["bundle_version"] is None


def test_bundle_version_carried_from_manifest(tmp_path):
    """When the index has a bundle_manifest, corpus_summary surfaces
    bundle_version so the LLM can cite the corpus version it's
    answering against."""
    fake = _make_index(
        papers={}, lexicon_to_papers={}, lexicon_mention_counts={},
        taxon_to_papers={}, taxon_mention_counts={}, taxon_display={},
        bundle_manifest={"bundle_version": "v0.5.0"},
    )
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    try:
        out = corpus_summary()
        assert out["bundle_version"] == "v0.5.0"
    finally:
        mcp_app.set_index(original)


# ---------------------------------------------------------------------------
# Empty-corpus edge cases
# ---------------------------------------------------------------------------


def test_empty_corpus_returns_zero_counts_and_null_year_range():
    fake = _make_index(
        papers={}, lexicon_to_papers={}, lexicon_mention_counts={},
        taxon_to_papers={}, taxon_mention_counts={}, taxon_display={},
    )
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    try:
        out = corpus_summary()
        assert out["n_papers"] == 0
        assert out["year_range"] is None
        assert out["by_decade"] == []
        assert out["lexicon_categories"] == []
        assert out["lexicon_coverage"] == {}
        assert out["top_taxa"] == []
        assert out["n_figures_total"] == 0
    finally:
        mcp_app.set_index(original)


def test_papers_without_year_skipped_from_year_aggregations():
    """A scanned old paper with no Grobid year shouldn't drag year_range
    or by_decade — the field is None and gets filtered out."""
    fake = _make_index(
        papers={
            "withhash00000": {"year": 2010, "n_figures": 0},
            "noyearhash000": {"year": None, "n_figures": 0},
        },
        lexicon_to_papers={}, lexicon_mention_counts={},
        taxon_to_papers={}, taxon_mention_counts={}, taxon_display={},
    )
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    try:
        out = corpus_summary()
        assert out["n_papers"] == 2
        assert out["year_range"] == {"min": 2010, "max": 2010}
        assert out["by_decade"] == [{"decade": 2010, "n_papers": 1}]
    finally:
        mcp_app.set_index(original)
