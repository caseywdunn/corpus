"""Tests for the workhorse dossier pair (#76, priority 2/6):

  - ``get_taxon_dossier(taxon, include=[...])``
  - ``get_chunks(paper_hash, chunk_ids=[...], with_text=...)``

These are the cache-friendly replacements for the
``search_taxon`` → ``get_papers_for_taxon`` → N× ``get_paper`` →
N× ``get_chunks_for_taxon`` → ``get_figures_for_taxon`` chain. The
contract these tests pin is the response shape — every dossier
section is opt-out via ``include``, every ordering is deterministic,
chunk_index returns IDs only (no full text — that's
``get_chunks``'s job).

Tests use a hand-built on-disk corpuscle so the per-paper JSON
reads exercise the real walk logic.
"""
from __future__ import annotations

import json
import types
from pathlib import Path
from typing import Dict, List

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools.chunks import get_chunks
from mcpsrv.tools.taxonomy import get_taxon_dossier


# ---------------------------------------------------------------------------
# Fixture: hand-built two-paper corpuscle on disk
# ---------------------------------------------------------------------------


def _write_paper(
    docs_dir: Path,
    hash_id: str,
    *,
    title: str,
    year: int,
    authors: List[Dict],
    chunks: List[Dict],
    taxa: List[Dict],
    mentions: List[Dict],
    figures: List[Dict],
    lexicons: Dict[str, List[Dict]] | None = None,
) -> Path:
    hash_dir = docs_dir / hash_id
    hash_dir.mkdir(parents=True)
    (hash_dir / "summary.json").write_text(json.dumps({
        "pdf_hash": hash_id, "processing_summary": {"stage_failures": []},
    }))
    (hash_dir / "metadata.json").write_text(json.dumps({
        "title": title, "year": year, "authors": authors,
    }))
    (hash_dir / "chunks.json").write_text(json.dumps({
        "chunks": chunks, "total_chunks": len(chunks),
    }))
    (hash_dir / "taxa.json").write_text(json.dumps({
        "taxa": taxa, "mentions": mentions, "unique_taxa": len(taxa),
    }))
    (hash_dir / "figures.json").write_text(json.dumps({
        "figures": figures, "total_figures": len(figures),
    }))
    for category, terms in (lexicons or {}).items():
        (hash_dir / f"{category}.json").write_text(json.dumps({
            "category": category, "terms": terms,
            "unique_terms": len(terms),
        }))
    return hash_dir


class _FakeTaxonomyDB:
    """Minimal stand-in matching the lookup() shape get_taxon_dossier needs."""

    def __init__(self, by_name: Dict[str, Dict]):
        self._by_name = {k.lower(): v for k, v in by_name.items()}

    def lookup(self, name: str):
        return self._by_name.get((name or "").strip().lower())


@pytest.fixture
def corpus(tmp_path: Path):
    """Two papers:
    - aaa (Dunn 2005): mentions Marrus (×3, in chunks c0+c1+c2), Agalma (×1, c2);
      one figure with 'Marrus' in caption; anatomy lexicon: pneumatophore, nectophore.
    - bbb (Pugh 2001): mentions Marrus (×1, c0); biogeography lexicon: pelagic.
    """
    output_dir = tmp_path / "output"
    docs = output_dir / "documents"
    docs.mkdir(parents=True)

    _write_paper(
        docs, "aaaaaaaaaaaa",
        title="Marrus claudanielis, new species",
        year=2005,
        authors=[{"surname": "Dunn", "forename": "Casey"}],
        chunks=[
            {"chunk_id": "c0", "section_class": "abstract",
             "headings": ["Abstract"],
             "text": "Marrus claudanielis is a deep-sea siphonophore.",
             "figure_refs": []},
            {"chunk_id": "c1", "section_class": "description",
             "headings": ["Description"],
             "text": "The nectosome of Marrus carries pneumatophores.",
             "figure_refs": ["fig1"]},
            {"chunk_id": "c2", "section_class": "discussion",
             "headings": ["Discussion"],
             "text": "Marrus and Agalma are siphonophores.",
             "figure_refs": []},
        ],
        taxa=[
            {"accepted_taxon_id": "t:marrus", "accepted_name": "Marrus",
             "rank": "genus", "mention_count": 3},
            {"accepted_taxon_id": "t:agalma", "accepted_name": "Agalma",
             "rank": "genus", "mention_count": 1},
        ],
        mentions=[
            {"accepted_taxon_id": "t:marrus", "chunk_id": "c0"},
            {"accepted_taxon_id": "t:marrus", "chunk_id": "c1"},
            {"accepted_taxon_id": "t:marrus", "chunk_id": "c2"},
            {"accepted_taxon_id": "t:agalma", "chunk_id": "c2"},
        ],
        figures=[
            {"figure_id": "fig1", "figure_type": "figure",
             "page": 4, "figure_number": 1,
             "caption_text": "Marrus claudanielis. Lateral view."},
            # An unclassified element that should NOT appear in figure_index.
            {"figure_id": "fig2", "figure_type": "graphical_element",
             "page": 1, "caption_text": ""},
        ],
        lexicons={
            "anatomy": [
                {"canonical": "pneumatophore", "mention_count": 3},
                {"canonical": "nectophore", "mention_count": 1},
            ],
        },
    )

    _write_paper(
        docs, "bbbbbbbbbbbb",
        title="The taxonomy of Erenna",
        year=2001,
        authors=[{"surname": "Pugh", "forename": "Philip"}],
        chunks=[
            {"chunk_id": "c0", "section_class": "introduction",
             "headings": ["Introduction"],
             "text": "Marrus is briefly discussed in the introduction.",
             "figure_refs": []},
        ],
        taxa=[
            {"accepted_taxon_id": "t:marrus", "accepted_name": "Marrus",
             "rank": "genus", "mention_count": 1},
        ],
        mentions=[
            {"accepted_taxon_id": "t:marrus", "chunk_id": "c0"},
        ],
        figures=[],
        lexicons={
            "biogeography": [
                {"canonical": "pelagic", "mention_count": 2},
            ],
        },
    )

    # In-memory index slice that mirrors what CorpusIndex.load() builds.
    papers = {
        "aaaaaaaaaaaa": {
            "hash": "aaaaaaaaaaaa",
            "title": "Marrus claudanielis, new species",
            "year": 2005,
            "authors": [{"surname": "Dunn", "forename": "Casey"}],
            "hash_dir": str(docs / "aaaaaaaaaaaa"),
        },
        "bbbbbbbbbbbb": {
            "hash": "bbbbbbbbbbbb",
            "title": "The taxonomy of Erenna",
            "year": 2001,
            "authors": [{"surname": "Pugh", "forename": "Philip"}],
            "hash_dir": str(docs / "bbbbbbbbbbbb"),
        },
    }
    taxon_to_papers = {
        "t:marrus": ["aaaaaaaaaaaa", "bbbbbbbbbbbb"],
        "t:agalma": ["aaaaaaaaaaaa"],
    }
    taxon_mention_counts = {
        "t:marrus": {"aaaaaaaaaaaa": 3, "bbbbbbbbbbbb": 1},
        "t:agalma": {"aaaaaaaaaaaa": 1},
    }
    taxon_display = {
        "t:marrus": {"accepted_name": "Marrus", "rank": "genus"},
        "t:agalma": {"accepted_name": "Agalma", "rank": "genus"},
    }
    taxonomy_db = _FakeTaxonomyDB({
        "Marrus": {"accepted_taxon_id": "t:marrus",
                   "accepted_name": "Marrus", "matched_name": "Marrus",
                   "rank": "genus", "authority": "Totton, 1954"},
    })

    fake = types.SimpleNamespace(
        papers=papers,
        taxon_to_papers=taxon_to_papers,
        taxon_mention_counts=taxon_mention_counts,
        taxon_display=taxon_display,
        taxonomy_db=taxonomy_db,
    )
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    yield fake
    mcp_app.set_index(original)


# ---------------------------------------------------------------------------
# get_chunks — drill-down pair
# ---------------------------------------------------------------------------


def test_get_chunks_returns_all_when_ids_none(corpus):
    rows = get_chunks(paper_hash="aaaaaaaaaaaa")
    assert [r["chunk_id"] for r in rows] == ["c0", "c1", "c2"]


def test_get_chunks_with_text_true_carries_full_text(corpus):
    rows = get_chunks(paper_hash="aaaaaaaaaaaa")
    assert "Marrus claudanielis is a deep-sea siphonophore." in rows[0]["text"]


def test_get_chunks_with_text_false_emits_len_chars_only(corpus):
    rows = get_chunks(paper_hash="aaaaaaaaaaaa", with_text=False)
    assert all("text" not in r for r in rows)
    assert all("len_chars" in r for r in rows)
    assert rows[0]["len_chars"] == len(
        "Marrus claudanielis is a deep-sea siphonophore.",
    )


def test_get_chunks_with_id_subset(corpus):
    """Caller drills down on specific chunk_ids from a dossier index."""
    rows = get_chunks(paper_hash="aaaaaaaaaaaa", chunk_ids=["c0", "c2"])
    assert [r["chunk_id"] for r in rows] == ["c0", "c2"]


def test_get_chunks_unknown_ids_silently_skipped(corpus):
    """A dossier might list ids that have rotated out — don't crash."""
    rows = get_chunks(paper_hash="aaaaaaaaaaaa",
                      chunk_ids=["c0", "does-not-exist"])
    assert [r["chunk_id"] for r in rows] == ["c0"]


def test_get_chunks_returns_paper_order_not_request_order(corpus):
    """chunk_ids is a set — the caller isn't requesting a sequence,
    they're requesting a slice. Paper-order is what makes the
    response coherent for an LLM reading multiple sections."""
    rows = get_chunks(paper_hash="aaaaaaaaaaaa", chunk_ids=["c2", "c0"])
    assert [r["chunk_id"] for r in rows] == ["c0", "c2"]


def test_get_chunks_carries_figure_refs(corpus):
    rows = get_chunks(paper_hash="aaaaaaaaaaaa", chunk_ids=["c1"])
    assert rows[0]["figure_refs"] == ["fig1"]


def test_get_chunks_unknown_paper_hash_errors(corpus):
    rows = get_chunks(paper_hash="ffffffffffff")
    assert "error" in rows[0]


# ---------------------------------------------------------------------------
# get_taxon_dossier — taxon metadata + structure
# ---------------------------------------------------------------------------


def test_dossier_returns_taxon_metadata(corpus):
    d = get_taxon_dossier("Marrus")
    assert d["taxon"]["taxon_id"] == "t:marrus"
    assert d["taxon"]["accepted_name"] == "Marrus"
    assert d["taxon"]["rank"] == "genus"
    assert d["n_papers_mentioning"] == 2


def test_dossier_taxon_not_found(corpus):
    d = get_taxon_dossier("nonexistent_genus")
    assert d["not_found"] is True
    assert d["queried"] == "nonexistent_genus"


def test_dossier_sparse_encodes_matched_name_when_synonym(corpus):
    """The fixture's lookup returns matched_name == accepted_name, so
    matched_name should be omitted from the output. The synonym path
    is covered by other tests if needed."""
    d = get_taxon_dossier("Marrus")
    assert "matched_name" not in d["taxon"]


# ---------------------------------------------------------------------------
# get_taxon_dossier — papers section
# ---------------------------------------------------------------------------


def test_dossier_papers_sorted_by_mention_count(corpus):
    d = get_taxon_dossier("Marrus")
    hashes = [p["hash"] for p in d["papers"]]
    assert hashes == ["aaaaaaaaaaaa", "bbbbbbbbbbbb"]  # 3 mentions, then 1


def test_dossier_papers_carries_first_author_year_n_mentions(corpus):
    d = get_taxon_dossier("Marrus")
    p0 = d["papers"][0]
    assert p0["first_author"] == "Dunn"
    assert p0["year"] == 2005
    assert p0["n_mentions"] == 3


# ---------------------------------------------------------------------------
# get_taxon_dossier — chunk_index (IDs + section + headings, no text)
# ---------------------------------------------------------------------------


def test_dossier_chunk_index_has_no_full_text(corpus):
    """The whole point of #76: the dossier returns IDs + headings, not
    text. Full text comes via get_chunks. A text key sneaking in
    would silently invalidate the cache-cost win."""
    d = get_taxon_dossier("Marrus")
    for entry in d["chunk_index"]:
        assert "text" not in entry


def test_dossier_chunk_index_lists_all_matching_chunks(corpus):
    d = get_taxon_dossier("Marrus")
    # Paper aaa: c0, c1, c2 all mention Marrus. Paper bbb: c0 only.
    expected = {
        ("aaaaaaaaaaaa", "c0"), ("aaaaaaaaaaaa", "c1"),
        ("aaaaaaaaaaaa", "c2"), ("bbbbbbbbbbbb", "c0"),
    }
    actual = {(c["paper_hash"], c["chunk_id"]) for c in d["chunk_index"]}
    assert actual == expected


def test_dossier_chunk_index_carries_section_and_headings(corpus):
    d = get_taxon_dossier("Marrus")
    c1 = next(c for c in d["chunk_index"]
              if c["paper_hash"] == "aaaaaaaaaaaa" and c["chunk_id"] == "c1")
    assert c1["section_class"] == "description"
    assert c1["headings"] == ["Description"]


def test_dossier_max_chunks_respected(corpus):
    d = get_taxon_dossier("Marrus", max_chunks=2)
    assert len(d["chunk_index"]) == 2


# ---------------------------------------------------------------------------
# get_taxon_dossier — figure_index
# ---------------------------------------------------------------------------


def test_dossier_figure_index_skips_graphical_elements(corpus):
    """fig2 is type=graphical_element and shouldn't surface — only
    figure / plate types belong in the index."""
    d = get_taxon_dossier("Marrus")
    fig_ids = {f["figure_id"] for f in d["figure_index"]}
    assert fig_ids == {"fig1"}


def test_dossier_figure_caption_match_flagged(corpus):
    d = get_taxon_dossier("Marrus")
    fig1 = d["figure_index"][0]
    assert fig1["caption_has_taxon"] is True


# ---------------------------------------------------------------------------
# get_taxon_dossier — lexicon_aggregated
# ---------------------------------------------------------------------------


def test_dossier_lexicon_aggregated_across_in_scope_papers(corpus):
    """The two in-scope papers ship distinct lexicon categories. Both
    should appear, with term-mention counts summed across papers."""
    d = get_taxon_dossier("Marrus")
    cats = sorted(d["lexicon_aggregated"].keys())
    assert cats == ["anatomy", "biogeography"]
    anatomy = d["lexicon_aggregated"]["anatomy"]
    pneumatophore = next(t for t in anatomy if t["term"] == "pneumatophore")
    assert pneumatophore["n_mentions"] == 3
    assert pneumatophore["n_papers"] == 1


# ---------------------------------------------------------------------------
# get_taxon_dossier — cooccurring_taxa
# ---------------------------------------------------------------------------


def test_dossier_cooccurring_taxa_includes_only_other_taxa(corpus):
    d = get_taxon_dossier("Marrus")
    ids = {c["taxon_id"] for c in d["cooccurring_taxa"]}
    assert "t:agalma" in ids
    assert "t:marrus" not in ids


def test_dossier_cooccurring_taxa_carries_shared_paper_count(corpus):
    d = get_taxon_dossier("Marrus")
    agalma = next(c for c in d["cooccurring_taxa"]
                  if c["taxon_id"] == "t:agalma")
    assert agalma["n_shared_papers"] == 1
    assert agalma["name"] == "Agalma"


# ---------------------------------------------------------------------------
# get_taxon_dossier — `include` filter
# ---------------------------------------------------------------------------


def test_include_filter_drops_other_sections(corpus):
    d = get_taxon_dossier("Marrus", include=["papers"])
    assert "papers" in d
    assert "chunk_index" not in d
    assert "figure_index" not in d
    assert "lexicon_aggregated" not in d
    assert "cooccurring_taxa" not in d


def test_include_filter_taxon_metadata_always_present(corpus):
    """include cannot drop the taxon-metadata + n_papers_mentioning
    headline — they're the cheap orientation fields."""
    d = get_taxon_dossier("Marrus", include=[])  # nothing requested
    assert "taxon" in d
    assert "n_papers_mentioning" in d
    assert "papers" not in d


def test_include_unknown_section_silently_ignored(corpus):
    d = get_taxon_dossier("Marrus", include=["papers", "made_up_section"])
    assert "papers" in d
    assert "made_up_section" not in d
