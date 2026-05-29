"""Tests for the figure-dossier pair (#76, priority 5/6):

  - ``get_figure_dossier_for_taxon(taxon)``
  - ``get_figure_dossier_for_term(category, term)``

The lever vs the existing ``get_figures_for_taxon`` /
``get_figures_for_lexicon_term`` is the ``linked_chunks`` field on
each figure: chunk_ids that reference the figure via
``chunks.json:figure_refs``. The LLM pulls one targeted
``get_chunks(paper_hash, chunk_ids=[...])`` instead of fishing
through ``get_chunks_for_taxon``.
"""
from __future__ import annotations

import json
import types
from pathlib import Path
from typing import Dict, List

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools.figures import (
    get_figure_dossier_for_taxon,
    get_figure_dossier_for_term,
)


def _write_paper(
    docs: Path,
    hash_id: str,
    *,
    chunks: List[Dict],
    figures: List[Dict],
    taxa_mentions: List[Dict] | None = None,
) -> None:
    hash_dir = docs / hash_id
    hash_dir.mkdir(parents=True)
    (hash_dir / "metadata.json").write_text(json.dumps({"title": hash_id}))
    (hash_dir / "chunks.json").write_text(json.dumps(
        {"chunks": chunks, "total_chunks": len(chunks)},
    ))
    (hash_dir / "figures.json").write_text(json.dumps(
        {"figures": figures, "total_figures": len(figures)},
    ))
    if taxa_mentions is not None:
        (hash_dir / "taxa.json").write_text(json.dumps({
            "mentions": taxa_mentions,
            "taxa": [], "unique_taxa": 0,
        }))


class _FakeTaxonomyDB:
    def __init__(self, by_name: Dict[str, Dict]):
        self._by_name = {k.lower(): v for k, v in by_name.items()}

    def lookup(self, name: str):
        return self._by_name.get((name or "").strip().lower())


@pytest.fixture
def corpus(tmp_path: Path):
    """Two papers:
    - aaa: figure fig1 (caption mentions Marrus), referenced by chunks
      c0 + c2. graphical_element fig2 (must NOT appear).
    - bbb: figure fig1 (caption doesn't mention Marrus), referenced by c0.
    """
    docs = tmp_path / "output" / "documents"
    docs.mkdir(parents=True)

    _write_paper(
        docs, "aaaaaaaaaaaa",
        chunks=[
            {"chunk_id": "c0", "section_class": "description",
             "headings": ["Description"], "figure_refs": ["fig1"]},
            {"chunk_id": "c1", "section_class": "discussion",
             "headings": [], "figure_refs": []},
            {"chunk_id": "c2", "section_class": "results",
             "headings": ["Results"], "figure_refs": ["fig1"]},
        ],
        figures=[
            {"figure_id": "fig1", "figure_type": "figure",
             "page": 4, "figure_number": 1,
             "filename": "fig1.png",
             "caption_text": "Marrus claudanielis. Lateral view of "
                             "pneumatophore and nectophores.",
             "panels_from_caption": [
                 {"label": "A", "description": "Whole animal"},
                 {"label": "B", "description": "Pneumatophore detail"},
             ],
             "panel_count_from_caption": 2,
             "rois": [{"label": "A", "bbox": [0, 0, 100, 100]}]},
            {"figure_id": "fig2", "figure_type": "graphical_element",
             "page": 1, "filename": "fig2.png", "caption_text": ""},
        ],
        taxa_mentions=[
            {"accepted_taxon_id": "t:marrus", "chunk_id": "c0"},
            {"accepted_taxon_id": "t:marrus", "chunk_id": "c1"},
        ],
    )

    _write_paper(
        docs, "bbbbbbbbbbbb",
        chunks=[
            {"chunk_id": "c0", "section_class": "abstract",
             "headings": ["Abstract"], "figure_refs": ["fig1"]},
        ],
        figures=[
            {"figure_id": "fig1", "figure_type": "figure",
             "page": 2, "figure_number": 1,
             "filename": "fig1.png",
             "caption_text": "An unrelated species illustration. "
                             "Pneumatophore tip shown."},
        ],
        taxa_mentions=[
            {"accepted_taxon_id": "t:marrus", "chunk_id": "c0"},
        ],
    )

    papers = {
        "aaaaaaaaaaaa": {"title": "Paper A", "year": 2005,
                         "hash_dir": str(docs / "aaaaaaaaaaaa")},
        "bbbbbbbbbbbb": {"title": "Paper B", "year": 2010,
                         "hash_dir": str(docs / "bbbbbbbbbbbb")},
    }
    taxon_to_papers = {"t:marrus": ["aaaaaaaaaaaa", "bbbbbbbbbbbb"]}
    taxon_mention_counts = {
        "t:marrus": {"aaaaaaaaaaaa": 2, "bbbbbbbbbbbb": 1},
    }
    lexicon_to_papers = {"anatomy": {}}  # category declared by the bundle
    taxonomy_db = _FakeTaxonomyDB({
        "Marrus": {"accepted_taxon_id": "t:marrus",
                   "accepted_name": "Marrus", "matched_name": "Marrus",
                   "rank": "genus"},
    })

    fake = types.SimpleNamespace(
        papers=papers,
        taxon_to_papers=taxon_to_papers,
        taxon_mention_counts=taxon_mention_counts,
        lexicon_to_papers=lexicon_to_papers,
        taxonomy_db=taxonomy_db,
    )
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    yield fake
    mcp_app.set_index(original)


# ---------------------------------------------------------------------------
# get_figure_dossier_for_taxon
# ---------------------------------------------------------------------------


def test_taxon_dossier_returns_real_figures_only(corpus):
    """graphical_element must NOT appear — only figure / plate."""
    out = get_figure_dossier_for_taxon("Marrus")
    types_seen = {f["figure_type"] for f in out["figures"]}
    assert "graphical_element" not in types_seen


def test_taxon_dossier_ranks_caption_match_above_mere_paper_mention(corpus):
    """fig1 of paper aaa has 'Marrus' in its caption (score 100+2);
    fig1 of paper bbb doesn't (score 0+1). aaa's figure should come
    first."""
    out = get_figure_dossier_for_taxon("Marrus")
    assert out["figures"][0]["paper_hash"] == "aaaaaaaaaaaa"


def test_taxon_dossier_linked_chunks_collected(corpus):
    """fig1 of paper aaa is referenced by c0 + c2 (not c1).
    linked_chunks must list both, with section + headings."""
    out = get_figure_dossier_for_taxon("Marrus")
    fig_a = next(f for f in out["figures"]
                 if f["paper_hash"] == "aaaaaaaaaaaa")
    chunk_ids = {c["chunk_id"] for c in fig_a["linked_chunks"]}
    assert chunk_ids == {"c0", "c2"}
    c0 = next(c for c in fig_a["linked_chunks"] if c["chunk_id"] == "c0")
    assert c0["section_class"] == "description"


def test_taxon_dossier_linked_chunks_carry_no_text(corpus):
    """The chunk drill-down pair (get_chunks) is responsible for
    full text. linked_chunks must surface IDs + section + headings
    only, otherwise the cache-cost lever is lost."""
    out = get_figure_dossier_for_taxon("Marrus")
    for f in out["figures"]:
        for c in f["linked_chunks"]:
            assert "text" not in c


def test_taxon_dossier_caption_has_taxon_flag(corpus):
    out = get_figure_dossier_for_taxon("Marrus")
    fig_a = next(f for f in out["figures"]
                 if f["paper_hash"] == "aaaaaaaaaaaa")
    fig_b = next(f for f in out["figures"]
                 if f["paper_hash"] == "bbbbbbbbbbbb")
    assert fig_a["caption_has_taxon"] is True
    assert fig_b["caption_has_taxon"] is False


def test_taxon_dossier_rois_summarized(corpus):
    """include_rois=True (default) emits panel count + labels +
    pixel-bbox count, not the full ROI list. Full ROIs stay
    accessible via list_figure_rois."""
    out = get_figure_dossier_for_taxon("Marrus")
    fig_a = next(f for f in out["figures"]
                 if f["paper_hash"] == "aaaaaaaaaaaa")
    assert fig_a["rois"]["panel_count_from_caption"] == 2
    assert fig_a["rois"]["panel_labels"] == ["A", "B"]
    assert fig_a["rois"]["n_rois_with_pixel_bbox"] == 1


def test_taxon_dossier_include_rois_false_omits(corpus):
    out = get_figure_dossier_for_taxon("Marrus", include_rois=False)
    for f in out["figures"]:
        assert "rois" not in f


def test_taxon_dossier_caption_preview_bounded(corpus):
    """Long captions should be trimmed to ~200 chars with an ellipsis."""
    # Re-write a paper with a long caption.
    p = corpus.papers["aaaaaaaaaaaa"]
    figs_path = Path(p["hash_dir"]) / "figures.json"
    data = json.loads(figs_path.read_text())
    data["figures"][0]["caption_text"] = "Marrus " + ("x" * 500)
    figs_path.write_text(json.dumps(data))

    out = get_figure_dossier_for_taxon("Marrus")
    fig_a = next(f for f in out["figures"]
                 if f["paper_hash"] == "aaaaaaaaaaaa")
    assert len(fig_a["caption_preview"]) <= 220
    assert fig_a["caption_preview"].endswith("…")


def test_taxon_dossier_not_found(corpus):
    out = get_figure_dossier_for_taxon("no_such_genus")
    assert out["not_found"] is True


def test_taxon_dossier_max_figures_respected(corpus):
    out = get_figure_dossier_for_taxon("Marrus", max_figures=1)
    assert len(out["figures"]) == 1


# ---------------------------------------------------------------------------
# get_figure_dossier_for_term
# ---------------------------------------------------------------------------


def test_term_dossier_returns_figures_with_caption_match(corpus):
    """pneumatophore appears in both captions; both figures should
    surface, with caption_match_count populated."""
    out = get_figure_dossier_for_term("anatomy", "pneumatophore")
    assert out["n_figures"] == 2
    for f in out["figures"]:
        assert f["caption_match_count"] >= 1


def test_term_dossier_ranks_by_caption_match_count(corpus):
    """Term that appears multiple times in a caption ranks above
    a single-mention caption (when fixtures provide it)."""
    p = corpus.papers["aaaaaaaaaaaa"]
    figs_path = Path(p["hash_dir"]) / "figures.json"
    data = json.loads(figs_path.read_text())
    # Make paper aaa's caption mention 'tentacle' three times.
    data["figures"][0]["caption_text"] = "Tentacle tentacle tentacle."
    figs_path.write_text(json.dumps(data))

    out = get_figure_dossier_for_term("anatomy", "tentacle")
    assert out["figures"][0]["paper_hash"] == "aaaaaaaaaaaa"
    assert out["figures"][0]["caption_match_count"] == 3


def test_term_dossier_linked_chunks_populated(corpus):
    out = get_figure_dossier_for_term("anatomy", "pneumatophore")
    fig_a = next(f for f in out["figures"]
                 if f["paper_hash"] == "aaaaaaaaaaaa")
    chunk_ids = {c["chunk_id"] for c in fig_a["linked_chunks"]}
    assert chunk_ids == {"c0", "c2"}


def test_term_dossier_unknown_category_error(corpus):
    out = get_figure_dossier_for_term("not_a_category", "pneumatophore")
    assert out["error"] == "unknown_category"
    assert "anatomy" in out["available"]


def test_term_dossier_empty_term_error(corpus):
    out = get_figure_dossier_for_term("anatomy", "")
    assert out["error"] == "empty_term"


def test_term_dossier_max_linked_chunks_respected(corpus):
    out = get_figure_dossier_for_term("anatomy", "pneumatophore",
                                       max_linked_chunks=1)
    fig_a = next(f for f in out["figures"]
                 if f["paper_hash"] == "aaaaaaaaaaaa")
    assert len(fig_a["linked_chunks"]) == 1
