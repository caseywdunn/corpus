"""Figure retrieval is synonym-aware (#143).

Ingestion always was: each mention in a per-paper `<category>.json`
records the `matched_text` found and the `canonical` term it belongs to,
so the stored data collapses synonyms correctly. Retrieval threw that
away and did a case-insensitive substring count of the *single string the
caller passed*, so:

* querying a canonical missed captions using a synonym (`wing` → `ala`),
* querying a synonym missed captions using the canonical (`ala` → `wing`),
* querying a synonym missed its siblings (`forewing` → `ala`),
* and the tool never even checked that the term was a known lexicon term.

All four directions are covered below. The expansion is built from
surface forms observed in the corpus's own text, because a distilled
bundle doesn't ship the lexicon YAML — so the last test pins the
documented degradation for a term nothing wrote.
"""
from __future__ import annotations

import json
import types
from pathlib import Path

from mcpsrv import app as mcp_app
from mcpsrv.tools.figures import (
    get_figure_dossier_for_term,
    get_figures_for_lexicon_term,
)

H1 = "aaaaaaaaaaaa"
H2 = "bbbbbbbbbbbb"


def _paper(tmp_path, h, *, captions, mentions):
    """Build one paper with figures whose captions use the given strings,
    and an anatomy.json recording the surface→canonical mapping that
    extraction would have produced."""
    hash_dir = tmp_path / "documents" / h
    (hash_dir / "figures").mkdir(parents=True)
    (hash_dir / "figures.json").write_text(json.dumps({"figures": [
        {
            "figure_id": f"docling_{i}",
            "filename": f"fig_{i}.png",
            "figure_type": "figure",
            "caption_text": cap,
        }
        for i, cap in enumerate(captions, start=1)
    ]}))
    (hash_dir / "chunks.json").write_text(json.dumps({"chunks": []}))
    (hash_dir / "anatomy.json").write_text(json.dumps({
        "category": "anatomy",
        "terms": sorted({m[1] for m in mentions}) and [
            {"canonical": c, "mention_count": 1}
            for c in sorted({m[1] for m in mentions})
        ],
        "mentions": [
            {"chunk_id": "chunk_0", "matched_text": surface, "canonical": canon}
            for surface, canon in mentions
        ],
    }))
    return hash_dir


def _make_index(tmp_path):
    """Two papers. One writes 'wing', the other only ever writes 'ala' —
    the exact shape that made the old literal match miss half the corpus."""
    d1 = _paper(
        tmp_path, H1,
        captions=["Figure 1. The wing venation in detail."],
        mentions=[("wing", "wing"), ("wings", "wing")],
    )
    d2 = _paper(
        tmp_path, H2,
        captions=["Figure 1. Ala and surrounding sclerites.",
                  "Figure 2. Unrelated habitus view."],
        mentions=[("ala", "wing"), ("forewing", "wing")],
    )
    idx = types.SimpleNamespace(
        papers={
            H1: {"hash_dir": str(d1), "title": "Paper one"},
            H2: {"hash_dir": str(d2), "title": "Paper two"},
        },
        lexicon_to_papers={"anatomy": {"wing": [H1, H2]}},
        lexicon_mention_counts={"anatomy": {"wing": {H1: 2, H2: 2}}},
        lexicon_surface_to_canonical={"anatomy": {
            "wing": "wing", "wings": "wing", "ala": "wing", "forewing": "wing",
        }},
        lexicon_canonical_surfaces={"anatomy": {
            "wing": {"wing", "wings", "ala", "forewing"},
        }},
        biblio_db=None,
        default_profile="report",
    )
    mcp_app.set_index(idx)
    return idx


def _hashes(rows):
    return {r["paper_hash"] for r in rows if "paper_hash" in r}


# --- get_figures_for_lexicon_term ---------------------------------------------


def test_canonical_query_finds_synonym_captions(tmp_path):
    """`wing` must find the paper whose caption only says `Ala`."""
    _make_index(tmp_path)
    rows = get_figures_for_lexicon_term("anatomy", "wing")
    assert H2 in _hashes(rows), (
        "querying the canonical missed a caption using a declared synonym"
    )
    assert H1 in _hashes(rows)


def test_synonym_query_finds_canonical_captions(tmp_path):
    """`ala` must find the paper whose caption says `wing`."""
    _make_index(tmp_path)
    rows = get_figures_for_lexicon_term("anatomy", "ala")
    assert H1 in _hashes(rows), (
        "querying a synonym missed a caption using the canonical term"
    )
    assert H2 in _hashes(rows)


def test_sibling_synonym_query(tmp_path):
    """`forewing` and `ala` are siblings — neither appears in the other's
    caption, but both resolve to `wing`."""
    _make_index(tmp_path)
    rows = get_figures_for_lexicon_term("anatomy", "forewing")
    assert {H1, H2} <= _hashes(rows)


def test_rows_report_which_surface_matched(tmp_path):
    _make_index(tmp_path)
    rows = get_figures_for_lexicon_term("anatomy", "wing")
    by_hash = {r["paper_hash"]: r for r in rows if "paper_hash" in r}
    assert by_hash[H1]["matched_surfaces"] == ["wing"]
    assert by_hash[H2]["matched_surfaces"] == ["ala"]
    # And the concept they resolved to.
    assert by_hash[H2]["canonical"] == "wing"


def test_offtopic_captions_are_not_matched(tmp_path):
    """Expansion must not become a firehose — the habitus figure has no
    surface form of the concept in its caption."""
    _make_index(tmp_path)
    rows = get_figures_for_lexicon_term("anatomy", "wing")
    ids = {(r["paper_hash"], r["figure_id"]) for r in rows if "paper_hash" in r}
    assert (H2, "docling_2") not in ids


def test_unknown_term_degrades_to_literal_and_says_so(tmp_path):
    """A term outside the lexicon still gets a substring search — but the
    caller is told no expansion happened, so an empty or thin result isn't
    mistaken for a synonym-aware one."""
    _make_index(tmp_path)
    rows = get_figures_for_lexicon_term("anatomy", "sclerites")
    notes = [r for r in rows if "note" in r]
    assert notes, rows
    assert notes[0]["resolved"] is False
    # The literal string does occur in paper two's caption.
    assert H2 in _hashes(rows)


def test_known_term_has_no_degradation_note(tmp_path):
    _make_index(tmp_path)
    rows = get_figures_for_lexicon_term("anatomy", "ala")
    assert not [r for r in rows if "note" in r]


# --- get_figure_dossier_for_term ----------------------------------------------


def test_dossier_expands_synonyms_too(tmp_path):
    """The dossier tool had the same literal-match bug at its own call
    site, so it needs its own coverage."""
    _make_index(tmp_path)
    out = get_figure_dossier_for_term("anatomy", "ala")
    assert out["canonical"] == "wing"
    assert out["resolved"] is True
    assert set(out["surfaces_searched"]) >= {"ala", "wing", "forewing"}
    assert out["n_figures"] >= 2, out
    hashes = {f["paper_hash"] for f in out["figures"]}
    assert H1 in hashes, "dossier missed the canonical-only paper"


def test_dossier_reports_matched_surfaces_per_figure(tmp_path):
    _make_index(tmp_path)
    out = get_figure_dossier_for_term("anatomy", "wing")
    by_hash = {f["paper_hash"]: f for f in out["figures"]}
    assert by_hash[H2]["matched_surfaces"] == ["ala"]


def test_dossier_unknown_term_reports_unresolved(tmp_path):
    _make_index(tmp_path)
    out = get_figure_dossier_for_term("anatomy", "sclerites")
    assert out["resolved"] is False
    assert out["canonical"] is None


def test_unknown_category_still_errors(tmp_path):
    _make_index(tmp_path)
    out = get_figure_dossier_for_term("biogeography", "wing")
    assert out["code"] == "invalid_argument"
    rows = get_figures_for_lexicon_term("biogeography", "wing")
    assert rows["code"] == "invalid_argument"
