"""Tests for #76 priority 6/6: get_papers + get_taxon_subtree_dossier.

The pair closes out the #76 priority list. get_papers is the
trivial batched lookup the suite has been imagining since priority
1; get_taxon_subtree_dossier is the monographic-prompt workhorse
that replaces N×get_papers_for_taxon enumerations under a clade.
"""
from __future__ import annotations

import json
import sqlite3
import types
from pathlib import Path
from typing import Dict, List

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools.papers import get_papers
from mcpsrv.tools.taxonomy import get_taxon_subtree_dossier


# ---------------------------------------------------------------------------
# get_papers fixtures + tests
# ---------------------------------------------------------------------------


def _write_paper_for_get_papers(docs: Path, hash_id: str, *, title: str,
                                year: int) -> None:
    hash_dir = docs / hash_id
    hash_dir.mkdir(parents=True)
    (hash_dir / "metadata.json").write_text(json.dumps({"title": title}))
    (hash_dir / "taxa.json").write_text(json.dumps({"taxa": []}))
    (hash_dir / "chunks.json").write_text(json.dumps({"chunks": []}))


@pytest.fixture
def papers_index(tmp_path: Path):
    docs = tmp_path / "output" / "documents"
    docs.mkdir(parents=True)
    _write_paper_for_get_papers(docs, "aaaaaaaaaaaa",
                                title="Paper A", year=2005)
    _write_paper_for_get_papers(docs, "bbbbbbbbbbbb",
                                title="Paper B", year=2010)
    papers = {
        "aaaaaaaaaaaa": {
            "hash": "aaaaaaaaaaaa", "title": "Paper A", "year": 2005,
            "authors": [{"surname": "Dunn", "forename": "Casey"}],
            "doi": "10.1/A", "n_chunks": 0, "n_figures": 0,
            "hash_dir": str(docs / "aaaaaaaaaaaa"),
        },
        "bbbbbbbbbbbb": {
            "hash": "bbbbbbbbbbbb", "title": "Paper B", "year": 2010,
            "authors": [],
            "doi": None, "n_chunks": 0, "n_figures": 0,
            "hash_dir": str(docs / "bbbbbbbbbbbb"),
        },
    }
    fake = types.SimpleNamespace(papers=papers)
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    yield fake
    mcp_app.set_index(original)


def test_get_papers_returns_input_order(papers_index):
    out = get_papers(hashes=["bbbbbbbbbbbb", "aaaaaaaaaaaa"])
    assert [r["hash"] for r in out] == ["bbbbbbbbbbbb", "aaaaaaaaaaaa"]


def test_get_papers_unknown_hash_surfaces_error_in_place(papers_index):
    out = get_papers(hashes=["aaaaaaaaaaaa", "fffffffffff0", "bbbbbbbbbbbb"])
    assert out[1] == {"hash": "fffffffffff0", "error": "not_found"}
    assert out[0]["title"] == "Paper A"
    assert out[2]["title"] == "Paper B"


def test_get_papers_field_whitelist_trims(papers_index):
    out = get_papers(
        hashes=["aaaaaaaaaaaa"], fields=["title", "year", "first_author"],
    )
    row = out[0]
    # Always-included hash plus requested fields only.
    assert set(row.keys()) == {"hash", "title", "year", "first_author"}
    assert row["first_author"] == "Dunn"


def test_get_papers_synthesizes_first_author(papers_index):
    """Caller convenience: get_papers exposes first_author so the
    LLM doesn't walk authors[0].surname every time."""
    out = get_papers(hashes=["aaaaaaaaaaaa"])
    assert out[0]["first_author"] == "Dunn"


def test_get_papers_first_author_none_when_no_authors(papers_index):
    out = get_papers(hashes=["bbbbbbbbbbbb"])
    assert out[0]["first_author"] is None


# ---------------------------------------------------------------------------
# get_taxon_subtree_dossier fixtures + tests
# ---------------------------------------------------------------------------


class _SqliteTaxonomyDB:
    """Wraps a sqlite conn with a `taxa` table matching the DwC
    schema get_taxon_subtree_dossier walks. Includes a lookup()
    method so the tool's resolver works."""

    def __init__(self, conn: sqlite3.Connection, lookups: Dict[str, Dict]):
        self.conn = conn
        self._lookups = {k.lower(): v for k, v in lookups.items()}

    def lookup(self, name: str):
        return self._lookups.get((name or "").strip().lower())


@pytest.fixture
def subtree_index(tmp_path: Path):
    """A two-level taxonomy: genus Marrus, two species claudanielis
    + orthocannoides. Plus an unrelated genus Erenna with one species,
    so we can verify the subtree filter doesn't leak."""
    db = sqlite3.connect(":memory:")
    db.executescript("""
        CREATE TABLE taxa (
          taxon_id TEXT PRIMARY KEY,
          scientific_name TEXT,
          scientific_name_authorship TEXT,
          taxon_rank TEXT,
          taxonomic_status TEXT,
          parent_name_usage_id TEXT
        );
    """)
    db.executemany(
        "INSERT INTO taxa VALUES (?, ?, ?, ?, ?, ?)",
        [
            ("t:marrus",        "Marrus",                None,
             "genus",   "accepted", None),
            ("t:claudanielis",  "Marrus claudanielis",   "Dunn, 2005",
             "species", "accepted", "t:marrus"),
            ("t:orthocannoides","Marrus orthocannoides", "Totton, 1954",
             "species", "accepted", "t:marrus"),
            ("t:erenna",        "Erenna",                None,
             "genus",   "accepted", None),
            ("t:richardi",      "Erenna richardi",       "Bedot, 1904",
             "species", "accepted", "t:erenna"),
            # A subspecies under claudanielis to confirm subspecies inclusion.
            ("t:cl-sub",        "Marrus claudanielis subsp. fake",
             "Smith, 2099", "subspecies", "accepted", "t:claudanielis"),
        ],
    )
    db.commit()

    papers = {
        "aaaaaaaaaaaa": {"title": "Paper A", "year": 2005,
                         "hash_dir": "/dev/null/aaa"},
        "bbbbbbbbbbbb": {"title": "Paper B", "year": 2010,
                         "hash_dir": "/dev/null/bbb"},
        "cccccccccccc": {"title": "Paper C (Erenna)", "year": 2012,
                         "hash_dir": "/dev/null/ccc"},
    }
    taxon_to_papers = {
        "t:claudanielis": ["aaaaaaaaaaaa", "bbbbbbbbbbbb"],
        "t:orthocannoides": ["aaaaaaaaaaaa"],
        "t:richardi": ["cccccccccccc"],
    }
    taxon_mention_counts = {
        "t:claudanielis":   {"aaaaaaaaaaaa": 5, "bbbbbbbbbbbb": 2},
        "t:orthocannoides": {"aaaaaaaaaaaa": 1},
        "t:richardi":       {"cccccccccccc": 3},
    }
    taxonomy_db = _SqliteTaxonomyDB(db, lookups={
        "Marrus": {"accepted_taxon_id": "t:marrus",
                   "accepted_name": "Marrus", "rank": "genus"},
        "Erenna": {"accepted_taxon_id": "t:erenna",
                   "accepted_name": "Erenna", "rank": "genus"},
        "EmptyClade": {"accepted_taxon_id": "t:empty",
                       "accepted_name": "EmptyClade", "rank": "genus"},
    })

    fake = types.SimpleNamespace(
        papers=papers,
        taxon_to_papers=taxon_to_papers,
        taxon_mention_counts=taxon_mention_counts,
        taxonomy_db=taxonomy_db,
        taxon_display={},
        lexicon_to_papers={},
    )
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    yield fake
    mcp_app.set_index(original)


def test_subtree_counts_descendants(subtree_index):
    """Marrus has two species + one subspecies — all 'accepted' rank
    in {species, subspecies}, all should count."""
    out = get_taxon_subtree_dossier("Marrus")
    assert out["n_species_in_subtree"] == 3


def test_subtree_species_only_with_corpus_coverage(subtree_index):
    """The subspecies has no papers; species_with_corpus_coverage
    should include only the two species that DO have papers."""
    out = get_taxon_subtree_dossier("Marrus")
    species_ids = {s["taxon_id"] for s in out["species_with_corpus_coverage"]}
    assert species_ids == {"t:claudanielis", "t:orthocannoides"}


def test_subtree_species_sorted_by_n_papers_desc(subtree_index):
    out = get_taxon_subtree_dossier("Marrus")
    paper_counts = [s["n_papers"]
                    for s in out["species_with_corpus_coverage"]]
    assert paper_counts == sorted(paper_counts, reverse=True)


def test_subtree_species_carries_authorship_when_present(subtree_index):
    out = get_taxon_subtree_dossier("Marrus")
    cla = next(s for s in out["species_with_corpus_coverage"]
               if s["taxon_id"] == "t:claudanielis")
    assert cla["authorship"] == "Dunn, 2005"
    assert cla["n_papers"] == 2
    assert cla["n_mentions"] == 7


def test_subtree_aggregate_papers_deduplicated(subtree_index):
    """Paper aaa mentions BOTH species — appears once in aggregate
    with n_species_covered=2."""
    out = get_taxon_subtree_dossier("Marrus")
    paper_a = next(p for p in out["aggregate_papers"]
                   if p["hash"] == "aaaaaaaaaaaa")
    assert paper_a["n_species_covered"] == 2
    # claudanielis (5 mentions) + orthocannoides (1) in paper aaa.
    assert paper_a["total_mentions"] == 6


def test_subtree_aggregate_excludes_papers_outside_subtree(subtree_index):
    """Paper ccc mentions Erenna richardi, NOT a Marrus species.
    Must not appear in Marrus's aggregate."""
    out = get_taxon_subtree_dossier("Marrus")
    hashes = {p["hash"] for p in out["aggregate_papers"]}
    assert "cccccccccccc" not in hashes


def test_subtree_empty_clade_returns_empty_lists(subtree_index):
    """Clade with no descendants should still return a clean shape."""
    out = get_taxon_subtree_dossier("EmptyClade")
    assert out["n_species_in_subtree"] == 0
    assert out["species_with_corpus_coverage"] == []
    assert out["aggregate_papers"] == []


def test_subtree_not_found(subtree_index):
    out = get_taxon_subtree_dossier("not_a_clade")
    assert out["not_found"] is True


def test_subtree_max_species_respected(subtree_index):
    out = get_taxon_subtree_dossier("Marrus", max_species=1)
    assert len(out["species_with_corpus_coverage"]) == 1


def test_subtree_max_aggregate_papers_respected(subtree_index):
    out = get_taxon_subtree_dossier("Marrus", max_aggregate_papers=1)
    assert len(out["aggregate_papers"]) == 1
