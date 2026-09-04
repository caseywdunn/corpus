"""Tests for the observation-backed missing-reference audit (#155)."""
from __future__ import annotations

import importlib.util
import json
import sqlite3
import sys
from pathlib import Path

import pytest

from bib import authority


_SCRIPT = Path(__file__).parent.parent / "tools" / "qc" / "reference_reconciliation.py"
_SPEC = importlib.util.spec_from_file_location("reference_reconciliation_qc", _SCRIPT)
qc = importlib.util.module_from_spec(_SPEC)
sys.modules["reference_reconciliation_qc"] = qc
_SPEC.loader.exec_module(qc)


def _write_paper(
    output_dir: Path,
    corpus_hash: str,
    title: str,
    year: int,
    surname: str,
    references: list[dict],
) -> None:
    doc_dir = output_dir / "documents" / corpus_hash
    doc_dir.mkdir(parents=True, exist_ok=True)
    (doc_dir / "metadata.json").write_text(json.dumps({
        "title": title,
        "year": year,
        "authors": [{"surname": surname, "forename": "A"}],
    }))
    (doc_dir / "references.json").write_text(json.dumps({
        "references": references,
        "total_references": len(references),
    }))


def _build_fixture(tmp_path: Path) -> Path:
    output_dir = tmp_path / "output"
    present_title = "A sufficiently distinctive monograph already in corpus"
    _write_paper(
        output_dir, "aaa", present_title, 1901, "Canonical", [],
    )
    # A second in-corpus row with the same title/year makes automatic
    # reconciliation correctly stay conservative; the QC report must expose
    # both candidates for review rather than silently choosing one.
    _write_paper(
        output_dir, "aad", present_title, 1901, "Other", [],
    )
    escaped = {
        "raw": "Wrong X. 1901. A sufficiently distinctive monograph...",
        "title": present_title,
        "year": 1901,
        "authors": ["X Wrong"],
        "xml_id": "b0",
    }
    genuine = {
        "raw": "Missing M. 1888. A genuinely absent foundational study.",
        "title": "A genuinely absent foundational study",
        "year": 1888,
        "authors": ["M Missing"],
        "xml_id": "b1",
    }
    _write_paper(
        output_dir, "bbb", "First citing paper", 2001, "Beta",
        [escaped, genuine],
    )
    _write_paper(
        output_dir, "ccc", "Second citing paper", 2002, "Gamma",
        [{**escaped, "xml_id": "b7"}],
    )

    db_path = output_dir / "biblio_authority.sqlite"
    with sqlite3.connect(db_path) as conn:
        authority.create_schema(conn)
        authority.phase1_corpus_papers(conn, output_dir)
        authority.phase2_references(conn, output_dir)
    return db_path


def test_report_explains_ranked_work_from_raw_observations(tmp_path: Path):
    report = qc.build_report(
        _build_fixture(tmp_path), min_citations=1, limit=None,
    )

    assert report["population"] == {
        "current_reference_observations": 3,
        "mapped_current_observations": 3,
        "unmapped_current_observations": 0,
        "unmapped_citing_documents": 0,
        "compatibility_citation_edges": 3,
        "observations_collapsed_by_compatibility_key": 0,
        "missing_works_at_threshold": 2,
        "missing_work_citation_edges": 3,
        "possible_in_corpus_duplicates": 1,
        "possible_duplicate_citation_edges": 2,
        "exact_identity_candidates": 0,
        "exact_identity_citation_edges": 0,
    }
    flagged = next(
        item for item in report["top_missing"]
        if item["possible_in_corpus_matches"]
    )
    assert flagged["cited_by_count"] == 2
    assert flagged["observation_count"] == 2
    assert flagged["citing_paper_count"] == 2
    assert flagged["review_reason"] == "same_normalized_title_and_year"
    assert {
        candidate["corpus_hash"]
        for candidate in flagged["possible_in_corpus_matches"]
    } == {"aaa", "aad"}
    assert not any(
        candidate["author_set_match"]
        for candidate in flagged["possible_in_corpus_matches"]
    )
    assert {d["method"] for d in flagged["mapping_decisions"]} == {
        "new_work", "alias_exact",
    }
    assert all(d["producer_version"] for d in flagged["mapping_decisions"])
    assert flagged["observation_examples"][0]["raw_citation"]
    assert flagged["observation_examples"][0]["authors"] == ["X Wrong"]
    assert len(flagged["observation_examples"]) == 1


def test_report_limit_does_not_change_population_counts(tmp_path: Path):
    db_path = _build_fixture(tmp_path)
    full = qc.build_report(db_path, min_citations=1, limit=None)
    summary_only = qc.build_report(db_path, min_citations=1, limit=0)

    assert summary_only["top_missing"] == []
    assert summary_only["population"] == full["population"]


def test_pre_observation_database_fails_with_migration_instruction(
    tmp_path: Path,
):
    db_path = tmp_path / "old.sqlite"
    with sqlite3.connect(db_path) as conn:
        authority.create_schema(conn)
        for table in (
            "observation_work",
            "reference_current_sets",
            "reference_observation_memberships",
            "reference_observation_sets",
            "reference_observations",
        ):
            conn.execute(f"DROP TABLE {table}")

    with pytest.raises(ValueError, match="re-run bib.authority"):
        qc.build_report(db_path)


def test_report_rejects_invalid_bounds(tmp_path: Path):
    db_path = _build_fixture(tmp_path)
    with pytest.raises(ValueError, match="min_citations"):
        qc.build_report(db_path, min_citations=0)
    with pytest.raises(ValueError, match="limit"):
        qc.build_report(db_path, limit=-1)
    with pytest.raises(ValueError, match="examples"):
        qc.build_report(db_path, examples=-1)
