"""Marker pagination counts and transports bounded evidence without silent gaps (#325)."""
from __future__ import annotations

import json
import sqlite3
import types

import pytest

from mcpsrv import app
from mcpsrv.tools import bibliography
from mcpsrv.tools.bibliography import get_excerpts_citing


@pytest.fixture
def served(tmp_path, monkeypatch):
    conn = sqlite3.connect(":memory:")
    conn.row_factory = sqlite3.Row
    conn.execute("CREATE TABLE citations (citing_corpus_hash, grobid_xml_id, cited_work_id)")
    papers = {}
    # Reverse insertion order, duplicate authority edges and multiple targets
    # for the same work must not affect source-marker order or duplicate rows.
    for paper_hash, n in [("bbbbbbbbbbbb", 3), ("aaaaaaaaaaaa", 4)]:
        directory = tmp_path / paper_hash
        directory.mkdir()
        citations = [{"target_xml_id": "#b1" if i % 2 else "#b0",
                      "surface": f"Author {i}", "section": "Discussion", "para_index": 0}
                     for i in range(n)]
        citations.insert(1, {"target_xml_id": None, "surface": "Unresolved", "para_index": 0})
        citations.append({"target_xml_id": "#other", "surface": "Other", "para_index": 1})
        data = {"paragraphs": ["Evidence with −2 °C and several citations.", "Other work"],
                "citations": citations}
        (directory / "intext_citations.json").write_text(json.dumps(data))
        papers[paper_hash] = {"hash_dir": str(directory), "title": "Citing paper"}
        for target in ("b1", "b0", "b0"):
            conn.execute("INSERT INTO citations VALUES (?, ?, 'work')", (paper_hash, target))
    conn.execute("INSERT INTO citations VALUES ('absent', 'b0', 'work')")
    conn.execute("INSERT INTO citations VALUES ('aaaaaaaaaaaa', NULL, 'work')")
    idx = types.SimpleNamespace(biblio_db=types.SimpleNamespace(conn=conn), papers=papers)
    monkeypatch.setattr(app, "_INDEX", idx)
    yield idx
    conn.close()


def test_counts_are_markers_not_unique_paragraphs(served):
    out = get_excerpts_citing("work", limit=2)
    assert out["excerpts_available"] == 7
    assert out["n_excerpts"] == out["excerpts_returned"] == 2
    assert out["counting_unit"] == "citation_markers"
    assert out["truncated"] is True
    assert out["truncated_reason"] == ["limit"]
    assert out["next_offset"] == 2
    assert out["excerpts"][0]["paragraph"] == out["excerpts"][1]["paragraph"]
    assert [r["citation_index"] for r in out["excerpts"]] == [0, 2]


@pytest.mark.parametrize("limit", [1, 2, 5, 50])
def test_pagination_recovers_every_marker_exactly_once(served, limit):
    offset = 0
    collected = []
    while offset is not None:
        out = get_excerpts_citing("work", limit=limit, offset=offset)
        assert out["offset"] == offset
        assert out["excerpts_available"] == 7
        assert out == get_excerpts_citing("work", limit=limit, offset=offset)
        collected.extend((r["citing_paper_hash"], r["citation_index"]) for r in out["excerpts"])
        assert out["next_offset"] is None or out["next_offset"] > offset
        offset = out["next_offset"]
    assert collected == [("aaaaaaaaaaaa", i) for i in (0, 2, 3, 4)] + [
        ("bbbbbbbbbbbb", i) for i in (0, 2, 3)]
    assert out["truncated"] is False
    assert out["truncated_reason"] == []


def test_byte_limited_pages_recover_all_markers(served, monkeypatch):
    monkeypatch.setattr(bibliography, "EXCERPTS_MAX_BYTES", 1300)
    offset = 0
    collected = []
    while offset is not None:
        out = get_excerpts_citing("work", limit=500, offset=offset)
        assert out["response_bytes"] == len(json.dumps(out).encode("utf-8"))
        assert out["response_bytes"] <= 1300
        collected.extend((r["citing_paper_hash"], r["citation_index"]) for r in out["excerpts"])
        assert all("truncated_fields" not in r for r in out["excerpts"])
        if out["next_offset"] is not None:
            assert out["next_offset"] > offset
            assert out["truncated_reason"] == ["response_bytes"]
        offset = out["next_offset"]
    assert len(collected) == len(set(collected)) == 7


def test_one_oversized_marker_preserves_identity_and_admits_preview(served, monkeypatch):
    monkeypatch.setattr(bibliography, "EXCERPTS_MAX_BYTES", 1600)
    paragraph = "𝛼−2 μm\n" * 500
    original_load = bibliography._load_json

    def load(path, default=None):
        data = original_load(path, default=default)
        return {**data, "paragraphs": [paragraph, "Other work"]}

    monkeypatch.setattr(bibliography, "_load_json", load)
    out = get_excerpts_citing("work", limit=50)
    assert out["n_excerpts"] == 1
    assert out["next_offset"] == 1
    assert out["response_bytes"] == len(json.dumps(out).encode("utf-8")) <= 1600
    assert out["truncated_reason"] == ["response_bytes"]
    row = out["excerpts"][0]
    assert row["truncated_fields"]["paragraph"] == len(paragraph)
    assert len(row["paragraph"]) < len(paragraph)
    assert paragraph.startswith(row["paragraph"])
    assert (row["citing_paper_hash"], row["citation_index"], row["para_index"], row["target_xml_id"]) == (
        "aaaaaaaaaaaa", 0, 0, "#b0")


@pytest.mark.parametrize("kwargs", [{"offset": -1}, {"limit": -1}])
def test_invalid_pagination_is_structured(served, kwargs):
    assert get_excerpts_citing("work", **kwargs)["code"] == "invalid_argument"


def test_zero_limit_is_explicit_counts_only(served):
    out = get_excerpts_citing("work", limit=0)
    assert out["n_excerpts"] == 0
    assert out["excerpts_available"] == 7
    assert out["next_offset"] == 0
    assert out["truncated_reason"] == ["limit"]


@pytest.mark.parametrize("work_id, offset, available", [("work", 7, 7), ("work", 99, 7), ("missing", 0, 0)])
def test_empty_or_past_end_page_is_complete(served, work_id, offset, available):
    out = get_excerpts_citing(work_id, offset=offset)
    assert out["n_excerpts"] == 0
    assert out["excerpts_available"] == available
    assert out["next_offset"] is None
    assert out["truncated"] is False
