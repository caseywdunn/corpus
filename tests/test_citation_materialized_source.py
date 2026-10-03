"""Captured TEI → real authority phases → served citation evidence (#309/#317)."""

from __future__ import annotations

import hashlib
import json
import sqlite3
from types import SimpleNamespace

import pytest

from bib.authority import create_schema, phase1_corpus_papers, phase2_references
from bib.documents import find_work
from bib.parser import bib_entry_to_metadata
from mcpsrv import app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import (
    get_citation_graph,
    get_excerpts_citing,
    get_intext_citations,
)
from pipeline.grobid_client import parse_tei_intext_citations, parse_tei_references
from tests.test_citation_source_spans import FIXTURES, _capture


@pytest.fixture
def source_materialization(tmp_path, monkeypatch):
    """Only citing headers are curated; every reference/edge comes from TEI."""
    headers = json.loads((FIXTURES / "citing_headers.json").read_text())
    cases = {}
    papers = {}
    for name, entry in headers["entries"].items():
        evidence = json.loads((FIXTURES / f"{name}.json").read_text())
        xml = (FIXTURES / f"{name}.tei.xml").read_bytes()
        assert hashlib.sha256(xml).hexdigest() == evidence["fragment_sha256"]
        sha = evidence["source_sha256"][:12]
        folder = tmp_path / "documents" / sha
        folder.mkdir(parents=True)
        (folder / "grobid.tei.xml").write_bytes(xml)
        metadata = bib_entry_to_metadata(entry, entry["file"])
        (folder / "metadata.json").write_text(json.dumps(metadata))
        references = parse_tei_references(xml.decode())
        (folder / "references.json").write_text(json.dumps({"references": references}))
        # Start with the actual parser's explicit TEI-only reading. The source
        # replay below replaces this artifact with coordinate-supported text.
        before = parse_tei_intext_citations(xml.decode())
        (folder / "intext_citations.json").write_text(json.dumps(before))
        cases[name] = {
            "hash": sha,
            "folder": folder,
            "references": references,
            "tei_only": before,
            "xml": xml,
        }
        papers[sha] = {"hash_dir": str(folder), "title": metadata["title"]}

    db = tmp_path / "biblio.sqlite"
    conn = sqlite3.connect(db)
    create_schema(conn)
    assert phase1_corpus_papers(conn, tmp_path) == 3
    assert phase2_references(conn, tmp_path) == (63, 63)
    raw_before = list(
        conn.execute("SELECT * FROM reference_observations ORDER BY observation_id")
    )
    for name, case in cases.items():
        observed = conn.execute(
            """SELECT grobid_xml_id,raw_citation,title,year,journal,doi,authors_json
               FROM reference_observations WHERE citing_corpus_hash=? ORDER BY ordinal""",
            (case["hash"],),
        ).fetchall()
        assert [(*row[:-1], json.loads(row[-1])) for row in observed] == [
            (
                r["xml_id"],
                r["raw"],
                r["title"],
                r["year"],
                r["journal"],
                r["doi"],
                r["authors"],
            )
            for r in case["references"]
        ]
        # _capture verifies the original fragment hash and substitutes only the
        # PDF character reader with pinned source geometry. The parser is real.
        data = _capture(name, monkeypatch)
        (case["folder"] / "intext_citations.json").write_text(json.dumps(data))
        case["data"] = data
        case["work_id"] = find_work(conn, case["hash"])
        case["targets"] = dict(
            conn.execute(
                "SELECT grobid_xml_id,cited_work_id FROM citations WHERE citing_corpus_hash=?",
                (case["hash"],),
            )
        )
        assert len(case["targets"]) == len(case["references"])
    assert (
        list(
            conn.execute("SELECT * FROM reference_observations ORDER BY observation_id")
        )
        == raw_before
    )
    conn.commit()
    biblio = BiblioAuthority(db)
    monkeypatch.setattr(app, "_INDEX", SimpleNamespace(biblio_db=biblio, papers=papers))
    try:
        yield SimpleNamespace(
            root=tmp_path, conn=conn, cases=cases, raw_before=raw_before
        )
    finally:
        biblio.conn.close()
        conn.close()


def _assert_graph(case):
    result = get_citation_graph(paper_hash=case["hash"], direction="cited_by")
    expected = set(case["targets"].values())
    assert result["root"]["work_id"] == case["work_id"]
    assert {row["work_id"] for row in result["cited_by"]} == expected
    assert (
        result["edges_available"]
        == result["edges_returned"]
        == {"cited_by": len(expected)}
    )
    assert not result["truncated"]
    return result


def _assert_served_marker(case, marker, *, author, year_text):
    marker_index = case["data"]["citations"].index(marker)
    paragraph = case["data"]["paragraphs"][marker["para_index"]]
    intext = get_intext_citations(case["hash"], offset=marker_index, limit=1)
    assert intext["total_citations"] == len(case["data"]["citations"])
    assert intext["citations"] == [dict(marker, para_index=0)]
    assert intext["paragraphs"] == [paragraph]
    assert paragraph[slice(*marker["author_span"])] == author
    assert paragraph[slice(*marker["year_span"])] == year_text
    assert marker["text_source"] == "pdf_coordinates"
    target = case["targets"][marker["target_xml_id"].lstrip("#")]
    result = get_excerpts_citing(target)
    excerpt = [
        row
        for row in result["excerpts"]
        if row["citing_paper_hash"] == case["hash"]
        and row["citation_index"] == marker_index
    ]
    assert len(excerpt) == 1
    assert excerpt[0]["paragraph"] == paragraph
    for field in (
        "surface",
        "target_xml_id",
        "citation_year",
        "author_span",
        "year_span",
        "validation_status",
        "text_source",
    ):
        assert excerpt[0][field] == marker[field]
    inbound = get_citation_graph(work_id=target, direction="citing")
    reference = next(
        ref
        for ref in case["references"]
        if "#" + ref["xml_id"] == marker["target_xml_id"]
    )
    assert inbound["root"]["title"] == reference["title"]
    assert inbound["root"]["year"] == reference["year"]
    assert case["work_id"] in {row["work_id"] for row in inbound["citing"]}
    return excerpt[0]


def test_oderberg_actual_reference_edges_exclude_wrong_synopsis_excerpt(
    source_materialization,
):
    case = source_materialization.cases["Oderberg2020"]
    graph = _assert_graph(case)
    markers = case["data"]["citations"]
    correct = next(row for row in markers if row.get("citation_year") == "1960")
    excerpt = _assert_served_marker(
        case, correct, author="Totton and Mackie", year_text="1960"
    )
    assert correct["target_xml_id"] == "#b31"
    assert correct["tei_observations"][0]["target_xml_id"] == "#b30"
    assert "cormidium contains" in excerpt["paragraph"]
    assert "(Totton and Mackie 1960, Bardi and Marques 2007)." in excerpt["paragraph"]
    bardi = next(row for row in markers if row.get("citation_year") == "2007")
    _assert_served_marker(case, bardi, author="Bardi and Marques", year_text="2007")
    assert bardi["target_xml_id"] == "#b0"
    wrong_excerpt = get_excerpts_citing(case["targets"]["b30"])
    assert wrong_excerpt["excerpts_available"] == wrong_excerpt["n_excerpts"] == 0
    assert all(row["target_xml_id"] != "#b30" for row in markers)
    # The bibliography genuinely lists both works. Its graph must retain the
    # 1965 reference even though this captured paragraph does not cite it.
    years = {row["work_id"]: row["year"] for row in graph["cited_by"]}
    assert years[case["targets"]["b30"]] == 1965
    assert years[case["targets"]["b31"]] == 1960


def test_fraser_source_paragraph_and_both_actual_targets_are_served(
    source_materialization,
):
    case = source_materialization.cases["Pugh1974"]
    _assert_graph(case)
    assert "Fraser ( , 1967) )" in case["tei_only"]["paragraphs"][0]
    for year, xml_id in [("1961", "#b27"), ("1967", "#b28")]:
        marker = next(
            row for row in case["data"]["citations"] if row.get("citation_year") == year
        )
        assert marker["target_xml_id"] == xml_id
        excerpt = _assert_served_marker(case, marker, author="Fraser", year_text=year)
        assert (
            "The findings of Fraser (1961, 1967) showed that the species"
            in excerpt["paragraph"]
        )
        assert "Fraser ( , 1967) )" not in excerpt["paragraph"]
    # Preserve the unrelated real upstream bibliography error as unresolved.
    unresolved = next(
        row for row in case["data"]["citations"] if row.get("citation_year") == "1955"
    )
    assert unresolved["target_xml_id"] is None
    assert unresolved["validation_status"] == "unresolved_author_year"
    assert get_excerpts_citing(case["targets"]["b39"])["excerpts_available"] == 0


def test_pugh_source_range_spans_and_ambiguous_target_survive_full_path(
    source_materialization,
):
    case = source_materialization.cases["Mapstone2009"]
    _assert_graph(case)
    paragraph = case["data"]["paragraphs"][0]
    printed = "Pugh (1992a-c, 1995, 1999a and b, 2001,2003,2005, 2006a and b)"
    assert printed in paragraph
    assert "Pugh (1992aPugh" not in paragraph
    start = paragraph.index(printed)
    markers = [
        row
        for row in case["data"]["citations"]
        if row.get("author_span") == [start, start + 4]
    ]
    expected = [
        ("1992a", "#b515", "1992a-c"),
        ("1992b", "#b516", "1992a-c"),
        ("1992c", "#b517", "1992a-c"),
        ("1995", "#b518", "1995"),
        ("1999a", "#b520", "1999a and b"),
        ("1999b", None, "1999a and b"),
        ("2001", "#b522", "2001"),
        ("2003", "#b524", "2003"),
        ("2005", "#b525", "2005"),
        ("2006a", "#b526", "2006a and b"),
        ("2006b", "#b527", "2006a and b"),
    ]
    assert [(row["citation_year"], row["target_xml_id"]) for row in markers] == [
        pair[:2] for pair in expected
    ]
    for marker, (_, target, year_text) in zip(markers, expected):
        if target:
            excerpt = _assert_served_marker(
                case, marker, author="Pugh", year_text=year_text
            )
            assert printed in excerpt["paragraph"]
    ambiguous = markers[5]
    assert paragraph[slice(*ambiguous["year_span"])] == "1999a and b"
    assert ambiguous["validation_status"] == "ambiguous_author_year"
    assert ambiguous["candidate_target_xml_ids"] == ["#b20", "#b521"]
    marker_index = case["data"]["citations"].index(ambiguous)
    served = get_intext_citations(case["hash"], offset=marker_index, limit=1)
    assert served["citations"] == [dict(ambiguous, para_index=0)]
    assert served["paragraphs"] == [paragraph]
    for candidate in ambiguous["candidate_target_xml_ids"]:
        excerpts = get_excerpts_citing(case["targets"][candidate.lstrip("#")])
        assert not any(
            row["citing_paper_hash"] == case["hash"]
            and row["citation_index"] == marker_index
            for row in excerpts["excerpts"]
        )


def test_source_regeneration_and_unchanged_authority_refresh_preserve_raw_history(
    source_materialization,
    monkeypatch,
):
    replay = source_materialization
    before = list(replay.conn.iterdump())
    changes = replay.conn.total_changes
    for name, case in replay.cases.items():
        path = case["folder"] / "intext_citations.json"
        old = path.read_bytes()
        regenerated = _capture(name, monkeypatch)
        path.write_text(json.dumps(regenerated))
        assert path.read_bytes() == old
        assert (case["folder"] / "grobid.tei.xml").read_bytes() == case["xml"]
    assert phase1_corpus_papers(replay.conn, replay.root) == 0
    assert phase2_references(replay.conn, replay.root) == (0, 0)
    assert replay.conn.total_changes == changes
    assert list(replay.conn.iterdump()) == before
    assert (
        list(
            replay.conn.execute(
                "SELECT * FROM reference_observations ORDER BY observation_id"
            )
        )
        == replay.raw_before
    )
