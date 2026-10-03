"""Source-captured regressions for clipped targets and TEI expansions (#309/#317)."""
from __future__ import annotations

import hashlib
import json
import sqlite3
import types
from pathlib import Path

from lxml import etree

from bib.authority import create_schema
from mcpsrv import app
from mcpsrv.indexes import BiblioAuthority
from mcpsrv.tools.bibliography import get_citation_graph, get_excerpts_citing, get_intext_citations
from pipeline import citation_spans
from pipeline.grobid_client import GrobidClient, parse_tei_intext_citations
from tests.test_format_citation_tool import _insert_work
from tests.test_intext_citations import _wrap

FIXTURES = Path(__file__).parent / "fixtures/citation_spans"


class CapturedSource(citation_spans.PdfCitationSource):
    """Captured PDF characters exercise the actual coordinate-to-text code."""
    def __init__(self, evidence):
        self.evidence = evidence

    def page_chars(self, page):
        return self.evidence["pdf_characters"][str(page)]

    def close(self):
        pass


def _capture(name, monkeypatch):
    evidence = json.loads((FIXTURES / f"{name}.json").read_text())
    xml = (FIXTURES / f"{name}.tei.xml").read_bytes()
    assert hashlib.sha256(xml).hexdigest() == evidence["fragment_sha256"]
    assert evidence["source_sha256"].startswith({
        "Pugh1974": "b5a7af6140ca", "Oderberg2020": "c59691fa845a", "Mapstone2009": "4efbb4af134b",
    }[name])
    monkeypatch.setattr(citation_spans, "PdfCitationSource", lambda _: CapturedSource(evidence))
    return parse_tei_intext_citations(xml.decode(), pdf_path=Path("captured-processed.pdf"))


def test_fraser_source_text_keeps_two_complete_year_mappings(monkeypatch):
    out = _capture("Pugh1974", monkeypatch)
    paragraph = out["paragraphs"][0]
    assert "The findings of Fraser (1961, 1967) showed that the species" in paragraph
    rows = [r for r in out["citations"] if r.get("citation_year") in {"1961", "1967"}]
    assert [r["target_xml_id"] for r in rows] == ["#b27", "#b28"]
    for row in rows:
        assert row["text_source"] == "pdf_coordinates"
        assert paragraph[slice(*row["author_span"])] == "Fraser"
        assert paragraph[slice(*row["year_span"])] == row["citation_year"]
    assert rows[1]["tei_observations"][1]["surface"] == "Fraser ( , 1967) )"


def test_complete_coauthors_and_year_override_clipped_wrong_target(monkeypatch):
    out = _capture("Oderberg2020", monkeypatch)
    paragraph = out["paragraphs"][0]
    assert "(Totton and Mackie 1960, Bardi and Marques 2007)." in paragraph
    totton = next(r for r in out["citations"] if r.get("citation_year") == "1960")
    assert totton["target_xml_id"] == "#b31"
    assert totton["surface"] == "Totton and Mackie 1960"
    assert paragraph[slice(*totton["author_span"])] == "Totton and Mackie"
    assert paragraph[slice(*totton["year_span"])] == "1960"
    assert totton["tei_observations"][0]["target_xml_id"] == "#b30"
    assert all(r["target_xml_id"] != "#b30" for r in out["citations"])
    bardi = next(r for r in out["citations"] if r.get("citation_year") == "2007")
    assert bardi["target_xml_id"] == "#b0"
    assert bardi["surface"] == "Bardi and Marques 2007"


def test_pugh_letter_ranges_and_shared_years_preserve_original_group(monkeypatch):
    out = _capture("Mapstone2009", monkeypatch)
    paragraph = out["paragraphs"][0]
    expected = "Pugh (1992a-c, 1995, 1999a and b, 2001,2003,2005, 2006a and b)"
    assert expected in paragraph
    assert "Pugh (1992aPugh" not in paragraph
    start = paragraph.index(expected)
    rows = [r for r in out["citations"] if r.get("author_span") == [start, start + 4]]
    assert [r["citation_year"] for r in rows] == [
        "1992a", "1992b", "1992c", "1995", "1999a", "1999b", "2001", "2003", "2005", "2006a", "2006b",
    ]
    assert [r["target_xml_id"] for r in rows[:3]] == ["#b515", "#b516", "#b517"]
    assert all(paragraph[slice(*r["year_span"])] == "1992a-c" for r in rows[:3])
    ambiguous = next(r for r in rows if r["citation_year"] == "1999b")
    assert ambiguous["validation_status"] == "ambiguous_author_year"
    assert ambiguous["target_xml_id"] is None
    assert ambiguous["candidate_target_xml_ids"] == ["#b20", "#b521"]
    # Two distinct reference entries with the same author/year remain a
    # reviewable ambiguity; this parser must not guess their work identity.
    assert rows[-1]["target_xml_id"] == "#b527"


def _author_year_tei(body):
    xml = _wrap(body)
    return xml.replace("<analytic><title>Some paper</title></analytic>", """
        <analytic><title>First work</title><author><persName><surname>Fraser</surname></persName></author></analytic>
        <monogr><imprint><date type="published" when="1961"/></imprint></monogr>
    """).replace("</listBibl>", """
        <biblStruct xml:id="b1"><analytic><title>Second work</title>
        <author><persName><surname>Fraser</surname></persName></author></analytic>
        <monogr><imprint><date type="published" when="1967"/></imprint></monogr></biblStruct>
        </listBibl>
    """)


def test_repeated_names_in_prose_are_not_deleted():
    body = '<p>Fraser discussed Fraser. <ref type="bibr" target="#b0">Fraser (1961)</ref> supported Fraser, while <ref type="bibr" target="#b1">Fraser (1967)</ref> disagreed with Fraser.</p>'
    out = parse_tei_intext_citations(_author_year_tei(body))
    assert out["paragraphs"] == ["Fraser discussed Fraser. Fraser (1961) supported Fraser, while Fraser (1967) disagreed with Fraser."]
    assert [r["target_xml_id"] for r in out["citations"]] == ["#b0", "#b1"]


def test_shared_author_years_and_unresolved_group_remain_complete():
    out = parse_tei_intext_citations(_author_year_tei(
        '<p>Compare <ref type="bibr" target="#b0">Fraser (1961, 1967, 1969)</ref>.</p>'))
    assert out["paragraphs"] == ["Compare Fraser (1961, 1967, 1969)."]
    assert [r["target_xml_id"] for r in out["citations"]] == ["#b0", "#b1", None]
    assert out["citations"][-1]["validation_status"] == "unresolved_author_year"


def test_numeric_ranges_are_preserved_without_inventing_individual_works():
    out = parse_tei_intext_citations(_author_year_tei(
        '<p>See <ref type="bibr" target="#b0">Fraser (1961–1967)</ref>.</p>'))
    assert out["paragraphs"] == ["See Fraser (1961–1967)."]
    assert len(out["citations"]) == 1
    assert out["citations"][0]["surface"] == "Fraser (1961–1967)"
    assert out["citations"][0]["validation_status"] == "unverified_tei"


def test_missing_source_does_not_invent_clean_text_or_keep_conflicting_target():
    xml = (FIXTURES / "Oderberg2020.tei.xml").read_text()
    out = parse_tei_intext_citations(xml)
    assert "Totton andMackie 1960" in out["paragraphs"][0]
    assert all(r["text_source"] == "tei" for r in out["citations"])
    assert all(r["target_xml_id"] != "#b30" for r in out["citations"])
    assert any(o["target_xml_id"] == "#b30" for r in out["citations"] for o in r["tei_observations"])


def test_wrong_pdf_coordinates_do_not_replace_tei_with_unrelated_text(monkeypatch):
    evidence = json.loads((FIXTURES / "Pugh1974.json").read_text())
    for characters in evidence["pdf_characters"].values():
        for character in characters:
            if character[0].isalpha():
                character[0] = "Z"
    source = CapturedSource(evidence)
    root = etree.parse(str(FIXTURES / "Pugh1974.tei.xml"))
    refs = root.findall(".//tei:body//tei:ref", citation_spans.NS)
    assert source.group_text(refs) is None


def test_rebuilt_excerpts_exclude_wrong_work_but_keep_valid_1965_evidence(tmp_path, monkeypatch):
    _capture("Oderberg2020", monkeypatch)
    # Legitimate occurrences of both years are independently retained. Their
    # TEI spans are already complete, so target validation must leave them.
    xml = (FIXTURES / "Oderberg2020.tei.xml").read_text()
    xml = xml.replace("</body>", '''<div><p>Valid source: <ref type="bibr" target="#b31">Totton and Mackie (1960)</ref>.</p>
        <p>Valid synopsis: <ref type="bibr" target="#b30">Totton (1965)</ref>.</p></div></body>''')
    out = parse_tei_intext_citations(xml, pdf_path=Path("captured.pdf"))
    document = tmp_path / "document"
    document.mkdir()
    (document / "intext_citations.json").write_text(json.dumps(out))
    db = tmp_path / "biblio.sqlite"
    conn = sqlite3.connect(db)
    create_schema(conn)
    for work, year in [("source", 2020), ("physalia", 1960), ("synopsis", 1965)]:
        _insert_work(conn, work, title=work, year=year, in_corpus=work == "source",
                     corpus_hash="c59691fa845a" if work == "source" else None)
    for xml_id, work in [("b31", "physalia"), ("b30", "synopsis")]:
        conn.execute("INSERT INTO citations (citing_work_id,cited_work_id,citing_corpus_hash,grobid_xml_id,match_method) VALUES ('source',?,'c59691fa845a',?,'source_fixture')", (work, xml_id))
    conn.commit()
    conn.close()
    biblio = BiblioAuthority(db)
    monkeypatch.setattr(app, "_INDEX", types.SimpleNamespace(
        biblio_db=biblio, papers={"c59691fa845a": {"hash_dir": str(document), "title": "Oderberg2020"}}))
    try:
        earlier = get_excerpts_citing("physalia")
        later = get_excerpts_citing("synopsis")
        assert earlier["n_excerpts"] == 2
        assert any("cormidium contains" in row["paragraph"] for row in earlier["excerpts"])
        assert later["n_excerpts"] == 1
        assert later["excerpts"][0]["paragraph"] == "Valid synopsis: Totton (1965)."
        intext = get_intext_citations("c59691fa845a", limit=1)
        assert intext["citations"][0]["target_xml_id"] == "#b31"
        assert "Totton and Mackie 1960" in intext["paragraphs"][0]
        graph = get_citation_graph(work_id="source", direction="cited_by")
        assert {r["work_id"] for r in graph["cited_by"]} == {"physalia", "synopsis"}
    finally:
        biblio.conn.close()


def test_grobid_requests_ref_coordinates_from_source(tmp_path, monkeypatch):
    pdf = tmp_path / "input.pdf"
    pdf.write_bytes(b"test input")
    requests = []

    def post(url, **kwargs):
        requests.append(kwargs["data"])
        return types.SimpleNamespace(text="<TEI/>", raise_for_status=lambda: None)

    monkeypatch.setattr("pipeline.grobid_client.requests.post", post)
    assert GrobidClient().process_fulltext(pdf) == "<TEI/>"
    assert requests[0]["teiCoordinates"] == ["ref"]


def test_prepared_pdf_characters_recover_group_text_and_spans(tmp_path):
    import pymupdf

    pdf = tmp_path / "prepared.pdf"
    with pymupdf.open() as doc:
        page = doc.new_page()
        page.insert_text((72, 72), "Fraser (1961, 1967)")
        words = page.get_text("words")
        first = pymupdf.Rect(words[0][:4]) | pymupdf.Rect(words[1][:4])
        second = pymupdf.Rect(words[2][:4])
        doc.save(pdf)

    def coords(rect):
        return f"1,{rect.x0},{rect.y0},{rect.width},{rect.height}"

    tei = _author_year_tei(
        '<p>The findings of '
        f'<ref type="bibr" target="#b0" coords="{coords(first)}">Fraser (1961</ref>'
        f'<ref type="bibr" target="#b1" coords="{coords(second)}">Fraser (, 1967))</ref>'
        ' are discussed.</p>')
    out = parse_tei_intext_citations(tei, pdf_path=pdf)
    assert out["paragraphs"] == ["The findings of Fraser (1961, 1967) are discussed."]
    assert [r["target_xml_id"] for r in out["citations"]] == ["#b0", "#b1"]
    assert all(r["text_source"] == "pdf_coordinates" for r in out["citations"])


def test_backfill_coordinates_require_matching_pdf_and_tei_receipt(tmp_path):
    from pipeline.intext_citations import _matching_source_pdf

    pdf = tmp_path / "processed.pdf"
    pdf.write_bytes(b"prepared PDF bytes")
    xml = "<TEI/>"
    receipt = tmp_path / "grobid.tei.xml.provenance.json"
    assert _matching_source_pdf(tmp_path, xml) is None
    receipt.write_text(json.dumps({"tei_sha256": hashlib.sha256(xml.encode()).hexdigest(),
                                   "inputs": {"pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest()}}))
    assert _matching_source_pdf(tmp_path, xml) == pdf
    assert _matching_source_pdf(tmp_path, xml + " ") is None
    pdf.write_bytes(b"different PDF")
    assert _matching_source_pdf(tmp_path, xml) is None
