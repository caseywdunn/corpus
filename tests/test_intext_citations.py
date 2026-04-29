"""Tests for parse_tei_intext_citations (issue #7).

Hand-crafted TEI fragments cover the cases that matter:
* a resolved citation (``target="#bN"`` set by Grobid)
* an unresolved citation (no target — keep surface text)
* multiple citations in one paragraph (paragraph dedup)
* citations across separate sections (heading capture)
* refs outside ``<text><body>`` (header / listBibl) ignored

Plus one round-trip against the on-disk siphonophore corpus when
``CORPUS_OUTPUT_DIR`` is set, as a sanity check on real Grobid output.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from grobid_client import parse_tei_intext_citations


def _wrap(body_xml: str) -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<TEI xmlns="http://www.tei-c.org/ns/1.0">
  <teiHeader/>
  <text>
    <body>{body_xml}</body>
    <back>
      <div>
        <listBibl>
          <biblStruct xml:id="b0"><analytic><title>Some paper</title></analytic></biblStruct>
        </listBibl>
      </div>
    </back>
  </text>
</TEI>"""


def test_resolved_citation():
    tei = _wrap("""
      <div><head>Introduction</head>
        <p>Earlier work <ref type="bibr" target="#b0">Smith 2010</ref> showed X.</p>
      </div>
    """)
    out = parse_tei_intext_citations(tei)
    assert len(out["citations"]) == 1
    c = out["citations"][0]
    assert c["target_xml_id"] == "#b0"
    assert c["surface"] == "Smith 2010"
    assert c["section"] == "Introduction"
    assert out["paragraphs"][c["para_index"]].startswith("Earlier work")


def test_unresolved_citation_preserves_surface():
    """When Grobid emits <ref> with no target= attribute, we still keep
    the surface text so a later fuzzy-resolution pass can recover it."""
    tei = _wrap("""
      <div><head>Discussion</head>
        <p>Older lit <ref type="bibr">Cervigon, 1958</ref> noted Y.</p>
      </div>
    """)
    out = parse_tei_intext_citations(tei)
    assert len(out["citations"]) == 1
    c = out["citations"][0]
    assert c["target_xml_id"] is None
    assert c["surface"] == "Cervigon, 1958"
    assert c["section"] == "Discussion"


def test_paragraph_dedup_for_multiple_cites_in_one_p():
    tei = _wrap("""
      <div><head>Methods</head>
        <p>We follow <ref type="bibr" target="#b0">A 1990</ref> and
           <ref type="bibr" target="#b1">B 1991</ref>.</p>
      </div>
    """)
    out = parse_tei_intext_citations(tei)
    assert len(out["citations"]) == 2
    assert len(out["paragraphs"]) == 1, "same paragraph should be deduped"
    assert out["citations"][0]["para_index"] == out["citations"][1]["para_index"] == 0


def test_separate_sections_get_separate_section_strings():
    tei = _wrap("""
      <div><head>Introduction</head>
        <p>See <ref type="bibr" target="#b0">A 1990</ref>.</p>
      </div>
      <div><head>Results</head>
        <p>We confirm <ref type="bibr" target="#b1">B 1991</ref>.</p>
      </div>
    """)
    out = parse_tei_intext_citations(tei)
    sections = [c["section"] for c in out["citations"]]
    assert sections == ["Introduction", "Results"]


def test_listBibl_refs_are_ignored():
    """Bibliography-internal cross-references (within <listBibl>) and
    other non-body refs must not show up as in-text citations."""
    tei = _wrap("""
      <div><head>Introduction</head>
        <p>Body cite: <ref type="bibr" target="#b0">A 1990</ref>.</p>
      </div>
    """)
    # Inject a bogus ref inside listBibl — would be picked up by a naive
    # //ref query but our parser restricts to text/body.
    tei = tei.replace(
        '<biblStruct xml:id="b0">',
        '<biblStruct xml:id="b0"><note><ref type="bibr" target="#bX">SHOULD NOT MATCH</ref></note>',
    )
    out = parse_tei_intext_citations(tei)
    assert len(out["citations"]) == 1
    assert out["citations"][0]["surface"] == "A 1990"


def test_empty_surface_skipped():
    tei = _wrap("""<div><head>X</head><p>Empty <ref type="bibr" target="#b0"/> ref.</p></div>""")
    out = parse_tei_intext_citations(tei)
    assert out["citations"] == []


def test_smoke_real_corpus():
    """If CORPUS_OUTPUT_DIR points at a real run, parse one TEI and
    check the shape — guards against schema drift in Grobid output."""
    output_dir = os.environ.get("CORPUS_OUTPUT_DIR")
    if not output_dir:
        pytest.skip("CORPUS_OUTPUT_DIR not set")
    docs = Path(output_dir) / "documents"
    if not docs.is_dir():
        pytest.skip(f"No documents/ in {output_dir}")
    teis = sorted(docs.glob("*/grobid.tei.xml"))
    if not teis:
        pytest.skip("No grobid.tei.xml files in corpus")
    out = parse_tei_intext_citations(teis[0].read_text(encoding="utf-8"))
    assert "paragraphs" in out and "citations" in out
    for c in out["citations"]:
        assert "target_xml_id" in c and "surface" in c
        assert isinstance(c["para_index"], int)
        assert 0 <= c["para_index"] < len(out["paragraphs"])
