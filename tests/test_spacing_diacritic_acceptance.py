"""Named spacing-diacritic source and explicit author adjudication (#316)."""

import copy
import hashlib
import json
import os
from pathlib import Path

import pytest

from bib.authority import phase2_references
from bib.documents import find_work
from pipeline.text_encoding import repair_region, recover_text_encoding
from tests.test_bibliographic_integrity import build

FIXTURE = Path(__file__).parent / "fixtures/spacing_diacritics"
CASE = json.loads((FIXTURE / "nino_source.json").read_text())


def test_source_nino_nina_preserve_printed_accents_and_word_boundaries():
    assert (
        hashlib.sha256((FIXTURE / "nino_source.png").read_bytes()).hexdigest()
        == CASE["crop_sha256"]
    )
    text, repairs, unresolved = repair_region(CASE["observed"], CASE["native_glyphs"])
    assert text == "El Niño/La Niña" == CASE["expected"]
    assert not unresolved
    assert len(repairs) == 2
    for repair in repairs:
        start, end = repair["charspan"]
        assert CASE["observed"][start:end] == repair["original"] == "n ˜"
        assert repair["replacement"] == "ñ"
        assert repair["evidence"] == "overlapping_native_accent_and_base_glyph"
    assert repair_region(text, CASE["native_glyphs"]) == (text, [], [])


def test_nino_nearby_but_nonoverlapping_tildes_do_not_authorize_repair():
    # Moving the actual source accents away from their base letters removes
    # the evidence; the familiar phrase alone must never authorize replacement.
    glyphs = copy.deepcopy(CASE["native_glyphs"])
    for glyph in glyphs:
        if glyph["c"] == "˜":
            glyph["bbox"][0] += 1000
            glyph["bbox"][2] += 1000
    assert repair_region(CASE["observed"], glyphs) == (CASE["observed"], [], [])


def test_optional_original_nino_pdf_replay_retains_observation_and_provenance():
    library = os.environ.get("CORPUS_LIBRARY_DIR")
    if not library:
        pytest.skip("Set CORPUS_LIBRARY_DIR for original Niño/Niña PDF replay")
    from docling_core.types.doc import (
        DoclingDocument,
        DocItemLabel,
        Size,
        BoundingBox,
        ProvenanceItem,
    )

    source = Path(library) / CASE["source"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == CASE["source_sha256"]
    doc = DoclingDocument(name="Niño/Niña source replay")
    width, height = CASE["page_size"]
    doc.add_page(page_no=1, size=Size(width=width, height=height))
    left, top, right, bottom = CASE["bbox"]
    item = doc.add_text(
        label=DocItemLabel.TEXT,
        text=CASE["observed"],
        orig=CASE["observed"],
        prov=ProvenanceItem(
            page_no=1,
            bbox=BoundingBox(l=left, t=top, r=right, b=bottom, coord_origin="TOPLEFT"),
            charspan=(0, len(CASE["observed"])),
        ),
    )
    report = recover_text_encoding(doc, source)
    assert item.text == CASE["expected"] and item.orig == CASE["observed"]
    assert len(report["repairs"]) == 2 and report["unresolved"] == []
    assert all(
        repair["page"] == 1 and repair["item_ref"] == item.self_ref
        for repair in report["repairs"]
    )
    assert item.meta.corpus__text_encoding[0]["repairs"]
    before = copy.deepcopy(doc.export_to_dict())
    assert recover_text_encoding(doc, source)["repairs"] == []
    assert doc.export_to_dict() == before


def test_adjudicated_complete_author_rebuilds_edge_and_retains_raw_history(tmp_path):
    # The curated entry's fields are from Kolliker1853 in the pinned source
    # library. The two parsed-author states are a controlled upstream correction
    # fixture, not a claim that we repaired an actual Grobid record automatically.
    title = "Die Schwimmpolypen oder Siphonophoren von Messina"
    target_entry = {
        "_key": "Kolliker1853",
        "author": "Kölliker, A.",
        "year": "1853",
        "title": title,
        "journal": "96 pp. Wilhelm Engelmann, Leipzig",
        "doi": "10.5962/bhl.title.12447",
    }
    conn, _ = build(
        tmp_path,
        {
            "article": target_entry,
            "citing": {
                "author": "Observer, B.",
                "year": "1900",
                "title": "A citing study",
            },
        },
    )
    target = find_work(conn, "article")
    source = find_work(conn, "citing")
    path = tmp_path / "documents/citing/references.json"
    observed = {
        "xml_id": "b0",
        "authors": ["A ¨lliker"],
        "year": 1853,
        "title": title,
        "raw": f"Kölliker, A. 1853. {title}. Leipzig.",
    }
    path.write_text(json.dumps({"references": [observed]}))
    phase2_references(conn, tmp_path)
    old_observation = tuple(
        conn.execute("SELECT * FROM reference_observations").fetchone()
    )
    old_id = conn.execute(
        "SELECT observation_id FROM reference_observations"
    ).fetchone()[0]
    old_work, method = conn.execute(
        "SELECT work_id,match_method FROM observation_work"
    ).fetchone()
    assert method == "unresolved_author" and old_work != target
    assert not conn.execute(
        "SELECT 1 FROM work_aliases WHERE work_id=?", (old_work,)
    ).fetchone()
    assert (
        conn.execute(
            "SELECT disposition FROM reference_observation_quality"
        ).fetchone()[0]
        == "review_needed"
    )
    assert (
        conn.execute(
            "SELECT cited_work_id FROM citations WHERE citing_work_id=?", (source,)
        ).fetchone()[0]
        == old_work
    )
    changes = conn.total_changes
    assert phase2_references(conn, tmp_path) == (0, 0)
    assert conn.total_changes == changes

    # Model a regenerated reference artifact after explicit source adjudication.
    # No source observation, served bundle, or guessed surname alias is edited.
    corrected = dict(observed, authors=["A Kölliker"])
    path.write_text(json.dumps({"references": [corrected]}))
    phase2_references(conn, tmp_path)
    assert (
        tuple(
            conn.execute(
                "SELECT * FROM reference_observations WHERE observation_id=?", (old_id,)
            ).fetchone()
        )
        == old_observation
    )
    rows = conn.execute(
        "SELECT authors_json,raw_citation FROM reference_observations"
    ).fetchall()
    assert {tuple(json.loads(authors)) for authors, _ in rows} == {
        ("A ¨lliker",),
        ("A Kölliker",),
    }
    assert all(raw == observed["raw"] for _, raw in rows)
    assert conn.execute("SELECT work_id FROM observation_work").fetchall() == [
        (target,)
    ]
    assert conn.execute(
        "SELECT citing_work_id,cited_work_id FROM citations"
    ).fetchall() == [(source, target)]
    assert not conn.execute(
        "SELECT 1 FROM observation_work WHERE observation_id=?", (old_id,)
    ).fetchone()
    assert not conn.execute(
        "SELECT 1 FROM work_aliases WHERE work_id=?", (old_work,)
    ).fetchone()
    assert (
        conn.execute(
            "SELECT disposition FROM reference_observation_quality"
        ).fetchone()[0]
        == "usable"
    )
    assert conn.execute(
        "SELECT surname,forename FROM work_authors WHERE work_id=? ORDER BY position",
        (target,),
    ).fetchall() == [("Kölliker", "A.")]
    changes = conn.total_changes
    assert phase2_references(conn, tmp_path) == (0, 0)
    assert conn.total_changes == changes
