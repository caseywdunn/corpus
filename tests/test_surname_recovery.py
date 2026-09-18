"""Graded source/name controls and immutable reference-edge rematerialization."""

import copy
import json
import os
from pathlib import Path

import pytest

from bib.authority import phase2_references
from bib.documents import find_work
from bib.surname_evidence import supported_reference
from pipeline.surname_recovery import (
    author_catalog,
    propose,
    surname_recovery_producer,
    adjudicate_crop,
    recover_citation_surnames,
)
from tests.test_bibliographic_integrity import build

FIXTURE = Path(__file__).parent / "fixtures/surname_recovery"
ENTRIES = json.loads((FIXTURE / "catalog.json").read_text())
CAPTURE = json.loads((FIXTURE / "sample.json").read_text())
SAMPLE = CAPTURE["rows"]
CATALOG = author_catalog(ENTRIES)


def live_producer():
    producer = surname_recovery_producer(CATALOG)
    required = {"spa", "deu", "fra", "lat"}
    if not required.issubset(producer["models"]):
        pytest.skip("Regional image verification requires Tesseract spa/deu/fra/lat")
    return producer


def test_detector_and_recorded_source_decisions_on_graded_sample():
    import hashlib

    counts = {"source_pdf": [0, 0, 0], "synthetic": [0, 0, 0]}
    for row in SAMPLE:
        assert (
            hashlib.sha256((FIXTURE / row["image"]).read_bytes()).hexdigest()
            == row["image_sha256"]
        )
        candidates = propose(row["observed"], CATALOG)
        assert bool(candidates) == row["candidate_expected"]
        group = "source_pdf" if row["origin"] == "source_pdf" else "synthetic"
        counts[group][0] += 1
        counts[group][1] += bool(candidates)
        if not candidates:
            assert not row["requires_repair"]
            continue
        verdict = row["capture_ocr"]
        if verdict["status"] == "verified":
            assert row["requires_repair"]
            assert verdict["replacement"] in row["printed"]
            counts[group][2] += 1
        elif not row["requires_repair"]:
            assert verdict["status"] == "source_supports_observed_spelling"
        else:
            assert verdict["status"] == "ocr_disagreement"
    # Deliberately selected cases: these are not corpus precision/recall rates.
    assert counts == {"source_pdf": [16, 14, 10], "synthetic": [4, 1, 0]}
    assert CAPTURE["capture_producer"]["catalog_sha256"] == CATALOG["sha256"]


def test_live_regional_ocr_replays_graded_source_and_printed_misspelling():
    producer = live_producer()
    for row in SAMPLE:
        candidates = propose(row["observed"], CATALOG)
        if not candidates:
            continue
        result = adjudicate_crop(
            (FIXTURE / row["image"]).read_bytes(), candidates[0], producer
        )
        assert result["status"] == row["capture_ocr"]["status"], row["id"]
        if result["status"] == "verified":
            assert row["requires_repair"] and result["replacement"] in row["printed"]


def test_catalog_policy_is_order_independent_and_language_selection_is_declared():
    assert CATALOG == author_catalog(list(reversed(ENTRIES)))
    casing = [
        {"author": "ALVARIÑO, A", "year": "1991", "_key": "Upper"},
        {"author": "Alvariño, A", "year": "1991", "_key": "Mixed"},
    ]
    assert author_catalog(casing) == author_catalog(list(reversed(casing)))
    candidate = propose("Alvarifio 1971", CATALOG)[0]["candidates"][0]
    assert candidate["languages"] == ["spa"]
    assert candidate["language_basis"] == "same_author_curated_publications"
    assert any(
        row["bib_key"] == "Alvarino1991" for row in candidate["language_sources"]
    )
    updated = copy.deepcopy(ENTRIES)
    updated[0]["author"] = "Different, A"
    assert author_catalog(updated)["sha256"] != CATALOG["sha256"]


def test_non_candidate_catalog_fields_do_not_invalidate_extraction():
    ascii_entry = {
        "author": "Pacifici, A",
        "year": "1991",
        "title": "Original title",
        "_key": "Pacifici1991",
        "ocrlang": "eng",
    }
    accented = {
        "author": "Alvariño, A",
        "year": "1991",
        "title": "Curated source title",
        "_key": "Alvarino1991",
        "ocrlang": "spa",
    }
    baseline = author_catalog([ascii_entry, accented])
    producer = surname_recovery_producer(baseline)
    for field, value in [("title", "Corrected title"), ("_key", "RevisedKey")]:
        edited = author_catalog([dict(ascii_entry, **{field: value}), accented])
        assert edited == baseline
        assert surname_recovery_producer(edited) == producer
    for field, value in [("author", "Alvarifio, A"), ("year", "1992")]:
        edited = author_catalog([dict(ascii_entry, **{field: value}), accented])
        assert edited["sha256"] != baseline["sha256"]
        assert (
            surname_recovery_producer(edited)["catalog_sha256"]
            != producer["catalog_sha256"]
        )
    # Known ASCII spellings still stop a near-name proposal, including a name
    # that otherwise resembles the accented author's spelling in the same year.
    known = author_catalog([dict(ascii_entry, author="Alvarifio, A"), accented])
    assert propose("Alvarifio 1991", known) == []
    assert propose("Alvarifio 1991", baseline)
    for field, value in [("title", "Corrected candidate title"), ("ocrlang", "fra")]:
        edited = author_catalog([ascii_entry, dict(accented, **{field: value})])
        assert edited["sha256"] != baseline["sha256"]
        assert (
            surname_recovery_producer(edited)["catalog_sha256"]
            != producer["catalog_sha256"]
        )
    candidate = propose("Alvarifio 1991", baseline)[0]["candidates"][0]
    assert candidate["sources"] == [
        {"bib_key": "Alvarino1991", "title": "Curated source title"}
    ]
    assert candidate["languages"] == ["spa"]


def test_near_names_without_unique_citation_context_are_never_automatic():
    assert propose("Alvarifio is a quoted spelling", CATALOG) == []
    assert propose("Alvarifio 1800", CATALOG) == []
    assert propose("Alvariño 1991", CATALOG) == []
    assert (
        propose("“Alvarifio 1991” [sic]", CATALOG)[0]["status"]
        == "quoted_or_sic_context"
    )
    extra = dict(ENTRIES[-1], author="Alvarífio, A", year="1991", _key="AnotherName")
    assert (
        len(
            propose("Alvarifio 1991", author_catalog(ENTRIES + [extra]))[0][
                "candidates"
            ]
        )
        == 2
    )


def source_report():
    producer = CAPTURE["capture_producer"]
    row = SAMPLE[0]
    proposal = propose(row["observed"], CATALOG)[0]
    verdict = row["capture_ocr"]
    assert verdict["status"] == "verified"
    decision = {
        **proposal,
        **verdict,
        "page": 47,
        "item_ref": "source-page47",
        "crop_sha256": row["image_sha256"],
    }
    return {
        "decisions": [decision],
        "policy": producer["policy"],
        "producer": producer,
        "catalog_sha256": CATALOG["sha256"],
        "source_pdf_sha256": row["pdf_sha256"],
    }


def test_actual_preserved_reference_edge_repairs_only_with_title_and_source_evidence(
    tmp_path,
):
    ref = json.loads((FIXTURE / "reference.json").read_text())["reference"]
    article = next(e for e in ENTRIES if e["_key"] == "Alvarino1971a")
    conn, _ = build(
        tmp_path,
        {
            "citing": {
                "author": "Mapstone, G",
                "title": "A citing work",
                "year": "2009",
            },
            "article": article,
        },
    )
    folder = tmp_path / "documents/citing"
    (folder / "references.json").write_text(json.dumps({"references": [ref]}))
    phase2_references(conn, tmp_path)
    target = find_work(conn, "article")
    assert conn.execute("SELECT work_id FROM observation_work").fetchone()[0] != target
    raw_before = list(conn.execute("SELECT * FROM reference_observations"))
    report = source_report()
    (folder / "text.json").write_text(
        json.dumps(
            {"text": "Source text", "source_text_integrity": {"surnames": report}}
        )
    )
    phase2_references(conn, tmp_path)
    assert list(conn.execute("SELECT * FROM reference_observations")) == raw_before
    assert json.loads(raw_before[0][9]) == ["A Alvarifio"]
    assert conn.execute("SELECT cited_work_id FROM citations").fetchone()[0] == target
    reason = json.loads(
        conn.execute(
            "SELECT reasons_json FROM reference_observation_quality"
        ).fetchone()[0]
    )[0]
    assert reason["code"] == "source_ocr_surname_supported"
    assert reason["derived_author"] == "A Alvariño"
    changes = conn.total_changes
    assert phase2_references(conn, tmp_path) == (0, 0)
    assert conn.total_changes == changes
    # Source evidence changes rederive unchanged parsed observations again.
    report["decisions"][0]["status"] = "ocr_disagreement"
    (folder / "text.json").write_text(
        json.dumps(
            {"text": "Source text", "source_text_integrity": {"surnames": report}}
        )
    )
    phase2_references(conn, tmp_path)
    assert conn.execute("SELECT cited_work_id FROM citations").fetchone()[0] != target
    assert list(conn.execute("SELECT * FROM reference_observations")) == raw_before


def test_missing_or_damaged_reference_titles_and_conflicting_candidates_stay_raw():
    ref = json.loads((FIXTURE / "reference.json").read_text())["reference"]
    report = source_report()
    for title in (
        "",
        "A different title entirely",
        "Siphonophores of the Paciñc with an uncertain reading",
    ):
        observed = dict(ref, title=title)
        derived, reasons = supported_reference(observed, report)
        assert derived == observed
        assert reasons[0]["requires_source_review"]
    conflict = copy.deepcopy(report["decisions"][0])
    conflict["replacement"] = "Another"
    report["decisions"].append(conflict)
    derived, reasons = supported_reference(ref, report)
    assert derived == ref and reasons[0]["requires_source_review"]


def test_unavailable_ocr_preserves_reviewable_candidate():
    proposal = propose("Alvarifio 1991", CATALOG)[0]
    assert (
        adjudicate_crop(b"", proposal, {"models": {}})["status"]
        == "ocr_language_unavailable"
    )


def test_optional_source_pdf_replay_preserves_original_and_reaches_real_chunker(
    tmp_path, monkeypatch
):
    library = os.environ.get("CORPUS_LIBRARY_DIR")
    if not library:
        pytest.skip("Set CORPUS_LIBRARY_DIR for original PDF and chunking replay")
    import fitz
    from docling_core.types.doc import (
        DoclingDocument,
        DocItemLabel,
        Size,
        BoundingBox,
        ProvenanceItem,
    )
    from pipeline.chunking import chunk_text

    producer = live_producer()
    doc = DoclingDocument(name="Source surname replay")
    pdf_path = Path(library) / "M/Mapstone2009.pdf"
    with fitz.open(pdf_path) as pdf:
        width, height = pdf[46].rect.width, pdf[46].rect.height
    doc.add_page(page_no=47, size=Size(width=width, height=height))
    for row in SAMPLE[:2]:
        l, t, r, b = row["bbox"]
        doc.add_text(
            label=DocItemLabel.TEXT,
            text=row["observed"],
            orig=row["observed"],
            prov=ProvenanceItem(
                page_no=47,
                bbox=BoundingBox(l=l, t=t, r=r, b=b, coord_origin="TOPLEFT"),
                charspan=(0, len(row["observed"])),
            ),
        )
    report = recover_citation_surnames(doc, pdf_path, CATALOG, producer=producer)
    assert [d["status"] for d in report["decisions"]] == ["verified", "verified"]
    assert all(
        "Alvariño" in item.text and "Alvarifio" in item.orig for item in doc.texts
    )
    before = copy.deepcopy(doc.export_to_dict())
    assert (
        recover_citation_surnames(doc, pdf_path, CATALOG, producer=producer)[
            "decisions"
        ]
        == []
    )
    assert doc.export_to_dict() == before
    doc.save_as_json(tmp_path / "docling_doc.json")
    (tmp_path / "text.json").write_text(json.dumps({"text": doc.export_to_markdown()}))
    monkeypatch.setenv("HF_HUB_OFFLINE", "1")
    chunk_text(tmp_path / "text.json", chunks_output=tmp_path / "chunks.json")
    chunks = json.loads((tmp_path / "chunks.json").read_text())
    assert chunks["chunker"] == "hybrid_chunker"
    assert any("Alvariño 1971" in c["text"] for c in chunks["chunks"])
    assert not any("Alvarifio" in c["text"] for c in chunks["chunks"])
    source_decisions = [
        decision
        for chunk in chunks["chunks"]
        for decision in chunk.get("text_integrity", [])
        if decision.get("policy") == report["policy"]
    ]
    assert {decision["year"] for decision in source_decisions} == {1971, 1991}
    assert all(decision["original"] == "Alvarifio" for decision in source_decisions)


@pytest.mark.parametrize("rotation", [0, 90, 180, 270])
def test_rotated_source_anchor_and_upright_crop_preserve_printed_alternative(rotation):
    import fitz
    from pipeline.surname_recovery import _source_crop, _render_source_crop

    proposal = propose("Alvarifio 1991", CATALOG)[0]
    producer = live_producer()
    with fitz.open() as pdf:
        page = pdf.new_page(width=360, height=180)
        page.insert_text((45, 70), "Alvarifio 1991", fontname="tiro", fontsize=12)
        native_box = page.search_for("Alvarifio 1991")[0]
        page.set_rotation(rotation)
        displayed_box = native_box * page.rotation_matrix
        crop = _source_crop(page, displayed_box, proposal)
        assert crop is not None
        assert not (crop & displayed_box).is_empty
        png = _render_source_crop(page, crop, 600)
        assert page.rotation == rotation
        image = fitz.Pixmap(png)
        assert image.width > image.height
        result = adjudicate_crop(png, proposal, producer)
        assert result["status"] == "source_supports_observed_spelling"


@pytest.mark.parametrize(
    "status",
    [
        "verified",
        "source_supports_observed_spelling",
        "quoted_or_sic_context",
        "ocr_disagreement",
        "ocr_language_unavailable",
    ],
)
def test_extraction_quality_gate_reports_only_unresolved_surname_decisions(
    tmp_path, monkeypatch, status
):
    import fitz
    from docling_core.types.doc import (
        DoclingDocument,
        DocItemLabel,
        Size,
        BoundingBox,
        ProvenanceItem,
    )
    from pipeline.stages import _run_quality_gates
    import pipeline.surname_recovery as recovery

    pdf_path = tmp_path / "source.pdf"
    with fitz.open() as pdf:
        page = pdf.new_page(width=360, height=180)
        page.insert_text((45, 70), "Alvarifio 1991", fontname="tiro", fontsize=12)
        bbox = page.search_for("Alvarifio 1991")[0]
        pdf.save(pdf_path)
    observed = "Alvarifio 1991"
    if status == "quoted_or_sic_context":
        observed = "“Alvarifio 1991” [sic]"
    doc = DoclingDocument(name="Source quality gate")
    doc.add_page(page_no=1, size=Size(width=360, height=180))
    doc.add_text(
        label=DocItemLabel.TEXT,
        text=observed,
        orig=observed,
        prov=ProvenanceItem(
            page_no=1,
            bbox=BoundingBox(
                l=bbox.x0, t=bbox.y0, r=bbox.x1, b=bbox.y1, coord_origin="TOPLEFT"
            ),
            charspan=(0, len(observed)),
        ),
    )
    monkeypatch.setattr(
        recovery,
        "adjudicate_crop",
        lambda *args: {"status": status, "replacement": "Alvariño", "ocr_outputs": []},
    )
    report = recovery.recover_citation_surnames(
        doc, pdf_path, CATALOG, producer=CAPTURE["capture_producer"]
    )
    assert len(report["decisions"]) == 1
    assert report["decisions"][0]["status"] == status
    needs_review = status in {"ocr_disagreement", "ocr_language_unavailable"}
    assert report["unresolved"] == (report["decisions"] if needs_review else [])
    (tmp_path / "text.json").write_text(
        json.dumps(
            {
                "text": "Recovered scientific text. " * 250,
                "pages": 1,
                "source_text_integrity": {"surnames": report},
            }
        )
    )
    flags = [
        flag
        for flag in _run_quality_gates(tmp_path)
        if flag["gate"] == "source_text_integrity"
    ]
    assert len(flags) == int(needs_review)
    if needs_review:
        assert flags[0]["severity"] == "warning"
        assert flags[0]["metric"] == 1
        assert flags[0]["detail"].startswith(
            "surnames: 1 source-text candidates need review"
        )
