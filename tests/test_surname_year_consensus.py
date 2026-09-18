"""Surname agreement and independently valid dates on recorded OCR evidence."""

import copy
import json
from types import SimpleNamespace

import pytest

import pipeline.surname_recovery as recovery
from tests.test_surname_recovery import CATALOG, CAPTURE, FIXTURE, SAMPLE

HOSTED = json.loads((FIXTURE / "hosted_ocr_553.json").read_text())


def replay(monkeypatch, readings, *, modes=(6, 13)):
    outputs = iter(readings)
    monkeypatch.setattr(
        recovery.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(stdout=next(outputs).encode()),
    )
    producer = copy.deepcopy(CAPTURE["capture_producer"])
    producer.update(executable="recorded-ocr", ocr_modes=list(modes))
    proposal = recovery.propose("Alvarifio 1981a", CATALOG)[0]
    return recovery.adjudicate_crop(b"recorded source crop", proposal, producer)


@pytest.mark.parametrize("reverse", [False, True])
def test_actual_hosted_dual_name_readings_need_only_one_valid_year(
    monkeypatch, reverse
):
    row = next(row for row in SAMPLE if row["id"] == HOSTED["case"])
    assert row["image_sha256"] == HOSTED["image_sha256"]
    readings = [row["text"] for row in HOSTED["ocr_outputs"]]
    if reverse:
        readings.reverse()
    result = replay(monkeypatch, readings)
    assert result["status"] == "verified"
    assert result["replacement"] == "Alvariño"
    assert [row["text"] for row in result["ocr_outputs"]] == readings
    assert [row["names"] for row in result["ocr_outputs"]] == [["Alvariño"]] * 2
    evidence = result["year_evidence"]
    assert evidence["source_anchor_year"] == 1981
    assert evidence["matching_readings"] == [
        {"language": "spa", "psm": 6 if reverse else 13}
    ]
    assert evidence["invalid_tokens"] == ["19814"]
    assert evidence["conflicting_valid_years"] == []
    invalid = next(row for row in result["ocr_outputs"] if "19814" in row["text"])
    assert invalid["year_tokens"] == ["19814"]
    assert invalid["valid_years"] == invalid["adjacent_years"] == []


@pytest.mark.parametrize(
    "readings, conflicting",
    [
        (["Alvariño 1982", "Alvariño 1981a"], [1982]),
        (["Alvariño 1981a", "Alvariño 1982"], [1982]),
        (["Alvariño1982", "Alvariño 1981a"], [1982]),
        (["Alvariño 1981a 1982", "Alvariño 1981a"], [1982]),
        (["Alvariño 19814", "Alvariño 1981al"], []),
        (["Alvariño 19814", "Alvariño 1981æ"], []),
        (["Alvariño 19814", "Alvariño 1981_"], []),
        (["Alvariño 19814", "Alvarifio 1981a"], []),
        (["Alvarifio 1981a", "Alvariño 19814"], []),
        (["Alvariñ 19814", "Alvariño 1981a"], []),
        (["Alvariño 1981a Alvariño 1981a", "Alvariño 1981a"], []),
        (["Alvariño 19814, note from 1981", "Alvariño 19814"], []),
    ],
)
def test_conflicts_invalid_dates_and_incomplete_names_never_repair(
    monkeypatch, readings, conflicting
):
    result = replay(monkeypatch, readings)
    assert result["status"] == "ocr_disagreement"
    assert "replacement" not in result
    assert result["year_evidence"]["conflicting_valid_years"] == conflicting
    assert [row["text"] for row in result["ocr_outputs"]] == readings


def test_one_segmentation_mode_cannot_verify(monkeypatch):
    result = replay(monkeypatch, ["Alvariño 1981a"], modes=(6,))
    assert result["status"] == "ocr_disagreement"
    assert "replacement" not in result


def test_dual_observed_spelling_with_valid_year_stays_unchanged(monkeypatch):
    result = replay(monkeypatch, ["Alvarifio 19814", "Alvarifio 1981a"])
    assert result["status"] == "source_supports_observed_spelling"
    assert "replacement" not in result


def test_surname_recovery_preserves_source_year_suffix_and_raw_ocr(
    tmp_path, monkeypatch
):
    import fitz
    from docling_core.types.doc import (
        DoclingDocument,
        DocItemLabel,
        Size,
        BoundingBox,
        ProvenanceItem,
    )

    source = tmp_path / "source.pdf"
    observed = "Alvarifio 1981a"
    with fitz.open() as pdf:
        page = pdf.new_page(width=360, height=180)
        page.insert_text((45, 70), observed, fontname="tiro", fontsize=12)
        box = page.search_for(observed)[0]
        pdf.save(source)
    document = DoclingDocument(name="Surname and date evidence")
    document.add_page(page_no=1, size=Size(width=360, height=180))
    document.add_text(
        label=DocItemLabel.TEXT,
        text=observed,
        orig=observed,
        prov=ProvenanceItem(
            page_no=1,
            bbox=BoundingBox(
                l=box.x0, t=box.y0, r=box.x1, b=box.y1, coord_origin="TOPLEFT"
            ),
            charspan=(0, len(observed)),
        ),
    )
    readings = iter(row["text"] for row in HOSTED["ocr_outputs"])
    monkeypatch.setattr(
        recovery.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(stdout=next(readings).encode()),
    )
    producer = dict(CAPTURE["capture_producer"], executable="recorded-ocr")
    report = recovery.recover_citation_surnames(
        document, source, CATALOG, producer=producer
    )
    assert document.texts[0].text == "Alvariño 1981a"
    assert document.texts[0].orig == observed
    assert report["unresolved"] == []
    decision = report["decisions"][0]
    assert decision["status"] == "verified"
    assert decision["ocr_outputs"][0]["text"] == "[Alvariño, 19814]."
    assert decision["year_evidence"]["invalid_tokens"] == ["19814"]
    assert document.texts[0].meta.corpus__surname_recovery == report["decisions"]
    # A malformed native anchor cannot borrow the OCR suffix interpretation.
    assert recovery.propose("Alvarifio 19814", CATALOG) == []


def test_v3_policy_invalidates_extraction_and_all_consumers_only():
    from pipeline.build_inputs import config_fingerprints

    current = recovery.surname_recovery_producer(CATALOG)
    assert current["policy"] == "curated-author-year-source-ocr-consensus-v3"
    previous = dict(current, policy="curated-author-year-source-ocr-consensus-v2")
    before = config_fingerprints({}, panel_mode="ocr", surname_producer=previous)
    after = config_fingerprints({}, panel_mode="ocr", surname_producer=current)
    changed = {stage for stage in before if before[stage] != after[stage]}
    assert changed == {
        "docling_extraction",
        "text_chunking",
        "taxa_and_lexicon_extraction",
        "figure_materialization",
        "figure_crossref",
    }
