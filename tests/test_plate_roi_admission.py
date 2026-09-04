"""Grouped plate captions feed numeric regions to Pass 3 (#203)."""
from __future__ import annotations

import json

from PIL import Image

from pipeline import figure_passes
from pipeline.figure_passes import (
    _pass25_annotate_figures,
    _pass3a_annotate_rois,
    _pass3b_annotate_rois,
)
from pipeline.figures import detect_figure_rois, detect_figure_rois_via_vision


def _shared_plate(tmp_path):
    image = tmp_path / "plate_31.png"
    Image.new("RGB", (200, 100), "white").save(image)
    return {
        "figures": [
            {
                "figure_id": "docling_18",
                # A legend-backed plate is eligible even if the conservative
                # classifier left its host in the review bucket.
                "figure_type": "unclassified",
                "figure_number": "31",
                "caption_text": "Fig. 31. Colony.",
                "filename": image.name,
                "file_path": str(image),
            },
            {
                "figure_id": "docling_18_fig32",
                "figure_type": "figure",
                "figure_number": "32",
                "caption_text": "Fig. 32. Eudoxid.",
                "filename": image.name,
                "file_path": str(image),
                "caption_source": "plate_legend",
                "shares_image_with": 18,
            },
        ],
    }


def _write_pass25_inputs(tmp_path):
    figures_file = tmp_path / "figures.json"
    text_file = tmp_path / "text.json"
    figures_file.write_text(json.dumps(_shared_plate(tmp_path)), encoding="utf-8")
    text_file.write_text(
        json.dumps({"text": "The eudoxid is illustrated in Fig. 32."}),
        encoding="utf-8",
    )
    return text_file, figures_file


def test_pass25_keeps_plate_numbers_distinct_from_panel_letters(tmp_path):
    text_file, figures_file = _write_pass25_inputs(tmp_path)

    _pass25_annotate_figures(text_file, figures_file)

    records = json.loads(figures_file.read_text(encoding="utf-8"))["figures"]
    host, sibling = records
    assert host["panels_from_caption"] == []
    assert sibling["panels_from_caption"] == []
    assert host["plate_roi_host"] is True
    assert [
        target["figure_number"]
        for target in host["plate_figures_from_caption"]
    ] == ["31", "32"]
    assert host["plate_missing_figure_crosscheck"] == ["32"]
    target_32 = host["plate_figures_from_caption"][1]
    assert target_32["kind"] == "figure"
    assert target_32["missing_figure_crosscheck"] is True


def test_pass3a_runs_shared_plate_once_and_distributes_rois(tmp_path, monkeypatch):
    text_file, figures_file = _write_pass25_inputs(tmp_path)
    _pass25_annotate_figures(text_file, figures_file)
    calls = []

    def fake_detect(image_path, targets, *, target_kind="panel"):
        calls.append((image_path, targets, target_kind))
        return {
            "rois": [
                {"type": "figure", "label": "31", "figure_number": "31",
                 "roi_px": [0, 0, 100, 100]},
                {"type": "figure", "label": "32", "figure_number": "32",
                 "roi_px": [100, 0, 200, 100]},
            ],
            "pass3_status": "completed",
            "ocr_token_count": 2,
            "image_size_px": [200, 100],
        }

    monkeypatch.setattr(figure_passes, "detect_figure_rois", fake_detect)
    _pass3a_annotate_rois(figures_file)

    assert len(calls) == 1
    assert calls[0][2] == "figure"
    records = json.loads(figures_file.read_text(encoding="utf-8"))["figures"]
    assert records[0]["rois"][0]["figure_number"] == "31"
    assert records[1]["rois"][0]["figure_number"] == "32"
    assert all(record["pass3_target_kind"] == "figure" for record in records)
    assert all(record["pass3_status"] == "completed" for record in records)


def test_ocr_numeric_targets_require_the_caption_allow_list(tmp_path, monkeypatch):
    image = tmp_path / "plate.png"
    Image.new("RGB", (200, 100), "white").save(image)
    monkeypatch.setattr(
        "pipeline.figures._ocr_figure_tokens",
        lambda _path: [
            {"text": "31.", "conf": 90, "bbox_px": [5, 5, 15, 15]},
            {"text": "999", "conf": 99, "bbox_px": [90, 5, 110, 15]},
        ],
    )

    result = detect_figure_rois(
        image,
        [{"label": "31"}, {"label": "32"}],
        target_kind="figure",
    )

    assert result["pass3_status"] == "partial_ocr"
    assert [roi["figure_number"] for roi in result["rois"]] == ["31", "32"]
    assert result["rois"][0]["roi_px"] is not None
    assert result["rois"][1]["roi_px"] is None


def test_vision_numeric_targets_use_embedded_figure_regions(tmp_path):
    image = tmp_path / "plate.png"
    Image.new("RGB", (200, 100), "white").save(image)

    class Backend:
        name = "vision:test"

        @staticmethod
        def detect_figure_panels(_path, _caption, labels):
            assert labels == ["31", "32"]
            return [
                {"type": "embedded_figure", "figure_number": "31",
                 "bbox_px": [0, 0, 100, 100], "source": "vision:test"},
                {"type": "embedded_figure", "figure_number": "32",
                 "bbox_px": [100, 0, 200, 100], "source": "vision:test"},
            ]

    result = detect_figure_rois_via_vision(
        image,
        [{"label": "31"}, {"label": "32"}],
        Backend(),
        target_kind="figure",
    )

    assert result["pass3_status"] == "completed"
    assert [roi["label"] for roi in result["rois"]] == ["31", "32"]


def test_bare_plate_discovery_keeps_only_distinct_high_confidence_numbers(tmp_path):
    image = tmp_path / "plate_xvi.png"
    Image.new("RGB", (300, 200), "white").save(image)

    class Backend:
        name = "vision:test"

        @staticmethod
        def detect_figure_panels(_path, _caption, labels):
            assert labels == []
            return [
                {"type": "embedded_figure", "figure_number": "1",
                 "bbox_px": [0, 0, 100, 100], "confidence": 0.95,
                 "source": "vision:test"},
                # Some models put a visible number in the panel collection.
                {"type": "panel", "label": "2",
                 "bbox_px": [100, 0, 200, 100], "confidence": 0.90,
                 "source": "vision:test"},
                {"type": "embedded_figure", "figure_number": "1965",
                 "bbox_px": [200, 0, 300, 100], "confidence": 0.99,
                 "source": "vision:test"},
                {"type": "embedded_figure", "figure_number": "3",
                 "bbox_px": [0, 100, 100, 200], "confidence": 0.40,
                 "source": "vision:test"},
                {"type": "embedded_figure", "figure_number": "1",
                 "bbox_px": [5, 5, 95, 95], "confidence": 0.85,
                 "source": "vision:test"},
            ]

    result = detect_figure_rois_via_vision(
        image, [], Backend(), caption_text="PLATE XVI",
        target_kind="figure_discovery",
    )

    assert result["pass3_status"] == "discovery_materialized"
    assert [roi["figure_number"] for roi in result["rois"]] == ["1", "2"]
    candidates = result["figure_number_candidates"]
    assert [c["figure_number"] for c in candidates if c["accepted"]] == ["1", "2"]
    assert candidates[2]["rejection_reason"] == "unsupported_number_format"
    assert candidates[3]["rejection_reason"] == "below_confidence_threshold"
    assert candidates[4]["rejection_reason"] == "lower_confidence_duplicate"


def test_bare_plate_discovery_requires_a_group(tmp_path):
    image = tmp_path / "plate_xvi.png"
    Image.new("RGB", (200, 100), "white").save(image)

    class Backend:
        name = "vision:test"

        @staticmethod
        def detect_figure_panels(_path, _caption, _labels):
            return [{
                "type": "embedded_figure", "figure_number": "1",
                "bbox_px": [0, 0, 100, 100], "confidence": 0.95,
                "source": "vision:test",
            }]

    result = detect_figure_rois_via_vision(
        image, [], Backend(), caption_text="PLATE XVI",
        target_kind="figure_discovery",
    )

    assert result["pass3_status"] == "discovery_insufficient_evidence"
    assert result["rois"] == []
    assert result["figure_number_candidates"][0]["rejection_reason"] == (
        "fewer_than_two_distinct_numbers"
    )


def test_pass3b_materializes_bare_plate_numbers_without_inventing_captions(tmp_path):
    image = tmp_path / "plate_16.png"
    Image.new("RGB", (200, 100), "white").save(image)
    figures_file = tmp_path / "figures.json"
    figures_file.write_text(json.dumps({"figures": [{
        "figure_id": "docling_174",
        "figure_type": "plate",
        "figure_number": "16",
        "caption_text": "PLATE XVI",
        "caption_kind": "bare_label",
        "caption_status": "bound",
        "caption_confidence": "high",
        "filename": image.name,
        "file_path": str(image),
        "page": 266,
        "bbox": [0, 0, 500, 700],
        "bbox_coord_system": "pdf_pts_bottom_left",
        "extraction_method": "docling",
    }]}), encoding="utf-8")

    class Backend:
        name = "vision:test"

        @staticmethod
        def detect_figure_panels(_path, caption, labels):
            assert caption == "PLATE XVI"
            assert labels == []
            return [
                {"type": "embedded_figure", "figure_number": "1",
                 "bbox_px": [0, 0, 100, 100], "confidence": 0.95,
                 "source": "vision:test"},
                {"type": "embedded_figure", "figure_number": "2",
                 "bbox_px": [100, 0, 200, 100], "confidence": 0.90,
                 "source": "vision:test"},
            ]

    _pass3b_annotate_rois(figures_file, Backend())
    # A refresh replaces the derived generation rather than appending it.
    _pass3b_annotate_rois(figures_file, Backend())

    records = json.loads(figures_file.read_text(encoding="utf-8"))["figures"]
    assert len(records) == 3
    host = records[0]
    assert host["pass3_status"] == "discovery_materialized"
    assert host["plate_number_discovery"]["accepted_count"] == 2
    discovered = records[1:]
    assert [record["figure_number"] for record in discovered] == ["1", "2"]
    assert all(record["figure_number_source"] == "vision_plate_discovery"
               for record in discovered)
    assert all(record["image_shared_with"] == "docling_174"
               for record in discovered)
    assert all(record["caption_text"] == "" for record in discovered)
    assert all(record["caption_status"] == "unbound" for record in discovered)
