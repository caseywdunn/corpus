import json

import fitz
import pytest

from pipeline.page_report import generate_page_report, parse_page_spec


def _write_fixture(hash_dir, page_count=1):
    hash_dir.mkdir()
    pdf = fitz.open()
    for page_number in range(page_count):
        page = pdf.new_page(width=200, height=300)
        page.insert_text(
            (20, 30),
            "Figure 1. Nanomia nectophore" if page_number == 0 else "Page 2",
        )
    pdf.save(hash_dir / "processed.pdf")
    pdf.close()

    text = "Figure 1. Nanomia nectophore"
    (hash_dir / "docling_doc.json").write_text(json.dumps({
        "version": "test-docling",
        "texts": [{
            "label": "caption",
            "text": text,
            "prov": [{
                "page_no": 1,
                "bbox": {
                    "l": 20, "b": 35, "r": 180, "t": 55,
                    "coord_origin": "BOTTOMLEFT",
                },
                "charspan": [0, len(text)],
            }],
        }],
    }), encoding="utf-8")
    (hash_dir / "metadata.json").write_text(json.dumps({
        "title": "Audit fixture", "year": 2026,
        "authors": [{"forename": "A.", "surname": "Author"}],
    }), encoding="utf-8")
    (hash_dir / "scan_detection.json").write_text(json.dumps({
        "file_type": "born_digital",
        "tesseract_packs": ["eng"],
        "keeppages_selected": [41],
    }), encoding="utf-8")
    (hash_dir / "figures.json").write_text(json.dumps({
        "pipeline_version": "1.3-test",
        "figures": [{
            "figure_id": "docling_1",
            "figure_number": "1",
            "page": 1,
            "bbox": [20, 70, 180, 270],
            "bbox_coord_system": "pdf_pts_bottom_left",
            "image_size_px": [160, 200],
            "caption_text": text,
            "caption_source": "heuristic_proximity",
            "caption_status": "uncertain",
            "caption_confidence": "low",
            "caption_kind": "prose_caption",
            "caption_page_distance": 0,
            "caption_candidates": [
                {
                    "caption_text": text,
                    "caption_page": 1,
                    "caption_bbox": [20, 35, 180, 55],
                    "caption_source": "heuristic_proximity",
                    "confidence": "low",
                    "chosen": True,
                    "distance_pts": 15,
                    "rejection_reason": None,
                },
                {
                    "caption_text": "Figure 2. Wrong owner",
                    "caption_page": 1,
                    "caption_bbox": [20, 5, 180, 20],
                    "caption_source": "heuristic_proximity",
                    "confidence": "low",
                    "chosen": False,
                    "distance_pts": 50,
                    "rejection_reason": "lower_ranked_candidate",
                },
            ],
            "rois": [{"type": "panel", "label": "A", "roi_px": [0, 0, 80, 100]}],
        }],
    }), encoding="utf-8")


def test_page_spec_parses_ranges_and_rejects_out_of_bounds():
    assert parse_page_spec("3,1-2,2", 4) == [1, 2, 3]
    with pytest.raises(ValueError, match="outside"):
        parse_page_spec("5", 4)


def test_page_report_joins_page_text_figures_rois_and_caption_evidence(tmp_path):
    hash_dir = tmp_path / "abc123def456"
    _write_fixture(hash_dir)

    output = generate_page_report(hash_dir, pages=[1])
    rendered = output.read_text(encoding="utf-8")

    assert output == hash_dir / "page_report.html"
    assert "data:image/jpeg;base64," in rendered
    assert "Processed page 1 · source page 41" in rendered
    assert "[caption] Figure 1. Nanomia nectophore" in rendered
    assert "docling_1 · Fig. 1" in rendered
    assert "status <strong>uncertain</strong>" in rendered
    assert "Figure 2. Wrong owner" in rendered
    assert "lower_ranked_candidate" in rendered
    assert '"figures":[[20.0,30.0,180.0,230.0,"1"]]' in rendered
    assert '"rois":[[20.0,30.0,100.0,130.0,"A"]]' in rendered
    assert '"chosen":[[20.0,245.0,180.0,265.0,"1"]]' in rendered
    assert '"rejected":[[20.0,280.0,180.0,295.0,"1"]]' in rendered


def test_page_report_is_rerunnable_and_honors_custom_output(tmp_path):
    hash_dir = tmp_path / "abc123def456"
    _write_fixture(hash_dir)
    output = tmp_path / "reports" / "one.html"

    first = generate_page_report(hash_dir, output, pages=[1])
    second = generate_page_report(hash_dir, output, pages=[1])

    assert first == second == output.resolve()
    assert output.is_file()


def test_safety_limit_requires_explicit_zero_override(tmp_path):
    hash_dir = tmp_path / "abc123def456"
    _write_fixture(hash_dir, page_count=2)

    with pytest.raises(ValueError, match="safety limit"):
        generate_page_report(hash_dir, pages=[1, 2], max_pages=1)
    assert generate_page_report(hash_dir, pages=[1, 2], max_pages=0).is_file()
