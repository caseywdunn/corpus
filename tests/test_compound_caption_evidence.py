import json

from pipeline.figures import resolve_compound_figures


def test_compound_split_marks_inferred_caption_as_uncertain(tmp_path):
    figures_file = tmp_path / "figures.json"
    figures_file.write_text(json.dumps({
        "figures": [{
            "figure_id": "docling_1",
            "figure_type": "figure",
            "figure_number": "3",
            "page": 8,
            "bbox": [10, 20, 300, 500],
            "bbox_coord_system": "pdf_pts_bottom_left",
            "pass3_status": "completed_compound",
            "panels_from_caption": [{"label": "A"}, {"label": "B"}],
            "rois": [
                {"type": "panel", "label": "A", "parent_figure_index": 0},
                {"type": "panel", "label": "B", "parent_figure_index": 0},
                {"type": "panel", "label": "C", "parent_figure_index": 1},
                {"type": "panel", "label": "D", "parent_figure_index": 1},
            ],
        }],
        "missing_figures": [{
            "figure_number": "4",
            "caption_text_candidate": (
                "Fig. 4. (A) colony; (B) bud; (C) young and (D) mature "
                "developing nectophores."
            ),
        }],
        "total_missing_figures": 1,
    }), encoding="utf-8")

    result = resolve_compound_figures(figures_file)
    data = json.loads(figures_file.read_text(encoding="utf-8"))
    recovered = data["figures"][1]

    assert result["resolved"] == 1
    assert recovered["figure_number"] == "4"
    assert recovered["figure_number_source"] == (
        "compound_split_from_missing_figures"
    )
    assert recovered["caption_status"] == "uncertain"
    assert recovered["caption_confidence"] == "low"
    assert recovered["caption_kind"] == "prose_caption"
    assert recovered["caption_page"] is None
    assert recovered["page"] == 8
    assert recovered["bbox"] == [10, 20, 300, 500]
    assert recovered["caption_candidates"][0]["chosen"] is True
    assert data["missing_figures"] == []
