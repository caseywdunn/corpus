"""The operator figure report makes caption decisions inspectable (#195)."""

import json

from pipeline.figures import generate_figures_report


def test_report_shows_caption_status_and_candidate_decisions(tmp_path):
    (tmp_path / "metadata.json").write_text(json.dumps({"title": "A paper"}))
    (tmp_path / "figures.json").write_text(json.dumps({
        "figures": [{
            "figure_id": "docling_1",
            "figure_type": "figure",
            "figure_number": "8",
            "page": 12,
            "caption_text": "FIGURE 8. A complete caption.",
            "caption_source": "heuristic_proximity",
            "caption_status": "bound",
            "caption_confidence": "medium",
            "caption_kind": "prose_caption",
            "caption_page_distance": 0,
            "caption_candidates": [
                {
                    "caption_text": "FIGURE 8. A complete caption.",
                    "caption_source": "heuristic_proximity",
                    "caption_page": 12,
                    "distance_pts": 8.5,
                    "confidence": "medium",
                    "chosen": True,
                    "rejection_reason": None,
                },
                {
                    "caption_text": "FIGURE 9. Wrong local owner.",
                    "caption_source": "heuristic_proximity",
                    "caption_page": 13,
                    "distance_pts": 20.0,
                    "confidence": "low",
                    "chosen": False,
                    "rejection_reason": "substantial_picture_on_caption_page",
                },
            ],
        }],
    }))

    report = generate_figures_report(tmp_path)
    rendered = report.read_text()

    assert "confidence: medium" in rendered
    assert "prose_caption" in rendered
    assert "caption evidence" in rendered
    assert "substantial_picture_on_caption_page" in rendered
