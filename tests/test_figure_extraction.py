"""Tests for figure extraction quality."""

import json

import pytest


def test_min_figure_count(paper_figures):
    """Pipeline finds at least the expected number of figures."""
    gt, fig_data, _doc_dir = paper_figures
    min_count = gt["figures"].get("min_count")
    if min_count is None:
        pytest.skip("No min_count specified")
    actual = fig_data.get("total_figures", len(fig_data.get("figures", [])))
    assert actual >= min_count, f"Too few figures: {actual} < {min_count}"


def test_expected_figures(paper_figures):
    """Specific expected figures are found with correct attributes."""
    gt, fig_data, _doc_dir = paper_figures
    expected = gt["figures"].get("expected", [])
    if not expected:
        pytest.skip("No specific figures to check")
    figures = fig_data.get("figures", [])
    for exp in expected:
        fig_num = exp.get("figure_number")
        # Figure numbers may be stored as int or str
        matches = [
            f for f in figures
            if str(f.get("figure_number", "")) == str(fig_num)
        ]
        assert matches, f"Figure {fig_num} not found in extraction"
        fig = matches[0]
        if "caption_contains" in exp:
            caption = (fig.get("caption_text") or "").lower()
            assert exp["caption_contains"].lower() in caption, (
                f"Figure {fig_num} caption missing {exp['caption_contains']!r}; "
                f"got: {caption[:120]!r}"
            )
        if "figure_type" in exp:
            assert fig.get("figure_type") == exp["figure_type"], (
                f"Figure {fig_num} type: expected {exp['figure_type']!r}, "
                f"got {fig.get('figure_type')!r}"
            )


def test_figure_files_exist(paper):
    """Every figure in figures.json has a corresponding file on disk."""
    _gt, doc_dir = paper
    fig_path = doc_dir / "figures.json"
    if not fig_path.exists():
        pytest.skip("No figures.json")
    with open(fig_path) as f:
        fig_data = json.load(f)
    for fig in fig_data.get("figures", []):
        fname = fig.get("filename")
        if fname:
            assert (doc_dir / "figures" / fname).exists(), (
                f"Figure file missing: {fname}"
            )
