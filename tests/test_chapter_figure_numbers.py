"""Printed chapter figure identity survives captions and body links (#343)."""

from pipeline.figures import (_is_bare_figure_label, caption_figure_entries,
                              detect_missing_figures, link_chunks_to_figures,
                              parse_figure_number)


def test_source_checked_simpson_caption_is_one_figure():
    caption = "Figure 1-2 Calcitic spicules (SEMs)."
    assert parse_figure_number(caption) == "1-2"
    assert caption_figure_entries(caption) == [
        {"figure_number": "1-2", "caption_text": caption}]
    chunks = [{"chunk_id": "body", "text": "See Fig. 1-2."}]
    figures = [{"figure_id": "docling_4", "figure_number": "1-2"}]
    link_chunks_to_figures(chunks, figures)
    assert chunks[0]["figure_refs"] == ["docling_4"]
    assert figures[0]["referenced_in_chunks"] == ["body"]
    assert detect_missing_figures("Figure 1-2 Calcitic spicules (SEMs).\nSee Fig. 1-2.",
                                  ["1-2"]) == []


def test_real_range_remains_multiple_figures():
    caption = "Figs. 58-63. Illustrations."
    assert [r["figure_number"] for r in caption_figure_entries(caption)] == [
        str(n) for n in range(58, 64)]
    chunks = [{"chunk_id": "body", "text": "See Figs. 58-63."}]
    figures = [{"figure_id": str(n), "figure_number": str(n)} for n in range(58, 64)]
    link_chunks_to_figures(chunks, figures)
    assert chunks[0]["figure_refs"] == [str(n) for n in range(58, 64)]


def test_compound_reference_prefers_exact_printed_identity():
    chunks = [{"chunk_id": "body", "text": "Fig. 1-2 and Fig. 1."}]
    figures = [{"figure_id": "chapter", "figure_number": "1-2"},
               {"figure_id": "simple", "figure_number": "1"},
               {"figure_id": "second", "figure_number": "2"}]
    link_chunks_to_figures(chunks, figures)
    assert chunks[0]["figure_refs"] == ["chapter", "simple"]


def test_single_label_with_dash_is_not_assumed_a_chapter_if_abbreviated():
    assert [r["figure_number"] for r in caption_figure_entries(
        "Fig. 1-2. Two illustrations.")] == ["1", "2"]


def test_source_chapter_namespace_identifies_missing_compound_figure():
    missing = detect_missing_figures(
        "Fig. 1-3 Acicular spicules.\nFig. 1-2 Calcitic spicules.",
        ["1-3"],
    )
    assert [row["figure_number"] for row in missing] == ["1-2"]


def test_plural_numeric_range_does_not_become_chapter_identity():
    chunks = [{"chunk_id": "body", "text": "See Figs. 1-2."}]
    figures = [{"figure_id": "compound", "figure_number": "1-2"},
               {"figure_id": "first", "figure_number": "1"},
               {"figure_id": "second", "figure_number": "2"}]
    link_chunks_to_figures(chunks, figures)
    assert chunks[0]["figure_refs"] == ["first", "second"]


def test_bare_chapter_label_stays_a_bare_label():
    assert _is_bare_figure_label("Figure 1-2.")
    assert not _is_bare_figure_label("Figure 1-2 Calcitic spicules (SEMs).")
