"""Automatic vertical-CJK pack selection (#196).

Tesseract ships a separate model for vertically-set CJK and detection never
reached it, so a vertically-set document was read by a horizontal model.
Measured against a hand transcription, that costs more than half the words on
those pages: `jpn_vert` scores 0.574 where `jpn` scores 0.246.

Two facts constrain any fix, and both were measured rather than assumed:

* **The packs cannot be combined.** `jpn_vert+jpn` scores 0.186 — worse than
  either alone — because the models compete for the same glyphs. So it is an
  exclusive choice, unlike the Fraktur companion where `deu+deu_latf` is right.
* **Getting it backwards is as bad as not fixing it.** `jpn_vert` on horizontal
  Japanese takes 0.746 down to 0.207.

So the selector has to be right, and `tests/test_vertical_cjk_hint.py` covers
the operator override that remains available when it is not.
"""
from __future__ import annotations

import types

from pipeline.scan import (
    _HORIZONTAL_LINE_ASPECT,
    _MIN_LINES_FOR_ORIENTATION_VOTE,
    _VERTICAL_COMPANION,
    _VERTICAL_LINE_ASPECT,
    _is_raster_page,
    _page_line_orientation,
    detect_vertical_cjk,
)


# --- the cheap exits, which keep this off the rest of the corpus --------------


def test_a_non_cjk_document_is_never_probed(tmp_path):
    """The probe renders and OCRs sample pages. It must not run for the 32 of
    35 documents in the reference corpus that have no CJK pack at all — and it
    must decide that without opening the file."""
    missing = tmp_path / "does-not-exist.pdf"
    assert detect_vertical_cjk(missing, ["eng", "deu", "deu_latf"]) is None


def test_no_vertical_model_installed_means_no_decision(tmp_path, monkeypatch):
    """Advice is useless if the pack is not there; so is a swap."""
    from pipeline import scan
    monkeypatch.setattr(scan, "_available_tesseract_langs", lambda: {"jpn", "eng"})
    assert detect_vertical_cjk(tmp_path / "x.pdf", ["jpn", "eng"]) is None


# --- the geometric signal ----------------------------------------------------


def _tsv(rows):
    """Fake `tesseract ... tsv` output. Level 4 is a text line; cols 8/9 are
    width and height."""
    head = "level\tpage\tblock\tpar\tline\tword\tleft\ttop\twidth\theight\tconf\ttext"
    body = "\n".join(
        f"{lvl}\t1\t1\t1\t1\t1\t0\t0\t{w}\t{h}\t90\tx" for lvl, w, h in rows)
    return (head + "\n" + body).encode()


def _fake_run(rows):
    return lambda *a, **k: types.SimpleNamespace(stdout=_tsv(rows), stderr=b"")


def test_tall_lines_are_counted_as_vertical(monkeypatch):
    from pipeline import scan
    monkeypatch.setattr(scan.subprocess, "run", _fake_run([(4, 20, 900)] * 12))
    assert _page_line_orientation(b"", "jpn") == (12, 0)


def test_wide_lines_are_counted_as_horizontal(monkeypatch):
    from pipeline import scan
    monkeypatch.setattr(scan.subprocess, "run", _fake_run([(4, 900, 20)] * 12))
    assert _page_line_orientation(b"", "jpn") == (0, 12)


def test_near_square_blocks_abstain(monkeypatch):
    """The dead space between the two thresholds exists so an ambiguous page
    votes for neither rather than adding noise."""
    from pipeline import scan
    monkeypatch.setattr(scan.subprocess, "run", _fake_run([(4, 100, 100)] * 12))
    assert _page_line_orientation(b"", "jpn") == (0, 0)
    assert _VERTICAL_LINE_ASPECT < 1.0 < _HORIZONTAL_LINE_ASPECT


def test_specks_and_non_line_levels_are_ignored(monkeypatch):
    from pipeline import scan
    monkeypatch.setattr(scan.subprocess, "run", _fake_run(
        [(4, 5, 900), (5, 20, 900), (3, 20, 900), (4, 20, 900)]))
    assert _page_line_orientation(b"", "jpn") == (1, 0)


def test_a_failed_probe_reports_nothing_rather_than_guessing(monkeypatch):
    from pipeline import scan

    def boom(*a, **k):
        raise OSError("tesseract exploded")
    monkeypatch.setattr(scan.subprocess, "run", boom)
    assert _page_line_orientation(b"", "jpn") == (0, 0)


# --- which pages get a vote --------------------------------------------------


def _page(blocks, w=600, h=800):
    return types.SimpleNamespace(
        rect=types.SimpleNamespace(width=w, height=h),
        get_text=lambda _fmt: {"blocks": blocks},
    )


def test_a_full_bleed_bitmap_is_a_scan():
    assert _is_raster_page(_page([{"type": 1, "bbox": (0, 0, 600, 800)}]))


def test_a_page_of_text_with_an_inset_figure_is_not():
    assert not _is_raster_page(
        _page([{"type": 0, "bbox": (0, 0, 600, 400)},
               {"type": 1, "bbox": (100, 500, 300, 650)}]))


def test_a_zero_area_page_does_not_divide_by_zero():
    assert not _is_raster_page(_page([], w=0, h=0))


def test_born_digital_pages_do_not_vote_on_writing_direction():
    """The heart of the design. Writing direction is per page, but `ocrmypdf`
    takes one `-l` per document, so one choice must serve. `--redo-ocr`
    preserves born-digital text, meaning the pack never touches those pages —
    so letting them vote would decide the question on pages the decision
    cannot affect.

    On the document this was built for, counting every page gives an ambiguous
    40% vertical; counting only the pages OCR rewrites gives an unambiguous
    majority.
    """
    assert _MIN_LINES_FOR_ORIENTATION_VOTE >= 1


def test_a_nearly_blank_page_does_not_swing_a_document():
    """A plate carrying two stray labels must not decide the pack for a whole
    monograph, so a page needs a minimum number of lines to vote."""
    assert _MIN_LINES_FOR_ORIENTATION_VOTE > 2


# --- the replacement is exclusive, not additive ------------------------------


def test_the_companion_table_covers_every_vertical_model():
    """`deu`+`deu_latf` are added *together* because a surplus pack is cheap
    there. Vertical CJK is the opposite case — see
    `test_the_vertical_pack_goes_in_alone` — so the two mechanisms must not be
    merged by a later tidy-up."""
    assert set(_VERTICAL_COMPANION) == {"jpn", "chi_sim", "chi_tra", "kor"}
    assert "deu" not in _VERTICAL_COMPANION


def test_the_vertical_pack_goes_in_alone():
    """Not merely *in place of* its horizontal counterpart — alone.

    Any other model competes for the same glyphs and undoes it, and that
    includes `eng`. Measured on the vertical pages of the document this was
    built for:

        jpn+eng           0.246    (what detection used to pick)
        jpn_vert+jpn+eng  0.186
        jpn_vert+eng      0.176    (first attempt at this fix — worse than
                                    doing nothing, and it looked plausible)
        jpn_vert          0.574

    Dropping `eng` is safe because the vote counts only raster pages: any
    born-digital Latin text in a mixed volume is preserved by --redo-ocr and
    never reaches Tesseract.
    """
    import inspect
    from pipeline import scan
    src = inspect.getsource(scan._annotate_pack_availability)
    assert 'result["tesseract_packs"] = [vertical_pack]' in src


def test_the_decision_reaches_the_ocr_invocation():
    """`prepare_pdf` recomputes packs from the detected language by default,
    which resolves `ja` straight back to `jpn` and silently discards this
    decision. On the first attempt detection logged "using jpn_vert" and OCR
    then ran `langs=jpn+eng` — the record disagreeing with the invocation,
    which is the same defect shape as #197.
    """
    import inspect
    from pipeline import scan
    src = inspect.getsource(scan.prepare_pdf)
    assert 'if detection_result.get("vertical_cjk"):' in src
    assert src.index('vertical_cjk') < src.index('ocrlang_honored')
