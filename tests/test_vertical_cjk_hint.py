"""Vertical-CJK operator hint (#196).

Tesseract ships a separate model for vertically-set text in each CJK script
and detection never selects one. Unlike the Fraktur companion — where
`deu`+`deu_latf` are added together and degrade gracefully — the vertical
models cannot be added automatically. Measured against a hand transcription
of the same pages, as the fraction of printed words recovered:

                       horizontal Japanese   vertical Japanese
    jpn                       0.75                 0.25
    jpn_vert                  0.21                 0.57
    jpn_vert+jpn               --                  0.19

The union is worse than either pack alone and each is catastrophic on the
other's direction, so there is no static rule to apply — and `ocrmypdf` takes
one pack list per document, so it cannot be decided per page either. The
supported answer is the `ocrlang` directive (#176); these tests cover the hint
that tells an operator when to reach for it.
"""
from __future__ import annotations

from pipeline.scan import _VERTICAL_COMPANION


def test_companion_table_covers_every_vertical_model_tesseract_ships():
    assert _VERTICAL_COMPANION == {
        "jpn": "jpn_vert",
        "chi_sim": "chi_sim_vert",
        "chi_tra": "chi_tra_vert",
        "kor": "kor_vert",
    }


def test_fraktur_companion_is_not_in_this_table():
    """`deu` → `deu_latf` is a different mechanism and must stay separate.

    Fraktur packs are added *together* by `_resolve_tesseract_packs` because a
    surplus pack is cheap there. Doing the same for vertical CJK measured
    worse than doing nothing, so the two must not be merged into one table by
    a later tidy-up.
    """
    assert "deu" not in _VERTICAL_COMPANION


def _hint(monkeypatch, langs, honored=False, available=None):
    """Call the real hint function, so these tests break if it changes."""
    from pipeline import scan

    monkeypatch.setattr(
        scan, "_available_tesseract_langs",
        lambda: available if available is not None
        else {"jpn", "jpn_vert", "eng", "deu", "chi_sim", "chi_sim_vert"},
    )
    return scan._vertical_cjk_hint(langs, {"ocrlang_honored": honored})


def test_hint_fires_for_a_cjk_pack_and_names_the_vertical_model(monkeypatch):
    msg = _hint(monkeypatch, ["jpn", "eng"])
    assert "ocrlang = {jpn_vert}" in msg
    # The warning against the union is the part that is easy to get wrong,
    # since the Fraktur precedent points the other way.
    assert "worse than either alone" in msg


def test_hint_is_silent_for_a_latin_document(monkeypatch):
    assert _hint(monkeypatch, ["deu", "deu_latf", "eng"]) is None


def test_hint_is_silent_when_the_operator_already_pinned_the_packs(monkeypatch):
    """They have made the call themselves; repeating it is noise."""
    assert _hint(monkeypatch, ["jpn_vert"], honored=True) is None


def test_hint_is_silent_when_the_vertical_model_is_not_installed(monkeypatch):
    """Advice to pin a pack that isn't there would send an operator in circles."""
    assert _hint(monkeypatch, ["jpn", "eng"], available={"jpn", "eng"}) is None


def test_hint_handles_chinese_too(monkeypatch):
    assert "ocrlang = {chi_sim_vert}" in _hint(monkeypatch, ["chi_sim"])


def test_prepare_pdf_actually_logs_the_hint():
    """Guard the wiring, not just the function: the string must reach the log
    line the operator reads, beside `langs=`."""
    import inspect
    from pipeline import scan
    src = inspect.getsource(scan.prepare_pdf)
    assert "_vertical_cjk_hint(langs, detection_result)" in src
