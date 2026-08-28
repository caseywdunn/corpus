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


# ---------------------------------------------------------------------------
# A pin that contradicts measured page geometry (#188 follow-up)
# ---------------------------------------------------------------------------
#
# The hint above goes silent once `ocrlang` is honored, on the reasoning that
# the operator has made the call themselves. That reasoning holds for a
# hand-written tag and fails completely for a derived one.
#
# An annotation pass deriving `ocrlang` from `doclang` (#214, #215) cannot get
# this right: BCP-47 describes language and script, and vertical setting is
# typesetting — deliberately out of `bcp47_to_tesseract`'s scope. So the
# derivation emits `jpn` for a vertically-set Japanese paper every time, the
# pin overrides the geometric verdict, and the one mechanism that would have
# said so is disabled *by the thing that caused it*. Found on Kawamura1911a,
# the document #196 was written for: `jpn_vert` 0.574 against `jpn` 0.246.

def _annotate(ocrlang, packs, vertical):
    from pipeline import scan
    result = {"filename": "Kawamura1911a.pdf", "needs_ocr": True,
              "tesseract_packs": list(packs)}
    if vertical:
        result["vertical_cjk"] = True
    return scan._annotate_pack_availability(ocrlang, result)


def test_a_horizontal_pin_over_a_vertical_verdict_is_recorded(monkeypatch):
    from pipeline import scan
    monkeypatch.setattr(scan, "_available_tesseract_langs",
                        lambda: frozenset({"jpn", "jpn_vert", "eng"}))
    out = _annotate("jpn+eng", ["jpn"], vertical=True)

    # The pin still wins — it is documented to beat every inferred signal,
    # and silently overriding an explicit instruction is worse than obeying
    # a bad one.
    assert out["tesseract_packs"] == ["jpn", "eng"]
    # But the conflict is on the record, not merely in a log line that
    # scrolled past during a 40-minute build.
    assert out["ocrlang_overrides_vertical_cjk"] == "jpn_vert"


def test_the_warning_names_the_pack_that_should_have_been_pinned(monkeypatch, caplog):
    import logging
    from pipeline import scan
    monkeypatch.setattr(scan, "_available_tesseract_langs",
                        lambda: frozenset({"jpn", "jpn_vert", "eng"}))
    with caplog.at_level(logging.WARNING, logger="pipeline.scan"):
        _annotate("jpn+eng", ["jpn"], vertical=True)
    text = caplog.text
    assert "jpn_vert" in text
    assert "2x" in text
    # And says why a derived tag gets it wrong, so the reader fixes the
    # generator rather than one bib entry.
    assert "BCP-47" in text


def test_no_warning_when_the_pin_agrees_with_the_geometry(monkeypatch):
    from pipeline import scan
    monkeypatch.setattr(scan, "_available_tesseract_langs",
                        lambda: frozenset({"jpn", "jpn_vert", "eng"}))
    out = _annotate("jpn_vert", ["jpn"], vertical=True)
    assert "ocrlang_overrides_vertical_cjk" not in out
    assert out["tesseract_packs"] == ["jpn_vert"]


def test_no_warning_on_a_horizontally_set_document(monkeypatch):
    """Yamamori2014's shape: Japanese, but set horizontally.

    The geometric probe returns None there, so `jpn+eng` is simply correct
    and must not be second-guessed.
    """
    from pipeline import scan
    monkeypatch.setattr(scan, "_available_tesseract_langs",
                        lambda: frozenset({"jpn", "jpn_vert", "eng"}))
    out = _annotate("jpn+eng", ["jpn"], vertical=False)
    assert "ocrlang_overrides_vertical_cjk" not in out


def test_a_non_cjk_pin_is_untouched(monkeypatch):
    from pipeline import scan
    monkeypatch.setattr(scan, "_available_tesseract_langs",
                        lambda: frozenset({"deu", "deu_latf", "jpn_vert"}))
    out = _annotate("deu_latf+deu", ["deu"], vertical=True)
    assert "ocrlang_overrides_vertical_cjk" not in out
