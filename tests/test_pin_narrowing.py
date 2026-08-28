"""A pin discards what detection resolved, and that should be visible (#245).

Pinning `ocrlang` on 31 of 35 gold documents moved corpus-wide prose coverage
0.9474 -> 0.9450 with every language correct. The loss came entirely from each
pin replacing a *union* detection had chosen with a single pack:

    Linnaeus1735   eng+deu+deu_latf+fra+lat+spa+por -> lat    -0.079
    Tilesius1814   swe+cat+fra+eng                  -> swe    -0.063
    Vanhoeffen1906 deu+deu_latf+eng                 -> deu    +0.019

Narrowing is not uniformly bad — two documents improved — so the pin still
wins. What was wrong is that a directive costing 0.05 looked identical to one
gaining it. See dev_docs/OCR_LANGUAGES.md §5.
"""
from __future__ import annotations

import logging

from pipeline import scan


def _annotate(ocrlang, detected, monkeypatch, available=None):
    """Annotate a result whose *targeted* resolution is ``detected``.

    Patching `_resolve_tesseract_packs` rather than seeding
    `tesseract_packs` directly, because that key is recomputed inside the
    function — and the distinction matters: `tesseract_packs` folds in the
    configured fallback union when detection found nothing, while targeted
    resolution is what the document's own evidence supports.
    """
    monkeypatch.setattr(scan, "_available_tesseract_langs",
                        lambda: frozenset(available or
                                          {"lat", "eng", "deu", "deu_latf",
                                           "fra", "swe", "por"}))
    monkeypatch.setattr(scan, "_resolve_tesseract_packs",
                        lambda iso, visual_script=None: list(detected))
    result = {"filename": "X.pdf", "needs_ocr": True,
              "detected_language": "xx"}
    return scan._annotate_pack_availability(ocrlang, result)


def test_pinning_over_a_bare_fallback_is_not_narrowing(monkeypatch):
    """Detection found nothing, so there is nothing to narrow.

    `tesseract_packs` still holds the configured default union in that
    case, and reporting it as displaced would fire on exactly the documents
    `ocrlang` exists for.
    """
    out = _annotate("por", [], monkeypatch)
    assert "ocrlang_narrowed_from" not in out


def test_a_narrowing_pin_records_what_it_displaced(monkeypatch):
    """Linnaeus1735's shape: seven packs resolved, one pinned."""
    out = _annotate("lat", ["eng", "deu", "deu_latf", "fra", "lat", "por"],
                    monkeypatch)
    assert out["tesseract_packs"] == ["lat"]          # the pin still wins
    assert out["ocrlang_narrowed_from"] == ["eng", "deu", "deu_latf", "fra", "por"]


def test_a_pin_that_adds_packs_is_not_narrowing(monkeypatch):
    """`deu` detected, `deu+deu_latf` pinned — nothing was discarded."""
    out = _annotate("deu+deu_latf", ["deu"], monkeypatch)
    assert out["tesseract_packs"] == ["deu", "deu_latf"]
    assert "ocrlang_narrowed_from" not in out


def test_a_pin_matching_detection_is_silent(monkeypatch):
    out = _annotate("por", ["por"], monkeypatch)
    assert "ocrlang_narrowed_from" not in out


def test_a_pin_that_swaps_records_the_displaced_pack(monkeypatch):
    """Quoy_Gaimard1834Plates' shape: `eng` detected, `fra` pinned.

    Not a narrowing by count, but detection's pack is still discarded and
    that is the thing worth seeing — it scored 0.062 worse.
    """
    out = _annotate("fra", ["eng"], monkeypatch)
    assert out["ocrlang_narrowed_from"] == ["eng"]


def test_the_log_line_names_both_sides_and_points_at_the_evidence(monkeypatch, caplog):
    with caplog.at_level(logging.INFO, logger="pipeline.scan"):
        _annotate("lat", ["eng", "lat", "fra"], monkeypatch)
    text = caplog.text
    assert "eng+lat+fra" in text          # what detection resolved
    assert "eng+fra" in text              # what the pin drops
    assert "OCR_LANGUAGES" in text        # where the measurement lives


def test_an_unhonoured_pin_records_nothing(monkeypatch):
    """A pin naming no installed pack falls back to detection untouched."""
    out = _annotate("xxx", ["eng", "fra"], monkeypatch, available={"eng", "fra"})
    assert out["tesseract_packs"] == ["eng", "fra"]
    assert "ocrlang_narrowed_from" not in out


def test_no_pin_records_nothing(monkeypatch):
    out = _annotate(None, ["eng", "fra"], monkeypatch)
    assert "ocrlang_narrowed_from" not in out
