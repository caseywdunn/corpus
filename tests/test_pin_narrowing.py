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


def test_pinning_over_a_fallback_union_is_still_narrowing(monkeypatch):
    """Detection had no targeted resolution, but packs were still displaced.

    `Linnaeus1735`'s shape, and the case the first version of this missed:
    seven packs would have been used, `lat` was pinned, and the document lost
    0.079. Reporting only targeted displacement flagged 4 of the 22 documents
    whose pack list the pin changed, and missed the largest regression in the
    set. What was displaced is the fact worth acting on; whether detection
    was confident about it is a qualifier, kept separately.
    """
    out = _annotate("lat", [], monkeypatch)
    assert out["ocrlang_narrowed_from"]                  # the fallback union
    assert "ocrlang_narrowed_from_targeted" not in out   # none of it targeted


def test_a_narrowing_pin_records_what_it_displaced(monkeypatch):
    """Linnaeus1735's shape: seven packs resolved, one pinned."""
    out = _annotate("lat", ["eng", "deu", "deu_latf", "fra", "lat", "por"],
                    monkeypatch)
    assert out["tesseract_packs"] == ["lat"]          # the pin still wins
    assert set(out["ocrlang_narrowed_from"]) >= {"eng", "deu", "deu_latf", "fra", "por"}
    # Detection resolved these from the document's own evidence, so they are
    # a stronger signal than packs the configured fallback contributed.
    assert set(out["ocrlang_narrowed_from_targeted"]) == {"eng", "deu", "deu_latf",
                                                          "fra", "por"}


def test_a_pin_naming_everything_that_would_have_run_is_silent(monkeypatch):
    """Nothing displaced, nothing reported."""
    out = _annotate("por+eng", ["por"], monkeypatch)
    assert out["tesseract_packs"] == ["por", "eng"]
    assert "ocrlang_narrowed_from" not in out


def test_dropping_the_automatic_eng_suffix_counts_as_narrowing(monkeypatch):
    """`Candeias1932`'s shape, and worth 0.048.

    Detection targets `por`, composition adds `eng`, and a pin of `por`
    alone silently drops it. That the dropped pack came from composition
    rather than from targeted resolution does not make it cost less — this
    document lost 0.048 — so it is reported, with the provenance kept as a
    separate qualifier.
    """
    out = _annotate("por", ["por"], monkeypatch)
    assert out["tesseract_packs"] == ["por"]
    assert out["ocrlang_narrowed_from"] == ["eng"]
    assert "ocrlang_narrowed_from_targeted" not in out


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
    assert "eng+lat+fra" in text          # what OCR would have used
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
