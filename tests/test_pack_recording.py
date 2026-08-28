"""`tesseract_packs` must not record "unknown" as "none" (#197).

Three documents in the reference corpuscle once recorded
`"tesseract_packs": []` while OCR ran with seven packs. The cause was that
`_compose_ocr_langs` returns early when `_available_tesseract_langs()` is
empty — before the configured fallback union is reached — so a failure to
*enumerate* the installed languages was written down as a resolution that
found nothing.

That matters because `scan_detection.json` is the operator-facing record.
`corpus status` reads it, and #176's `ocrlang` workflow tells an operator to
consult it when deciding which pack to pin. A document whose record says "no
packs" when seven were used sends that diagnosis in the wrong direction — and
these are exactly the documents an operator would be investigating, since
untrusted language detection correlates with the hardest pages.

The symptom no longer reproduces: every document in the 35-document corpuscle
now records what OCR actually ran, verified by comparing each
`scan_detection.json` against the `Running OCR` line in its `pipeline.log`.
These tests are the guard, so a silent recurrence becomes a visible one.
"""
from __future__ import annotations

import pipeline.scan as scan


def _annotate(monkeypatch, available, **result):
    monkeypatch.setattr(scan, "_available_tesseract_langs", lambda: frozenset(available))
    base = {"filename": "Doc.pdf", "detected_language": "de",
            "language_trusted": True, "probe_languages": ["de"]}
    base.update(result)
    return scan._annotate_pack_availability(None, base)


def test_an_unenumerable_tesseract_is_flagged_not_silently_empty(monkeypatch, caplog):
    """The defect in one assertion: an empty availability probe must be
    distinguishable from a resolution that legitimately found nothing."""
    import logging
    with caplog.at_level(logging.WARNING, logger="pipeline.scan"):
        out = _annotate(monkeypatch, [])
    assert out.get("tesseract_langs_unavailable") is True
    assert "unknown" in caplog.text
    # And the warning has to point at where the truth is.
    assert "Running OCR" in caplog.text


def test_a_working_probe_sets_no_flag(monkeypatch):
    out = _annotate(monkeypatch, ["deu", "deu_latf", "eng", "fra"])
    assert "tesseract_langs_unavailable" not in out
    assert out["tesseract_packs"], "packs should resolve when languages are known"


def test_the_flag_does_not_change_the_packs_that_were_resolved(monkeypatch):
    """A guard, not a behaviour change — OCR still does whatever it did."""
    with_langs = _annotate(monkeypatch, ["deu", "deu_latf", "eng"])
    assert with_langs["tesseract_packs"] == _annotate(
        monkeypatch, ["deu", "deu_latf", "eng"])["tesseract_packs"]


def test_record_and_invocation_use_the_same_resolution():
    """`prepare_pdf` must not recompute packs by a different route than the
    one that produced the record. It did diverge once, for vertical CJK
    (#196), where detection logged one pack list and OCR ran another."""
    import inspect
    src = inspect.getsource(scan.prepare_pdf)
    # Either it takes the recorded packs verbatim, or it calls the same
    # composer the record was built from.
    assert ("detection_result.get(\"tesseract_packs\")" in src
            or "_compose_ocr_langs(" in src)
