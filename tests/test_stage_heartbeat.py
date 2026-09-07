"""A long stage says it is still working (#170).

During docling layout analysis the run log went silent for minutes —
measured 3m20s on a 27-page scan and far longer on the 314-page Totton
monograph — and the last line before the gap is a docling banner. A
reader tailing `run.log` could not tell working from hung, and the
README's "a few dozen PDFs run on a laptop in a couple of hours" sets no
expectation for per-document latency.

More pressing since v1.0 re-OCRs scans rather than trusting their text
layers: 31 of 35 papers in the smoke corpus now OCR, where 4 did before.

The beat lives in `_stage`, which is the one place every stage passes
through, so it covers OCR and the vision pass too rather than only
docling.
"""
from __future__ import annotations

import logging
import time

import pytest

from pipeline.runner import _pages_note
from pipeline.stages import _HEARTBEAT_SECONDS, _heartbeat, _stage


def test_a_long_stage_reports_itself(caplog):
    log = logging.getLogger("test.heartbeat")
    with caplog.at_level("INFO"):
        with _heartbeat("docling_extraction", log, interval=0.05):
            time.sleep(0.28)
    beats = [r for r in caplog.records if "still running" in r.getMessage()]
    assert len(beats) >= 3
    assert "docling_extraction still running after" in beats[0].getMessage()


def test_a_fast_stage_stays_silent(caplog):
    """The first beat lands one interval in, so no ordinary stage emits
    one at the 60 s default. A heartbeat on every stage of every document
    would be its own noise problem."""
    log = logging.getLogger("test.heartbeat")
    with caplog.at_level("INFO"):
        with _heartbeat("text_chunking", log, interval=5.0):
            pass
    assert not [r for r in caplog.records if "still running" in r.getMessage()]


def test_zero_disables_it(caplog):
    log = logging.getLogger("test.heartbeat")
    with caplog.at_level("INFO"):
        with _heartbeat("x", log, interval=0):
            time.sleep(0.15)
    assert not [r for r in caplog.records if "still running" in r.getMessage()]


def test_the_elapsed_time_is_readable(caplog):
    log = logging.getLogger("test.heartbeat")
    with caplog.at_level("INFO"):
        with _heartbeat("s", log, interval=0.05):
            time.sleep(0.12)
    msg = [r.getMessage() for r in caplog.records if "still running" in r.getMessage()][0]
    # Never "after 0s" — a sub-second first beat still reads as elapsed time.
    assert "after 0s" not in msg
    assert "after 1s" in msg


def test_an_exception_still_stops_the_beat(caplog):
    """A stage that raises must not leave a thread logging behind it."""
    log = logging.getLogger("test.heartbeat")
    with pytest.raises(ValueError):
        with _heartbeat("s", log, interval=0.05):
            raise ValueError("boom")
    with caplog.at_level("INFO"):
        time.sleep(0.15)
    assert not [r for r in caplog.records if "still running" in r.getMessage()]


def test_the_thread_is_a_daemon():
    """An interpreter exit must never be held up waiting for a beat."""
    log = logging.getLogger("test.heartbeat")
    import threading
    before = {t.name for t in threading.enumerate()}
    with _heartbeat("s", log, interval=30):
        new = [t for t in threading.enumerate() if t.name not in before]
        assert new and all(t.daemon for t in new)


# --- what the line says -------------------------------------------------


def test_the_page_count_explains_the_wait():
    """Page count *is* the answer to "why is this slow", and it is
    recorded by huge_document_check, which runs first."""
    assert _pages_note({"page_count": 314}) == " (314 pages)"
    assert _pages_note({"page_count": 1}) == " (1 pages)"


def test_an_unknown_page_count_says_nothing():
    assert _pages_note({}) == ""
    assert _pages_note({"page_count": None}) == ""
    assert _pages_note({"page_count": 0}) == ""


# --- wiring -------------------------------------------------------------


def test_stage_emits_the_beat_with_the_paper_logger(caplog):
    """`_stage` is the integration point, and the per-paper adapter is
    what makes a line attributable in an interleaved multi-paper stream."""
    from pipeline.runner import _PaperLogAdapter

    plog = _PaperLogAdapter(logging.getLogger("test.stage"),
                            {"pdf": "Totton1965"})
    summary = {}
    with caplog.at_level("INFO"):
        with _stage(summary, "docling_extraction", logger_=plog,
                    heartbeat_detail=" (314 pages)"):
            time.sleep(0.05)
    # The default interval is 60 s, so nothing fired — the wiring is what
    # is asserted here, and the interval behaviour above.
    assert summary["stage_timings"][0]["stage"] == "docling_extraction"


def test_every_stage_call_site_passes_the_paper_logger():
    """Miss one and that stage's silence is the only one left unexplained."""
    import inspect
    import re

    from pipeline import runner
    src = inspect.getsource(runner)
    calls = re.findall(r"with _stage\(processing_summary, \"([a-z_0-9]+)\","
                       r"(?: logger_=plog,)?", src)
    without = re.findall(r"with _stage\(processing_summary, \"([a-z_0-9]+)\","
                         r"(?! logger_=plog,)", src)
    assert calls, "no _stage call sites found; this guard would pass vacuously"
    assert not without, f"stages without the paper logger: {without}"


def test_the_interval_is_configurable():
    from pipeline.config_schema import LoggingConfig

    assert LoggingConfig().heartbeat_seconds == _HEARTBEAT_SECONDS
    assert LoggingConfig(heartbeat_seconds=0).heartbeat_seconds == 0
    with pytest.raises(Exception):
        LoggingConfig(heartbeat_seconds=-1)
