"""Tests for the shared rich console layer (#63).

Asserts the ASCII-fallback path (non-TTY) — which is what every
SLURM ``.out`` file, CI log, and ``journalctl`` reader sees. The
TTY/rich path is harder to assert in pytest (escape sequences,
terminal width detection); a separate ``demo.py`` exercise is
sufficient for visual verification of that path.
"""
from __future__ import annotations

import io
from contextlib import redirect_stdout

from pipeline import console as console_mod
from pipeline.console import print_status, progress


def test_console_is_not_terminal_under_pytest():
    # Pytest runs without a controlling TTY, so the ASCII path must engage.
    assert console_mod.console.is_terminal is False


def test_print_status_ascii_fallback_uses_bracket_tags():
    buf = io.StringIO()
    with redirect_stdout(buf):
        print_status("packaging done", status="ok")
        print_status("vision skipped", status="warn")
        print_status("config error", status="fail")
        print_status("starting embed batch", status="info")
    out = buf.getvalue()
    # Plain ASCII tags survive — rich markup is bypassed off-TTY.
    assert "[ok] packaging done" in out
    assert "[warn] vision skipped" in out
    assert "[fail] config error" in out
    assert "[info] starting embed batch" in out


def test_progress_disabled_off_tty_completes_cleanly():
    # Disabled bars don't render but the context manager still works
    # so call sites stay uniform.
    with progress("hashing pdfs") as p:
        assert p.disable is True
    with progress("papers", total=5) as p:
        for tid in p.task_ids:
            for _ in range(5):
                p.update(tid, advance=1)
        assert p.disable is True
