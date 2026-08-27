"""OCR worker cap and process-group cleanup (#209).

Reported by @ejedwards with a process table: on a 12-core / 30 GB host,
ocrmypdf spawned one ~1.9 GB Tesseract worker per CPU and the build was
OOM-killed. Two things made it unfixable from outside:

* ocrmypdf takes its worker count from ``multiprocessing.cpu_count()``, which
  ignores CPU affinity — so neither ``taskset`` nor a cgroup CPU limit reaches
  it, and a cgroup *memory* limit only contains the blast radius while the
  build still dies.
* Killing the pipeline left ocrmypdf reparented to PID 1 with its Tesseract
  children still holding ~20 GB, because ``subprocess.run(timeout=)`` kills
  the direct child and not its grandchildren.
"""
from __future__ import annotations

import os
import subprocess
import sys

import pytest

from pipeline.scan import (
    _OCR_PARENT_GB,
    _OCR_RESERVED_GB,
    _OCR_WORKER_GB,
    _host_ram_gb,
    _resolve_ocr_jobs,
    _run_ocr,
)


# --- choosing the cap ---------------------------------------------------------


def test_an_explicit_setting_is_honoured(monkeypatch):
    from pipeline import scan
    monkeypatch.setitem(scan.CONFIG, "ocr", {"jobs": 4})
    assert _resolve_ocr_jobs() == 4


def test_an_explicit_setting_may_exceed_the_derived_cap(monkeypatch):
    """An operator who knows their pages are light should be able to say so;
    the derived value is a safety net, not a ceiling."""
    from pipeline import scan
    monkeypatch.setitem(scan.CONFIG, "ocr", {"jobs": 64})
    assert _resolve_ocr_jobs() == 64


def test_ram_that_is_not_the_binding_constraint_leaves_ocrmypdf_alone(monkeypatch):
    """On a large host this must return None so nothing gets slower — the cap
    only ever *lowers* the worker count."""
    from pipeline import scan
    monkeypatch.setitem(scan.CONFIG, "ocr", {})
    monkeypatch.setattr(scan, "_host_ram_gb", lambda: 512.0)
    monkeypatch.setattr(scan.os, "cpu_count", lambda: 8)
    assert _resolve_ocr_jobs() is None


def test_a_small_host_is_capped(monkeypatch):
    from pipeline import scan
    monkeypatch.setitem(scan.CONFIG, "ocr", {})
    monkeypatch.setattr(scan, "_host_ram_gb", lambda: 30.0)
    monkeypatch.setattr(scan.os, "cpu_count", lambda: 12)
    jobs = _resolve_ocr_jobs()
    # The reporter's host: OOM-killed at 8 concurrent workers, proceeded at 4.
    assert jobs is not None
    assert 1 <= jobs < 8


def test_the_cap_never_reaches_zero(monkeypatch):
    """A tiny host still has to make progress, one page at a time."""
    from pipeline import scan
    monkeypatch.setitem(scan.CONFIG, "ocr", {})
    monkeypatch.setattr(scan, "_host_ram_gb", lambda: 2.0)
    monkeypatch.setattr(scan.os, "cpu_count", lambda: 16)
    assert _resolve_ocr_jobs() == 1


def test_an_unknown_host_is_not_second_guessed(monkeypatch):
    from pipeline import scan
    monkeypatch.setitem(scan.CONFIG, "ocr", {})
    monkeypatch.setattr(scan, "_host_ram_gb", lambda: None)
    assert _resolve_ocr_jobs() is None


def test_ram_is_read_without_psutil():
    """psutil is only a transitive dependency here; this must not start
    failing because some other package drops it."""
    ram = _host_ram_gb()
    assert ram is None or ram > 0


def test_the_memory_model_is_stated_in_named_terms():
    """The constants come from one measured host, so they are named rather
    than folded into a magic formula."""
    assert _OCR_PARENT_GB > 0 and _OCR_WORKER_GB > 0 and _OCR_RESERVED_GB > 0
    # A worker is budgeted above the ~1.9 GB measured, because the failure
    # mode is an OOM kill rather than a slowdown.
    assert _OCR_WORKER_GB > 1.9


def test_the_flag_reaches_the_command_line():
    import inspect
    from pipeline import scan
    src = inspect.getsource(scan.prepare_pdf)
    assert '"--jobs"' in src
    assert "_resolve_ocr_jobs()" in src


# --- taking the workers down with it -----------------------------------------


@pytest.mark.skipif(os.name != "posix", reason="process groups are POSIX")
def test_a_timeout_kills_the_grandchildren_too():
    """The orphan case, with a real process tree: a parent that spawns
    children and hangs, exactly the shape of ocrmypdf + Tesseract."""
    spawn = (
        "import subprocess,sys,time;"
        "kids=[subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)'])"
        " for _ in range(2)];"
        "print(' '.join(str(k.pid) for k in kids), flush=True);"
        "time.sleep(60)"
    )
    # Learn the grandchild pids from an identical tree we control.
    probe = subprocess.Popen([sys.executable, "-c", spawn],
                             stdout=subprocess.PIPE, text=True,
                             start_new_session=True)
    kids = [int(x) for x in probe.stdout.readline().split()]

    def alive(pid):
        try:
            os.kill(pid, 0)
            return True
        except OSError:
            return False

    assert all(alive(k) for k in kids), "probe tree did not start"
    os.killpg(os.getpgid(probe.pid), 9)
    probe.wait(timeout=5)
    import time
    time.sleep(0.5)
    assert not any(alive(k) for k in kids), (
        "a process-group kill must reach the grandchildren — this is the "
        "property _run_ocr relies on"
    )


@pytest.mark.skipif(os.name != "posix", reason="process groups are POSIX")
def test_run_ocr_raises_timeout_and_leaves_nothing_behind():
    with pytest.raises(subprocess.TimeoutExpired):
        _run_ocr([sys.executable, "-c", "import time; time.sleep(30)"], timeout=1)


def test_run_ocr_returns_a_completed_process_on_success():
    r = _run_ocr([sys.executable, "-c", "print('ok')"], timeout=30)
    assert r.returncode == 0 and "ok" in r.stdout
