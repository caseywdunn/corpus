"""Tests for the macOS arm64 Python arch gate in `corpus check` (#77).

The check is implemented as ``pipeline.cli._check_python_arch``, a
small helper extracted from ``_cmd_check`` so the contract is
unit-testable without standing up the full config + grobid + disk
environment ``_cmd_check`` otherwise requires.

The helper returns ``(arch_label, failure_message)`` per its
docstring; these tests pin the three branches: not darwin (skip),
darwin/arm64 (ok), darwin/anything-else (fail loud).
"""
from __future__ import annotations

from unittest import mock

from pipeline.cli import _check_python_arch


def test_skip_on_linux():
    with mock.patch("pipeline.cli.sys.platform", "linux"):
        arch, failure = _check_python_arch()
    assert arch == ""
    assert failure is None


def test_skip_on_windows():
    with mock.patch("pipeline.cli.sys.platform", "win32"):
        arch, failure = _check_python_arch()
    assert arch == ""
    assert failure is None


def test_pass_on_darwin_arm64():
    with mock.patch("pipeline.cli.sys.platform", "darwin"), \
         mock.patch("platform.machine", return_value="arm64"):
        arch, failure = _check_python_arch()
    assert arch == "arm64"
    assert failure is None


def test_fail_on_darwin_x86_64():
    with mock.patch("pipeline.cli.sys.platform", "darwin"), \
         mock.patch("platform.machine", return_value="x86_64"):
        arch, failure = _check_python_arch()
    assert arch == "x86_64"
    assert failure is not None
    assert "Rosetta" in failure
    assert "arm64-native" in failure
    assert "README" in failure


def test_fail_message_includes_actual_arch():
    """If the arch is some unexpected non-arm64 string (e.g., a future
    Apple ISA), the failure message still surfaces it so the operator
    can see what they got."""
    with mock.patch("pipeline.cli.sys.platform", "darwin"), \
         mock.patch("platform.machine", return_value="ppc64"):
        arch, failure = _check_python_arch()
    assert arch == "ppc64"
    assert failure is not None
    assert "ppc64" in failure
