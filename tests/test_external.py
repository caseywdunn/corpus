"""Unit tests for the shared external-service primitives (#37)."""
from __future__ import annotations

import time

import pytest
import requests

import external
from external import (
    CircuitBreaker,
    CircuitOpenError,
    is_transient,
    retry_with_backoff,
    set_strict_network,
)


# ---------------------------------------------------------------------------
# is_transient
# ---------------------------------------------------------------------------


def test_is_transient_for_connection_and_timeout_errors():
    assert is_transient(requests.ConnectionError("x"))
    assert is_transient(requests.Timeout("x"))


@pytest.mark.parametrize("status_code, transient", [
    (408, True),
    (429, True),
    (500, True),
    (502, True),
    (503, True),
    (504, True),
    (400, False),
    (401, False),
    (403, False),
    (404, False),
    (200, False),  # would never reach here in practice
])
def test_is_transient_classifies_http_status(status_code, transient):
    class _Resp:
        pass
    resp = _Resp()
    resp.status_code = status_code
    err = requests.HTTPError("x")
    err.response = resp
    assert is_transient(err) is transient


def test_is_transient_false_for_non_network():
    assert not is_transient(ValueError("nope"))
    assert not is_transient(RuntimeError("nope"))


# ---------------------------------------------------------------------------
# retry_with_backoff
# ---------------------------------------------------------------------------


def test_retry_succeeds_after_transient_failures():
    calls = []

    def flaky():
        calls.append(None)
        if len(calls) < 3:
            raise requests.ConnectionError("transient")
        return "ok"

    result = retry_with_backoff(flaky, max_attempts=5, sleep=lambda _: None)
    assert result == "ok"
    assert len(calls) == 3


def test_retry_short_circuits_on_permanent_error():
    calls = []

    def perm():
        calls.append(None)
        raise ValueError("permanent — not network")

    with pytest.raises(ValueError):
        retry_with_backoff(perm, max_attempts=5, sleep=lambda _: None)
    assert len(calls) == 1


def test_retry_exhausts_attempts_on_persistent_transient():
    calls = []

    def always_flaky():
        calls.append(None)
        raise requests.Timeout("never recovers")

    with pytest.raises(requests.Timeout):
        retry_with_backoff(always_flaky, max_attempts=4, sleep=lambda _: None)
    assert len(calls) == 4


def test_retry_respects_strict_network():
    calls = []

    def flaky():
        calls.append(None)
        raise requests.ConnectionError("transient")

    set_strict_network(True)
    try:
        with pytest.raises(requests.ConnectionError):
            retry_with_backoff(flaky, max_attempts=5, sleep=lambda _: None)
        assert len(calls) == 1, "STRICT_NETWORK should fail on first attempt"
    finally:
        set_strict_network(False)


# ---------------------------------------------------------------------------
# CircuitBreaker
# ---------------------------------------------------------------------------


def test_breaker_opens_after_threshold_failures():
    b = CircuitBreaker("svc", threshold=3, cooldown_s=10)
    assert not b.is_open()
    b.record_failure()
    b.record_failure()
    assert not b.is_open()
    b.record_failure()
    assert b.is_open()
    with pytest.raises(CircuitOpenError):
        b.check()


def test_breaker_resets_after_success():
    b = CircuitBreaker("svc", threshold=3, cooldown_s=10)
    b.record_failure()
    b.record_failure()
    b.record_success()
    # Counter reset
    b.record_failure()
    b.record_failure()
    assert not b.is_open()  # still under threshold


def test_breaker_reopens_after_cooldown():
    b = CircuitBreaker("svc", threshold=2, cooldown_s=0.05)
    b.record_failure()
    b.record_failure()
    assert b.is_open()
    time.sleep(0.1)
    assert not b.is_open()


def test_breaker_check_passes_when_closed():
    b = CircuitBreaker("svc", threshold=10, cooldown_s=10)
    b.check()  # no raise
