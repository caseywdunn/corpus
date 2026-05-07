"""Shared primitives for talking to flaky external services (#37).

The pipeline depends on a small set of network services — Grobid (local
or HPC), BHL, and the future CrossRef / OpenAlex paths. v0.1 had bespoke
retry loops on some calls and none on others; a single network blip
mid-batch could silently degrade the corpus output. This module
centralizes:

* :func:`retry_with_backoff` — exponential backoff with jitter; classifies
  transient vs. permanent failures via :func:`is_transient`.
* :class:`CircuitBreaker` — after ``threshold`` consecutive failures,
  short-circuits subsequent calls until ``cooldown_s`` elapses. Stops a
  down service from soaking up retry budget on every paper in an array
  job.
* :data:`STRICT_NETWORK` — process-wide toggle. When True, the retry
  helper raises immediately on the first transient failure instead of
  burning attempts. Use for release-build runs where silent partial
  data is worse than aborting.

Usage:

    from pipeline.external import retry_with_backoff, CircuitBreaker, is_transient

    breaker = CircuitBreaker("grobid", threshold=10, cooldown_s=300)

    def do():
        breaker.check()
        try:
            r = requests.post(url, ..., timeout=30)
            r.raise_for_status()
            breaker.record_success()
            return r
        except Exception as e:
            if is_transient(e):
                breaker.record_failure()
            raise

    return retry_with_backoff(do, max_attempts=3)

Future services (CrossRef, OpenAlex) should grow their own breaker
instances and use the same retry helper. Don't pre-build harnesses for
services that aren't yet wired into the pipeline.
"""
from __future__ import annotations

import logging
import os
import random
import time
from typing import Any, Callable, Optional, Tuple, Type

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# STRICT_NETWORK
# ---------------------------------------------------------------------------
# Process-wide toggle. Set via environment (CORPUS_STRICT_NETWORK=1) or
# by callers (e.g. process_corpus.py main() sets this from --strict-network).
# Off by default — production batch runs prefer "skip the bad paper, keep
# going" semantics.
STRICT_NETWORK: bool = os.environ.get("CORPUS_STRICT_NETWORK", "").lower() in ("1", "true", "yes")


def set_strict_network(value: bool) -> None:
    """Public setter so callers don't import-and-assign a module global."""
    global STRICT_NETWORK
    STRICT_NETWORK = bool(value)


# ---------------------------------------------------------------------------
# Transient-failure classification
# ---------------------------------------------------------------------------


def is_transient(exc: BaseException) -> bool:
    """True if ``exc`` looks worth retrying.

    Conservative: connect errors, read timeouts, and 5xx / 429 from
    HTTPError. Everything else (4xx other than 429, JSONDecodeError, etc.)
    is treated as permanent so we don't burn retry budget on a programming
    error or a genuine server-side rejection.
    """
    try:
        import requests
    except ImportError:
        return False
    if isinstance(exc, (requests.ConnectionError, requests.Timeout)):
        return True
    if isinstance(exc, requests.HTTPError):
        resp = getattr(exc, "response", None)
        if resp is not None:
            return resp.status_code in (408, 429, 500, 502, 503, 504)
    return False


# ---------------------------------------------------------------------------
# Retry with exponential backoff
# ---------------------------------------------------------------------------


def retry_with_backoff(
    fn: Callable[[], Any],
    *,
    max_attempts: int = 3,
    base_delay: float = 1.0,
    factor: float = 2.0,
    jitter: bool = True,
    transient_predicate: Callable[[BaseException], bool] = is_transient,
    sleep: Callable[[float], None] = time.sleep,
) -> Any:
    """Call ``fn()``; retry transient failures with exponential backoff.

    On a non-transient exception or exhaustion of attempts, the original
    exception is re-raised so the caller sees the underlying type
    (e.g. ``GrobidUnavailableError``). Permanent errors fail immediately
    without burning attempts.

    When :data:`STRICT_NETWORK` is set, the first transient failure raises
    immediately — useful for release-build runs that prefer aborting over
    silently degrading data.

    Args:
        fn: nullary callable returning a result
        max_attempts: total attempts including the first
        base_delay: seconds before retry 1; doubles each attempt by default
        factor: multiplicative growth (use 1.0 for fixed delay)
        jitter: scale each delay by uniform(0.5, 1.5) to avoid thundering herd
        transient_predicate: classifier; defaults to :func:`is_transient`
        sleep: injectable sleep (tests pass a no-op)
    """
    last_exc: Optional[BaseException] = None
    for attempt in range(1, max_attempts + 1):
        try:
            return fn()
        except Exception as e:
            last_exc = e
            if not transient_predicate(e):
                raise
            if STRICT_NETWORK:
                logger.warning("STRICT_NETWORK: aborting on transient %s", type(e).__name__)
                raise
            if attempt >= max_attempts:
                raise
            delay = base_delay * (factor ** (attempt - 1))
            if jitter:
                delay *= 0.5 + random.random()
            logger.warning(
                "Transient %s (attempt %d/%d); retrying in %.1fs",
                type(e).__name__, attempt, max_attempts, delay,
            )
            sleep(delay)
    # Unreachable in practice — defensive
    if last_exc is not None:
        raise last_exc
    return None


# ---------------------------------------------------------------------------
# CircuitBreaker
# ---------------------------------------------------------------------------


class CircuitOpenError(Exception):
    """Raised when a CircuitBreaker is open (too many recent failures)."""


class CircuitBreaker:
    """Per-service breaker. Opens after ``threshold`` consecutive failures
    and stays open for ``cooldown_s``. Reset on first success after cooldown.

    State is in-process — appropriate for a single CLI run. For long-lived
    services (the MCP server) consider per-request reset semantics.
    """

    def __init__(self, name: str, threshold: int = 10, cooldown_s: float = 300.0):
        self.name = name
        self.threshold = threshold
        self.cooldown_s = cooldown_s
        self._failures = 0
        self._opened_at: Optional[float] = None

    def is_open(self) -> bool:
        if self._opened_at is None:
            return False
        elapsed = time.monotonic() - self._opened_at
        if elapsed >= self.cooldown_s:
            # Half-open: allow next call through to probe the service.
            return False
        return True

    def check(self) -> None:
        """Raise :class:`CircuitOpenError` if the breaker is currently open."""
        if self.is_open():
            raise CircuitOpenError(
                f"{self.name} circuit open ({self._failures} consecutive failures; "
                f"retry after cooldown)"
            )

    def record_success(self) -> None:
        if self._failures or self._opened_at:
            logger.info("%s circuit reset (success after %d failures)", self.name, self._failures)
        self._failures = 0
        self._opened_at = None

    def record_failure(self) -> None:
        self._failures += 1
        if self._failures >= self.threshold and self._opened_at is None:
            logger.warning(
                "%s circuit opening: %d consecutive failures, cooldown %.0fs",
                self.name, self._failures, self.cooldown_s,
            )
            self._opened_at = time.monotonic()


__all__ = [
    "STRICT_NETWORK",
    "set_strict_network",
    "is_transient",
    "retry_with_backoff",
    "CircuitBreaker",
    "CircuitOpenError",
]
