"""Logging setup for the pipeline.

Two helpers:

* :func:`setup_root_logging` — single stderr stream handler on the
  root logger, idempotent across re-invocations.
* :func:`per_pdf_file_log` — context manager that attaches a
  ``FileHandler`` writing to ``<hash_dir>/pipeline.log`` for the
  duration of one document's processing. Output ends up in both the
  per-paper log file and the root stream handler.
"""
from __future__ import annotations

import logging
from contextlib import contextmanager
from pathlib import Path


def setup_root_logging(level: int = logging.INFO) -> None:
    """Configure the root logger with a single stderr stream handler.

    Per-PDF file handlers are added/removed around each document by
    ``per_pdf_file_log``; they coexist with this stream handler so output
    appears both on the terminal and in the per-paper pipeline.log.
    """
    root = logging.getLogger()
    root.setLevel(level)
    # Avoid duplicate stream handlers if called twice.
    for h in root.handlers:
        if isinstance(h, logging.StreamHandler) and getattr(h, "_corpus_stream", False):
            return
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter("%(levelname)s %(name)s: %(message)s"))
    sh._corpus_stream = True  # marker so we don't re-add on repeated calls
    root.addHandler(sh)


@contextmanager
def per_pdf_file_log(hash_dir: Path):
    """Attach a FileHandler writing to ``<hash_dir>/pipeline.log`` for the
    duration of the context. All ``logging`` calls made inside the block are
    captured to that file in addition to the root stream handler.
    """
    hash_dir.mkdir(parents=True, exist_ok=True)
    log_path = hash_dir / "pipeline.log"
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
    )
    root = logging.getLogger()
    root.addHandler(fh)
    try:
        yield log_path
    finally:
        root.removeHandler(fh)
        fh.close()
