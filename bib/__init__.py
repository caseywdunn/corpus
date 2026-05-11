"""Bibliographic round-trip + parser package (#43).

Three sibling modules that already share state, now grouped:

  bib.parser    — BibTeX parser + BibIndex.
  bib.export    — biblio_authority.sqlite → BibTeX (formerly ``bib_export.py``).
  bib.importer  — hand-edited BibTeX → biblio_authority.sqlite
                  (formerly ``bib_import.py``; the file is named ``importer``
                  rather than ``import`` to avoid shadowing the keyword).

Top-level CLIs (``bib_export.py``, ``bib_import.py``) stay at the repo
root as thin ``__main__`` shims that delegate here.
"""
# Honor CORPUS_LOG_LEVEL env var before any submodule's basicConfig
# runs (see pipeline/__init__.py for the rationale).
import logging as _logging  # noqa: E402
import os as _os  # noqa: E402
_log_level = _os.environ.get("CORPUS_LOG_LEVEL", "").upper()
if _log_level in {"WARNING", "INFO", "DEBUG"}:
    _logging.basicConfig(
        level=getattr(_logging, _log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

from .parser import (  # noqa: E402
    BibIndex,
    bib_entry_to_metadata,
    parse_bibtex,
)
from .export import export_bibtex  # noqa: E402
from .importer import import_bibtex  # noqa: E402

__all__ = [
    "BibIndex",
    "bib_entry_to_metadata",
    "parse_bibtex",
    "export_bibtex",
    "import_bibtex",
]
