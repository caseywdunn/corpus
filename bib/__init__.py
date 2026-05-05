"""Bibliographic round-trip + parser package (#43).

Three sibling modules that already share state, now grouped:

  bib.parser    — BibTeX parser + BibIndex (formerly ``bib_metadata.py``).
  bib.export    — biblio_authority.sqlite → BibTeX (formerly ``bib_export.py``).
  bib.importer  — hand-edited BibTeX → biblio_authority.sqlite
                  (formerly ``bib_import.py``; the file is named ``importer``
                  rather than ``import`` to avoid shadowing the keyword).

Top-level CLIs (``bib_export.py``, ``bib_import.py``) stay at the repo
root as thin ``__main__`` shims that delegate here. The legacy
``bib_metadata`` module is kept as a re-export shim during a
deprecation window — rip out when no caller imports it directly.
"""
from .parser import (
    BibIndex,
    bib_entry_to_metadata,
    parse_bibtex,
)
from .export import export_bibtex
from .importer import import_bibtex

__all__ = [
    "BibIndex",
    "bib_entry_to_metadata",
    "parse_bibtex",
    "export_bibtex",
    "import_bibtex",
]
