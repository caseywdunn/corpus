"""Backwards-compat shim — `bib_metadata` moved into the `bib` package (#43).

New code should import from ``bib`` (or ``bib.parser`` for private helpers):

    from bib import BibIndex, parse_bibtex, bib_entry_to_metadata

This shim re-exports the public surface and keeps existing
``from bib_metadata import ...`` lines working during the
deprecation window. Slated for removal once no caller depends on it.
"""
from bib.parser import (  # noqa: F401  (re-export)
    BibIndex,
    bib_entry_to_metadata,
    parse_bibtex,
    _split_authors,
    _strip_outer_braces,
)
