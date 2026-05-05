#!/usr/bin/env python3
"""CLI shim for ``bib.importer`` (#43).

The implementation moved to ``bib/importer.py``. This file stays at
the repo root so existing ``python bib_import.py …`` invocations keep
working unchanged. The internal module is named ``importer`` rather
than ``import`` to avoid shadowing the keyword.
"""
import sys

from bib.importer import main

if __name__ == "__main__":
    sys.exit(main())
