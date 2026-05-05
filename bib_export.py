#!/usr/bin/env python3
"""CLI shim for ``bib.export`` (#43).

The implementation moved to ``bib/export.py`` along with the rest of
the bib round-trip code. This file stays at the repo root so existing
``python bib_export.py …`` invocations keep working unchanged.
"""
import sys

from bib.export import main

if __name__ == "__main__":
    sys.exit(main())
