"""Authority tests use the production schema instead of drifting copies."""
from __future__ import annotations

import re
from pathlib import Path


def test_no_test_redeclares_the_works_table():
    """A schema column addition must reach every fixture automatically (#237)."""
    tests_dir = Path(__file__).parent
    declaration = re.compile(
        "CREATE" + r"\s+TABLE(?:\s+IF\s+NOT\s+EXISTS)?\s+works\b",
        re.IGNORECASE,
    )
    offenders = []
    for path in sorted(tests_dir.glob("test_*.py")):
        if declaration.search(path.read_text(encoding="utf-8")):
            offenders.append(path.name)
    assert offenders == [], (
        "test fixtures must call bib.authority.create_schema; hand-written "
        f"works schemas drifted in five places before #237: {offenders}"
    )

