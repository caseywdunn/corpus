"""Dependencies point one way: `pipeline/` never imports `tools/` or `skills/`.

A *library* is the upstream collection of PDFs and a `.bib`; a *corpuscle* is
what corpus builds from one. Code that helps curate a library pulls toward
living in the library repo that needs it — and then a second library carries a
second copy, and the copies drift. The siphonophore library already carries a
temporary local copy of the BCP-47 → Tesseract table (#215) for that reason.

So library-facing work lives here, in tiers (see AGENTS.md, "Where
library-facing work goes"): generic PDF/OCR knowledge in `pipeline/`,
read-only inspection in `tools/`, collection-specific judgment in `skills/`.

This test is what keeps the tiers from collapsing. Without it the natural
drift is for a `pipeline/` module to reach into `tools/` for something
convenient, at which point the product depends on its own measurement
tooling and the wall is gone.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
PRODUCT_DIRS = ("pipeline", "mcpsrv", "bib")
FORBIDDEN = ("tools", "skills")


def _imported_top_levels(path: Path):
    """Top-level module names imported by one file, absolute imports only."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except SyntaxError:  # pragma: no cover - a syntax error is its own failure
        return set()
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            # level > 0 is a relative import and cannot escape the package.
            if not node.level and node.module:
                names.add(node.module.split(".")[0])
    return names


@pytest.mark.parametrize("product_dir", PRODUCT_DIRS)
def test_the_product_does_not_import_its_tooling(product_dir):
    offenders = []
    for py in sorted((REPO_ROOT / product_dir).rglob("*.py")):
        bad = _imported_top_levels(py) & set(FORBIDDEN)
        if bad:
            offenders.append(f"{py.relative_to(REPO_ROOT)} imports {sorted(bad)}")
    assert not offenders, (
        "pipeline/mcpsrv/bib must not depend on tools/ or skills/. Dependencies "
        "point one way — a skill may import from pipeline/, never the reverse. "
        "See AGENTS.md, 'Where library-facing work goes'.\n\n" + "\n".join(offenders)
    )


def test_tooling_may_import_the_product():
    """The permitted direction, asserted so the rule is not read as symmetric.

    `tools/qc/` loads `pipeline.figures` for `parse_figure_number`, which is
    the point: the scorer measures the shipped parser rather than a copy of it.
    """
    src = (REPO_ROOT / "tools" / "qc" / "caption_binding.py").read_text(encoding="utf-8")
    assert "from pipeline.figures import" in src
