"""Caption preview default on the older figure-listing tools (#85).

get_figures_for_taxon / get_figures_for_lexicon_term return a trimmed
caption preview by default; full_caption=True restores the verbatim
caption. Pins the token-cost contract for figure listings.
"""
from __future__ import annotations

import json
import types
from pathlib import Path

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools.figures import _FIGURE_CAPTION_PREVIEW_CHARS, get_figures_for_lexicon_term

_LONG_CAPTION = "Figure 1. A nectophore " + ("blah " * 80) + "end."


@pytest.fixture
def index(tmp_path: Path):
    h = "aaaaaaaaaaaa"
    hash_dir = tmp_path / "documents" / h
    hash_dir.mkdir(parents=True)
    (hash_dir / "figures.json").write_text(json.dumps({"figures": [
        {"figure_id": "docling_1", "filename": "fig_1.png",
         "figure_type": "figure", "caption_text": _LONG_CAPTION},
    ]}))
    fake = types.SimpleNamespace(
        papers={h: {"hash_dir": str(hash_dir), "title": "A paper"}},
        lexicon_to_papers={"anatomy": {"nectophore": [h]}},
    )
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    yield fake
    mcp_app.set_index(original)


def test_caption_is_previewed_by_default(index):
    rows = get_figures_for_lexicon_term("anatomy", "nectophore")
    cap = rows[0]["caption_text"]
    assert cap != _LONG_CAPTION                  # trimmed
    assert len(cap) <= _FIGURE_CAPTION_PREVIEW_CHARS + 1   # +1 for the ellipsis
    assert cap.endswith("…")


def test_full_caption_opt_in(index):
    rows = get_figures_for_lexicon_term("anatomy", "nectophore", full_caption=True)
    assert rows[0]["caption_text"] == _LONG_CAPTION
