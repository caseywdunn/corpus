"""`lexicon_matrix(detail=True)` is bounded and says so (#83, #88).

#83 proposed a column-store row shape to cut per-row key repetition, and
was explicitly held pending evidence that it matters. Measured, it does
not — not in the way proposed:

* the default `detail=False` view is **469-1,606 bytes** on the reference
  corpora, where key repetition is irrelevant;
* the opt-in `detail=True` grid is **382 kB over 1,775 rows**
  (siphonophore) and 143-175 kB over 699 (viburnum), and a column store
  saves a measured **16.4-20.0%** there — which does not make a 382 kB
  payload deliverable, it makes an undeliverable one 19% smaller.

What the measurement did expose is that #88 made the grid opt-in because
it was a multi-MB runaway and never bounded it: no cap, no flag, so a
caller could not tell a complete grid from one its transport dropped.
Same defect class as #166, so the same treatment.
"""
from __future__ import annotations

import json
import types
from pathlib import Path

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools import lexicon as lexicon_tools
from mcpsrv.tools.lexicon import lexicon_matrix


@pytest.fixture
def served(tmp_path, request):
    """A corpus of `n` papers, each with `terms` lexicon terms."""
    n_papers, n_terms, pad = getattr(request, "param", (5, 3, ""))
    terms = [f"term_{i}{pad}" for i in range(n_terms)]
    papers, lex_to_papers, counts = {}, {}, {}
    for i in range(n_papers):
        h = f"{i:012d}"
        hd = tmp_path / h
        hd.mkdir()
        (hd / "anatomy.json").write_text(json.dumps({
            "category": "anatomy",
            "terms": [{"canonical": t, "mention_count": i + 1} for t in terms],
        }))
        papers[h] = {"hash": h, "title": f"Paper {i}" + pad, "year": 1900 + i,
                     "hash_dir": str(hd)}
    for t in terms:
        lex_to_papers[t] = list(papers)
        counts[t] = {h: 3 for h in papers}
    idx = types.SimpleNamespace(
        papers=papers,
        lexicon_to_papers={"anatomy": lex_to_papers},
        lexicon_mention_counts={"anatomy": counts},
        lexicon_surface_to_canonical={"anatomy": {}},
    )
    original = mcp_app._INDEX
    mcp_app.set_index(idx)
    yield idx
    mcp_app.set_index(original)


# --- the default view is not the problem and must not change ------------


@pytest.mark.parametrize("served", [(200, 20, "")], indirect=True)
def test_the_default_view_stays_compact_and_ungrown(served):
    """469-1,606 bytes on the real corpora. Nothing here applies to it."""
    out = lexicon_matrix(category="anatomy")
    assert out["detail"] is False
    assert "rows" not in out
    assert "truncated" not in out
    assert len(json.dumps(out)) < 4_000


# --- the grid is bounded ------------------------------------------------


@pytest.mark.parametrize("served", [(400, 20, "x" * 120)], indirect=True)
def test_a_large_grid_is_bounded_and_admits_it(served, monkeypatch):
    monkeypatch.setattr(lexicon_tools, "LEXICON_MATRIX_MAX_BYTES", 30_000)
    out = lexicon_matrix(category="anatomy", detail=True)
    assert out["truncated"] is True
    assert out["truncated_reason"] == "response_bytes"
    assert out["rows_available"] == 400
    assert out["rows_returned"] < 400
    assert len(out["rows"]) == out["rows_returned"]


@pytest.mark.parametrize("served", [(400, 20, "x" * 120)], indirect=True)
def test_the_reported_size_is_the_real_size(served, monkeypatch):
    """The trap #166 hit and this repeated: a ceiling that does not count
    the fields it adds is not a ceiling. That first cut came back 51 bytes
    over its own limit."""
    monkeypatch.setattr(lexicon_tools, "LEXICON_MATRIX_MAX_BYTES", 30_000)
    out = lexicon_matrix(category="anatomy", detail=True)
    real = len(json.dumps(out, default=str).encode("utf-8"))
    assert out["response_bytes"] == real
    assert real <= 30_000


@pytest.mark.parametrize("served", [(4, 3, "")], indirect=True)
def test_a_small_grid_is_complete_and_says_so(served):
    out = lexicon_matrix(category="anatomy", detail=True)
    assert out["truncated"] is False
    assert "truncated_reason" not in out
    assert out["rows_available"] == out["rows_returned"] == 4


@pytest.mark.parametrize("served", [(400, 20, "x" * 400)], indirect=True)
def test_at_least_one_row_survives(served, monkeypatch):
    """A ceiling so low that no row fits should still return the shape,
    not an empty grid that reads as "no papers have these terms"."""
    monkeypatch.setattr(lexicon_tools, "LEXICON_MATRIX_MAX_BYTES", 200)
    out = lexicon_matrix(category="anatomy", detail=True)
    assert out["rows_returned"] == 1
    assert out["truncated"] is True
    assert out["rows_available"] == 400


@pytest.mark.parametrize("served", [(400, 20, "")], indirect=True)
def test_narrowing_the_paper_set_is_the_way_to_a_usable_grid(served):
    """What the docstring tells a caller to do instead of asking for
    1,775 rows."""
    hashes = [f"{i:012d}" for i in range(5)]
    out = lexicon_matrix(category="anatomy", detail=True, paper_hashes=hashes)
    assert out["rows_available"] == 5
    assert out["truncated"] is False


def test_the_ceiling_is_operator_overridable():
    import inspect
    assert "CORPUS_LEXICON_MATRIX_MAX_BYTES" in inspect.getsource(lexicon_tools)


# --- the declined proposal ---------------------------------------------


def test_the_input_schema_is_unchanged():
    """No column-store parameter, so the 1.0 input freeze holds and no
    caller has to change a call. The decision and its numbers live in the
    tool's own comments, not in PLAN.md, which turns over each release."""
    import inspect
    sig = inspect.signature(lexicon_matrix.fn if hasattr(lexicon_matrix, "fn")
                            else lexicon_matrix)
    assert list(sig.parameters) == [
        "category", "terms", "top_n", "paper_hashes",
        "year_from", "year_to", "detail",
    ]
    src = inspect.getsource(lexicon_tools.lexicon_matrix.fn
                            if hasattr(lexicon_tools.lexicon_matrix, "fn")
                            else lexicon_tools.lexicon_matrix)
    assert "row_schema" not in src


@pytest.mark.parametrize("served", [(400, 20, "x" * 120)], indirect=True)
def test_the_reported_count_matches_the_rows_beside_it(served, monkeypatch):
    """A count that disagrees with its own payload is the defect this
    change is about, so it must not be reintroduced by the trim loop —
    testing the loop condition with the stamping helper left
    `rows_returned` one ahead of the rows actually returned."""
    for ceiling in (200, 5_000, 30_000, 200_000):
        monkeypatch.setattr(lexicon_tools, "LEXICON_MATRIX_MAX_BYTES", ceiling)
        out = lexicon_matrix(category="anatomy", detail=True)
        assert out["rows_returned"] == len(out["rows"]), ceiling
        assert out["response_bytes"] == len(
            json.dumps(out, default=str).encode("utf-8")), ceiling
