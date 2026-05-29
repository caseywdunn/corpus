"""1.0 API freeze-contract meta-test (#88 §3 — the Phase-3 freeze gate).

Pins the MCP tool surface that 1.0 freezes, so any drift is a deliberate,
reviewed edit rather than an accident:

1. The exact set of registered tool names (39).
2. The pagination-naming convention — single-knob tools use ``limit``;
   only the documented multi-cap dossier/graph tools use ``max_*``.
3. The uniform error payload — every tool error return carries a human
   ``error`` message plus a machine ``code`` from a frozen vocabulary,
   built through :func:`mcpsrv.app.error`; no site uses a bare
   ``{"error": ...}`` dict.
"""
from __future__ import annotations

import inspect
import pathlib
import re

import pytest

import mcpsrv.main  # noqa: F401 — triggers @mcp.tool() registration
import mcpsrv.tools
from mcpsrv.app import ERROR_CODES, error, mcp


# The frozen 1.0 tool surface. Editing this set is editing the public
# API — pair any change with a CHANGELOG migration note.
EXPECTED_TOOLS = frozenset({
    "bundle_info", "corpus_summary", "format_citations", "get_active_profile",
    "get_bibliography", "get_chunks", "get_chunks_by_section",
    "get_chunks_for_taxon", "get_chunks_for_topic", "get_citation_graph",
    "get_excerpts_citing", "get_figure", "get_figure_dossier_for_taxon",
    "get_figure_dossier_for_term", "get_figure_image", "get_figure_roi_image",
    "get_figure_url", "get_figures_for_lexicon_term", "get_figures_for_taxon",
    "get_intext_citations", "get_lexicon_term_dossier",
    "get_missing_references", "get_original_description", "get_papers",
    "get_papers_by_author", "get_papers_for_taxon", "get_taxon_dossier",
    "get_taxon_lexicon_slice", "get_taxon_mentions",
    "get_taxon_subtree_dossier", "get_works_by_author", "lexicon_matrix",
    "list_figure_rois", "list_output_profiles", "list_papers",
    "list_valid_species_under", "resolve_reference", "search_taxon",
    "translate_chunk",
})

# The documented exception to "``limit`` everywhere": tools where one call
# returns several independently-bounded sections, so a single ``limit``
# can't express the caps. These keep their descriptive ``max_*`` names.
MULTI_CAP_TOOLS = frozenset({
    "get_citation_graph",
    "get_figure_dossier_for_taxon",
    "get_figure_dossier_for_term",
    "get_lexicon_term_dossier",
    "get_taxon_dossier",
    "get_taxon_lexicon_slice",
    "get_taxon_subtree_dossier",
})


def _registered():
    return {t.name: t for t in mcp._tool_manager.list_tools()}


def _fn(tool):
    return getattr(tool, "fn", None) or getattr(tool, "func", None)


# --- 1. tool surface ----------------------------------------------------------


def test_registered_tool_surface_is_frozen():
    names = set(_registered())
    assert names == set(EXPECTED_TOOLS), (
        "MCP tool surface changed — this set is what 1.0 freezes. Update "
        "EXPECTED_TOOLS deliberately + add a CHANGELOG migration note.\n"
        f"added={sorted(names - EXPECTED_TOOLS)} "
        f"removed={sorted(EXPECTED_TOOLS - names)}"
    )


def test_tool_count_is_39():
    assert len(_registered()) == 39


def test_removed_singular_tools_stay_removed():
    """Phase-2 (#88 §2.3) removals must not creep back in."""
    names = set(_registered())
    assert {"format_citation", "get_paper", "get_chunk"}.isdisjoint(names)


# --- 2. pagination naming -----------------------------------------------------


def test_pagination_naming_convention():
    """`limit` is the single-knob pagination name everywhere; only the
    documented multi-cap tools may expose `max_*` params."""
    offenders = {}
    for name, tool in _registered().items():
        fn = _fn(tool)
        if fn is None:
            continue
        max_params = sorted(
            p for p in inspect.signature(fn).parameters if p.startswith("max_")
        )
        if max_params and name not in MULTI_CAP_TOOLS:
            offenders[name] = max_params
    assert not offenders, (
        f"tools expose max_* without being in the documented dossier "
        f"exception: {offenders}. Rename the knob to `limit`, or (if the "
        f"tool genuinely has several independent caps) add it to "
        f"MULTI_CAP_TOOLS."
    )


def test_multi_cap_allowlist_is_live():
    """Every allowlisted tool is registered and really has >1 cap — keeps
    the exception from rotting into a rubber stamp."""
    reg = _registered()
    for name in MULTI_CAP_TOOLS:
        assert name in reg, f"stale MULTI_CAP_TOOLS entry: {name}"
        caps = [
            p for p in inspect.signature(_fn(reg[name])).parameters
            if p.startswith("max_")
        ]
        assert caps, f"{name} is allowlisted but exposes no max_* caps"


# --- 3. uniform error shape ---------------------------------------------------


def test_error_helper_shape():
    e = error("no such work", "not_found", queried="x")
    assert e["error"] == "no such work"
    assert e["code"] == "not_found"
    assert e["queried"] == "x"  # extra context preserved
    assert e["code"] in ERROR_CODES


def test_error_codes_are_a_closed_vocabulary():
    # The frozen set 1.0 clients branch on. Additions are fine; this just
    # documents the baseline so a typo'd code is caught in review.
    assert {
        "not_found", "ambiguous", "invalid_argument", "not_configured",
        "no_results", "unavailable", "empty_item", "forbidden",
    } <= set(ERROR_CODES)


def test_no_raw_error_dict_returns_in_tools():
    """Every tool error return goes through app.error() — no bare
    ``return {"error": ...}`` / ``return [{"error": ...}]`` that would
    ship a payload without the frozen ``code``."""
    tools_dir = pathlib.Path(mcpsrv.tools.__file__).parent
    raw = re.compile(r'return \[?\{\s*"error"')
    offenders = []
    for path in sorted(tools_dir.glob("*.py")):
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if raw.search(line):
                offenders.append(f"{path.name}:{i}: {line.strip()}")
    assert not offenders, (
        "raw error-dict returns bypass app.error() (no machine `code`):\n"
        + "\n".join(offenders)
    )
