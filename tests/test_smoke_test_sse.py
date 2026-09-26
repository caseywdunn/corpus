"""Focused tests for corpus-agnostic MCP smoke-test result handling."""

import asyncio
import inspect
import json
from contextlib import asynccontextmanager
from types import SimpleNamespace

import pytest

from tools.smoke_test_sse import (
    _parse_tool_result, _result_items, layer2_mcp_client, layer3_tool_coverage,
)


def _result(*texts):
    return SimpleNamespace(
        content=[SimpleNamespace(text=text) for text in texts]
    )


def test_empty_list_tool_result_is_a_valid_empty_collection():
    parsed = _parse_tool_result(_result())
    assert parsed is None
    assert _result_items(parsed) == []


def test_singleton_list_tool_result_becomes_one_discoverable_item():
    parsed = _parse_tool_result(
        _result('{"hash": "abc", "first_author": "Example"}')
    )
    assert _result_items(parsed) == [
        {"hash": "abc", "first_author": "Example"}
    ]


def test_multi_item_list_tool_result_is_reassembled():
    parsed = _parse_tool_result(
        _result('{"name": "Taxon one"}', '{"name": "Taxon two"}')
    )
    assert _result_items(parsed) == [
        {"name": "Taxon one"},
        {"name": "Taxon two"},
    ]


@pytest.mark.parametrize("value", ["unparseable", [None], {"error": "failed"}, [{"error": "failed"}]])
def test_errors_and_malformed_rows_are_not_empty_results(value):
    with pytest.raises(ValueError):
        _result_items(value)


@pytest.mark.parametrize("flag", ["is_error", "isError"])
def test_sdk_error_without_content_is_not_an_empty_collection(flag):
    result = _result()
    setattr(result, flag, True)
    with pytest.raises(ValueError, match="MCP tool returned an error"):
        _parse_tool_result(result)


def _install_session(monkeypatch, *, count=0, taxonomy=True, author=True, overrides=None):
    """Exercise the real coverage routine with SDK-shaped result blocks.

    The values deliberately describe a different collection from the reference
    corpus. Bind each request against the real tool signature so the fake
    session cannot hide unsupported arguments.
    """
    import mcp
    import mcp.client.sse

    from mcpsrv.tools import bibliography, chunks, figures, lexicon, papers, taxonomy as taxa
    from mcpsrv.tools import profiles

    modules = [bibliography, chunks, figures, lexicon, papers, taxa, profiles]
    calls = []
    paper_hash = "a1b2c3d4e5f6"
    surname = "del Río"
    taxon_name = "Quercus robur"
    responses = {
        "bundle_info": {"bundle_version": "test-fixture"},
        "list_papers": [{"hash": paper_hash, "first_author": "M. del Río" if author else ""}],
        "get_papers": [{"hash": paper_hash, "first_author": surname}],
        "corpus_summary": {
            "n_papers": 7,
            "top_taxa": [{"name": taxon_name}] if taxonomy else [],
            "lexicon_categories": ["leaf_shape"] if taxonomy else [],
        },
        "search_taxon": {"matched_taxon_id": "oak-1", "accepted_name": taxon_name},
        "list_valid_species_under": [{"accepted_taxon_id": f"oak-{i}"} for i in range(count)],
        "get_papers_for_taxon": [{"hash": paper_hash} for _ in range(count)],
        "get_taxon_dossier": {"taxon": {"name": taxon_name}, "papers": []},
        "get_chunks_for_topic": [{"chunk_id": f"c{i}"} for i in range(count)],
        "get_chunks_by_section": [{"chunk_id": f"c{i}"} for i in range(count)],
        "get_chunks": [{"chunk_id": f"c{i}"} for i in range(count)],
        "get_bibliography": [{"work_id": f"w{i}"} for i in range(count)],
        "get_missing_references": [{"work_id": f"w{i}"} for i in range(count)],
        "get_works_by_author": [{"work_id": f"w{i}"} for i in range(count)],
        "lexicon_matrix": {"term_totals": []},
        "get_figures_for_taxon": [{"figure_id": f"f{i}"} for i in range(count)],
        "list_output_profiles": {"profiles": ["report"]},
    }
    responses.update(overrides or {})

    class Session:
        def __init__(self, *_streams):
            pass

        async def __aenter__(self):
            return self

        async def __aexit__(self, *_exc):
            pass

        async def initialize(self):
            pass

        async def list_tools(self):
            return SimpleNamespace(tools=[SimpleNamespace(name=name) for name in responses])

        async def call_tool(self, name, arguments):
            calls.append((name, arguments))
            function = next(getattr(module, name) for module in modules if hasattr(module, name))
            inspect.signature(function).bind(**arguments)
            value = responses[name]
            if isinstance(value, SimpleNamespace):
                return value
            if isinstance(value, list):
                return _result(*(json.dumps(row) for row in value))
            return _result(json.dumps(value))

    @asynccontextmanager
    async def transport(*_args, **_kwargs):
        yield None, None

    monkeypatch.setattr(mcp, "ClientSession", Session)
    monkeypatch.setattr(mcp.client.sse, "sse_client", transport)
    return calls


@pytest.mark.parametrize("count", [0, 1, 3])
def test_layer3_discovers_foreign_collection_and_accepts_list_encodings(monkeypatch, count):
    calls = _install_session(monkeypatch, count=count)
    assert asyncio.run(layer3_tool_coverage("localhost", 1, "test-token")) == 0
    requests = dict(calls)
    assert requests["search_taxon"] == {"name": "Quercus robur"}
    assert requests["list_valid_species_under"] == {"parent_taxon_name": "Quercus robur", "limit": 3}
    assert requests["get_works_by_author"] == {"surname": "del Río", "limit": 3}
    assert requests["lexicon_matrix"] == {"category": "leaf_shape", "top_n": 3}
    assert requests["get_chunks_by_section"]["limit"] == 3
    assert requests["get_chunks"]["chunk_ids"] == [f"c{i}" for i in range(count)]
    assert requests["get_chunks"]["paper_hash"] == "a1b2c3d4e5f6"


def test_layer3_runs_paper_checks_without_taxonomy_lexicon_or_authors(monkeypatch):
    calls = _install_session(monkeypatch, taxonomy=False, author=False)
    assert asyncio.run(layer3_tool_coverage("localhost", 1, "test-token")) == 0
    names = {name for name, _ in calls}
    assert {"get_chunks", "get_bibliography", "list_output_profiles"} <= names
    assert not names & {"search_taxon", "list_valid_species_under", "lexicon_matrix", "get_works_by_author"}


@pytest.mark.parametrize("tool,value", [
    ("list_papers", []),
    ("list_papers", [{"error": "database broke"}]),
    ("get_chunks_by_section", [{"unexpected": "missing identifier"}]),
    ("get_chunks_for_topic", {"error": "embedding failed", "code": "unavailable"}),
    ("get_papers_for_taxon", [{"error": "lookup failed"}]),
    ("get_figures_for_taxon", SimpleNamespace(content=[], isError=True)),
    ("list_output_profiles", {"error": "profiles failed", "code": "unavailable"}),
    ("list_output_profiles", {"unexpected": "wrong shape"}),
])
def test_layer3_fails_discovery_and_tool_errors(monkeypatch, tool, value):
    _install_session(monkeypatch, overrides={tool: value})
    assert asyncio.run(layer3_tool_coverage("localhost", 1, "test-token")) != 0


def test_layer3_skips_explicitly_unconfigured_optional_features(monkeypatch, capsys):
    _install_session(monkeypatch, overrides={
        name: {"error": "not configured", "code": "not_configured"}
        for name in ["get_chunks_for_topic", "get_missing_references", "get_works_by_author"]
    })
    assert asyncio.run(layer3_tool_coverage("localhost", 1, "test-token")) == 0
    assert "skipping semantic search" in capsys.readouterr().out


@pytest.mark.parametrize("value,expected", [
    ([], 0),
    ([{"hash": "a1b2c3d4e5f6"}], 0),
    ({"error": "database failed", "code": "unavailable"}, 1),
    ("not a list", 1),
])
def test_layer2_validates_paper_results_without_layer3(monkeypatch, value, expected):
    _install_session(monkeypatch, overrides={"list_papers": value})
    assert asyncio.run(layer2_mcp_client("localhost", 1, "test-token")) == expected
