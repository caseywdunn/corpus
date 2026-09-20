"""Focused tests for corpus-agnostic MCP smoke-test result handling."""

from types import SimpleNamespace

from tools.smoke_test_sse import _parse_tool_result, _result_items


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
