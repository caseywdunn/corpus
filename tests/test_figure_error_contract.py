"""Check figure errors through actual MCP result conversion, not exception text."""
import asyncio
import json

import pytest

from mcp.server.mcpserver import Context
from mcpsrv.app import mcp
from mcpsrv.tools.figures import get_figure_url
from mcpsrv.tools.profiles import list_output_profiles
from tests.test_figure_licensing_states import _make_index, _CLEARED, _UNCLEARED, HASH


def call_image(**kwargs):
    ctx = Context(mcp_server=mcp, subscriptions=mcp._subscriptions)
    return asyncio.run(mcp._tool_manager.call_tool(
        "get_figure_image", {"paper_hash": HASH, "figure_id": "docling_1", **kwargs},
        ctx, convert_result=True,
    ))


def test_matching_refusal_has_structured_reason_and_transport_error(tmp_path):
    idx = _make_index(tmp_path, work=_UNCLEARED)
    idx.figure_url_base = "https://corpus.example.test"
    url = get_figure_url(HASH, "docling_1", profile="manuscript")
    image = call_image(profile="manuscript")
    assert image.is_error is True
    payload = image.structured_content
    for key in ("code", "profile", "publication_clearance", "license_source"):
        assert payload[key] == url[key]
    assert json.loads(image.content[0].text) == payload
    assert payload["code"] == "forbidden"
    assert all(block.type != "image" for block in image.content)


@pytest.mark.parametrize("kwargs,code", [
    ({"profile": "unknown"}, "invalid_argument"),
    ({"paper_hash": "ffffffffffff"}, "not_found"),
    ({"figure_id": "absent"}, "not_found"),
])
def test_other_refusals_are_machine_readable(tmp_path, kwargs, code):
    _make_index(tmp_path, work=_CLEARED)
    result = call_image(**kwargs)
    assert result.is_error is True
    assert result.structured_content["code"] == code


def test_successful_inline_image_contract_is_unchanged(tmp_path):
    _make_index(tmp_path, work=_CLEARED)
    result = call_image(profile="manuscript", label="A")
    assert not result.is_error
    assert result.structured_content is None
    assert len(result.content) == 1
    assert result.content[0].type == "image"
    assert result.content[0].mime_type == "image/png"


@pytest.mark.parametrize("requested_profile", [None, "manuscript", "presentation"])
def test_client_selects_in_chat_alternative_from_structured_profiles(tmp_path, requested_profile):
    """A client can change its output purpose without parsing refusal prose."""
    idx = _make_index(tmp_path, work=_UNCLEARED, default_profile="manuscript")
    idx.figure_url_base = "https://corpus.example.test"
    refused = call_image(profile=requested_profile)
    assert refused.is_error
    reason = refused.structured_content
    assert reason["code"] == "forbidden"
    assert reason["publication_clearance"] == "no_record"

    ctx = Context(mcp_server=mcp, subscriptions=mcp._subscriptions)
    discovered = asyncio.run(mcp._tool_manager.call_tool(
        list_output_profiles.__name__, {}, ctx, convert_result=True,
    ))
    assert not discovered.is_error
    # The SDK wraps generic dict return annotations in a `result` object.
    wire = discovered.structured_content
    catalog = wire.get("result", wire)
    policies = {profile["name"]: profile for profile in catalog["profiles"]}
    assert policies[reason["profile"]]["figure_licensing"] == "strict"
    # The client explicitly wants in-chat display now, not publication output.
    # Discover the permitted policy instead of extracting a name from `error`.
    selected = next(profile["name"] for profile in catalog["profiles"]
                    if profile["figure_licensing"] == "permissive")
    assert selected != reason["profile"]
    image = call_image(profile=selected)
    assert not image.is_error
    assert image.structured_content is None
    assert len(image.content) == 1 and image.content[0].type == "image"
    assert image.content[0].mime_type == "image/png"
    url = get_figure_url(HASH, "docling_1", profile=selected)
    assert url["profile"] == selected and "url" in url and "error" not in url

    # Per-call recovery neither changes the server default nor clears rights.
    assert idx.default_profile == catalog["server_default"] == "manuscript"
    assert call_image(profile=requested_profile).structured_content == reason
    strict_url = get_figure_url(HASH, "docling_1", profile=requested_profile)
    for field in ("code", "profile", "publication_clearance", "license_source"):
        assert strict_url[field] == reason[field]
