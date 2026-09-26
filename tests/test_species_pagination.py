"""Bounded taxonomy lists preserve list payloads and a usable MCP session (#338)."""
import asyncio
import hashlib
import json
import sqlite3
from types import SimpleNamespace

import anyio
import pytest
from mcp import ClientSession
from mcp.server.mcpserver import Context
from mcp.shared.memory import create_client_server_memory_streams

from mcpsrv import app
from mcpsrv.tools import taxonomy
from pipeline.taxa import TaxonomyDB
from pipeline.taxonomy_ingest import create_schema


@pytest.fixture
def snapshot(tmp_path, monkeypatch):
    opened = []

    def make(count=8, authorship="Mańko, 2026", cycle=False):
        path = tmp_path / f"taxonomy-{len(opened)}.sqlite"
        conn = sqlite3.connect(path)
        create_schema(conn)
        rows = [("root", "Porifera", "phylum", "accepted", "genus" if cycle else None, None),
                ("genus", "Example", "genus", "accepted", "root", None),
                ("alias", "Oldphylum", "phylum", "unaccepted", None, "root"),
                ("empty", "Empty", "genus", "accepted", "root", None)]
        rows += [(f"s{i:04}", f"Example species{i:04}", "species", "accepted", "genus", None)
                 for i in reversed(range(count))]
        rows += [("sub", "Example species0000 minor", "SUBSPECIES", "accepted", "s0000", None),
                 ("z", "Tied species", "species", "accepted", "genus", None),
                 ("a", "Tied species", "species", "accepted", "genus", None),
                 ("syn", "Old species", "species", "synonym", "genus", "s0000"),
                 ("var", "Example variety", "variety", "accepted", "genus", None)]
        conn.executemany("""INSERT INTO taxa (
            taxon_id,scientific_name,taxon_rank,taxonomic_status,parent_name_usage_id,
            accepted_name_usage_id,scientific_name_authorship) VALUES (?,?,?,?,?,?,?)""",
                         [(*row, authorship) for row in rows])
        conn.executemany("INSERT INTO names VALUES (?,?,?,?)",
                         [(name, name.lower(), tid, "accepted" if status == "accepted" else "synonym")
                          for tid, name, _, status, _, _ in rows])
        conn.commit()
        conn.close()
        db = TaxonomyDB(path)
        opened.append(db)
        idx = SimpleNamespace(taxonomy_db=db, taxon_to_papers={"s0000": ["p1", "p2"]})
        monkeypatch.setattr(app, "_INDEX", idx)
        return idx

    yield make
    for db in opened:
        db.close()


def dispatch(**arguments):
    ctx = Context(mcp_server=app.mcp, subscriptions=app.mcp._subscriptions)
    return asyncio.run(app.mcp._tool_manager.call_tool(
        "list_valid_species_under", {"parent_taxon_name": "Porifera", **arguments},
        ctx, convert_result=True,
    ))


def rows(result):
    return result.structured_content["result"]


def page(result):
    return result.model_dump(mode="json", by_alias=True, exclude_none=True)["_meta"]["pagination"]


def assert_bounded(result):
    size = len(result.model_dump_json(by_alias=True, exclude_none=True).encode())
    assert size <= taxonomy.SPECIES_MAX_RESULT_BYTES
    if result.meta:
        assert page(result)["result_bytes"] == size
    assert [json.loads(block.text) for block in result.content] == rows(result)
    event = "event: message\ndata: " + json.dumps({
        "jsonrpc": "2.0", "id": 1,
        "result": result.model_dump(mode="json", by_alias=True, exclude_none=True),
    }, ensure_ascii=False, separators=(",", ":")) + "\n\n"
    assert len(event.encode()) < 1_048_576


def test_legacy_small_list_shape_rank_synonym_and_tied_order(snapshot):
    snapshot()
    legacy = taxonomy.list_valid_species_under("Porifera")
    assert isinstance(legacy, list) and len(legacy) == 11
    assert set(legacy[0]) == {"accepted_taxon_id", "accepted_name", "authorship", "rank", "mentioning_paper_count"}
    assert legacy[0]["mentioning_paper_count"] == 2
    assert legacy[-2]["accepted_taxon_id"] == "a" and legacy[-1]["accepted_taxon_id"] == "z"
    assert {r["rank"] for r in legacy} == {"species", "SUBSPECIES"}
    assert taxonomy.list_valid_species_under("Oldphylum") == legacy
    assert [r["accepted_taxon_id"] for r in taxonomy.list_valid_species_under("Example species0000")] == ["sub"]
    assert taxonomy.list_valid_species_under("Empty") == []
    assert taxonomy.list_valid_species_under("missing") == []
    assert rows(dispatch()) == legacy


@pytest.mark.parametrize("cycle", [False, True])
def test_pagination_is_deterministic_gap_free_and_cycle_safe(snapshot, cycle):
    snapshot(cycle=cycle)
    expected = taxonomy.list_valid_species_under("Porifera")
    got, offset = [], 0
    while True:
        result = dispatch(limit=3, offset=offset)
        assert_bounded(result)
        assert result == dispatch(parent_taxon_name="Oldphylum", limit=3, offset=offset)
        info = page(result)
        assert info["returned"] == len(rows(result)) and info["total_available"] == len(expected)
        got.extend(rows(result))
        if info["next_offset"] is None:
            assert not info["truncated"]
            break
        assert info["next_offset"] == offset + len(rows(result)) > offset
        offset = info["next_offset"]
    assert got == expected
    assert len({r["accepted_taxon_id"] for r in got}) == len(expected)


def test_unicode_byte_shortened_pages_preserve_every_complete_row(snapshot, monkeypatch):
    snapshot(count=12, authorship="作者𝛼é" * 100)
    expected = taxonomy.list_valid_species_under("Porifera")
    monkeypatch.setattr(taxonomy, "SPECIES_MAX_RESULT_BYTES", 6000)
    got, offset = [], 0
    while True:
        result = dispatch(limit=10, offset=offset)
        assert_bounded(result)
        info = page(result)
        assert not result.is_error
        got.extend(rows(result))
        if info["next_offset"] is None:
            break
        assert "response_bytes" in info["truncated_reason"]
        assert 0 < info["returned"] < 10
        assert info["next_offset"] > offset
        offset = info["next_offset"]
    assert got == expected
    assert len({r["accepted_taxon_id"] for r in got}) == len(expected)


@pytest.mark.parametrize("arguments", [{"limit": 0}, {"limit": -1}, {"offset": 1},
                                       {"limit": 3, "offset": -1}, {"limit": 3, "offset": 2**63}])
def test_invalid_limits_and_offsets_are_bounded_errors(snapshot, arguments):
    snapshot()
    result = dispatch(**arguments)
    assert rows(result)[0]["code"] == "invalid_argument"
    assert_bounded(result)


def test_cap_empty_unknown_and_missing_taxonomy(snapshot):
    idx = snapshot(count=600)
    capped = dispatch(limit=100_000)
    assert_bounded(capped)
    assert page(capped)["limit"] == app.MAX_LIMIT
    assert len(rows(capped)) == app.MAX_LIMIT
    assert page(capped)["next_offset"] == app.MAX_LIMIT
    assert page(capped)["total_available"] == 603
    for arguments in ({"parent_taxon_name": "Empty"}, {"parent_taxon_name": "missing"}, {"offset": 1000}):
        result = dispatch(limit=3, **arguments)
        assert rows(result) == []
        assert page(result)["returned"] == 0 and page(result)["next_offset"] is None
        assert_bounded(result)
    idx.taxonomy_db = None
    for arguments in ({}, {"limit": 3}):
        result = dispatch(**arguments)
        assert rows(result)[0]["code"] == "not_configured"
        assert_bounded(result)


def test_single_oversized_row_errors_without_clipping_skipping_or_looping(snapshot, monkeypatch):
    snapshot(authorship="作者" * 10_000)
    monkeypatch.setattr(taxonomy, "SPECIES_MAX_RESULT_BYTES", 6000)
    result = dispatch(limit=3)
    assert result.is_error
    assert rows(result)[0]["reason"] == "row_exceeds_response_budget"
    assert rows(result)[0]["blocked_offset"] == 0
    assert page(result)["returned"] == 0 and page(result)["next_offset"] is None
    assert page(result)["truncated"]
    assert_bounded(result)
    assert rows(dispatch())[0]["reason"] == "pagination_required"


def test_large_root_does_not_break_actual_mcp_client_session(snapshot):
    idx = snapshot(count=80, authorship="作者" * 4096)
    before = hashlib.sha256(idx.taxonomy_db.db_path.read_bytes()).hexdigest()
    # The old unrestricted result alone would exceed the reported 1 MiB event.
    author = "作者" * 4096
    old_rows = [{"accepted_taxon_id": f"s{i:04}", "accepted_name": f"Example species{i:04}",
                 "authorship": author, "rank": "species", "mentioning_paper_count": 0}
                for i in range(80)]
    assert taxonomy._species_result_bytes(taxonomy._species_result(old_rows)) > 1_048_576

    async def run():
        async with create_client_server_memory_streams() as (client, server):
            async with anyio.create_task_group() as tasks:
                low = app.mcp._lowlevel_server
                tasks.start_soon(low.run, *server, low.create_initialization_options())
                async with ClientSession(*client) as session:
                    await session.initialize()
                    catalog = await session.list_tools()
                    assert any(t.name == "list_valid_species_under" for t in catalog.tools)
                    rejected = await session.call_tool("list_valid_species_under", {"parent_taxon_name": "Porifera"})
                    assert rows(rejected)[0]["reason"] == "pagination_required"
                    assert_bounded(rejected)
                    first = await session.call_tool("list_valid_species_under", {"parent_taxon_name": "Porifera", "limit": 3})
                    assert_bounded(first)
                    assert len(rows(first)) == 3 and page(first)["next_offset"] == 3
                    second = await session.call_tool("list_valid_species_under", {
                        "parent_taxon_name": "Porifera", "limit": 3, "offset": page(first)["next_offset"]})
                    assert_bounded(second)
                    assert not ({r["accepted_taxon_id"] for r in rows(first)} &
                                {r["accepted_taxon_id"] for r in rows(second)})
                    following = await session.call_tool("search_taxon", {"name": "Porifera"})
                    assert not following.is_error
                    await session.send_ping()
                tasks.cancel_scope.cancel()

    asyncio.run(run())
    assert hashlib.sha256(idx.taxonomy_db.db_path.read_bytes()).hexdigest() == before
