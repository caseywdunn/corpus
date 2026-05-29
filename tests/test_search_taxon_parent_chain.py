"""search_taxon parent_chain= ancestry walk (#88).

parent_chain=True walks the DwC parent_name_usage_id tree upward and
attaches the ancestor list (immediate parent → root). Additive: the
field is absent unless requested. The walk is cycle- and depth-guarded
so a malformed snapshot can't hang the tool.
"""
from __future__ import annotations

import sqlite3
import types

import pytest

from mcpsrv import app as mcp_app
from mcpsrv.tools.taxonomy import search_taxon


def _make_taxonomy_db(conn: sqlite3.Connection, lookups):
    return types.SimpleNamespace(
        conn=conn,
        lookup=lambda name: lookups.get(name),
    )


@pytest.fixture
def index():
    conn = sqlite3.connect(":memory:")
    conn.execute(
        """CREATE TABLE taxa (
             taxon_id TEXT PRIMARY KEY,
             scientific_name TEXT,
             taxon_rank TEXT,
             parent_name_usage_id TEXT
           )"""
    )
    # Apolemia uvaria → Apolemia → Apolemiidae → Siphonophorae (root).
    rows = [
        ("t:sipho", "Siphonophorae", "order", None),
        ("t:apolfam", "Apolemiidae", "family", "t:sipho"),
        ("t:apolemia", "Apolemia", "genus", "t:apolfam"),
        ("t:uvaria", "Apolemia uvaria", "species", "t:apolemia"),
        # A self-referential cycle, to exercise the guard.
        ("t:cycle", "Cyclus badus", "species", "t:cycle"),
    ]
    conn.executemany("INSERT INTO taxa VALUES (?, ?, ?, ?)", rows)
    conn.commit()

    lookups = {
        "Apolemia uvaria": {"accepted_taxon_id": "t:uvaria",
                            "scientific_name": "Apolemia uvaria",
                            "rank": "species", "name_type": "accepted"},
        "Siphonophorae": {"accepted_taxon_id": "t:sipho",
                          "scientific_name": "Siphonophorae",
                          "rank": "order", "name_type": "accepted"},
        "Cyclus badus": {"accepted_taxon_id": "t:cycle",
                         "scientific_name": "Cyclus badus",
                         "rank": "species", "name_type": "accepted"},
    }
    fake = types.SimpleNamespace(
        taxonomy_db=_make_taxonomy_db(conn, lookups),
        taxon_to_papers={"t:uvaria": ["h1"]},
    )
    original = mcp_app._INDEX
    mcp_app.set_index(fake)
    yield fake
    mcp_app.set_index(original)


def test_parent_chain_omitted_by_default(index):
    out = search_taxon("Apolemia uvaria")
    assert "parent_chain" not in out


def test_parent_chain_walks_immediate_parent_to_root(index):
    out = search_taxon("Apolemia uvaria", parent_chain=True)
    names = [a["scientific_name"] for a in out["parent_chain"]]
    assert names == ["Apolemia", "Apolemiidae", "Siphonophorae"]
    assert out["parent_chain"][0] == {
        "taxon_id": "t:apolemia",
        "scientific_name": "Apolemia",
        "rank": "genus",
    }


def test_parent_chain_empty_at_root(index):
    out = search_taxon("Siphonophorae", parent_chain=True)
    assert out["parent_chain"] == []


def test_parent_chain_cycle_is_guarded(index):
    # Self-referential parent must not loop forever.
    out = search_taxon("Cyclus badus", parent_chain=True)
    assert out["parent_chain"] == []
