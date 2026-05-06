"""Regression test for issue #49.

The taxon_mentions SQLite column is named ``taxon_rank`` (declared in
``build_taxon_mentions.create_schema``). Prior to the fix, the MCP
``get_taxon_mentions`` fast path read ``r["rank"]``, which raises
``KeyError`` on every fast-path query. The fallback path
(per-paper ``taxa.json`` scan) was fine because that JSON shape uses
``rank`` as the dict key.

This test seeds the schema by hand, exercises the wrapper, and asserts
that the row dict returned to the tool path can be indexed by
``taxon_rank`` — the column name the schema actually declares.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from build_taxon_mentions import create_schema
from mcpsrv.indexes import TaxonMentionDB


def _seed(db_path: Path) -> None:
    conn = sqlite3.connect(db_path)
    create_schema(conn)
    conn.execute(
        """INSERT INTO taxon_mentions
           (taxon_id, matched_name, accepted_name, taxon_rank,
            corpus_hash, chunk_id, chunk_index,
            char_start, char_end, mention_text, name_type, method)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        ("1371", "Nanomia", "Nanomia",
         "Genus", "abc123def456", "c0", 0, 10, 17, "Nanomia",
         "binomial", "regex_taxonomy"),
    )
    conn.commit()
    conn.close()


def test_mentions_for_taxon_returns_taxon_rank_key(tmp_path):
    db_path = tmp_path / "taxon_mentions.sqlite"
    _seed(db_path)

    db = TaxonMentionDB(db_path)
    rows = db.mentions_for_taxon("1371")

    assert rows, "expected one seeded row"
    row = rows[0]

    # The MCP fast-path used to do row["rank"] here and crashed.
    # The schema column is taxon_rank — the API output dict can
    # still be keyed by "rank" downstream, but the SQL row exposes
    # only the column names that exist in the schema.
    assert row["taxon_rank"] == "Genus"
    assert "rank" not in row, (
        "If the schema is ever extended to publish a separate `rank` "
        "column, update mcpsrv/tools/taxonomy.py accordingly."
    )


def test_mentions_for_taxon_filters_by_corpus_hash(tmp_path):
    db_path = tmp_path / "taxon_mentions.sqlite"
    _seed(db_path)
    db = TaxonMentionDB(db_path)

    assert db.mentions_for_taxon("1371", corpus_hash="abc123def456")
    assert not db.mentions_for_taxon("1371", corpus_hash="missinghash")
