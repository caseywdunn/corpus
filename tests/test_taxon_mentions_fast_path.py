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

import json
import os
import sqlite3
from pathlib import Path

import pytest

from pipeline.stages import _file_sha256
from pipeline.taxon_mentions import build, create_schema
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


# ---------------------------------------------------------------------------
# #95 — fingerprint-aware freshness
# ---------------------------------------------------------------------------
#
# The mtime gate alone misses the cross-component case: a taxonomy
# refresh + per-paper re-run regenerates taxa.json against a new
# backbone, but on an HPC filesystem the regenerated file's mtime can
# be <= the mtime we last recorded (clock skew across array-task
# nodes), so an mtime-only gate skips it forever. The fix gates on the
# taxonomy sha256 recorded in each paper's taxa-stage fingerprint.


_HASH = "abc123def456"


def _write_paper(output_dir: Path, taxa_mtime: float, taxonomy_sha: str | None) -> Path:
    """Write one paper's taxa.json + pipeline_state.json under
    documents/<hash>/, stamping the taxa stage with ``taxonomy_sha`` and
    forcing taxa.json's mtime to ``taxa_mtime``."""
    hash_dir = output_dir / "documents" / _HASH
    hash_dir.mkdir(parents=True, exist_ok=True)

    taxa_path = hash_dir / "taxa.json"
    taxa_path.write_text(json.dumps({
        "unique_taxa": 1,
        "mentions": [{
            "chunk_id": "chunk_0", "text_span": [10, 17],
            "accepted_taxon_id": "1371", "matched_text": "Nanomia",
            "accepted_name": "Nanomia", "rank": "Genus", "name_type": "binomial",
        }],
    }), encoding="utf-8")
    os.utime(taxa_path, (taxa_mtime, taxa_mtime))

    state = {"stages": {"taxa_and_lexicon_extraction": {
        "input_fingerprint": (
            {"taxonomy": {"sha256": taxonomy_sha}} if taxonomy_sha else {}
        ),
    }}}
    (hash_dir / "pipeline_state.json").write_text(json.dumps(state), encoding="utf-8")
    return taxa_path


def _write_taxonomy(output_dir: Path, content: bytes) -> str:
    """Write a stand-in taxonomy.sqlite (build only hashes it) and
    return its sha256."""
    p = output_dir / "taxonomy.sqlite"
    p.write_bytes(content)
    return _file_sha256(p)


def _build(output_dir: Path) -> dict:
    conn = sqlite3.connect(output_dir / "taxon_mentions.sqlite")
    create_schema(conn)
    stats = build(conn, output_dir)
    conn.close()
    return stats


def test_taxonomy_rebuild_reingests_despite_stale_mtime(tmp_path):
    """The exact #95 repro: rebuild taxonomy.sqlite, re-run the per-paper
    taxa stage (new recorded backbone), but with taxa.json's mtime NOT
    advanced. An mtime-only gate would skip; the fingerprint gate must
    re-ingest."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    sha_v1 = _write_taxonomy(output_dir, b"taxonomy version one")
    _write_paper(output_dir, taxa_mtime=1000.0, taxonomy_sha=sha_v1)
    stats = _build(output_dir)
    assert stats["papers"] == 1 and stats["mentions"] == 1

    # Rebuild the taxonomy and re-run the taxa stage against it, but
    # leave taxa.json's mtime at (in fact below) the recorded value.
    sha_v2 = _write_taxonomy(output_dir, b"taxonomy version TWO (different)")
    assert sha_v2 != sha_v1
    _write_paper(output_dir, taxa_mtime=999.0, taxonomy_sha=sha_v2)

    stats2 = _build(output_dir)
    assert stats2["refreshed"] == 1, "changed backbone must force a re-ingest"
    assert stats2["skipped"] == 0
    assert stats2["mentions"] == 1  # re-ingested, not duplicated

    # And the mention rows are present exactly once.
    conn = sqlite3.connect(output_dir / "taxon_mentions.sqlite")
    (n,) = conn.execute("SELECT COUNT(*) FROM taxon_mentions").fetchone()
    (stored_sha,) = conn.execute(
        "SELECT taxonomy_sha FROM papers_processed WHERE corpus_hash = ?",
        (_HASH,),
    ).fetchone()
    conn.close()
    assert n == 1
    assert stored_sha == sha_v2


def test_unchanged_backbone_skips_on_rerun(tmp_path):
    """Idempotency: same backbone + same mtime → the second build skips,
    no pointless re-ingest."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    sha = _write_taxonomy(output_dir, b"stable taxonomy")
    _write_paper(output_dir, taxa_mtime=1000.0, taxonomy_sha=sha)

    assert _build(output_dir)["papers"] == 1
    stats2 = _build(output_dir)
    assert stats2["skipped"] == 1
    assert stats2["refreshed"] == 0 and stats2["papers"] == 0


def test_lagging_backbone_flagged_without_spurious_reingest(tmp_path):
    """taxonomy.sqlite advanced but the per-paper taxa stage hasn't
    re-run (taxa.json still records the old backbone): the paper is
    flagged as lagging, and re-ingesting the still-stale taxa.json is
    NOT triggered (its recorded backbone is unchanged since our ingest)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    sha_v1 = _write_taxonomy(output_dir, b"taxonomy v1")
    _write_paper(output_dir, taxa_mtime=1000.0, taxonomy_sha=sha_v1)
    _build(output_dir)

    # Operator rebuilt taxonomy.sqlite but did NOT re-run `corpus run`,
    # so taxa.json + its fingerprint still point at v1.
    _write_taxonomy(output_dir, b"taxonomy v2 not yet propagated")
    stats2 = _build(output_dir)
    assert stats2["lagging"] == 1
    assert stats2["skipped"] == 1  # recorded backbone unchanged → no re-ingest
    assert stats2["refreshed"] == 0
