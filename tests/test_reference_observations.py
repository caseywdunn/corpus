"""Reference evidence and deterministic materialization tests (#240)."""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest

from bib import authority


def _write_paper(
    output_dir: Path,
    corpus_hash: str,
    *,
    title: str,
    year: int,
    surname: str,
    references: list[dict],
    doi: str = "",
) -> None:
    doc_dir = output_dir / "documents" / corpus_hash
    doc_dir.mkdir(parents=True, exist_ok=True)
    (doc_dir / "metadata.json").write_text(json.dumps({
        "title": title,
        "year": year,
        "journal": "Test Journal",
        "doi": doi,
        "authors": [{"surname": surname, "forename": "A"}],
    }))
    (doc_dir / "references.json").write_text(json.dumps({
        "references": references,
        "total_references": len(references),
    }))


def _reference(raw: str, title: str, xml_id: str) -> dict:
    return {
        "raw": raw,
        "title": title,
        "year": 1888,
        "journal": "Challenger Reports",
        "authors": ["E Haeckel"],
        "xml_id": xml_id,
    }


def _open_and_build(output_dir: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(":memory:")
    authority.create_schema(conn)
    authority.phase1_corpus_papers(conn, output_dir)
    authority.phase2_references(conn, output_dir)
    return conn


def _current_mappings(conn: sqlite3.Connection) -> list[tuple]:
    return conn.execute(
        """SELECT ro.observation_id, ow.work_id, ow.match_method,
                  ow.match_score, ow.producer_version
           FROM reference_current_sets current
           JOIN reference_observation_memberships member
             ON member.corpus_hash = current.corpus_hash
            AND member.source_fingerprint = current.source_fingerprint
           JOIN reference_observations ro
             ON ro.observation_id = member.observation_id
           JOIN observation_work ow
             ON ow.observation_id = ro.observation_id
           ORDER BY ro.observation_id"""
    ).fetchall()


def _snapshot(conn: sqlite3.Connection) -> dict[str, list[tuple]]:
    tables = (
        "works",
        "work_authors",
        "work_aliases",
        "citations",
        "reference_observations",
        "reference_observation_sets",
        "reference_observation_memberships",
        "reference_current_sets",
        "observation_work",
    )
    return {
        table: conn.execute(f"SELECT * FROM {table} ORDER BY 1").fetchall()
        for table in tables
    }


@pytest.mark.parametrize("artifact,broken", [
    ("references", "{"), ("references", "[]"),
    ("references", '{"references": [null]}'),
    ("references", '{"references": {}}'),
    ("references", '{"references": [{"authors": "Author"}]}'),
    ("metadata", "{"), ("metadata", "[]"),
    ("metadata", '{"authors": ["Author"]}'),
    ("metadata", '{"authors": {}}'),
])
def test_malformed_source_aborts_before_any_current_evidence_changes(tmp_path, artifact, broken):
    for sha in ("aaa", "zzz"):
        _write_paper(tmp_path, sha, title="A paper", year=2000, surname=sha,
                     references=[_reference("Raw", "Title", "b0")])
    conn = _open_and_build(tmp_path)
    before = _snapshot(conn)
    # An earlier valid changed source must not be partially ingested before
    # discovering the later invalid artifact.
    _write_paper(tmp_path, "aaa", title="Edited paper", year=2001, surname="Edited",
                 references=[_reference("New raw", "New title", "b0")])
    (tmp_path / "documents/zzz" / f"{artifact}.json").write_text(broken)
    phase = authority.phase1_corpus_papers if artifact == "metadata" else authority.phase2_references
    with pytest.raises(ValueError, match="Cannot ingest"):
        phase(conn, tmp_path)
    assert _snapshot(conn) == before
    conn.close()


def test_unreadable_reference_is_not_a_successful_cached_noop(tmp_path, monkeypatch):
    _write_paper(tmp_path, "aaa", title="Paper", year=2000, surname="Author", references=[])
    conn = _open_and_build(tmp_path)
    before = _snapshot(conn)
    original = Path.read_bytes
    def read(path):
        if path.name == "references.json":
            raise PermissionError("unreadable input")
        return original(path)
    monkeypatch.setattr(Path, "read_bytes", read)
    with pytest.raises(ValueError, match="unreadable input"):
        authority.phase2_references(conn, tmp_path)
    assert _snapshot(conn) == before
    conn.close()


def test_enabling_bhl_on_unchanged_sources_applies_the_requested_policy(tmp_path, monkeypatch):
    _write_paper(tmp_path, "aaa", title="Paper", year=2000, surname="Author",
                 references=[_reference("Raw", "An old publication", "b0")])
    conn = _open_and_build(tmp_path)
    calls = []
    def lookup(*args, **kwargs):
        calls.append(kwargs)
        return "found", ("bhl:part/123", "", "123"), None
    monkeypatch.setattr(authority, "_bhl_lookup", lookup)
    assert authority.phase2_references(conn, tmp_path, enrich_bhl=True, bhl_api_key="test-key")[0] == 1
    assert len(calls) == 1
    assert conn.execute("SELECT work_id FROM observation_work").fetchone()[0] == "bhl:part/123"
    before = conn.total_changes
    assert authority.phase2_references(conn, tmp_path, enrich_bhl=True, bhl_api_key="test-key") == (0, 0)
    assert conn.total_changes == before and len(calls) == 1
    # Disabling future enrichment does not erase already acquired BHL evidence.
    authority.phase2_references(conn, tmp_path, enrich_bhl=False)
    assert len(calls) == 1
    assert conn.execute("SELECT work_id FROM observation_work").fetchone()[0] == "bhl:part/123"
    conn.close()


def test_bhl_year_window_change_revisits_previously_ineligible_references(tmp_path, monkeypatch):
    _write_paper(tmp_path, "aaa", title="Paper", year=2000, surname="Author",
                 references=[_reference("Raw", "An old publication", "b0")])
    conn = _open_and_build(tmp_path)
    calls = []
    def lookup(*args, **kwargs):
        calls.append(kwargs["max_year"])
        return "skipped", None, "outside year window"
    monkeypatch.setattr(authority, "_bhl_lookup", lookup)
    authority.phase2_references(conn, tmp_path, enrich_bhl=True, bhl_api_key="test-key", bhl_max_year=1800)
    authority.phase2_references(conn, tmp_path, enrich_bhl=True, bhl_api_key="test-key", bhl_max_year=1900)
    assert calls == [1800, 1900]
    conn.close()


def test_bhl_key_becoming_available_invalidates_without_storing_the_key(tmp_path, monkeypatch):
    _write_paper(tmp_path, "aaa", title="Paper", year=2000, surname="Author",
                 references=[_reference("Raw", "An old publication", "b0")])
    conn = _open_and_build(tmp_path)
    calls = []
    def lookup(*args, **kwargs):
        calls.append(bool(kwargs["api_key"]))
        return "skipped", None, "unavailable test service"
    monkeypatch.setattr(authority, "_bhl_lookup", lookup)
    authority.phase2_references(conn, tmp_path, enrich_bhl=True)
    authority.phase2_references(conn, tmp_path, enrich_bhl=True, bhl_api_key="never-store-this-secret")
    assert calls == [False, True]
    before = conn.total_changes
    authority.phase2_references(conn, tmp_path, enrich_bhl=True, bhl_api_key="a-replacement-valid-key")
    assert conn.total_changes == before and calls == [False, True]
    assert "never-store-this-secret" not in "\n".join(conn.iterdump())
    conn.close()


def test_reference_evidence_is_separate_auditable_and_noop(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "output"
    ref = _reference(
        "Haeckel E. 1888. Report on the Siphonophorae.",
        "Report on the Siphonophorae",
        "b0",
    )
    _write_paper(
        output_dir, "aaa", title="A citing paper", year=2001,
        surname="Alpha", references=[ref],
    )
    conn = _open_and_build(output_dir)

    observation = conn.execute(
        """SELECT citing_corpus_hash, ordinal, grobid_xml_id, raw_citation,
                  title, year, journal, authors_json
           FROM reference_observations"""
    ).fetchone()
    assert observation == (
        "aaa", 0, "b0", ref["raw"], ref["title"], 1888,
        "Challenger Reports", '["E Haeckel"]',
    )
    mapping = conn.execute(
        """SELECT ow.match_method, ow.match_score, ow.producer_version,
                  ro.raw_citation
           FROM observation_work ow
           JOIN reference_observations ro USING (observation_id)"""
    ).fetchone()
    assert mapping[0] == "new_work"
    assert mapping[1] == 1.0
    assert mapping[2] == authority.REFERENCE_MAPPING_PRODUCER
    assert mapping[3] == ref["raw"]
    assert conn.execute("SELECT COUNT(*) FROM citations").fetchone()[0] == 1

    before = _snapshot(conn)
    assert authority.phase1_corpus_papers(conn, output_dir) == 0
    assert authority.phase2_references(conn, output_dir) == (0, 0)
    assert _snapshot(conn) == before
    conn.close()


def test_changed_source_keeps_superseded_observation_evidence(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "output"
    old_ref = _reference("Old raw citation", "Old parsed title", "b0")
    _write_paper(
        output_dir, "aaa", title="A citing paper", year=2001,
        surname="Alpha", references=[old_ref],
    )
    conn = _open_and_build(output_dir)
    old_observation_id = conn.execute(
        "SELECT observation_id FROM reference_observations"
    ).fetchone()[0]

    new_ref = _reference("Corrected raw citation", "Corrected title", "b0")
    _write_paper(
        output_dir, "aaa", title="A citing paper", year=2001,
        surname="Alpha", references=[new_ref],
    )
    mapped, _new_works = authority.phase2_references(conn, output_dir)
    assert mapped == 1

    # Both source observations and both source-set memberships survive. Only
    # the selected set and its observation have a current mapping.
    assert conn.execute(
        "SELECT COUNT(*) FROM reference_observations"
    ).fetchone()[0] == 2
    assert conn.execute(
        "SELECT raw_citation FROM reference_observations "
        "WHERE observation_id = ?", (old_observation_id,),
    ).fetchone()[0] == "Old raw citation"
    assert conn.execute(
        "SELECT COUNT(*) FROM reference_observation_sets"
    ).fetchone()[0] == 2
    assert conn.execute(
        "SELECT COUNT(*) FROM reference_observation_memberships"
    ).fetchone()[0] == 2
    assert len(_current_mappings(conn)) == 1
    current_raw = conn.execute(
        """SELECT ro.raw_citation
           FROM observation_work ow
           JOIN reference_observations ro USING (observation_id)"""
    ).fetchone()[0]
    assert current_raw == "Corrected raw citation"
    conn.close()


def test_clean_and_incremental_builds_have_identical_current_mapping(
    tmp_path: Path,
) -> None:
    long_ref = _reference(
        "Haeckel 1888. Full report.",
        "Report on the Siphonophorae collected by HMS Challenger",
        "b0",
    )
    short_ref = _reference(
        "Haeckel 1888. Report.", "Report on the Siphonophorae", "b7",
    )

    clean_dir = tmp_path / "clean"
    _write_paper(
        clean_dir, "aaa", title="First paper", year=2001,
        surname="Alpha", references=[long_ref],
    )
    _write_paper(
        clean_dir, "bbb", title="Second paper", year=2002,
        surname="Beta", references=[short_ref],
    )
    clean = _open_and_build(clean_dir)

    incremental_dir = tmp_path / "incremental"
    _write_paper(
        incremental_dir, "aaa", title="First paper", year=2001,
        surname="Alpha", references=[long_ref],
    )
    incremental = _open_and_build(incremental_dir)
    _write_paper(
        incremental_dir, "bbb", title="Second paper", year=2002,
        surname="Beta", references=[short_ref],
    )
    assert authority.phase1_corpus_papers(incremental, incremental_dir) == 1
    assert authority.phase2_references(incremental, incremental_dir)[0] == 2

    assert _current_mappings(incremental) == _current_mappings(clean)
    assert incremental.execute(
        """SELECT citing_work_id, cited_work_id, citing_corpus_hash,
                  grobid_xml_id, raw_citation, match_method, match_score
           FROM citations ORDER BY 1, 2, 3"""
    ).fetchall() == clean.execute(
        """SELECT citing_work_id, cited_work_id, citing_corpus_hash,
                  grobid_xml_id, raw_citation, match_method, match_score
           FROM citations ORDER BY 1, 2, 3"""
    ).fetchall()
    clean.close()
    incremental.close()


def test_existing_artifact_stamp_does_not_skip_observation_migration(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "output"
    _write_paper(
        output_dir, "aaa", title="A citing paper", year=2001,
        surname="Alpha", references=[
            _reference("Raw", "Report on the Siphonophorae", "b0")
        ],
    )
    conn = sqlite3.connect(":memory:")
    authority.create_schema(conn)
    authority.phase1_corpus_papers(conn, output_dir)
    refs_path = output_dir / "documents" / "aaa" / "references.json"
    authority._record_artifact(
        conn, "aaa", "references", refs_path.stat().st_mtime,
    )
    conn.commit()

    mapped, _new_works = authority.phase2_references(conn, output_dir)
    assert mapped == 1
    assert conn.execute(
        "SELECT COUNT(*) FROM reference_observations"
    ).fetchone()[0] == 1
    conn.close()


def test_empty_references_get_distinct_deterministic_work_ids(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "output"
    for corpus_hash, surname in (("aaa", "Alpha"), ("bbb", "Beta")):
        _write_paper(
            output_dir, corpus_hash, title=f"Paper {surname}", year=2001,
            surname=surname, references=[{"xml_id": "b0"}],
        )
    first = _open_and_build(output_dir)
    first_mapping = _current_mappings(first)
    assert len({row[1] for row in first_mapping}) == 2

    second = _open_and_build(output_dir)
    assert _current_mappings(second) == first_mapping
    first.close()
    second.close()


def test_old_mapping_producer_forces_complete_rematerialization(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "output"
    _write_paper(
        output_dir, "aaa", title="A citing paper", year=2001,
        surname="Alpha", references=[
            _reference("Raw", "Report on the Siphonophorae", "b0")
        ],
    )
    conn = _open_and_build(output_dir)
    conn.execute("UPDATE observation_work SET producer_version = 'obsolete'")
    conn.commit()

    mapped, _new_works = authority.phase2_references(conn, output_dir)

    assert mapped == 1
    assert conn.execute(
        "SELECT DISTINCT producer_version FROM observation_work"
    ).fetchall() == [(authority.REFERENCE_MAPPING_PRODUCER,)]
    conn.close()


def test_multiple_documents_share_a_work_without_losing_observations(
    tmp_path: Path,
) -> None:
    """Both volumes map references; repeating the mapping is a true no-op."""
    output_dir = tmp_path / "output"
    shared_doi = "10.6/bhl.shared-volume"
    for corpus_hash, title in (("aaa", "Volume one"), ("bbb", "Volume two")):
        _write_paper(
            output_dir, corpus_hash, title=title, year=1841,
            surname="Author", doi=shared_doi,
            references=[
                _reference("Raw", "Report on the Siphonophorae", "b0")
            ],
        )
    conn = sqlite3.connect(":memory:")
    authority.create_schema(conn)
    assert authority.phase1_corpus_papers(conn, output_dir) == 1
    assert authority.phase2_references(conn, output_dir)[0] == 2
    assert conn.execute(
        "SELECT COUNT(*) FROM reference_observations"
    ).fetchone()[0] == 2
    before = _snapshot(conn)

    assert authority.phase2_references(conn, output_dir) == (0, 0)
    assert _snapshot(conn) == before
    conn.close()
