"""Real materialization → bundle → MCP formatter acceptance (#296, #301).

The prepared boundary is metadata from a pinned library BibTeX fixture, not
PDF extraction. Everything from authority schema creation through the stdio
tool call is production code; no models, network service or stage mocks.
"""
from __future__ import annotations

import asyncio
from contextlib import closing
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import sys

from mcp import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client

from bib.authority import create_schema, phase1_corpus_papers
from bib.documents import document_metadata, find_work
from bib.export import export_bibtex, render_entry
from bib.fields import LOCATOR_FIELDS
from bib.importer import import_bibtex
from bib.parser import bib_entry_to_metadata, parse_bibtex
from mcpsrv.bundle import package


FIXTURE = Path(__file__).parent / "fixtures/bibliographic_integrity"
REPO = Path(__file__).resolve().parents[1]


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _inventory(root):
    return {str(path.relative_to(root)): _sha(path)
            for path in root.rglob("*") if path.is_file()}


def _decode(result):
    assert not result.is_error, result
    texts = [block.text for block in result.content if block.type == "text"]
    assert len(texts) == 1, result
    payload = json.loads(texts[0])
    assert "error" not in payload, payload
    return payload


async def _query_bundle(bundle, hashes, stderr):
    # A separate production server process must discover the packaged DB;
    # passing the build DB directly to BiblioAuthority would miss this boundary.
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "mcpsrv.main", str(bundle), "--transport", "stdio"],
        cwd=REPO,
        env={**os.environ, "HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1"},
    )
    async with stdio_client(params, errlog=stderr) as streams:
        async with ClientSession(*streams, read_timeout_seconds=30) as session:
            initialized = await session.initialize()
            tools = await session.list_tools()
            assert "format_citations" in {tool.name for tool in tools.tools}
            by_hash = _decode(await session.call_tool(
                "format_citations", {"paper_hashes": hashes, "style": "author-year"},
            ))
            assert by_hash["count"] == len(hashes)
            assert all("error" not in row for row in by_hash["citations"])
            by_work = _decode(await session.call_tool(
                "format_citations", {"work_ids": [row["work_id"] for row in by_hash["citations"]]},
            ))
            assert by_work == by_hash
            return {"server_name": initialized.server_info.name, "by_hash": by_hash,
                    "by_work": by_work}


def _package_and_query(build, bundle, entries, label, receipts):
    manifest = package(build, bundle, version="acceptance-301", include_pdfs=False, dry_run=False)
    assert manifest["paper_count"] == len(entries)
    assert manifest["embedding_model"] is None
    assert manifest["chunk_count"] == manifest["figure_count"] == 0
    assert _sha(build / "biblio_authority.sqlite") == _sha(bundle / "biblio_authority.sqlite")
    before = _inventory(bundle)
    with (bundle.parent / f"{label}-server.stderr").open("w", encoding="utf-8") as stderr:
        response = asyncio.run(asyncio.wait_for(
            _query_bundle(bundle, list(entries), stderr), timeout=60,
        ))
    assert _inventory(bundle) == before, "Formatting must not modify the served bundle"
    receipts[label] = {"manifest": manifest, "bundle_files_sha256": before, **response}
    return dict(zip(entries, response["by_hash"]["citations"]))


def _assert_source_values(citations, entries, provenance, bundle):
    for sha, entry in entries.items():
        citation = citations[sha]
        assert citation["corpus_hash"] == sha
        assert citation["provenance"] == "bib"
        assert not citation["warning"]
        fields = citation["fields"]
        assert fields["title"] == entry["title"]
        for key in LOCATOR_FIELDS:
            assert fields[key] == entry.get(key), (entry["_key"], key, fields[key])
        assert "keeppages" not in fields
        metadata = json.loads((bundle / "documents" / sha / "metadata.json").read_text())
        assert metadata["bib_key"] == entry["_key"]
        assert metadata["extraction_method"] == "bib"
        for key in LOCATOR_FIELDS:
            assert metadata[key] == entry.get(key)

    def item(key):
        return citations[provenance["entries"][key]["pdf_hash"]]

    church = item("Churchetal2015")
    assert "*324*(5), 435–449" in church["formatted"]
    # Preserve the supplied catalogue value; this is not a source-page
    # adjudication of the library's known Siebert pagination discrepancy.
    assert "*3702*(3), 202–232" in item("Siebertetal2013")["formatted"]
    article = item("Ahujaetal2024")
    assert "*16*(3), evae048" in article["formatted"]
    assert article["fields"]["eid"] is article["fields"]["articleno"] is None
    parts = [item(key) for key in provenance["entries"] if key.startswith("delleChiaje")]
    assert len(parts) == len({part["work_id"] for part in parts}) == 7
    assert all(part["fields"]["volume"] is None for part in parts)
    assert all(part["fields"]["title"] in part["formatted"] for part in parts)
    chun = item("Chun1898b")
    assert "309–313" in chun["formatted"]
    assert "2–6" not in chun["formatted"]
    chun_meta = json.loads((bundle / "documents" / chun["corpus_hash"] / "metadata.json").read_text())
    assert chun_meta["keeppages"] == "2--6"
    manko = item("MankoPugh2018")
    assert [(a["surname"], a["forename"]) for a in manko["fields"]["authors"]] == [
        ("Mańko", "M.K."), ("Pugh", "P.R."),
    ]
    assert manko["formatted"].startswith("Mańko, M. K., & Pugh, P. R.")


def _write_entry(path, entry):
    path.write_text(render_entry(entry["_key"], entry.get("_type", "article"),
                                {key: value for key, value in entry.items() if not key.startswith("_")}),
                    encoding="utf-8")


def test_source_locators_survive_real_bundle_live_formatter_and_unchanged_refresh(tmp_path):
    provenance = json.loads((FIXTURE / "provenance.json").read_text())
    source = {entry["_key"]: entry for entry in parse_bibtex((FIXTURE / "source.bib").read_text())}
    assert source.keys() == provenance["entries"].keys()
    entries = {provenance["entries"][key]["pdf_hash"]: entry for key, entry in source.items()}
    build, bundle = tmp_path / "build", tmp_path / "served"
    for sha, entry in entries.items():
        directory = build / "documents" / sha
        directory.mkdir(parents=True)
        metadata = bib_entry_to_metadata(entry, provenance["entries"][entry["_key"]]["filename"])
        (directory / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False), encoding="utf-8")
    db = build / "biblio_authority.sqlite"
    with closing(sqlite3.connect(db)) as conn, conn:
        create_schema(conn)
        assert phase1_corpus_papers(conn, build) == len(entries)
        assert set(LOCATOR_FIELDS) <= {row[1] for row in conn.execute("PRAGMA table_info(works)")}
        for sha, entry in entries.items():
            metadata = document_metadata(conn, sha)
            assert all(metadata[key] == entry.get(key) for key in LOCATOR_FIELDS)
        manko_id = find_work(conn, provenance["entries"]["MankoPugh2018"]["pdf_hash"])
        assert conn.execute("SELECT bib_source FROM works WHERE work_id=?", (manko_id,)).fetchone() == ("metadata",)
        assert conn.execute("SELECT DISTINCT origin FROM work_bib_sources WHERE work_id=?", (manko_id,)).fetchall() == [("metadata",)]
    receipts = {"fixture_provenance": provenance, "fixture_sha256": _sha(FIXTURE / "source.bib"),
                "scope": "Prepared BibTeX metadata; real authority, bundle, export/import and live stdio formatter. No PDF extraction or deployed-bundle history."}
    fresh = _package_and_query(build, bundle, entries, "fresh", receipts)
    _assert_source_values(fresh, entries, provenance, bundle)
    assert fresh[provenance["entries"]["MankoPugh2018"]["pdf_hash"]]["bib_key"] == "MankoPugh2018"

    # Deliberately wrong values are a synthetic precedence control, never
    # attributed to the real Church source. An explicit import first beats
    # metadata; importing the actual source then wins by stable source ID.
    church_sha = provenance["entries"]["Churchetal2015"]["pdf_hash"]
    conflicting = dict(source["Churchetal2015"], _key="ZDeliberateLocatorConflict",
                       volume="999", number="99", pages="1--2")
    path = tmp_path / "synthetic-conflict.bib"
    _write_entry(path, conflicting)
    receipts["synthetic_conflict_import"] = import_bibtex(db, path)
    with closing(sqlite3.connect(db)) as conn, conn:
        church_id = find_work(conn, church_sha)
        assert conn.execute("SELECT volume,number,pages FROM works WHERE work_id=?", (church_id,)).fetchone() == ("999", "99", "1--2")
    path = tmp_path / "restored-source.bib"
    _write_entry(path, source["Churchetal2015"])
    receipts["source_restoring_import"] = import_bibtex(db, path)
    restored = _package_and_query(build, bundle, entries, "source_precedence", receipts)
    _assert_source_values(restored, entries, provenance, bundle)
    conflicts = {row["field"]: row for row in restored[church_sha]["bibliographic_conflicts"]}
    for field in ("volume", "number", "pages"):
        assert conflicts[field]["policy"] == "explicit_import_then_source_id"
        assert conflicts[field]["selected_source"] == "import:Churchetal2015"
        assert {row["source_id"]: row["value"] for row in conflicts[field]["sources"]} == {
            "import:Churchetal2015": source["Churchetal2015"][field],
            "import:ZDeliberateLocatorConflict": conflicting[field],
        }

    exported = tmp_path / "exported.bib"
    exported.write_text(export_bibtex(db), encoding="utf-8")
    exported_by_hash = {entry["corpus_hash"]: entry for entry in parse_bibtex(exported.read_text())}
    assert exported_by_hash.keys() == entries.keys()
    for sha, entry in entries.items():
        assert all(exported_by_hash[sha].get(key) == entry.get(key) for key in LOCATOR_FIELDS)
    assert exported_by_hash[provenance["entries"]["Chun1898b"]["pdf_hash"]]["keeppages"] == "2--6"
    receipts["roundtrip_imports"] = [import_bibtex(db, exported), import_bibtex(db, exported)]
    assert all(result["no_changes"] == len(entries) and result["no_match"] == 0
               for result in receipts["roundtrip_imports"])
    with closing(sqlite3.connect(db)) as conn, conn:
        before = list(conn.iterdump())
        assert phase1_corpus_papers(conn, build) == 0
        assert list(conn.iterdump()) == before
    unchanged = _package_and_query(build, bundle, entries, "unchanged_refresh", receipts)
    _assert_source_values(unchanged, entries, provenance, bundle)
    # Export intentionally generates cite keys. Explicitly reimporting those
    # keys can select a different key/source ID without changing field values;
    # the packaged document metadata still retains its original library key.
    manko_sha = provenance["entries"]["MankoPugh2018"]["pdf_hash"]
    assert unchanged[manko_sha]["bib_key"] == exported_by_hash[manko_sha]["_key"]
    for sha in entries:
        for key in ("work_id", "fields", "formatted", "inline", "provenance", "warning"):
            assert unchanged[sha][key] == fresh[sha][key]
    unchanged_conflicts = unchanged[church_sha]["bibliographic_conflicts"]
    assert {row["field"] for row in unchanged_conflicts} == {"volume", "number", "pages"}
    for conflict in unchanged_conflicts:
        field = conflict["field"]
        assert field in ("volume", "number", "pages")
        values = {row["source_id"]: row["value"] for row in conflict["sources"]}
        assert values["import:Churchetal2015"] == source["Churchetal2015"][field]
        assert values["import:ZDeliberateLocatorConflict"] == conflicting[field]
        assert values[conflict["selected_source"]] == source["Churchetal2015"][field]
    (tmp_path / "roundtrip-acceptance.json").write_text(
        json.dumps(receipts, ensure_ascii=False, indent=2) + "\n", encoding="utf-8",
    )
