"""Source retirement is recoverable, strict and reflected in derived indexes."""
import json
import os
import sqlite3
from types import SimpleNamespace

import lancedb
import pytest

from bib import authority
from bib.documents import work_map
from mcpsrv.bundle import package
from pipeline import embed, io, taxon_mentions


def source(root, name, content):
    root.mkdir(exist_ok=True, parents=True)
    pdf = root / name
    pdf.write_bytes(content)
    return pdf, io.short_hash(io.calculate_pdf_hash(pdf))


def artifact(root, sha):
    hd = root / "documents" / sha
    hd.mkdir(parents=True, exist_ok=True)
    for name, value in {
        "summary": {"relative_paths": [sha + ".pdf"]},
        "metadata": {"title": sha, "year": 2000, "authors": [{"surname": "Author"}]},
        "references": {"references": [{"xml_id": "b0", "title": "Shared source", "year": 1990}]},
        "chunks": {"chunks": [{"chunk_id": "c0", "text": "Taxon evidence"}]},
        "taxa": {"unique_taxa": 1, "mentions": [{"accepted_taxon_id": "1", "chunk_id": "c0", "text_span": [0, 5], "matched_text": "Taxon"}]},
    }.items():
        (hd / f"{name}.json").write_text(json.dumps(value))
    return hd


def indexes(root):
    b = sqlite3.connect(root / "biblio_authority.sqlite")
    authority.create_schema(b)
    authority.phase1_corpus_papers(b, root)
    authority.phase2_references(b, root)
    t = sqlite3.connect(root / "taxon_mentions.sqlite")
    taxon_mentions.create_schema(t)
    stats = taxon_mentions.build(t, root)
    return b, t, stats


def test_retirement_matches_clean_indexes_and_vectors(tmp_path):
    library, root = tmp_path / "library", tmp_path / "build"
    a, sha = source(library, "a.pdf", b"source a")
    _, other = source(library, "b.pdf", b"source b")
    first, second = artifact(root, sha), artifact(root, other)
    b, t, _ = indexes(root)
    b.close()
    t.close()
    model = embed.make_chunk_model(2)
    backend = SimpleNamespace(model_name="test", dim=2, embed=lambda texts: [[1., 2.] for _ in texts])
    db = lancedb.connect(str(root / "vector_db" / "lancedb"))
    table = db.create_table("document_chunks", schema=model.to_arrow_schema())
    for hd in (first, second):
        embed.embed_document(hd, table, model, backend)
    served = root / "_serve"
    package(root, served, "test", False, False)
    a.unlink()
    before = set(root.rglob("*"))
    assert io.prune_orphans(library, root, force=True, dry_run=True)["doc_pruned"] == 1
    assert set(root.rglob("*")) == before
    result = io.prune_orphans(library, root, force=True)
    assert result["doc_pruned"] == result["vec_pruned"] == 1
    assert not first.exists()
    assert len(list((root / ".retired").glob(f"documents-*/{sha}/references.json"))) == 1
    b, t, stats = indexes(root)
    assert stats["retired"] == 1
    assert set(work_map(b)) == {other}
    assert b.execute("SELECT COUNT(*) FROM reference_observations").fetchone()[0] == 2
    clean = tmp_path / "clean"
    clean_hd = artifact(clean, other)
    cb, ct, _ = indexes(clean)
    for left, right, sql in (
        (b, cb, "SELECT citing_corpus_hash,citing_work_id,cited_work_id FROM citations ORDER BY 1,2,3"),
        (b, cb, "SELECT observation_id,work_id,match_method FROM observation_work ORDER BY 1"),
        (t, ct, "SELECT corpus_hash,taxon_id,chunk_id,mention_text FROM taxon_mentions ORDER BY 1,2,3"),
    ):
        assert left.execute(sql).fetchall() == right.execute(sql).fetchall()
    clean_table = lancedb.connect(str(clean / "vector_db" / "lancedb")).create_table("document_chunks", schema=model.to_arrow_schema())
    embed.embed_document(clean_hd, clean_table, model, backend)
    table.checkout_latest()
    def logical_rows(tbl):
        return [{k: v for k, v in row.items() if k != "embedding_generation"} for row in tbl.to_arrow().to_pylist()]
    assert logical_rows(table) == logical_rows(clean_table)
    package(root, served, "test", False, False)
    clean_served = clean / "_serve"
    package(clean, clean_served, "test", False, False)
    assert {p.name for p in (served / "documents").iterdir()} == {other}
    assert {p.relative_to(served / "documents") for p in (served / "documents").rglob("*")} == {
        p.relative_to(clean_served / "documents") for p in (clean_served / "documents").rglob("*")}
    for conn in (b, t, cb, ct):
        conn.close()


def test_nested_output_is_not_part_of_source_inventory(tmp_path):
    _, sha = source(tmp_path, "a.pdf", b"original")
    root = tmp_path / "output"
    hd = artifact(root, sha)
    (hd / "processed.pdf").write_bytes(b"original")
    (tmp_path / "a.pdf").unlink()
    assert io.prune_orphans(tmp_path, root, force=True)["doc_pruned"] == 1
    assert not io.find_all_pdfs(tmp_path, exclude_under=root, strict=True)


def test_incomplete_inventory_never_retires_even_when_forced(tmp_path, monkeypatch):
    library = tmp_path / "library"
    _, sha = source(library, "a.pdf", b"a")
    hd = artifact(tmp_path / "build", sha)
    def fail(path):
        raise PermissionError("unreadable source")
    monkeypatch.setattr(io, "calculate_pdf_hash", fail)
    with pytest.raises(RuntimeError, match="Cannot hash"):
        io.prune_orphans(library, tmp_path / "build", force=True)
    assert hd.exists()
    with pytest.raises(RuntimeError, match="missing input"):
        io.prune_orphans(tmp_path / "missing", tmp_path / "build", force=True)


def test_failed_vector_prune_aborts_before_retirement(tmp_path, monkeypatch):
    library = tmp_path / "library"
    library.mkdir()
    root = tmp_path / "build"
    hd = artifact(root, "abc")
    (root / "vector_db" / "lancedb").mkdir(parents=True)
    def fail(*args, **kwargs):
        raise RuntimeError("storage unavailable")
    monkeypatch.setattr(lancedb, "connect", fail)
    with pytest.raises(RuntimeError, match="update aborted"):
        io.prune_orphans(library, root, force=True)
    assert hd.exists()


def test_directory_scan_error_is_not_an_empty_inventory(tmp_path, monkeypatch):
    def unreadable(path, *, onerror):
        onerror(PermissionError("unreadable subtree"))
        return iter(())
    monkeypatch.setattr(io.os, "walk", unreadable)
    with pytest.raises(RuntimeError, match="Cannot inventory all"):
        io.find_all_pdfs(tmp_path, strict=True)


def test_taxon_content_change_with_same_mtime_and_corrupt_replacement(tmp_path):
    hd = artifact(tmp_path, "abc")
    b, t, _ = indexes(tmp_path)
    path = hd / "taxa.json"
    stamp = path.stat().st_mtime_ns
    path.write_text(json.dumps({"unique_taxa": 0, "mentions": []}))
    os.utime(path, ns=(stamp, stamp))
    assert taxon_mentions.build(t, tmp_path)["refreshed"] == 1
    assert t.execute("SELECT COUNT(*) FROM taxon_mentions").fetchone()[0] == 0
    assert taxon_mentions.build(t, tmp_path)["skipped"] == 1
    path.write_text("{")
    before = t.total_changes
    assert taxon_mentions.build(t, tmp_path)["errors"] == 1
    assert t.total_changes == before
    b.close()
    t.close()


def test_removed_legacy_scalar_membership_is_retired(tmp_path):
    (tmp_path / "documents").mkdir()
    conn = sqlite3.connect(":memory:")
    authority.create_schema(conn)
    authority.insert_work(conn, "legacy", "corpus_key", "Old", 2000, "", None, "gone", True, "corpus_paper")
    authority.phase1_corpus_papers(conn, tmp_path)
    assert not work_map(conn)
    conn.close()


def test_bundle_replacement_removes_old_optional_files_and_retains_previous(tmp_path):
    root = tmp_path / "build"
    hd = artifact(root, "abc")
    (hd / "processed.pdf").write_bytes(b"PDF")
    (root / "instructions.md").write_text("Previous guidance")
    served = root / "_serve"
    package(root, served, "test", True, False)
    previous_manifest = (served / "bundle_manifest.json").read_bytes()
    (root / "instructions.md").unlink()
    package(root, served, "test", False, False)
    assert not (served / "instructions.md").exists()
    assert not (served / "documents" / "abc" / "processed.pdf").exists()
    backups = list(root.glob("._serve-previous-*/bundle/bundle_manifest.json"))
    assert len(backups) == 1 and backups[0].read_bytes() == previous_manifest


def test_malformed_served_json_cannot_bypass_audit_or_replace_bundle(tmp_path):
    root = tmp_path / "build"
    hd = artifact(root, "abc")
    served = root / "_serve"
    package(root, served, "good", False, False)
    previous = (served / "bundle_manifest.json").read_bytes()
    (hd / "references.json").write_text('{"private_path":"/home/private/source.pdf",')
    with pytest.raises(ValueError, match="Cannot audit served JSON"):
        package(root, served, "broken", False, False)
    assert (served / "bundle_manifest.json").read_bytes() == previous
    assert json.loads((served / "documents/abc/references.json").read_text())["references"]


def test_bundle_copy_failure_preserves_previous_generation(tmp_path, monkeypatch):
    from mcpsrv import bundle
    root = tmp_path / "build"
    artifact(root, "abc")
    served = root / "_serve"
    package(root, served, "test", False, False)
    before = {p.relative_to(served): p.read_bytes() for p in served.rglob("*") if p.is_file()}
    def fail(*args, **kwargs):
        raise OSError("disk full")
    monkeypatch.setattr(bundle, "_copy_file", fail)
    with pytest.raises(OSError, match="disk full"):
        package(root, served, "test", False, False)
    assert {p.relative_to(served): p.read_bytes() for p in served.rglob("*") if p.is_file()} == before
