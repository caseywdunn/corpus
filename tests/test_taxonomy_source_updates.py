"""Taxonomy snapshots own source receipts and replace, rather than accumulate."""
import argparse
import sqlite3
import zipfile

import pytest

from pipeline import taxonomy_ingest as ti
from pipeline import orchestrator


def source(path, rows):
    path.write_text("taxonID\tscientificName\tparentNameUsageID\n" + "\n".join(rows) + "\n")


def run(monkeypatch, root, src, *extra):
    monkeypatch.setattr("sys.argv", ["ingest", str(root), "--source", "dwc", "--input", str(src), *extra])
    return ti.main()


def logical(path):
    conn = sqlite3.connect(path)
    result = (conn.execute("SELECT taxon_id,scientific_name,parent_name_usage_id FROM taxa ORDER BY taxon_id").fetchall(),
              conn.execute("SELECT name,name_lowercase,taxon_id,name_type FROM names ORDER BY 1,2,3,4").fetchall())
    conn.close()
    return result


def test_taxonomy_edit_removes_old_taxa_and_names_and_matches_clean(tmp_path, monkeypatch):
    src = tmp_path / "Taxon.tsv"
    root = tmp_path / "build"
    source(src, ["1\tOld name\t", "2\tRemoved name\t1"])
    assert run(monkeypatch, root, src) == 0
    old = (root / "taxonomy.sqlite").read_bytes()
    assert run(monkeypatch, root, src) == 0
    assert (root / "taxonomy.sqlite").read_bytes() == old
    source(src, ["1\tCorrect name\t", "3\tNew name\t1"])
    assert run(monkeypatch, root, src) == 0
    clean = tmp_path / "clean"
    assert run(monkeypatch, clean, src) == 0
    assert logical(root / "taxonomy.sqlite") == logical(clean / "taxonomy.sqlite")
    backups = list((root / ".retired").glob("taxonomy-*.sqlite"))
    assert len(backups) == 1 and backups[0].read_bytes() == old


def test_failed_or_empty_replacement_preserves_previous_snapshot(tmp_path, monkeypatch):
    src = tmp_path / "Taxon.tsv"
    source(src, ["1\tValid taxon\t"])
    run(monkeypatch, tmp_path / "build", src)
    db = tmp_path / "build" / "taxonomy.sqlite"
    before = db.read_bytes()
    src.write_text("not valid headers")
    with pytest.raises(ValueError):
        run(monkeypatch, tmp_path / "build", src, "--rebuild")
    assert db.read_bytes() == before
    source(src, [])
    assert run(monkeypatch, tmp_path / "build", src, "--rebuild") == 2
    assert db.read_bytes() == before


def test_root_change_and_removal_replaces_the_snapshot(tmp_path, monkeypatch):
    src = tmp_path / "Taxon.tsv"
    source(src, ["1\tRoot one\t", "2\tChild one\t1", "3\tRoot two\t"])
    root = tmp_path / "build"
    run(monkeypatch, root, src, "--root-id", "1")
    assert len(logical(root / "taxonomy.sqlite")[0]) == 2
    run(monkeypatch, root, src, "--root-id", "3")
    assert [r[0] for r in logical(root / "taxonomy.sqlite")[0]] == ["3"]
    run(monkeypatch, root, src)
    assert len(logical(root / "taxonomy.sqlite")[0]) == 3


def test_orchestrator_checks_receipts_not_existence(tmp_path, monkeypatch):
    src = tmp_path / "Taxon.tsv"
    source(src, ["1\tTaxon name\t"])
    root = tmp_path / "build"
    run(monkeypatch, root, src)
    args = argparse.Namespace(taxonomy_source="dwc", taxonomy_root_id=None, taxonomy_path=src,
                              output_dir=root, force_rebuild=False, force_rebuild_taxonomy=False)
    assert orchestrator._should_skip_taxonomy_ingest(args)
    source(src, ["2\tReplacement\t"])
    assert not orchestrator._should_skip_taxonomy_ingest(args)
    assert "pre-build" in orchestrator._check_taxonomy_available(args, [orchestrator.Step("extract", "", "")])


def test_archive_fingerprint_ignores_unconsumed_files(tmp_path):
    path = tmp_path / "tax.zip"
    def archive(extra):
        with zipfile.ZipFile(path, "w") as z:
            z.writestr("Taxon.tsv", "taxonID\tscientificName\n1\tTaxon name\n")
            z.writestr("README.txt", extra)
    archive("first")
    fp = ti.source_fingerprint("dwca", None, path)
    archive("changed but not consumed")
    assert ti.source_fingerprint("dwca", None, path) == fp


def test_worms_noop_is_pinned_until_explicit_rebuild(tmp_path, monkeypatch):
    calls = []
    def walk(*args, **kwargs):
        calls.append(True)
        return [ti.make_record(taxon_id="1", scientific_name="Taxon name")]
    monkeypatch.setattr(ti, "iter_worms_walk", walk)
    argv = ["ingest", str(tmp_path), "--source", "worms", "--root-id", "1"]
    monkeypatch.setattr("sys.argv", argv)
    assert ti.main() == 0
    assert ti.main() == 0
    assert len(calls) == 1
    monkeypatch.setattr("sys.argv", [*argv, "--rebuild"])
    assert ti.main() == 0
    assert len(calls) == 2


def test_snapshot_receipt_is_stable_across_reingests_of_identical_input(tmp_path, monkeypatch):
    """#278: the fingerprint must identify the source, not the file's bytes.

    ``taxonomy.sqlite`` embeds per-row ``fetched_at`` and ``meta.last_ingest_ts``,
    so re-ingesting byte-identical input produces a byte-different file. When the
    stage fingerprint hashed that file, every document's ``taxa.json`` changed on
    every rebuild and no two builds could ever compare equal.
    """
    src = tmp_path / "Taxon.tsv"
    source(src, ["1\tAlpha\t", "2\tBeta\t1"])

    a = tmp_path / "a"
    b = tmp_path / "b"
    assert run(monkeypatch, a, src) == 0
    assert run(monkeypatch, b, src) == 0
    db_a, db_b = a / "taxonomy.sqlite", b / "taxonomy.sqlite"

    # The premise: same logical content, different bytes.
    assert logical(db_a) == logical(db_b)
    assert db_a.read_bytes() != db_b.read_bytes()

    # The property: the recorded receipt is identical anyway.
    receipt_a = ti.snapshot_receipt(db_a)
    receipt_b = ti.snapshot_receipt(db_b)
    assert receipt_a is not None
    assert receipt_a == receipt_b

    # And a re-ingest over an existing snapshot keeps it stable.
    assert run(monkeypatch, a, src) == 0
    assert ti.snapshot_receipt(db_a) == receipt_a


def test_snapshot_receipt_is_none_for_a_legacy_snapshot(tmp_path, monkeypatch):
    """Snapshots predating receipts fall back to the file hash, not a crash."""
    src = tmp_path / "Taxon.tsv"
    source(src, ["1\tAlpha\t"])
    root = tmp_path / "legacy"
    assert run(monkeypatch, root, src) == 0
    db = root / "taxonomy.sqlite"
    conn = sqlite3.connect(db)
    conn.execute("DELETE FROM meta WHERE key='input_fingerprint'")
    conn.commit()
    conn.close()
    assert ti.snapshot_receipt(db) is None
    assert ti.snapshot_receipt(tmp_path / "absent.sqlite") is None
