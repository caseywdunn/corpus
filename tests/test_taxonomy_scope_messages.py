"""#298: honest single-root replacement messages through the actual CLI route."""

import logging
import signal
import sqlite3

import pytest

from pipeline import taxonomy_ingest as ti


@pytest.fixture
def ingest(tmp_path, monkeypatch):
    source = tmp_path / "Taxon.tsv"
    source.write_text(
        "taxonID\tscientificName\tparentNameUsageID\n"
        "1\tCommon ancestor\t\n10\tFirst clade\t1\n"
        "11\tFirst child\t10\n20\tOther clade\t1\n"
    )
    output = tmp_path / "build"
    old_handler = signal.getsignal(signal.SIGINT)

    def run(root=None, *extra):
        argv = ["ingest", str(output), "--source", "dwc", "--input", str(source)]
        if root is not None:
            argv += ["--root-id", root]
        monkeypatch.setattr("sys.argv", argv + list(extra))
        return ti.main()

    yield run, source, output
    signal.signal(signal.SIGINT, old_handler)


def _ids(output):
    with sqlite3.connect(output / "taxonomy.sqlite") as conn:
        return {row[0] for row in conn.execute("SELECT taxon_id FROM taxa")}


def _warnings(caplog):
    return "\n".join(record.message for record in caplog.records if record.levelno >= logging.WARNING)


@pytest.mark.parametrize("before,after", [(None, "10"), ("10", None)])
def test_scope_changes_with_unrestricted_snapshot_warn_and_retain_backup(ingest, caplog, before, after):
    run, _, output = ingest
    assert run(before) == 0
    previous = (output / "taxonomy.sqlite").read_bytes()
    caplog.clear()
    assert run(after) == 0
    warning = _warnings(caplog)
    assert "Replacing the snapshot" in warning
    assert "--root-id 10" in warning and "no root restriction" in warning
    assert "replace rather than merge" in warning and ".retired/" in warning
    assert _ids(output) == ({"10", "11"} if after else {"1", "10", "11", "20"})
    retired = list((output / ".retired").glob("taxonomy-*.sqlite"))
    assert len(retired) == 1 and retired[0].read_bytes() == previous


def test_ancestor_replacement_does_not_claim_previous_root_disappears(ingest, caplog):
    run, _, output = ingest
    assert run("10") == 0
    caplog.clear()
    assert run("1") == 0
    warning = _warnings(caplog)
    assert "--root-id 10" in warning and "--root-id 1" in warning
    assert "NOT be in the result" not in warning
    assert "previous taxa are not carried forward automatically" in warning
    assert {"1", "10", "11", "20"} == _ids(output)


@pytest.mark.parametrize("root", [None, "10"])
def test_fresh_and_unchanged_scope_are_quiet_and_noop_is_byte_identical(ingest, caplog, root):
    run, _, output = ingest
    assert run(root) == 0
    assert not _warnings(caplog)
    previous = (output / "taxonomy.sqlite").read_bytes()
    receipt = ti.snapshot_receipt(output / "taxonomy.sqlite")
    assert run(root) == 0
    assert not _warnings(caplog)
    assert (output / "taxonomy.sqlite").read_bytes() == previous
    assert not (output / ".retired").exists()
    assert run(root, "--rebuild") == 0
    assert not _warnings(caplog)
    assert ti.snapshot_receipt(output / "taxonomy.sqlite") == receipt


def test_dry_run_describes_proposal_without_replacing_or_retiring(ingest, caplog):
    run, _, output = ingest
    assert run("10") == 0
    previous = (output / "taxonomy.sqlite").read_bytes()
    caplog.clear()
    assert run("20", "--dry-run") == 0
    warning = _warnings(caplog)
    assert "Would replace the snapshot" in warning
    assert "Dry-run leaves the snapshot unchanged" in warning
    assert "Replacing the snapshot" not in warning
    assert "A successful replacement keeps" in warning
    assert (output / "taxonomy.sqlite").read_bytes() == previous
    assert not (output / ".retired").exists()


@pytest.mark.parametrize("keep_receipt", [True, False])
def test_legacy_scope_uses_receipt_or_reports_unknown_without_guessing(ingest, caplog, keep_receipt):
    run, _, output = ingest
    assert run("10") == 0
    with sqlite3.connect(output / "taxonomy.sqlite") as conn:
        conn.execute("DELETE FROM meta WHERE key='root_id'")
        if not keep_receipt:
            conn.execute("DELETE FROM meta WHERE key='input_fingerprint'")
    previous = (output / "taxonomy.sqlite").read_bytes()
    caplog.clear()
    assert run("20", "--dry-run") == 0
    warning = _warnings(caplog)
    assert ("--root-id 10" if keep_receipt else "previous root scope unknown") in warning
    assert "no root restriction" not in warning
    assert (output / "taxonomy.sqlite").read_bytes() == previous


def test_warning_precedes_source_read_failure_and_preserves_previous_snapshot(ingest, caplog, monkeypatch):
    run, _, output = ingest
    assert run("10") == 0
    previous = (output / "taxonomy.sqlite").read_bytes()
    caplog.clear()

    def fail_source_read(*args):
        assert "--root-id 10" in _warnings(caplog)
        assert "--root-id 20" in _warnings(caplog)
        raise OSError("source unavailable")

    monkeypatch.setattr(ti, "source_fingerprint", fail_source_read)
    with pytest.raises(OSError, match="source unavailable"):
        run("20")
    assert (output / "taxonomy.sqlite").read_bytes() == previous
    assert not (output / ".retired").exists()


def test_worms_dry_run_honestly_reports_no_walk(ingest, monkeypatch, caplog):
    run, _, output = ingest
    assert run("10") == 0
    previous = (output / "taxonomy.sqlite").read_bytes()

    def forbidden_walk(*args, **kwargs):
        pytest.fail("A WoRMS dry-run must not walk the API")

    monkeypatch.setattr(ti, "iter_worms_walk", forbidden_walk)
    monkeypatch.setattr("sys.argv", ["ingest", str(output), "--source", "worms",
                                     "--root-id", "20", "--dry-run"])
    caplog.clear()
    with caplog.at_level(logging.INFO):
        assert ti.main() == 0
    assert "would walk WoRMS from --root-id 20; no API requests" in caplog.text
    assert (output / "taxonomy.sqlite").read_bytes() == previous
    assert not (output / ".retired").exists()


def test_help_explains_force_refresh_and_single_root_without_merge_or_acknowledgment(monkeypatch, capsys):
    monkeypatch.setattr("sys.argv", ["ingest", "--help"])
    with pytest.raises(SystemExit) as result:
        ti.main()
    assert result.value.code == 0
    help_text = " ".join(capsys.readouterr().out.split())
    assert "One root restriction per snapshot; sequential ingests replace, never merge" in help_text
    assert "Force ingestion even when the source/config receipt is unchanged" in help_text
    assert "without walking the API or counting records" in help_text
    assert "REPLACE-merges" not in help_text and "acknowledges" not in help_text
