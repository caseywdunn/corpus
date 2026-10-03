"""A root change replaces the snapshot — and now says so (#298).

Replacing rather than accumulating is deliberate: see
`tests/test_taxonomy_source_updates.py`, whose docstring reads "Taxonomy
snapshots own source receipts and replace, rather than accumulate", and which
pins the behaviour directly. This file does not challenge that.

What was wrong was the silence, and a help string that pointed the other way.
`--rebuild` advertised "otherwise REPLACE-merges", which reads as *a second
ingest merges into the existing snapshot*. It does not — the merge is
`INSERT OR REPLACE` within one run's staging database, never across runs.

That mattered because a library scoped to several clades needs several roots,
and `--root-id` takes one. The obvious way to express "six families" is six
sequential ingests into one output, which keeps the sixth, exits 0, and logs
nothing about the five it discarded. The snapshot then covers a sixth of the
library's literature, and nothing downstream can detect it: a name that does not
resolve is simply dropped from `taxa.json`.

So the ingest now warns, names the root being replaced, and points at
`.retired/` where the previous snapshot survives.
"""
from __future__ import annotations

import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DEMO_ARCHIVE = REPO / "demo" / "taxonomy.zip"

# Two disjoint genera in the demo archive, asserted in the fixture so a changed
# fixture fails loudly rather than making these tests vacuous.
ROOT_A = "135365"
ROOT_B = "135362"

pytestmark = pytest.mark.skipif(
    not DEMO_ARCHIVE.is_file(), reason="demo taxonomy archive not present"
)


def _ingest(out_dir: Path, root_id: str, *extra: str):
    return subprocess.run(
        [sys.executable, "-m", "pipeline.taxonomy_ingest",
         str(out_dir), "--source", "dwca", "--input", str(DEMO_ARCHIVE),
         "--root-id", root_id, *extra],
        cwd=REPO, capture_output=True, text=True,
    )


def _taxa(out_dir: Path) -> set[str]:
    db = out_dir / "taxonomy.sqlite"
    if not db.is_file():
        return set()
    con = sqlite3.connect(db)
    try:
        return {r[0] for r in con.execute("SELECT taxon_id FROM taxa")}
    finally:
        con.close()


@pytest.fixture
def seeded(tmp_path: Path) -> Path:
    result = _ingest(tmp_path, ROOT_A)
    assert result.returncode == 0, result.stderr
    assert ROOT_A in _taxa(tmp_path), "fixture changed: ROOT_A not in the archive"
    return tmp_path


def test_a_root_change_warns_and_names_what_it_drops(seeded: Path):
    """The whole point: the loss is announced rather than silent."""
    result = _ingest(seeded, ROOT_B)

    assert result.returncode == 0, "replacing is allowed by design"
    assert ROOT_A in result.stderr, "the discarded root must be named"
    assert ROOT_B in result.stderr
    assert ".retired" in result.stderr, "say where the old snapshot went"


def test_the_warning_says_snapshots_do_not_merge(seeded: Path):
    """The misleading help text is what made sequential ingests look right."""
    stderr = _ingest(seeded, ROOT_B).stderr
    assert "merge" in stderr.lower()


def test_the_replacement_still_happens(seeded: Path):
    """Warning, not refusing — the documented behaviour is preserved."""
    assert _ingest(seeded, ROOT_B).returncode == 0
    ids = _taxa(seeded)
    assert ROOT_B in ids
    assert ROOT_A not in ids, "snapshots replace rather than accumulate"


def test_the_previous_snapshot_is_recoverable(seeded: Path):
    """The warning points at .retired/; that had better be true."""
    _ingest(seeded, ROOT_B)
    retired = list((seeded / ".retired").glob("taxonomy-*.sqlite"))
    assert retired, "no retired snapshot despite the warning naming one"

    con = sqlite3.connect(retired[0])
    try:
        ids = {r[0] for r in con.execute("SELECT taxon_id FROM taxa")}
    finally:
        con.close()
    assert ROOT_A in ids, "the retired copy should hold the replaced clade"


def test_no_warning_when_the_root_is_unchanged(seeded: Path):
    """A refresh of the same clade is ordinary and should stay quiet."""
    stderr = _ingest(seeded, ROOT_A, "--rebuild").stderr
    assert "Replacing the snapshot" not in stderr


def test_no_warning_on_a_fresh_output(tmp_path: Path):
    """Nothing to replace means nothing to say."""
    result = _ingest(tmp_path, ROOT_B)
    assert result.returncode == 0
    assert "Replacing the snapshot" not in result.stderr
