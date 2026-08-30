"""`corpus taxonomy ingest` must not grow the names table on a re-run (#262).

`names` shipped with no PRIMARY KEY and no UNIQUE constraint, and both writes
were plain `INSERT`. The only dedup was a `set` built inside
`insert_records`, which knows nothing about rows already on disk — so
re-ingesting appended a complete duplicate set. Measured on the siphonophore
dwca: 801 names became 1,602, then 2,403.

Latent for as long as the code existed, because a re-ingest was an operator
choice. v1.2.1 made it fire on every launch: the #251 fix added an
unconditional pre-build to `slurm/batch_pipeline.sh`, justified by the claim
— asserted without checking a row count — that the ingest no-ops.

Lookups stayed correct throughout (`name_set()` uses SELECT DISTINCT, and
`lookup()` picks among identical rows), so nothing surfaced this except the
file growing. That is why the test is about counts rather than results.
"""
from __future__ import annotations

import sqlite3

import pytest

from pipeline.taxonomy_ingest import create_schema, insert_records, make_record


def _records():
    return [
        make_record(taxon_id="1", scientific_name="Nanomia bijuga",
                    taxonomic_status="accepted"),
        make_record(taxon_id="2", scientific_name="Physalia physalis",
                    taxonomic_status="accepted", extra_names=["Physalia arethusa"]),
    ]


def _counts(conn):
    return (conn.execute("select count(*) from taxa").fetchone()[0],
            conn.execute("select count(*) from names").fetchone()[0])


@pytest.fixture
def conn():
    c = sqlite3.connect(":memory:")
    create_schema(c)
    return c


def test_a_single_ingest_writes_each_name_once(conn):
    insert_records(conn, _records())
    assert _counts(conn) == (2, 3)      # 2 primary + 1 synonym


def test_re_ingesting_does_not_grow_the_names_table(conn):
    """The defect, in one line: this was (2, 6) and then (2, 9)."""
    insert_records(conn, _records())
    first = _counts(conn)
    for _ in range(3):
        insert_records(conn, _records())
    assert _counts(conn) == first


def test_the_reported_count_describes_the_table_not_the_attempt(conn):
    """`n_names` counted attempted inserts, so the log said "801 names"
    against a table holding 1,602. It must count rows actually written."""
    first = insert_records(conn, _records())
    assert first["names"] == 3
    again = insert_records(conn, _records())
    assert again["names"] == 0, "nothing new was written; the log must say so"


def test_new_names_still_land_on_a_second_ingest(conn):
    """Idempotence must not become 'ignores everything after the first run'."""
    insert_records(conn, _records())
    extra = [make_record(taxon_id="3", scientific_name="Agalma elegans",
                         taxonomic_status="accepted")]
    assert insert_records(conn, extra)["names"] == 1
    assert _counts(conn) == (3, 4)


# --- repairing a database damaged by v1.2.1 ----------------------------------


def test_an_already_doubled_table_is_deduplicated_on_open():
    """A corpuscle built with 1.2.1 must be fixed in place rather than
    needing a rebuild — and CREATE UNIQUE INDEX would fail on it, so the
    dedup has to come first."""
    c = sqlite3.connect(":memory:")
    # Reconstruct the pre-#262 state from the real schema, so this exercises
    # the shipped table definitions rather than a hand-rolled approximation:
    # build it, then drop the one thing 1.2.2 added.
    create_schema(c)
    c.execute("DROP INDEX idx_names_unique")
    rows = [("Nanomia bijuga", "nanomia bijuga", "1", "accepted"),
            ("Physalia physalis", "physalia physalis", "2", "accepted")]
    for _ in range(3):                       # three launches under 1.2.1
        c.executemany("INSERT INTO names VALUES (?,?,?,?)", rows)
    assert c.execute("select count(*) from names").fetchone()[0] == 6

    create_schema(c)                          # what the next ingest does

    assert c.execute("select count(*) from names").fetchone()[0] == 2
    dupes = c.execute(
        "select count(*) from (select taxon_id, name_lowercase, count(*) n "
        "from names group by taxon_id, name_lowercase having n > 1)"
    ).fetchone()[0]
    assert dupes == 0


def test_the_constraint_makes_a_repeat_insert_a_no_op(conn):
    insert_records(conn, _records())
    before = _counts(conn)[1]
    conn.execute("INSERT OR IGNORE INTO names VALUES (?,?,?,?)",
                 ("Nanomia bijuga", "nanomia bijuga", "1", "accepted"))
    assert _counts(conn)[1] == before


# --- and the claims that were made about it ----------------------------------


def test_no_doc_still_claims_the_old_behaviour_without_qualification():
    """Three sites asserted idempotence before it was true. If the claim is
    made, it must name what changed."""
    from pathlib import Path
    repo = Path(__file__).resolve().parent.parent
    for rel in ("slurm/batch_pipeline.sh", "pipeline/config.template.yaml"):
        text = (repo / rel).read_text()
        if "idempotent" in text.lower() or "no-op" in text.lower():
            assert "#262" in text, f"{rel} claims idempotence without the caveat"
