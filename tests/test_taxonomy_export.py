"""Unit + round-trip tests for pipeline.taxonomy_export.

Two layers of coverage:

1. Unit: export_to_dwca() against a hand-crafted 4-row taxonomy.sqlite
   produces a valid DwC-A ZIP with meta.xml + taxon.tsv (+ optional
   vernacularname.tsv) and the expected DwC column order.

2. Round-trip: export -> ingest produces a SQLite with the same
   row count + same taxonomic_status distribution as the source.
   The corpus pipeline never sees provenance fields (source,
   fetched_at) on the re-ingest side, so we don't assert on those.
"""
from __future__ import annotations

import sqlite3
import zipfile
from pathlib import Path

import pytest

from pipeline.taxonomy_export import (
    _TAXON_COLUMNS,
    _VERNACULAR_COLUMNS,
    export_to_dwca,
)
from pipeline.taxonomy_ingest import create_schema, iter_dwca


def _make_minimal_taxonomy(path: Path, *, with_vernaculars: bool = False) -> None:
    """Build a synthetic taxonomy.sqlite with 4 taxa (1 accepted +
    1 synonym + 2 child species) so the round-trip can check the
    accepted/synonym distribution preserves correctly."""
    conn = sqlite3.connect(path)
    create_schema(conn)
    conn.executemany(
        """
        INSERT INTO taxa (
            taxon_id, scientific_name, scientific_name_authorship,
            taxon_rank, taxonomic_status,
            parent_name_usage_id, accepted_name_usage_id,
            accepted_name, kingdom, phylum, class, "order", family, genus,
            source, citation
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            ("100", "Apolemia",         "Eschscholtz, 1829", "Genus",
             "accepted", None, None, None,
             "Animalia", "Cnidaria", "Hydrozoa", "Siphonophorae", "Apolemiidae", None,
             "test", None),
            ("101", "Apolemia uvaria",  "(Lesueur, 1815)",   "Species",
             "accepted", "100", None, None,
             "Animalia", "Cnidaria", "Hydrozoa", "Siphonophorae", "Apolemiidae", "Apolemia",
             "test", "Lesueur 1815"),
            ("102", "Stephanomia uvaria", "Lesueur, 1815",   "Species",
             "synonym", "100", "101", "Apolemia uvaria",
             "Animalia", "Cnidaria", "Hydrozoa", "Siphonophorae", "Apolemiidae", "Stephanomia",
             "test", None),
            ("103", "Apolemia rubriversa", "Siebert et al., 2013", "Species",
             "accepted", "100", None, None,
             "Animalia", "Cnidaria", "Hydrozoa", "Siphonophorae", "Apolemiidae", "Apolemia",
             "test", None),
        ],
    )
    # `names` rows: one accepted name per taxon + the synonym name
    # linked back to the accepted taxon. Mirrors what taxonomy_ingest
    # writes when ingesting a real DwC-A.
    conn.executemany(
        "INSERT INTO names (name, name_lowercase, taxon_id, name_type) "
        "VALUES (?, ?, ?, ?)",
        [
            ("Apolemia",          "apolemia",          "100", "accepted"),
            ("Apolemia uvaria",   "apolemia uvaria",   "101", "accepted"),
            ("Stephanomia uvaria", "stephanomia uvaria", "102", "accepted"),
            ("Apolemia rubriversa", "apolemia rubriversa", "103", "accepted"),
        ],
    )
    if with_vernaculars:
        conn.executemany(
            "INSERT INTO names (name, name_lowercase, taxon_id, name_type) "
            "VALUES (?, ?, ?, ?)",
            [
                ("string siphonophore", "string siphonophore", "101", "vernacular"),
                ("Strangenkette",       "strangenkette",       "101", "vernacular"),
            ],
        )
    conn.commit()
    conn.close()


def test_export_emits_zip_with_expected_layout(tmp_path: Path):
    sqlite_path = tmp_path / "taxonomy.sqlite"
    _make_minimal_taxonomy(sqlite_path)
    out_zip = tmp_path / "out.zip"

    stats = export_to_dwca(sqlite_path, out_zip)

    assert out_zip.exists()
    assert stats["taxa_count"] == 4
    assert stats["vernacular_count"] == 0

    with zipfile.ZipFile(out_zip) as zf:
        names = set(zf.namelist())
        # No vernaculars in this fixture, so vernacularname.tsv is NOT
        # in the archive — we want minimal archives when possible.
        assert names == {"meta.xml", "taxon.tsv"}

        # Header row uses the DwC term names from _TAXON_COLUMNS, in
        # order. The ingester recognizes bare camelCase terms, so this
        # is the form it round-trips with no aliasing.
        with zf.open("taxon.tsv") as f:
            header = f.readline().decode("utf-8").rstrip("\n").split("\t")
        expected_header = [term for term, _ in _TAXON_COLUMNS]
        assert header == expected_header

        # meta.xml mentions the Taxon rowType + taxon.tsv as the core
        # location. Don't deeply parse — just sanity-check the shape.
        meta = zf.read("meta.xml").decode("utf-8")
        assert "rowType=\"http://rs.tdwg.org/dwc/terms/Taxon\"" in meta
        assert "<location>taxon.tsv</location>" in meta


def test_export_includes_vernacular_extension_when_names_present(tmp_path: Path):
    sqlite_path = tmp_path / "taxonomy.sqlite"
    _make_minimal_taxonomy(sqlite_path, with_vernaculars=True)
    out_zip = tmp_path / "out.zip"

    stats = export_to_dwca(sqlite_path, out_zip)
    assert stats["vernacular_count"] == 2

    with zipfile.ZipFile(out_zip) as zf:
        assert "vernacularname.tsv" in zf.namelist()
        with zf.open("vernacularname.tsv") as f:
            header = f.readline().decode("utf-8").rstrip("\n").split("\t")
        assert header == [term for term, _ in _VERNACULAR_COLUMNS]
        # meta.xml should now also declare the VernacularName extension.
        meta = zf.read("meta.xml").decode("utf-8")
        assert "VernacularName" in meta


def test_export_raises_when_sqlite_missing(tmp_path: Path):
    """The exporter must fail loudly with a useful hint when the
    operator points it at a corpus output dir that hasn't been
    built yet — silently writing an empty DwC-A would be much worse."""
    with pytest.raises(FileNotFoundError, match="not found"):
        export_to_dwca(tmp_path / "nonexistent.sqlite", tmp_path / "out.zip")


def test_export_then_ingest_preserves_taxa_count_and_status_distribution(
    tmp_path: Path,
):
    """The round-trip property: ingest(export(sqlite)) recovers the
    same row count + accepted/synonym distribution. This is the
    central guarantee that lets CI fixtures live as DwC-A ZIPs in
    the repo and stay in lockstep with the live ingest pipeline.
    """
    src = tmp_path / "src.sqlite"
    _make_minimal_taxonomy(src, with_vernaculars=True)
    zip_path = tmp_path / "trip.zip"

    export_to_dwca(src, zip_path)

    # Re-ingest via iter_dwca → load into a fresh SQLite. We exercise
    # the same code path the `corpus taxonomy ingest --source dwca`
    # CLI uses, just inline so the test stays in-process.
    rebuilt = tmp_path / "rebuilt.sqlite"
    conn = sqlite3.connect(rebuilt)
    create_schema(conn)
    for rec in iter_dwca(zip_path):
        conn.execute(
            """INSERT INTO taxa (
                taxon_id, scientific_name, scientific_name_authorship,
                taxon_rank, taxonomic_status,
                parent_name_usage_id, accepted_name_usage_id,
                accepted_name, kingdom, phylum, class, "order",
                family, genus, source, citation
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                rec["taxon_id"], rec["scientific_name"],
                rec.get("scientific_name_authorship"),
                rec.get("taxon_rank"), rec.get("taxonomic_status"),
                rec.get("parent_name_usage_id"),
                rec.get("accepted_name_usage_id"),
                rec.get("accepted_name"),
                rec.get("kingdom"), rec.get("phylum"), rec.get("class"),
                rec.get("order"), rec.get("family"), rec.get("genus"),
                rec.get("source"), rec.get("citation"),
            ),
        )
    conn.commit()

    def _summary(c):
        return {
            "total": c.execute("SELECT COUNT(*) FROM taxa").fetchone()[0],
            "accepted": c.execute(
                "SELECT COUNT(*) FROM taxa WHERE taxonomic_status='accepted'"
            ).fetchone()[0],
            "synonym": c.execute(
                "SELECT COUNT(*) FROM taxa WHERE taxonomic_status='synonym'"
            ).fetchone()[0],
        }

    src_conn = sqlite3.connect(src)
    assert _summary(conn) == _summary(src_conn), (
        "Round-trip changed the taxa row count or accepted/synonym "
        "distribution; the exporter and ingester have drifted."
    )
    src_conn.close()
    conn.close()
