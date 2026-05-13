"""Export a corpus ``taxonomy.sqlite`` to a Darwin Core Archive (.zip).

Round-trip companion to :mod:`pipeline.taxonomy_ingest`. An export
produced here can be re-ingested via
``corpus taxonomy ingest --source dwca <path.zip>``. This unblocks
two concrete use cases:

  - **Sharing a taxonomy snapshot.** A corpus operator who built
    ``taxonomy.sqlite`` from WoRMS (or curated it by hand) can ship
    the DwC-A to a colleague without forcing them to walk the
    WoRMS API. The recipient ingests it with the offline ``dwca``
    source.

  - **CI test fixtures.** A small DwC-A committed to the repo lets
    the ``dwca`` ingest code path run on every CI push without
    external network calls — see ``tests/fixtures/ci_corpus/`` for
    the canonical example.

DwC-A layout produced::

    taxonomy.zip
      ├── meta.xml        schema descriptor (TDWG DwC-A 1.1)
      └── taxon.tsv       one row per taxon; UTF-8, tab-separated,
                          single header row

Round-trip property: every column the corpus ingester reads (see
``CANONICAL_FIELD_MAP`` in :mod:`pipeline.taxonomy_ingest`) is
emitted here using its canonical Darwin Core term, so
ingest(export(s)) recovers the same ``taxa`` row set as ``s`` — same
row count, same accepted/synonym distribution, same parent /
accepted-name links. Provenance fields (``source``,
``source_record_id``, ``fetched_at``) are intentionally NOT in DwC,
so they don't round-trip; ``source`` is reset to ``"dwca"`` on
re-ingest.

The ``names`` table — used at runtime for fast name → taxon lookup
during corpus extraction — is not exported separately; it is rebuilt
by the ingester from ``scientificName`` + synonym rows. Vernacular
names recorded in ``names`` (``name_type='vernacular'``) are also
emitted in a sidecar ``vernacularname.tsv`` so external DwC tools see
them, but the corpus ingester re-derives its ``names`` index from
``taxa`` alone.

Usage::

    # Via the unified CLI (preferred):
    corpus taxonomy export -o taxonomy.zip

    # As a module:
    python -m pipeline.taxonomy_export <output_dir> -o taxonomy.zip
"""
from __future__ import annotations

import argparse
import csv
import logging
import sqlite3
import sys
import zipfile
from pathlib import Path
from typing import List, Tuple

logger = logging.getLogger("taxonomy_export")


# Column order in taxon.tsv. The ID column comes first (DwC-A meta.xml
# uses index=0 for the row identifier); the rest follow the canonical
# field map in pipeline.taxonomy_ingest, written under their bare
# Darwin Core term names so the ingester's bare-name aliases pick
# them up.
_TAXON_COLUMNS: List[Tuple[str, str]] = [
    # (DwC term, SQLite column)
    ("taxonID",                    "taxon_id"),
    ("scientificName",             "scientific_name"),
    ("scientificNameAuthorship",   "scientific_name_authorship"),
    ("taxonRank",                  "taxon_rank"),
    ("taxonomicStatus",            "taxonomic_status"),
    ("parentNameUsageID",          "parent_name_usage_id"),
    ("acceptedNameUsageID",        "accepted_name_usage_id"),
    ("acceptedNameUsage",          "accepted_name"),
    ("kingdom",                    "kingdom"),
    ("phylum",                     "phylum"),
    ("class",                      "class"),
    ("order",                      '"order"'),  # 'order' is a SQL reserved word
    ("family",                     "family"),
    ("genus",                      "genus"),
    ("bibliographicCitation",      "citation"),
]


# Optional second core: vernacular names. One row per (taxon_id, name)
# from the ``names`` table where name_type='vernacular'. Lets external
# DwC-A consumers (GBIF, ChecklistBank) see common names; not used by
# the corpus ingester itself.
_VERNACULAR_COLUMNS: List[Tuple[str, str]] = [
    ("taxonID",       "taxon_id"),
    ("vernacularName", "name"),
]


_META_XML_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<archive xmlns="http://rs.tdwg.org/dwc/text/">
  <core encoding="UTF-8"
        fieldsTerminatedBy="\\t"
        linesTerminatedBy="\\n"
        fieldsEnclosedBy=""
        ignoreHeaderLines="1"
        rowType="http://rs.tdwg.org/dwc/terms/Taxon">
    <files><location>taxon.tsv</location></files>
    <id index="0" />
{taxon_fields}  </core>
{vernacular_block}</archive>
"""

_VERNACULAR_EXTENSION_TEMPLATE = """  <extension encoding="UTF-8"
             fieldsTerminatedBy="\\t"
             linesTerminatedBy="\\n"
             fieldsEnclosedBy=""
             ignoreHeaderLines="1"
             rowType="http://rs.tdwg.org/dwc/terms/VernacularName">
    <files><location>vernacularname.tsv</location></files>
    <coreid index="0" />
    <field index="1" term="http://rs.tdwg.org/dwc/terms/vernacularName" />
  </extension>
"""


def _render_field_lines(columns: List[Tuple[str, str]], start_index: int) -> str:
    """Render ``<field index="N" term="..."/>`` lines for meta.xml."""
    lines = []
    for i, (term, _col) in enumerate(columns[1:], start=start_index):
        lines.append(
            f'    <field index="{i}" '
            f'term="http://rs.tdwg.org/dwc/terms/{term}" />'
        )
    return "\n".join(lines) + ("\n" if lines else "")


def _write_taxon_tsv(conn: sqlite3.Connection, out_path: Path) -> int:
    """Dump the ``taxa`` table to ``out_path``. Returns row count."""
    select_cols = ", ".join(col for _term, col in _TAXON_COLUMNS)
    cur = conn.execute(f"SELECT {select_cols} FROM taxa ORDER BY taxon_id")
    n = 0
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n",
                       quoting=csv.QUOTE_MINIMAL)
        w.writerow([term for term, _col in _TAXON_COLUMNS])
        for row in cur:
            # SQLite returns None for NULL; csv writes those as empty
            # strings, which is what DwC expects for missing values.
            w.writerow(["" if v is None else str(v) for v in row])
            n += 1
    return n


def _write_vernacular_tsv(conn: sqlite3.Connection, out_path: Path) -> int:
    """Dump vernacular names from ``names`` table. Returns row count.

    Skipped (returns 0, file not created) when no vernacular names
    exist — keeps the archive minimal for taxonomies that don't carry
    them. The caller omits the extension block from meta.xml in that
    case.
    """
    try:
        cur = conn.execute(
            "SELECT taxon_id, name FROM names "
            "WHERE name_type = 'vernacular' ORDER BY taxon_id, name"
        )
    except sqlite3.OperationalError:
        # Older corpora may not have the names table populated.
        return 0
    rows = cur.fetchall()
    if not rows:
        return 0
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n",
                       quoting=csv.QUOTE_MINIMAL)
        w.writerow([term for term, _col in _VERNACULAR_COLUMNS])
        for tid, name in rows:
            w.writerow([str(tid), str(name)])
    return len(rows)


def export_to_dwca(sqlite_path: Path, out_zip: Path) -> dict:
    """Read ``sqlite_path`` and write a DwC-A to ``out_zip``.

    Returns a dict with ``{taxa_count, vernacular_count, out_zip}``
    for caller logging. Raises FileNotFoundError if ``sqlite_path``
    doesn't exist.
    """
    if not sqlite_path.exists():
        raise FileNotFoundError(
            f"{sqlite_path} not found — run `corpus run` (or "
            f"`corpus taxonomy ingest`) first to build the SQLite."
        )

    out_zip = out_zip.resolve()
    out_zip.parent.mkdir(parents=True, exist_ok=True)

    # Stage the archive contents in a temp dir; ZIP-and-cleanup.
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        taxon_path = tdp / "taxon.tsv"
        vernacular_path = tdp / "vernacularname.tsv"

        conn = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
        try:
            n_taxa = _write_taxon_tsv(conn, taxon_path)
            n_vern = _write_vernacular_tsv(conn, vernacular_path)
        finally:
            conn.close()

        # meta.xml — describes whichever files we actually wrote.
        taxon_fields = _render_field_lines(_TAXON_COLUMNS, start_index=1)
        vernacular_block = ""
        if n_vern > 0:
            vernacular_block = _VERNACULAR_EXTENSION_TEMPLATE
        meta_xml = _META_XML_TEMPLATE.format(
            taxon_fields=taxon_fields,
            vernacular_block=vernacular_block,
        )
        (tdp / "meta.xml").write_text(meta_xml, encoding="utf-8")

        with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.write(tdp / "meta.xml", "meta.xml")
            zf.write(taxon_path, "taxon.tsv")
            if n_vern > 0:
                zf.write(vernacular_path, "vernacularname.tsv")

    return {
        "taxa_count": n_taxa,
        "vernacular_count": n_vern,
        "out_zip": str(out_zip),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "output_dir", type=Path,
        help="Corpus output directory; reads <output_dir>/taxonomy.sqlite "
             "unless --db overrides.",
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True,
        help="Output DwC-A path (typically a .zip).",
    )
    parser.add_argument(
        "--db", type=Path, default=None,
        help="Override the source SQLite path (default: "
             "<output_dir>/taxonomy.sqlite).",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    sqlite_path = args.db or (args.output_dir / "taxonomy.sqlite")
    try:
        stats = export_to_dwca(sqlite_path, args.output)
    except FileNotFoundError as e:
        logger.error("%s", e)
        return 1

    logger.info(
        "Exported %d taxa%s → %s",
        stats["taxa_count"],
        f" + {stats['vernacular_count']} vernacular names"
            if stats["vernacular_count"] else "",
        stats["out_zip"],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
