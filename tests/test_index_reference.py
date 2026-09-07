"""Snapshots detect logical index drift without conflating build history."""
import json
import sqlite3

from PIL import Image, PngImagePlugin

from tools.qc.index_reference import index_snapshot, figure_pixels, vector_snapshot
from tests.test_embedding_updates import build  # noqa: F401 — pytest fixture


def test_sqlite_order_timestamps_and_raw_history_are_not_current_content(tmp_path):
    from bib.authority import create_schema
    path = tmp_path / "biblio_authority.sqlite"
    conn = sqlite3.connect(path)
    create_schema(conn)
    conn.execute("INSERT INTO works(work_id,guid_type,source,title,created_at,updated_at) VALUES ('w','doi','corpus_paper','Title',1,1)")
    conn.execute("INSERT INTO observation_work VALUES ('o','w','doi_exact',1,'rule-v1',1)")
    conn.commit()
    before = index_snapshot(tmp_path)
    conn.execute("UPDATE works SET updated_at=2")
    conn.execute("UPDATE observation_work SET mapped_at=2")
    conn.execute("INSERT INTO reference_observations(observation_id,citing_corpus_hash,ordinal,authors_json,first_seen_at) VALUES ('old','h',0,'[]',1)")
    conn.commit()
    assert before == index_snapshot(tmp_path)
    conn.execute("UPDATE observation_work SET match_method='author_year_only',match_score=0.5")
    conn.commit()
    after = index_snapshot(tmp_path)
    assert before["bibliography"]["observation_work"] != after["bibliography"]["observation_work"]
    # Human curation is semantic even though the time of curation is not.
    conn.execute("UPDATE works SET bib_imported_at=3")
    conn.commit()
    assert after["bibliography"]["works"] != index_snapshot(tmp_path)["bibliography"]["works"]
    conn.close()


def test_mention_ids_are_not_evidence_but_spans_are(tmp_path):
    from pipeline.taxon_mentions import create_schema
    conn = sqlite3.connect(tmp_path / "taxon_mentions.sqlite")
    create_schema(conn)
    conn.execute("INSERT INTO taxon_mentions(matched_name,corpus_hash,chunk_id,chunk_index,char_start,char_end,mention_text) VALUES ('Taxon','h','c',0,0,5,'Taxon')")
    conn.commit()
    before = index_snapshot(tmp_path)
    conn.execute("UPDATE taxon_mentions SET mention_id=99")
    conn.commit()
    assert before == index_snapshot(tmp_path)
    conn.execute("UPDATE taxon_mentions SET char_end=6")
    conn.commit()
    assert before != index_snapshot(tmp_path)
    conn.close()


def test_vector_order_generation_and_duplicate_or_pixel_drift(build):
    hd = build.document()
    build.run(hd)
    path = build.root / "vector_db/lancedb"
    before = vector_snapshot(path)
    rows = build.table.to_arrow().to_pylist()
    for row in rows:
        row["embedding_generation"] = "another-transaction"
    build.table.delete("true")
    build.table.add(list(reversed(rows)))
    assert before == vector_snapshot(path)
    build.table.add([rows[0]])
    assert before != vector_snapshot(path)
    build.table.delete("true")
    rows[0]["vector"][0] += 0.01
    build.table.add(rows)
    assert before != vector_snapshot(path)


def test_missing_indexes_stay_missing_and_reads_do_not_write(build):
    hd = build.document()
    build.run(hd)
    paths_before = {p.relative_to(build.root): (p.stat().st_size, p.stat().st_mtime_ns)
                    for p in build.root.rglob("*") if p.is_file()}
    result = index_snapshot(build.root)
    assert result["bibliography"] is None and result["taxonomy"] is None
    paths_after = {p.relative_to(build.root): (p.stat().st_size, p.stat().st_mtime_ns)
                   for p in build.root.rglob("*") if p.is_file()}
    assert paths_before == paths_after


def test_pixel_evidence_ignores_png_metadata_but_not_pixel_changes(tmp_path):
    path = tmp_path / "figure.png"
    image = Image.new("RGB", (5, 5), "white")
    image.save(path)
    before = figure_pixels(tmp_path)
    metadata = PngImagePlugin.PngInfo()
    metadata.add_text("creation_time", "different encoder metadata")
    image.save(path, pnginfo=metadata, compress_level=1)
    assert before == figure_pixels(tmp_path)
    image.putpixel((0, 0), (0, 0, 0))
    image.save(path)
    assert before != figure_pixels(tmp_path)


def test_metadata_json_compares_as_data_not_serialization(tmp_path):
    from bib.authority import create_schema
    conn = sqlite3.connect(tmp_path / "biblio_authority.sqlite")
    create_schema(conn)
    conn.execute("INSERT INTO work_documents VALUES ('h','w',?, 'source')", (json.dumps({"title": "T", "year": 2000}),))
    conn.commit()
    before = index_snapshot(tmp_path)
    conn.execute("UPDATE work_documents SET metadata_json=?", ('{"year":2000,"title":"T"}',))
    conn.commit()
    assert before == index_snapshot(tmp_path)
    conn.close()


def test_taxonomy_reserved_column_names_and_fetch_timestamps(tmp_path):
    from pipeline.taxonomy_ingest import create_schema
    conn = sqlite3.connect(tmp_path / "taxonomy.sqlite")
    create_schema(conn)
    conn.execute('INSERT INTO taxa(taxon_id,scientific_name,"order",fetched_at) VALUES (\'1\',\'Taxon\',\'Example\',1)')
    conn.commit()
    before = index_snapshot(tmp_path)
    conn.execute("UPDATE taxa SET fetched_at=2")
    conn.commit()
    assert before == index_snapshot(tmp_path)
    conn.execute('UPDATE taxa SET "order"=\'Revised\'')
    conn.commit()
    assert before != index_snapshot(tmp_path)
    conn.close()
