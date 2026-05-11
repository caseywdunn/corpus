"""Idempotency tests for the four post-pipeline scripts (#30).

Verifies the two cases from the issue:
  1. Re-run with no changes → near-zero-work no-op.
  2. Add a paper, re-run → new paper picked up; existing rows untouched.

Cases 3 (input change) and 4 (orphan removal) are out of scope here —
input-change re-annotation belongs to #29 / #33; orphan handling lives
under #31's --audit-orphans.
"""
from __future__ import annotations

import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


# v0.3 (#60) clean break: root-level scripts moved into packages.
# Map the legacy script-name strings the tests use to the new module paths.
_SCRIPT_TO_MODULE = {
    "backfill_intext_citations.py": "pipeline.intext_citations",
    "build_taxon_mentions.py":      "pipeline.taxon_mentions",
    "build_biblio_authority.py":    "bib.authority",
    "reconcile_corpus_to_biblio.py": "bib.reconcile",
}


def _run(script: str, *argv: str) -> subprocess.CompletedProcess:
    """Invoke a moved CLI as ``python -m <module>``.

    Tests still pass the legacy script-name string for readability;
    the helper translates to the post-#60 module path.
    """
    module = _SCRIPT_TO_MODULE.get(script)
    if module is None:
        raise KeyError(f"no module mapping for legacy script {script!r}")
    return subprocess.run(
        [sys.executable, "-m", module, *argv],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )


def _row_counts(db_path: Path) -> dict:
    """Snapshot row counts across every user table — the comparison
    primitive for "did this re-run change anything?"
    """
    counts: dict = {}
    with sqlite3.connect(db_path) as conn:
        for (name,) in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"
        ).fetchall():
            counts[name] = conn.execute(f"SELECT COUNT(*) FROM {name}").fetchone()[0]
    return counts


# ---------------------------------------------------------------------------
# backfill_intext_citations
# ---------------------------------------------------------------------------


_MINIMAL_TEI = """\
<?xml version="1.0" encoding="UTF-8"?>
<TEI xmlns="http://www.tei-c.org/ns/1.0">
  <text>
    <body>
      <div>
        <head>Introduction</head>
        <p>This paper extends earlier work by
           <ref type="bibr" target="#b0">Smith (1995)</ref>
           and the survey of
           <ref type="bibr" target="#b1">Jones &amp; Lee (2001)</ref>.
        </p>
      </div>
    </body>
  </text>
</TEI>
"""


@pytest.fixture
def backfill_corpus(tmp_path: Path) -> Path:
    out = tmp_path / "out"
    docs = out / "documents"
    docs.mkdir(parents=True)
    # paper aaa: has TEI → script will produce intext_citations.json
    (docs / "aaa").mkdir()
    (docs / "aaa" / "grobid.tei.xml").write_text(_MINIMAL_TEI)
    # paper bbb: no TEI → script leaves alone
    (docs / "bbb").mkdir()
    return out


def test_backfill_intext_idempotent(backfill_corpus: Path):
    out = backfill_corpus
    intext = out / "documents" / "aaa" / "intext_citations.json"

    # 1st run — produces the file
    r1 = _run("backfill_intext_citations.py", str(out))
    assert intext.exists()
    payload_before = intext.read_text()
    assert "Done: 2 papers (1 parsed" in r1.stderr or "1 parsed" in r1.stderr

    # 2nd run — no work, file content stable
    r2 = _run("backfill_intext_citations.py", str(out))
    assert intext.read_text() == payload_before
    assert "1 skipped" in r2.stderr


def test_backfill_intext_picks_up_added_paper(backfill_corpus: Path):
    out = backfill_corpus

    # Run once — only aaa parsed
    _run("backfill_intext_citations.py", str(out))
    aaa_payload_before = (out / "documents" / "aaa" / "intext_citations.json").read_text()

    # Add a third paper; re-run
    (out / "documents" / "ccc").mkdir()
    (out / "documents" / "ccc" / "grobid.tei.xml").write_text(_MINIMAL_TEI)
    r = _run("backfill_intext_citations.py", str(out))

    # ccc got parsed; aaa unchanged
    assert (out / "documents" / "ccc" / "intext_citations.json").exists()
    assert (out / "documents" / "aaa" / "intext_citations.json").read_text() == aaa_payload_before
    assert "1 skipped" in r.stderr  # aaa skipped
    assert "1 parsed" in r.stderr   # ccc parsed


def test_backfill_intext_dry_run_writes_nothing(backfill_corpus: Path):
    out = backfill_corpus
    r = _run("backfill_intext_citations.py", str(out), "--dry-run")
    assert "No files written" in r.stderr
    assert not (out / "documents" / "aaa" / "intext_citations.json").exists()


# ---------------------------------------------------------------------------
# build_taxon_mentions
# ---------------------------------------------------------------------------


def _write_taxa(hd: Path, mentions: list) -> None:
    hd.mkdir(parents=True, exist_ok=True)
    (hd / "taxa.json").write_text(json.dumps({
        "total_mentions": len(mentions),
        "unique_taxa": len({m["accepted_taxon_id"] for m in mentions}),
        "mentions": mentions,
    }))


@pytest.fixture
def taxon_corpus(tmp_path: Path) -> Path:
    out = tmp_path / "out"
    docs = out / "documents"
    docs.mkdir(parents=True)

    _write_taxa(docs / "aaa", [
        {"accepted_taxon_id": "T1", "matched_text": "Agalma elegans",
         "accepted_name": "Agalma elegans", "rank": "species",
         "chunk_id": "chunk_0", "text_span": [10, 24], "name_type": "accepted"},
        {"accepted_taxon_id": "T2", "matched_text": "Hippopodius",
         "accepted_name": "Hippopodius", "rank": "genus",
         "chunk_id": "chunk_1", "text_span": [5, 16], "name_type": "accepted"},
    ])
    _write_taxa(docs / "bbb", [
        {"accepted_taxon_id": "T1", "matched_text": "A. elegans",
         "accepted_name": "Agalma elegans", "rank": "species",
         "chunk_id": "chunk_0", "text_span": [3, 13], "name_type": "synonym"},
    ])
    return out


def test_build_taxon_mentions_idempotent(taxon_corpus: Path):
    out = taxon_corpus
    db = out / "taxon_mentions.sqlite"

    _run("build_taxon_mentions.py", str(out))
    counts1 = _row_counts(db)
    assert counts1["taxon_mentions"] == 3   # 2 from aaa + 1 from bbb
    assert counts1["papers_processed"] == 2

    _run("build_taxon_mentions.py", str(out))
    counts2 = _row_counts(db)
    assert counts2 == counts1, "re-run produced new rows; expected no-op"


def test_build_taxon_mentions_picks_up_added_paper(taxon_corpus: Path):
    out = taxon_corpus
    db = out / "taxon_mentions.sqlite"

    _run("build_taxon_mentions.py", str(out))
    baseline = _row_counts(db)

    # Add a third paper
    _write_taxa(out / "documents" / "ccc", [
        {"accepted_taxon_id": "T3", "matched_text": "Forskalia",
         "accepted_name": "Forskalia", "rank": "genus",
         "chunk_id": "chunk_0", "text_span": [0, 9], "name_type": "accepted"},
    ])
    _run("build_taxon_mentions.py", str(out))
    after = _row_counts(db)

    assert after["taxon_mentions"] == baseline["taxon_mentions"] + 1
    assert after["papers_processed"] == baseline["papers_processed"] + 1


def test_build_taxon_mentions_dry_run_no_db(taxon_corpus: Path):
    out = taxon_corpus
    r = _run("build_taxon_mentions.py", str(out), "--dry-run")
    assert "No SQLite writes" in r.stderr
    assert not (out / "taxon_mentions.sqlite").exists()


# ---------------------------------------------------------------------------
# build_biblio_authority — smoke-level idempotency
# ---------------------------------------------------------------------------


@pytest.fixture
def biblio_corpus(tmp_path: Path) -> Path:
    out = tmp_path / "out"
    docs = out / "documents"
    docs.mkdir(parents=True)

    # paper aaa — minimal Grobid-like metadata + one cited reference
    a = docs / "aaa"
    a.mkdir()
    (a / "metadata.json").write_text(json.dumps({
        "title": "On siphonophores",
        "year": 2010,
        "journal": "Test J",
        "doi": "10.1/aaa",
        "authors": [{"surname": "Dunn", "forename": "C"}],
    }))
    # references.json authors are plain author-strings ("J Smith"),
    # not dicts — see grobid_client.parse_tei_references output shape.
    (a / "references.json").write_text(json.dumps({
        "references": [{
            "raw": "Smith J. 1995. Earlier work.",
            "title": "Earlier work",
            "year": 1995,
            "authors": ["J Smith"],
            "xml_id": "b0",
        }],
        "total_references": 1,
    }))
    return out


def test_build_biblio_authority_idempotent(biblio_corpus: Path):
    out = biblio_corpus
    db = out / "biblio_authority.sqlite"

    _run("build_biblio_authority.py", str(out))
    counts1 = _row_counts(db)
    assert counts1.get("works", 0) >= 1
    assert counts1.get("work_authors", 0) >= 1

    _run("build_biblio_authority.py", str(out))
    counts2 = _row_counts(db)
    # works/citations must be stable across re-runs (no double-inserts)
    for table in ("works", "work_authors", "work_aliases",
                  "citations", "taxon_work_links"):
        if table in counts1:
            assert counts2[table] == counts1[table], \
                f"{table} grew on re-run: {counts1[table]} → {counts2[table]}"
