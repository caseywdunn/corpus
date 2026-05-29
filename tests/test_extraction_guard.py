"""Unit tests for the zero-yield extraction guard (#99).

Exercises ``pipeline.main._audit_corpus_chunks`` — the helper behind the
silent-failure guard that makes ``corpus run`` exit non-zero when an
extraction produces no chunks instead of packaging an empty bundle.
"""

import json

from pipeline.main import _audit_corpus_chunks


def _make_doc(documents_dir, name, *, summary=True, chunks=None):
    """Create a per-document dir. ``chunks=None`` omits chunks.json;
    ``chunks=N`` writes N placeholder chunks; ``summary=False`` omits
    summary.json (the dir then looks un-attempted)."""
    d = documents_dir / name
    d.mkdir(parents=True)
    if summary:
        (d / "summary.json").write_text(json.dumps({"status": "success"}))
    if chunks is not None:
        (d / "chunks.json").write_text(
            json.dumps({"chunks": [{"id": i} for i in range(chunks)]})
        )
    return d


def test_missing_documents_dir_is_empty(tmp_path):
    attempted, total, zeros = _audit_corpus_chunks(tmp_path / "nope")
    assert (attempted, total, zeros) == (0, 0, [])


def test_all_documents_empty(tmp_path):
    docs = tmp_path / "documents"
    docs.mkdir()
    _make_doc(docs, "aaaa", chunks=0)
    _make_doc(docs, "bbbb", chunks=None)  # no chunks.json at all
    attempted, total, zeros = _audit_corpus_chunks(docs)
    assert attempted == 2
    assert total == 0
    assert sorted(zeros) == ["aaaa", "bbbb"]


def test_mixed_some_empty(tmp_path):
    docs = tmp_path / "documents"
    docs.mkdir()
    _make_doc(docs, "aaaa", chunks=5)
    _make_doc(docs, "bbbb", chunks=0)
    attempted, total, zeros = _audit_corpus_chunks(docs)
    assert attempted == 2
    assert total == 5
    assert zeros == ["bbbb"]


def test_all_populated(tmp_path):
    docs = tmp_path / "documents"
    docs.mkdir()
    _make_doc(docs, "aaaa", chunks=3)
    _make_doc(docs, "bbbb", chunks=7)
    attempted, total, zeros = _audit_corpus_chunks(docs)
    assert attempted == 2
    assert total == 10
    assert zeros == []


def test_unattempted_dirs_ignored(tmp_path):
    """A dir without summary.json hasn't been extracted — it must not
    count toward ``attempted`` (so a pure-resume no-op doesn't trip the
    guard) and a stray non-dir entry is skipped."""
    docs = tmp_path / "documents"
    docs.mkdir()
    _make_doc(docs, "aaaa", chunks=4)
    _make_doc(docs, "bbbb", summary=False, chunks=0)  # no summary -> ignored
    (docs / "loose_file.txt").write_text("not a doc dir")
    attempted, total, zeros = _audit_corpus_chunks(docs)
    assert attempted == 1
    assert total == 4
    assert zeros == []


def test_malformed_chunks_json_counts_as_zero(tmp_path):
    docs = tmp_path / "documents"
    docs.mkdir()
    d = _make_doc(docs, "aaaa", chunks=None)
    (d / "chunks.json").write_text("{ this is not valid json")
    attempted, total, zeros = _audit_corpus_chunks(docs)
    assert attempted == 1
    assert total == 0
    assert zeros == ["aaaa"]
