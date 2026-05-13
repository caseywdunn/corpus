"""Behavior tests for :func:`pipeline.io.find_all_pdfs`.

The ``exclude_under`` parameter was added so that re-runs against an
``input_pdfs`` that's an ancestor of ``output_dir`` (the demo's
``input_pdfs: .`` with ``output_dir: ./output``) don't pick up the
per-paper ``documents/<HASH>/processed.pdf`` artifacts from prior
runs and treat them as new documents.
"""
from __future__ import annotations

from pathlib import Path

from pipeline.io import find_all_pdfs


# Minimal valid PDF byte content for tests. Each variant differs by a
# trailing comment so we get distinct SHA-256 digests without invoking
# pdfgen libraries — find_all_pdfs only checks the suffix and reads
# bytes for hashing, never parses the PDF.
def _write_pdf(path: Path, marker: bytes) -> None:
    path.write_bytes(b"%PDF-1.0\n% " + marker + b"\n%%EOF\n")


def test_finds_all_top_level(tmp_path):
    _write_pdf(tmp_path / "a.pdf", b"a")
    _write_pdf(tmp_path / "b.pdf", b"b")
    result = find_all_pdfs(tmp_path)
    assert sum(len(v) for v in result.values()) == 2
    assert len(result) == 2  # two distinct hashes


def test_finds_nested(tmp_path):
    (tmp_path / "sub").mkdir()
    _write_pdf(tmp_path / "a.pdf", b"a")
    _write_pdf(tmp_path / "sub" / "b.pdf", b"b")
    result = find_all_pdfs(tmp_path)
    assert len(result) == 2


def test_groups_duplicates_by_hash(tmp_path):
    _write_pdf(tmp_path / "a.pdf", b"same")
    _write_pdf(tmp_path / "a_copy.pdf", b"same")
    result = find_all_pdfs(tmp_path)
    assert len(result) == 1
    (paths,) = result.values()
    assert sorted(p.name for p in paths) == ["a.pdf", "a_copy.pdf"]


def test_exclude_under_skips_output_subtree(tmp_path):
    """Demo-shape: input_dir is the parent, output_dir is a child."""
    output = tmp_path / "output"
    (output / "documents" / "abcd").mkdir(parents=True)
    _write_pdf(tmp_path / "real.pdf", b"real")
    _write_pdf(output / "documents" / "abcd" / "processed.pdf", b"processed")

    # Without exclusion: both are discovered.
    assert len(find_all_pdfs(tmp_path)) == 2

    # With exclusion: only the input-side PDF is discovered.
    result = find_all_pdfs(tmp_path, exclude_under=output)
    assert len(result) == 1
    (paths,) = result.values()
    assert paths[0].name == "real.pdf"


def test_exclude_under_handles_relative_paths(tmp_path, monkeypatch):
    """exclude_under should resolve before comparing, so a relative path
    pointing at the same physical directory still excludes correctly."""
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "output"
    output.mkdir()
    _write_pdf(tmp_path / "real.pdf", b"real")
    _write_pdf(output / "processed.pdf", b"processed")

    result = find_all_pdfs(Path("."), exclude_under=Path("./output"))
    assert len(result) == 1
    (paths,) = result.values()
    assert paths[0].name == "real.pdf"


def test_exclude_under_disjoint_is_noop(tmp_path):
    """If exclude_under is not under input_dir, nothing is filtered."""
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    here = tmp_path / "here"
    here.mkdir()
    _write_pdf(here / "a.pdf", b"a")
    _write_pdf(here / "b.pdf", b"b")
    result = find_all_pdfs(here, exclude_under=elsewhere)
    assert len(result) == 2


def test_exclude_under_none_is_no_op(tmp_path):
    """Passing ``exclude_under=None`` is equivalent to omitting it."""
    _write_pdf(tmp_path / "a.pdf", b"a")
    _write_pdf(tmp_path / "b.pdf", b"b")
    assert find_all_pdfs(tmp_path) == find_all_pdfs(tmp_path, exclude_under=None)
