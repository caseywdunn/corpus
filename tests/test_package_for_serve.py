"""Unit + integration tests for package_for_serve.py."""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "package_for_serve.py"

_spec = importlib.util.spec_from_file_location("package_for_serve", SCRIPT)
pkg = importlib.util.module_from_spec(_spec)
sys.modules["package_for_serve"] = pkg
_spec.loader.exec_module(pkg)


# ── _should_copy: mtime-based skip ──────────────────────────────────

def test_should_copy_when_dst_missing(tmp_path: Path):
    src = tmp_path / "src.txt"
    src.write_text("hello")
    assert pkg._should_copy(src, tmp_path / "dst.txt") is True


def test_should_not_copy_when_dst_fresher(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    src.write_text("old")
    time.sleep(0.01)
    dst.write_text("new")
    assert pkg._should_copy(src, dst) is False


def test_should_copy_when_src_fresher(tmp_path: Path):
    src = tmp_path / "src.txt"
    dst = tmp_path / "dst.txt"
    dst.write_text("old")
    time.sleep(0.01)
    src.write_text("new")
    assert pkg._should_copy(src, dst) is True


# ── package(): end-to-end on a minimal fake output tree ─────────────

def _make_fake_output(root: Path, paper_hashes=("abc", "def")) -> Path:
    """Build a minimal output/ tree with per-paper JSONs + figures +
    LanceDB and resource stubs.  Returns root."""
    for h in paper_hashes:
        hd = root / "documents" / h
        (hd / "figures").mkdir(parents=True)
        (hd / "summary.json").write_text(json.dumps({"pdf_hash": h}))
        (hd / "metadata.json").write_text(json.dumps({"title": f"Paper {h}"}))
        (hd / "references.json").write_text(json.dumps({"references": []}))
        (hd / "intext_citations.json").write_text(json.dumps(
            {"paragraphs": [], "citations": []}
        ))
        (hd / "text.json").write_text(json.dumps({"text": "body"}))
        (hd / "chunks.json").write_text(json.dumps(
            {"chunks": [{"id": 0, "text": "hello"}]}
        ))
        (hd / "figures.json").write_text(json.dumps(
            {"figures": [{"id": "f1"}]}
        ))
        (hd / "taxa.json").write_text(json.dumps({"taxa": []}))
        (hd / "anatomy.json").write_text(json.dumps({"terms": []}))
        (hd / "figures" / "fig1.png").write_bytes(b"\x89PNG...stub...")
        # Build-only: should NOT be copied unless --include-pdfs
        (hd / "processed.pdf").write_bytes(b"%PDF-1.5 stub")
        (hd / "docling_doc.json").write_text("{}")
        (hd / "pipeline.log").write_text("log text")
        (hd / "visualizations").mkdir()
        (hd / "visualizations" / "page_1.png").write_bytes(b"viz stub")
    # LanceDB dir
    vdb = root / "vector_db" / "lancedb"
    vdb.mkdir(parents=True)
    (vdb / "manifest.txt").write_bytes(b"lance manifest")
    (vdb / "data.bin").write_bytes(b"lance data blob")
    # Per-hash embedding marker (for manifest reading)
    (root / "vector_db" / "abc_embedded.done").write_text(json.dumps(
        {"embedding_model": "BAAI/bge-m3", "embedding_dim": 1024}
    ))
    return root


def test_package_copies_whitelisted_excludes_build_only(tmp_path: Path):
    src = _make_fake_output(tmp_path / "out")
    dst = tmp_path / "serve"
    manifest = pkg.package(
        output_dir=src, serve_dir=dst,
        version="v1.0.0", include_pdfs=False, dry_run=False,
    )

    # Per-paper whitelisted files present
    for h in ("abc", "def"):
        for f in pkg.PER_PAPER_FILES:
            assert (dst / "documents" / h / f).is_file(), f"missing {h}/{f}"
        assert (dst / "documents" / h / "figures" / "fig1.png").is_file()
        # Build-only files absent
        assert not (dst / "documents" / h / "processed.pdf").exists()
        assert not (dst / "documents" / h / "docling_doc.json").exists()
        assert not (dst / "documents" / h / "pipeline.log").exists()
        assert not (dst / "documents" / h / "visualizations").exists()

    # LanceDB copied in full
    assert (dst / "vector_db" / "lancedb" / "manifest.txt").is_file()
    assert (dst / "vector_db" / "lancedb" / "data.bin").is_file()

    # Manifest is correct
    assert manifest["bundle_version"] == "v1.0.0"
    assert manifest["paper_count"] == 2
    assert manifest["chunk_count"] == 2      # 1 chunk per paper, 2 papers
    assert manifest["figure_count"] == 2     # 1 figure per paper, 2 papers
    assert manifest["embedding_model"] == "BAAI/bge-m3"
    assert manifest["embedding_dim"] == 1024
    assert manifest["includes_pdfs"] is False
    assert manifest["created_at"].endswith("Z")

    # Manifest written to disk
    on_disk = json.loads((dst / "bundle_manifest.json").read_text())
    assert on_disk["bundle_version"] == "v1.0.0"


def test_package_include_pdfs_copies_processed_pdf(tmp_path: Path):
    src = _make_fake_output(tmp_path / "out")
    dst = tmp_path / "serve"
    pkg.package(
        output_dir=src, serve_dir=dst,
        version="v1.0.0", include_pdfs=True, dry_run=False,
    )
    for h in ("abc", "def"):
        assert (dst / "documents" / h / "processed.pdf").is_file()


def test_package_dry_run_writes_nothing(tmp_path: Path):
    src = _make_fake_output(tmp_path / "out")
    dst = tmp_path / "serve"
    manifest = pkg.package(
        output_dir=src, serve_dir=dst,
        version="v1.0.0", include_pdfs=False, dry_run=True,
    )
    # Dry run reports work but doesn't create anything
    assert manifest["_stats"]["dry_run"] is True
    assert manifest["_stats"]["n_files_copied"] > 0
    assert not dst.exists()


def test_package_is_idempotent(tmp_path: Path):
    """Second run should be a near no-op: the mtime check skips
    unchanged files."""
    src = _make_fake_output(tmp_path / "out")
    dst = tmp_path / "serve"
    first = pkg.package(
        output_dir=src, serve_dir=dst,
        version="v1.0.0", include_pdfs=False, dry_run=False,
    )
    n1 = first["_stats"]["n_files_copied"]
    assert n1 > 0

    second = pkg.package(
        output_dir=src, serve_dir=dst,
        version="v1.0.0", include_pdfs=False, dry_run=False,
    )
    n2 = second["_stats"]["n_files_copied"]
    # Manifest always rewrites; the rest should be skipped.  Allow a
    # small slop (<=1 file = manifest) for implementations that touch
    # the manifest.  Everything else should be mtime-skipped.
    assert n2 <= 1, f"idempotent run copied {n2} files"


def test_package_raises_when_documents_missing(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        pkg.package(
            output_dir=tmp_path / "no-such-dir",
            serve_dir=tmp_path / "serve",
            version="v0", include_pdfs=False, dry_run=False,
        )
