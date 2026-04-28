#!/usr/bin/env python3
"""Distill a Bouchet build bundle into a served bundle for AWS.

PLAN.md §10 separates the "build bundle" (everything `process_corpus.py`
emits — includes processed.pdf, raw docling dumps, QC visualizations,
per-paper logs) from the "served bundle" (just what MCP tools read
plus the precompiled indices).  For 2000 papers the served bundle is
~3 GB vs. ~10 GB for the build — small enough for a cheap EBS
volume, and fast enough to `aws s3 sync` on every release.

This script is the distiller.  Run on Bouchet after the pipeline
completes; the output goes to S3 and is then pulled onto the EC2
MCP host.

Files in the served-bundle whitelist (the contract):

  Per-paper (``documents/<HASH>/``):
    summary.json, metadata.json, references.json,
    text.json, chunks.json, figures.json,
    taxa.json, anatomy.json,
    figures/*.png                    (all extracted figure images)

  Top-level:
    vector_db/lancedb/               (the embedded chunks)
    resources/taxonomy.sqlite        (copy from repo's resources/)
    resources/biblio_authority.sqlite
    resources/taxon_mentions.sqlite
    resources/anatomy_lexicon.yaml
    bundle_manifest.json             (written by this script)

With ``--include-pdfs`` also: ``processed.pdf`` per document.

Usage:
    python package_for_serve.py /path/to/output /path/to/serve_bundle \\
        --version v1.0.0
    python package_for_serve.py /path/to/output /path/to/serve_bundle \\
        --version v1.0.0 --include-pdfs
    python package_for_serve.py /path/to/output /path/to/serve_bundle --dry-run
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

logger = logging.getLogger("package_for_serve")

# The per-paper whitelist.  Top-level files only — figures/ is handled
# separately as a directory copy.
PER_PAPER_FILES = (
    "summary.json",
    "metadata.json",
    "references.json",
    "text.json",
    "chunks.json",
    "figures.json",
    "taxa.json",
    "anatomy.json",
)

REPO_ROOT = Path(__file__).resolve().parent


# ── File copy with mtime skip ───────────────────────────────────────

def _should_copy(src: Path, dst: Path) -> bool:
    """Skip if destination exists and is at least as fresh as source."""
    if not dst.exists():
        return True
    try:
        return src.stat().st_mtime > dst.stat().st_mtime
    except OSError:
        return True


def _copy_file(src: Path, dst: Path, dry_run: bool) -> int:
    """Return bytes copied (0 if skipped)."""
    if not src.exists():
        return 0
    if not _should_copy(src, dst):
        return 0
    size = src.stat().st_size
    if dry_run:
        return size
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return size


def _copy_tree(src_dir: Path, dst_dir: Path, dry_run: bool,
               glob: str = "*") -> Tuple[int, int]:
    """Copy everything under src_dir/glob into dst_dir, preserving
    relative layout.  Returns (n_files_copied, total_bytes)."""
    n = 0
    total = 0
    if not src_dir.is_dir():
        return 0, 0
    for src in sorted(src_dir.rglob(glob)):
        if not src.is_file():
            continue
        rel = src.relative_to(src_dir)
        dst = dst_dir / rel
        bytes_written = _copy_file(src, dst, dry_run)
        if bytes_written:
            n += 1
            total += bytes_written
    return n, total


# ── Manifest helpers ────────────────────────────────────────────────

def _git_sha() -> Optional[str]:
    """Best-effort git HEAD SHA for the repo this script lives in.
    Returns None if git isn't available or the repo is missing."""
    try:
        out = subprocess.check_output(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
        )
        return out.decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def _read_one_embedding_marker(output_dir: Path) -> Tuple[Optional[str], Optional[int]]:
    """Sample a single ``<HASH>_embedded.done`` marker to learn the
    embedding model + dim used for the vector index.  Returns
    ``(model, dim)`` or ``(None, None)`` if no marker exists."""
    vdb = output_dir / "vector_db"
    if not vdb.is_dir():
        return None, None
    for marker in vdb.glob("*_embedded.done"):
        try:
            data = json.loads(marker.read_text())
        except Exception:
            continue
        return data.get("embedding_model"), data.get("embedding_dim")
    return None, None


def _count_figures_and_chunks(documents_dir: Path) -> Tuple[int, int]:
    """Sum figure_count and chunk_count across all per-paper files.
    Best-effort — missing/malformed JSONs are silently skipped."""
    fig_total = 0
    chunk_total = 0
    for hash_dir in documents_dir.iterdir():
        if not hash_dir.is_dir():
            continue
        fig_path = hash_dir / "figures.json"
        if fig_path.exists():
            try:
                fig_total += len(json.loads(fig_path.read_text()).get("figures", []))
            except Exception:
                pass
        chunks_path = hash_dir / "chunks.json"
        if chunks_path.exists():
            try:
                chunk_total += len(json.loads(chunks_path.read_text()).get("chunks", []))
            except Exception:
                pass
    return fig_total, chunk_total


def _taxonomy_snapshot_date(taxonomy_path: Path) -> Optional[str]:
    """Read the build timestamp from the taxonomy SQLite. ingest_taxonomy
    writes a ``meta(key='last_ingest_ts', value=...)`` row on each build;
    fall back to the file's mtime if that table is absent."""
    if not taxonomy_path.exists():
        return None
    try:
        conn = sqlite3.connect(f"file:{taxonomy_path}?mode=ro", uri=True)
        try:
            cur = conn.execute(
                "SELECT value FROM meta WHERE key IN "
                "('built_at', 'last_ingest_ts') LIMIT 1"
            )
            row = cur.fetchone()
            if row and row[0]:
                # last_ingest_ts is a Unix timestamp; rendered as ISO date.
                try:
                    ts = float(row[0])
                    return dt.datetime.fromtimestamp(ts, tz=dt.timezone.utc).strftime("%Y-%m-%d")
                except (TypeError, ValueError):
                    return str(row[0])
        except sqlite3.OperationalError:
            pass
        finally:
            conn.close()
    except sqlite3.Error:
        pass
    # Fallback: file mtime as ISO date
    try:
        ts = taxonomy_path.stat().st_mtime
        return dt.datetime.fromtimestamp(ts, tz=dt.timezone.utc).strftime("%Y-%m-%d")
    except OSError:
        return None


def _iso_now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ── Main packager ───────────────────────────────────────────────────

def package(output_dir: Path, serve_dir: Path, version: str,
            include_pdfs: bool, dry_run: bool) -> Dict:
    """Copy the served bundle and return the manifest dict."""
    documents_dir = output_dir / "documents"
    if not documents_dir.is_dir():
        raise FileNotFoundError(
            f"{documents_dir} not found — point --output-dir at a processed corpus"
        )

    if not dry_run:
        serve_dir.mkdir(parents=True, exist_ok=True)

    # Per-paper files + figures/
    n_papers = 0
    n_files = 0
    total_bytes = 0
    for hash_dir in sorted(documents_dir.iterdir()):
        if not hash_dir.is_dir():
            continue
        n_papers += 1
        dest_hash_dir = serve_dir / "documents" / hash_dir.name

        for fname in PER_PAPER_FILES:
            src = hash_dir / fname
            dst = dest_hash_dir / fname
            bw = _copy_file(src, dst, dry_run)
            if bw:
                n_files += 1
                total_bytes += bw

        # figures/ directory
        nf, nb = _copy_tree(hash_dir / "figures", dest_hash_dir / "figures",
                            dry_run, glob="*.png")
        n_files += nf
        total_bytes += nb

        if include_pdfs:
            bw = _copy_file(
                hash_dir / "processed.pdf",
                dest_hash_dir / "processed.pdf",
                dry_run,
            )
            if bw:
                n_files += 1
                total_bytes += bw

        if n_papers % 200 == 0:
            logger.info("Copied %d papers, %d files so far", n_papers, n_files)

    # LanceDB — directory tree
    vdb_src = output_dir / "vector_db" / "lancedb"
    vdb_dst = serve_dir / "vector_db" / "lancedb"
    nf, nb = _copy_tree(vdb_src, vdb_dst, dry_run)
    n_files += nf
    total_bytes += nb

    # Resources — copy from the repo's resources/, not from output/, so
    # we don't require the user to symlink them into their output dir.
    resources_src = REPO_ROOT / "resources"
    resources_dst = serve_dir / "resources"
    for fname in ("taxonomy.sqlite", "biblio_authority.sqlite",
                  "taxon_mentions.sqlite", "anatomy_lexicon.yaml"):
        bw = _copy_file(resources_src / fname, resources_dst / fname, dry_run)
        if bw:
            n_files += 1
            total_bytes += bw

    # Manifest
    model, dim = _read_one_embedding_marker(output_dir)
    fig_count, chunk_count = _count_figures_and_chunks(documents_dir)
    manifest = {
        "bundle_version": version,
        "created_at": _iso_now(),
        "pipeline_git_sha": _git_sha(),
        "embedding_model": model,
        "embedding_dim": dim,
        "taxonomy_snapshot_date": _taxonomy_snapshot_date(resources_src / "taxonomy.sqlite"),
        "paper_count": n_papers,
        "figure_count": fig_count,
        "chunk_count": chunk_count,
        "includes_pdfs": include_pdfs,
    }
    manifest_path = serve_dir / "bundle_manifest.json"
    if not dry_run:
        manifest_path.write_text(json.dumps(manifest, indent=2))

    manifest["_stats"] = {
        "n_files_copied": n_files,
        "total_bytes": total_bytes,
        "dry_run": dry_run,
    }
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("output_dir", type=Path,
                        help="Pipeline output (contains documents/ and vector_db/)")
    parser.add_argument("serve_dir", type=Path,
                        help="Destination for the served bundle")
    parser.add_argument("--version", default="v0.0.0-dev",
                        help="Bundle version string (e.g., v1.0.0)")
    parser.add_argument("--include-pdfs", action="store_true",
                        help="Also copy processed.pdf per document (~3 GB extra at 2000 papers)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report what would be copied without writing anything")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    try:
        manifest = package(
            output_dir=args.output_dir.resolve(),
            serve_dir=args.serve_dir.resolve(),
            version=args.version,
            include_pdfs=args.include_pdfs,
            dry_run=args.dry_run,
        )
    except FileNotFoundError as e:
        logger.error("%s", e)
        return 1

    stats = manifest.pop("_stats")
    logger.info("═══ package_for_serve %s ═══",
                "dry-run" if args.dry_run else "complete")
    logger.info("  papers:      %d", manifest["paper_count"])
    logger.info("  chunks:      %d", manifest["chunk_count"])
    logger.info("  figures:     %d", manifest["figure_count"])
    logger.info("  files copied: %d", stats["n_files_copied"])
    logger.info("  bytes:       %.2f MB", stats["total_bytes"] / 1e6)
    logger.info("  version:     %s", manifest["bundle_version"])
    logger.info("  git sha:     %s", manifest["pipeline_git_sha"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
