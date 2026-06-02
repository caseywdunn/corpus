#!/usr/bin/env python3
"""Distill a Bouchet build bundle into a served bundle for AWS.

dev_docs/PLAN.md §10 separates the "build bundle" (everything `process_corpus.py`
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
    intext_citations.json,
    text.json, chunks.json, figures.json,
    taxa.json,
    <category>.json                  (one per category in --lexicon —
                                      anatomy, biogeography, methods, …
                                      — discovered at copy time)
    figures/*.png                    (all extracted figure images)

  Top-level (corpuscle layout — flat at the bundle root):
    vector_db/lancedb/               (the embedded chunks)
    taxonomy.sqlite                  (copied from <output_dir>/)
    biblio_authority.sqlite
    taxon_mentions.sqlite
    instructions.md                  (per-corpus client guidance,
                                      surfaced via MCP InitializeResult.instructions)
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
import re
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from pipeline.version import __version__

logger = logging.getLogger("package_for_serve")

# A string value that *is* an absolute filesystem path needs to be
# scrubbed before this bundle leaves the build host. The pattern is
# anchored at the start of a string value (not a substring search) so
# body-text JSONs (chunks.json, text.json) — which routinely contain
# URLs and natural-language with "/home/..." mid-string — don't trip it.
# The audit also screens out values that match URL-ish shapes (``scheme:``
# prefix; see ``_walk_strings`` filter) before the regex sees them.
#
# Coverage is "starts with /" and ends with a path-shaped second
# component, which catches every common filesystem root: /home, /Users,
# /tmp, /private/tmp (macOS), /var, /opt, /mnt, /nfs, /srv, /Volumes,
# /workspace, /scratch, /data, /export, etc. Earlier revisions allow-
# listed a short set of roots — that left /private, /tmp, /Volumes
# leaking through the audit.
#
# Two false-positive classes excluded via regex + caller guard:
#   1. PDF glyph codes (/G03/G38/, /G3/G71/, /c81f5/c8167/...): first
#      component is one letter + one or more hex digits — excluded by
#      (?![A-Za-z][0-9A-Fa-f]+/). Real path roots (nfs, home, tmp, var,
#      …) always contain at least one non-hex character after the first
#      letter (e.g. 'n'+'f's', not all-hex), so the lookahead never fires
#      on them. Purely-hex names like /cafe/ would be silently skipped —
#      acceptable on HPC systems where roots are human-readable words.
#   2. DOI/URL fragments (/10.1016/j.ocemod.2012.10.007) and glyph+text
#      strings (/journal/rsos R. Soc., /g20/g17 This is...): DOIs start
#      with a digit, excluded by requiring [A-Za-z] as the first char.
#      Strings with spaces are excluded by the ``' ' not in s`` guard in
#      the caller (real filesystem paths never contain spaces on HPC).
_ABS_PATH_RE = re.compile(r'^/(?![A-Za-z][0-9A-Fa-f]+/)[A-Za-z][A-Za-z0-9_.+-]+/')

# The per-paper whitelist.  Top-level files only — figures/ is handled
# separately as a directory copy.
PER_PAPER_FILES = (
    "summary.json",
    "metadata.json",
    "references.json",
    "intext_citations.json",
    "text.json",
    "chunks.json",
    "figures.json",
    "taxa.json",
)


def _per_paper_lexicon_outputs(hash_dir: Path) -> List[str]:
    """Return the filenames of every lexicon-output ``<category>.json``
    in ``hash_dir``.

    The pipeline writes one ``<category>.json`` per category configured
    via ``--lexicon`` (anatomy, biogeography, …). The set is corpus-
    specific, so the served-bundle whitelist has to be discovered at
    copy time rather than hardcoded. Detection: any ``*.json`` whose
    root carries a ``category`` field — exactly what
    ``taxa.extract_lexicon_mentions`` stamps.
    """
    static = set(PER_PAPER_FILES) | {"scan_detection.json", "docling_doc.json",
                                     "pipeline_state.json"}
    out: List[str] = []
    for child in hash_dir.glob("*.json"):
        if child.name in static:
            continue
        try:
            with child.open(encoding="utf-8") as f:
                payload = json.load(f)
        except Exception:
            continue
        if isinstance(payload, dict) and "category" in payload:
            out.append(child.name)
    return out

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


def _count_figures_and_chunks(
    documents_dir: Path,
    excluded_hashes: Optional[Iterable[str]] = None,
) -> Tuple[int, int]:
    """Sum figure_count and chunk_count across the per-paper files we
    actually plan to serve. ``excluded_hashes`` is the set of hashes
    flagged ``works.serve = 0`` (#54) — they live in the build bundle
    but are filtered out of the distilled bundle, so the manifest must
    not count them. Best-effort — missing/malformed JSONs are silently
    skipped."""
    skip = set(excluded_hashes or ())
    fig_total = 0
    chunk_total = 0
    for hash_dir in documents_dir.iterdir():
        if not hash_dir.is_dir():
            continue
        if hash_dir.name in skip:
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


def _drop_chunks_for_hashes(vdb_path: Path, hashes: set) -> int:
    """Delete LanceDB rows whose ``metadata.pdf_hash`` is in ``hashes``.

    Returns the number of rows deleted. Used during distillation to
    keep skipped papers (#54) from leaking through semantic search.
    Falls back gracefully when lancedb isn't importable or the table
    doesn't exist.
    """
    if not hashes:
        return 0
    try:
        import lancedb  # type: ignore
        from pipeline.embeddings import lancedb_table_names
    except ImportError:
        logger.info("lancedb not importable; skipping LanceDB filter")
        return 0
    try:
        db = lancedb.connect(str(vdb_path))
        if "document_chunks" not in lancedb_table_names(db):
            return 0
        table = db.open_table("document_chunks")
        # The schema's hash column is nested as ``metadata.pdf_hash``
        # (see pipeline/embed.py:ChunkMetadata).
        before = table.count_rows()
        in_clause = ", ".join(f"'{h}'" for h in hashes)
        table.delete(f"metadata.pdf_hash IN ({in_clause})")
        return before - table.count_rows()
    except Exception as e:
        logger.warning("Could not filter LanceDB by skipped hashes: %s", e)
        return 0


def _load_skipped_hashes(biblio_path: Path) -> set:
    """Return the set of corpus_hash values whose works.serve = 0 (#54).

    Empty set if biblio_authority.sqlite is missing, lacks the v0.3
    serve column, or contains no skipped works. Read-only.
    """
    if not biblio_path.is_file():
        return set()
    try:
        conn = sqlite3.connect(f"file:{biblio_path}?mode=ro", uri=True)
        try:
            cols = {r[1] for r in conn.execute("PRAGMA table_info(works)")}
            if "serve" not in cols:
                return set()
            rows = conn.execute(
                "SELECT corpus_hash FROM works "
                "WHERE serve = 0 AND corpus_hash IS NOT NULL"
            ).fetchall()
            return {r[0] for r in rows}
        finally:
            conn.close()
    except sqlite3.Error:
        return set()


def _citation_for_manifest() -> Optional[dict]:
    """Read CITATION.cff (packaged into pipeline/CITATION.cff via
    package_data; falls back to repo root for editable installs that
    pre-date the package_data move) + return the resolved preferred
    citation so it can be stamped into bundle_manifest.json (#61).
    Returns None if CITATION.cff is missing.
    """
    cff = None
    try:
        from importlib import resources
        import yaml
        src = resources.files("pipeline").joinpath("CITATION.cff")
        with resources.as_file(src) as cff_path:
            if cff_path.exists():
                cff = yaml.safe_load(cff_path.read_text(encoding="utf-8")) or {}
    except (ModuleNotFoundError, FileNotFoundError, Exception):
        cff = None
    if cff is None:
        # Editable-install fallback.
        repo_root = Path(__file__).resolve().parent.parent
        cff_path = repo_root / "CITATION.cff"
        if not cff_path.exists():
            return None
        try:
            import yaml
            cff = yaml.safe_load(cff_path.read_text(encoding="utf-8")) or {}
        except Exception:
            return None
    pref = cff.get("preferred-citation") or cff
    return {
        "title": pref.get("title"),
        "year": pref.get("year"),
        "doi": pref.get("doi"),
        "url": pref.get("url"),
        "authors": [
            {"family-names": a.get("family-names"), "given-names": a.get("given-names")}
            for a in (pref.get("authors") or [])
        ],
    }


# ── Path scrubbing (dev_docs/PLAN.md §10) ────────────────────────────────────
#
# `process_corpus.py` writes summary.json with absolute Bouchet paths in
# input_dir / output_directory / processing_summary.original_pdf /
# processing_summary.files_created, and figures.json with absolute
# file_path / figures_directory. None of those are useful (or correct)
# once the bundle is served from a different machine. Rewrite them to
# corpus-root-relative form (e.g. ``documents/<HASH>/figures/fig_3.png``)
# at distillation time so the served bundle is portable, and assert
# absence of any remaining absolute path before declaring success.

def _to_corpus_relative(s: str, output_root: Path) -> Optional[str]:
    """Map an absolute path string to corpus-root-relative.

    Returns:
    - ``None`` if ``s`` is not an absolute path (already relative — leave alone).
    - The relative path under ``output_root`` if ``s`` is under that root
      (the common case during a same-machine distillation).
    - A ``documents/<HASH>/...`` slice if the path contains that segment
      but isn't under ``output_root`` (e.g. distilling a build copied
      from a different mount point — the corpus shape is preserved
      even though the absolute prefix differs).
    - The basename as a last resort (e.g. an input-PDF path with no
      ``documents/`` segment — the absolute prefix is provenance, not
      something the served bundle should carry, so keep the filename).
    """
    p = Path(s)
    if not p.is_absolute():
        return None
    try:
        return str(p.relative_to(output_root))
    except ValueError:
        parts = p.parts
        if "documents" in parts:
            i = parts.index("documents")
            return str(Path(*parts[i:]))
        return p.name


def _scrub_summary(serve_path: Path, output_root: Path) -> bool:
    """Rewrite absolute-path fields in a copied summary.json.

    Drops ``input_dir`` and ``output_directory`` (machine-specific
    provenance). Replaces ``processing_summary.original_pdf`` with the
    basename. Rewrites ``processing_summary.files_created`` entries to
    corpus-root-relative paths.

    Returns True if any change was made.
    """
    try:
        data = json.loads(serve_path.read_text())
    except Exception:
        return False
    changed = False
    for key in ("input_dir", "output_directory"):
        if key in data:
            data.pop(key)
            changed = True
    ps = data.get("processing_summary")
    if isinstance(ps, dict):
        op = ps.get("original_pdf")
        if isinstance(op, str) and Path(op).is_absolute():
            ps["original_pdf"] = Path(op).name
            changed = True
        fc = ps.get("files_created")
        if isinstance(fc, list):
            new_list: List[str] = []
            for entry in fc:
                if isinstance(entry, str):
                    rel = _to_corpus_relative(entry, output_root)
                    new_list.append(rel if rel is not None else entry)
                else:
                    new_list.append(entry)
            if new_list != fc:
                ps["files_created"] = new_list
                changed = True
    if changed:
        serve_path.write_text(json.dumps(data, indent=2, ensure_ascii=False))
    return changed


def _scrub_input_fingerprint_path(serve_path: Path) -> bool:
    """Rewrite the absolute-path field in any per-paper artifact that
    carries an ``input_fingerprint`` block (taxa.json + every lexicon
    ``<category>.json`` written by pipeline.annotate).

    ``input_fingerprint.path`` records the build-host path to the source
    file the artifact was derived from — e.g. taxonomy.sqlite for
    taxa.json, or lexicon.yaml for anatomy.json / biogeography.json /
    methods.json (#70). That path is not meaningful on the serve host;
    replace it with the basename so the portability audit passes while
    preserving the other fingerprint fields (sha256, size) for
    cache-invalidation checks.

    Returns True if any change was made.
    """
    try:
        data = json.loads(serve_path.read_text())
    except Exception:
        return False
    fp = data.get("input_fingerprint")
    if not isinstance(fp, dict):
        return False
    p = fp.get("path")
    if not isinstance(p, str) or not Path(p).is_absolute():
        return False
    fp["path"] = Path(p).name
    serve_path.write_text(json.dumps(data, indent=2, ensure_ascii=False))
    return True


def _scrub_figures(serve_path: Path, output_root: Path) -> bool:
    """Rewrite absolute-path fields in a copied figures.json.

    Rewrites ``figures[].file_path`` and the top-level
    ``figures_directory`` to corpus-root-relative form. Returns True
    if any change was made.
    """
    try:
        data = json.loads(serve_path.read_text())
    except Exception:
        return False
    changed = False
    for fig in data.get("figures") or []:
        if not isinstance(fig, dict):
            continue
        fp = fig.get("file_path")
        if isinstance(fp, str):
            rel = _to_corpus_relative(fp, output_root)
            if rel is not None and rel != fp:
                fig["file_path"] = rel
                changed = True
    fd = data.get("figures_directory")
    if isinstance(fd, str):
        rel = _to_corpus_relative(fd, output_root)
        if rel is not None and rel != fd:
            data["figures_directory"] = rel
            changed = True
    if changed:
        serve_path.write_text(json.dumps(data, indent=2, ensure_ascii=False))
    return changed


def _walk_strings(node) -> Iterable[str]:
    """Yield every string leaf in a parsed JSON tree (depth-first)."""
    if isinstance(node, str):
        yield node
    elif isinstance(node, dict):
        for v in node.values():
            yield from _walk_strings(v)
    elif isinstance(node, list):
        for v in node:
            yield from _walk_strings(v)


def _audit_no_absolute_paths(serve_dir: Path) -> List[Tuple[str, str]]:
    """Walk every JSON in the served bundle and flag string values that
    *are* absolute filesystem paths. Returns a list of
    ``(relative_path, offending_value)`` — empty when the bundle is clean.

    This is the §10 "no absolute paths in any served JSON" audit. Run
    after scrubbing; raises in the caller on non-empty result. Body-text
    fields (chunks.json text, intext_citations.json paragraph excerpts)
    can contain URLs with ``/home/...``-shaped path components — those
    don't trigger because the regex is anchored at start-of-string and
    the surrounding body text isn't itself a path.
    """
    offenders: List[Tuple[str, str]] = []
    for jp in sorted(serve_dir.rglob("*.json")):
        try:
            data = json.loads(jp.read_text())
        except Exception:
            continue
        for s in _walk_strings(data):
            if ' ' not in s and _ABS_PATH_RE.match(s):
                offenders.append((str(jp.relative_to(serve_dir)), s[:120]))
                break  # one offender per file is enough to flag it
    return offenders


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

    # #54 — load the set of hashes flagged works.serve = 0 in
    # biblio_authority.sqlite. These papers stay in the build bundle
    # (operators can still see them in `corpus status`) but get
    # excluded from the served bundle the MCP server reads.
    skipped_hashes = _load_skipped_hashes(output_dir / "biblio_authority.sqlite")
    if skipped_hashes:
        logger.info(
            "Skip flag: %d paper(s) marked works.serve = 0 will be "
            "excluded from the served bundle (#54)", len(skipped_hashes),
        )

    # #54 follow-up: re-distillation must prune any per-paper directory
    # that was copied previously but is now in skipped_hashes. Without
    # this, flipping `serve = false` and re-running `corpus run` would
    # leave the prior copy in serve_dir.
    if skipped_hashes and not dry_run:
        serve_documents_dir = serve_dir / "documents"
        if serve_documents_dir.is_dir():
            n_pruned = 0
            import shutil
            for h in skipped_hashes:
                stale = serve_documents_dir / h
                if stale.is_dir():
                    shutil.rmtree(stale)
                    n_pruned += 1
            if n_pruned:
                logger.info(
                    "Pruned %d stale per-paper dir(s) from serve_dir "
                    "(now flagged works.serve = 0; #54)", n_pruned,
                )

    # Per-paper files + figures/
    n_papers = 0
    n_skipped = 0
    n_files = 0
    total_bytes = 0
    for hash_dir in sorted(documents_dir.iterdir()):
        if not hash_dir.is_dir():
            continue
        if hash_dir.name in skipped_hashes:
            n_skipped += 1
            continue
        n_papers += 1
        dest_hash_dir = serve_dir / "documents" / hash_dir.name

        per_paper = list(PER_PAPER_FILES) + _per_paper_lexicon_outputs(hash_dir)
        for fname in per_paper:
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

    if n_skipped:
        logger.info(
            "Excluded %d paper(s) from the served bundle "
            "(works.serve = 0; #54)", n_skipped,
        )

    # LanceDB — directory tree
    vdb_src = output_dir / "vector_db" / "lancedb"
    vdb_dst = serve_dir / "vector_db" / "lancedb"
    nf, nb = _copy_tree(vdb_src, vdb_dst, dry_run)
    n_files += nf
    total_bytes += nb

    # #54 follow-up: filter LanceDB to drop chunks for skipped papers.
    # Without this, get_chunks_for_topic semantic search would still
    # surface excluded papers even though their per-paper directories
    # are gone.
    if skipped_hashes and not dry_run and vdb_dst.is_dir():
        n_dropped = _drop_chunks_for_hashes(vdb_dst, skipped_hashes)
        if n_dropped:
            logger.info(
                "Dropped %d LanceDB chunk(s) for skipped papers (#54)",
                n_dropped,
            )

    # Top-level corpuscle files — flat at the root of the output dir
    # alongside documents/, mirrored into the bundle root. All entries
    # here are expected; missing files emit a warning so the operator
    # notices before deploy. The earlier silent-skip behavior let v0.2
    # ship without biblio_authority.sqlite and taxon_mentions.sqlite.
    # The lexicon YAML is not bundled — extraction happens at process
    # time and lands in each paper's <category>.json artifacts.
    for fname in ("taxonomy.sqlite", "biblio_authority.sqlite",
                  "taxon_mentions.sqlite", "instructions.md"):
        src = output_dir / fname
        if not src.exists():
            logger.warning(
                "Expected top-level file %s missing from %s — "
                "bundle will not include it", fname, output_dir,
            )
            continue
        bw = _copy_file(src, serve_dir / fname, dry_run)
        if bw:
            n_files += 1
            total_bytes += bw

    # Path scrubbing (§10): rewrite absolute paths in copied summary.json
    # and figures.json to corpus-root-relative form, then audit the whole
    # served bundle to confirm no absolute paths leaked through. Skipped
    # in dry-run since the destination files don't exist.
    n_scrubbed = 0
    if not dry_run:
        for hash_dir in sorted((serve_dir / "documents").iterdir()):
            if not hash_dir.is_dir():
                continue
            if _scrub_summary(hash_dir / "summary.json", output_dir):
                n_scrubbed += 1
            if _scrub_figures(hash_dir / "figures.json", output_dir):
                n_scrubbed += 1
            # taxa.json + every lexicon <category>.json carry an
            # input_fingerprint.path → strip absolute prefix (#70).
            fingerprint_files = ["taxa.json"] + _per_paper_lexicon_outputs(hash_dir)
            for fname in fingerprint_files:
                if _scrub_input_fingerprint_path(hash_dir / fname):
                    n_scrubbed += 1
        offenders = _audit_no_absolute_paths(serve_dir)
        if offenders:
            logger.error("Absolute paths leaked into served bundle:")
            for rel, snippet in offenders[:10]:
                logger.error("  %s: …%s…", rel, snippet)
            raise RuntimeError(
                f"Absolute-path audit failed: {len(offenders)} file(s) "
                "still contain absolute paths. Update _scrub_summary / "
                "_scrub_figures (or add a new scrubber for the affected "
                "JSON) before re-running."
            )
        logger.info("Path scrub: rewrote %d files; audit clean.", n_scrubbed)

    # Manifest
    model, dim = _read_one_embedding_marker(output_dir)
    fig_count, chunk_count = _count_figures_and_chunks(
        documents_dir, excluded_hashes=skipped_hashes,
    )
    manifest = {
        "bundle_version": version,
        "created_at": _iso_now(),
        "pipeline_git_sha": _git_sha(),
        "embedding_model": model,
        "embedding_dim": dim,
        "taxonomy_snapshot_date": _taxonomy_snapshot_date(output_dir / "taxonomy.sqlite"),
        "paper_count": n_papers,
        "figure_count": fig_count,
        "chunk_count": chunk_count,
        "includes_pdfs": include_pdfs,
        # #61: stamp the resolved citation alongside version + git_sha so
        # downstream LLM clients (via the MCP `bundle_info` tool) can
        # cite the tool that gave them the literature without the user
        # having to dig for it. Single source: CITATION.cff at repo root.
        "citation": _citation_for_manifest(),
    }
    manifest_path = serve_dir / "bundle_manifest.json"
    if not dry_run:
        manifest_path.write_text(json.dumps(manifest, indent=2))

    manifest["_stats"] = {
        "n_files_copied": n_files,
        "total_bytes": total_bytes,
        "n_files_scrubbed": n_scrubbed if not dry_run else 0,
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
    parser.add_argument("--version", default=f"v{__version__}",
                        help="Bundle version string (e.g., v1.0.0). "
                             "Defaults to the pipeline's __version__.")
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
