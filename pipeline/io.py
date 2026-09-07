"""Hashing, orphan audit, summary writing, and per-document directory layout.

* :func:`calculate_pdf_hash` / :func:`short_hash` / :data:`HASH_PREFIX_LEN`
  — content-addressed identification.
* :func:`find_all_pdfs` — recursive discovery + deduplication by hash.
* :func:`audit_orphans` — read-only audit of stale ``documents/<HASH>/``
  directories whose source PDF is no longer in the input set (#31).
* :func:`create_output_structure` — make documents/ + vector_db/.
* :func:`get_relative_paths` / :func:`create_summary_json` —
  per-paper provenance writeback.
* :func:`_verify_or_raise_collision` — resume-time check that a hash
  prefix collision isn't silently overwriting a different PDF.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import time
from pathlib import Path
from typing import Dict, List, Optional

from . import stamp_artifact

logger = logging.getLogger(__name__)


HASH_PREFIX_LEN = 12  # 48 bits; collision-safe up to ~1e7 documents


def calculate_pdf_hash(pdf_path: Path) -> str:
    """Calculate the full SHA-256 hex digest of a PDF file.

    The per-document directory uses a 12-char lowercase prefix of this digest
    (see :data:`HASH_PREFIX_LEN`); the full digest is recorded in each
    ``summary.json`` so we can verify on resume that a re-encountered prefix
    really identifies the same PDF (rather than a hash-prefix collision).
    """
    hasher = hashlib.sha256()
    with open(pdf_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def short_hash(full_hash: str) -> str:
    """Return the 12-char lowercase prefix used as the per-document dir name."""
    return full_hash[:HASH_PREFIX_LEN].lower()


def find_all_pdfs(
    input_dir: Path,
    exclude_under: Optional[Path] = None,
    *,
    strict: bool = False,
) -> Dict[str, List[Path]]:
    """Recursively find all PDFs under ``input_dir`` and group by full SHA-256.

    Returns a dict mapping full hex digest → list of paths that share it
    (duplicates). The caller derives the short directory name via
    :func:`short_hash`.

    ``exclude_under`` skips PDFs whose resolved path is under the given
    directory. Required when ``input_dir`` is an ancestor of ``output_dir``
    (e.g., the demo's ``input_pdfs: .`` with ``output_dir: ./output``):
    each per-paper ``documents/<HASH>/processed.pdf`` is a different PDF
    (OCR adds a text layer) and would otherwise be ingested as a new
    document on every re-run, doubling the corpus.

    ``strict`` raises on missing/unreadable directories or PDFs. Retirement
    must use it: a partial inventory cannot prove a source was removed.
    """
    excl = exclude_under.resolve() if exclude_under is not None else None
    if strict and not input_dir.is_dir():
        raise RuntimeError(f"Cannot inventory missing input directory: {input_dir}")

    def scan_error(error):
        if strict:
            raise RuntimeError(f"Cannot inventory all source PDFs: {error}") from error
        logger.warning("Cannot scan source directory: %s", error)

    def candidates():
        for root, dirs, files in os.walk(input_dir, onerror=scan_error):
            if excl is not None:
                dirs[:] = [d for d in dirs if not (Path(root) / d).resolve().is_relative_to(excl)]
            dirs.sort()
            for name in sorted(files):
                if name.endswith(".pdf"):
                    yield Path(root) / name

    pdf_map: Dict[str, List[Path]] = {}
    for pdf_path in candidates():
        if not pdf_path.is_file():
            if strict:
                raise RuntimeError(f"Cannot inventory source PDF: {pdf_path}")
            continue
        if excl is not None:
            try:
                pdf_path.resolve().relative_to(excl)
                continue  # under the exclusion dir — skip
            except ValueError:
                pass
        try:
            full_hash = calculate_pdf_hash(pdf_path)
            pdf_map.setdefault(full_hash, []).append(pdf_path)
        except Exception as e:
            if strict:
                raise RuntimeError(f"Cannot hash source PDF {pdf_path}: {e}") from e
            logger.warning("Could not hash %s: %s", pdf_path, e)
    return pdf_map


def prune_orphans(
    input_dir: Path,
    output_dir: Path,
    *,
    dry_run: bool = False,
    force: bool = False,
    safety_pct: float = 0.25,
) -> dict:
    """Retire orphaned document directories and remove vector rows (#265).

    An orphan is a hash directory whose source PDF is no longer in the
    configured input set, or a vector-index row whose hash no longer
    has a backing dir.

    Safety rail: refuses to prune when more than ``safety_pct`` of the
    hash directories would be removed (likely a config mistake or
    unmounted volume). ``force=True`` bypasses the rail.

    Returns a dict ``{doc_pruned, vec_pruned, doc_total, would_remove}``.
    With ``dry_run=True`` nothing is moved; the same dict reports the
    counts that *would* have been removed.
    """
    documents_dir = output_dir / "documents"
    if not documents_dir.is_dir():
        return {"doc_pruned": 0, "vec_pruned": 0, "doc_total": 0, "would_remove": []}

    input_pdf_map = find_all_pdfs(input_dir, exclude_under=output_dir, strict=True)
    input_hashes = {short_hash(h) for h in input_pdf_map}
    doc_hashes = {p.name for p in sorted(documents_dir.iterdir()) if p.is_dir()}
    doc_orphans = sorted(doc_hashes - input_hashes)
    doc_total = len(doc_hashes)

    # Safety rail
    if doc_total > 0 and not force:
        pct = len(doc_orphans) / doc_total
        if pct > safety_pct:
            raise RuntimeError(
                f"orphan-prune safety rail: {len(doc_orphans)} of {doc_total} "
                f"hash dirs ({pct:.0%}) would be removed (threshold "
                f"{safety_pct:.0%}). Likely a config mistake or unmounted "
                f"volume. Re-run with --force-prune to override."
            )

    n_vec_pruned = 0
    n_doc_pruned = 0
    if not dry_run:
        # First drop stale vectors. Failure aborts before retiring documents;
        # never publish an apparently successful run with stale search rows.
        # Schema: pipeline/embed.py:ChunkMetadata stores the hash as
        # ``metadata.pdf_hash`` (nested), not a flat ``hash`` column.
        # The LanceDB tables live under ``vector_db/lancedb/`` —
        # ``pipeline.embed`` writes there at line 312.
        vector_db_path = output_dir / "vector_db" / "lancedb"
        if vector_db_path.is_dir():
            try:
                import lancedb  # type: ignore
                db = lancedb.connect(str(vector_db_path))
                from .embeddings import lancedb_table_names
                if "document_chunks" in lancedb_table_names(db):
                    table = db.open_table("document_chunks")
                    surviving_hashes = doc_hashes - set(doc_orphans)
                    if surviving_hashes:
                        before = table.count_rows()
                        in_clause = ", ".join("'" + h.replace("'", "''") + "'" for h in surviving_hashes)
                        table.delete(f"metadata.pdf_hash NOT IN ({in_clause})")
                        n_vec_pruned = before - table.count_rows()
                    else:
                        # No surviving docs — drop the whole table.
                        n_vec_pruned = table.count_rows()
                        db.drop_table("document_chunks")
            except Exception as e:
                raise RuntimeError(f"Could not prune LanceDB; update aborted: {e}") from e
        # Keep superseded extraction evidence recoverable, outside the active
        # documents tree. The output exclusion above also excludes this archive.
        if doc_orphans:
            import tempfile
            retired_root = output_dir / ".retired"
            retired_root.mkdir(exist_ok=True)
            archive = Path(tempfile.mkdtemp(prefix="documents-", dir=retired_root))
            for h in doc_orphans:
                (documents_dir / h).rename(archive / h)
                marker = output_dir / "vector_db" / f"{h}_embedded.done"
                if marker.exists():
                    marker.rename(archive / f"{h}_embedded.done")
                n_doc_pruned += 1
            logger.warning("Retired %d documents to %s (recoverable)", n_doc_pruned, archive)
    else:
        n_doc_pruned = len(doc_orphans)

    return {
        "doc_pruned": n_doc_pruned,
        "vec_pruned": n_vec_pruned,
        "doc_total": doc_total,
        "would_remove": doc_orphans,
    }


def audit_orphans(input_dir: Path, output_dir: Path) -> int:
    """Read-only orphan audit (#31). Returns count of orphans found.

    Two orphan classes:

    1. Document orphans — ``documents/<HASH>/`` whose hash is no longer
       present in the input set.
    2. Vector-index orphans — LanceDB rows whose ``hash`` column has
       no corresponding ``documents/<HASH>/`` directory (which implies
       the document was deleted but its embeddings linger).

    Re-hashes input PDFs to make the check path-independent — moving or
    renaming the input dir does not produce false orphans.
    """
    documents_dir = output_dir / "documents"
    if not documents_dir.is_dir():
        logger.error("No documents/ directory at %s", documents_dir)
        return 0

    logger.info("Hashing input PDFs under %s …", input_dir)
    input_pdf_map = find_all_pdfs(input_dir, exclude_under=output_dir, strict=True)
    input_hashes = {short_hash(h) for h in input_pdf_map}
    logger.info("Found %d unique input PDFs", len(input_hashes))

    doc_hashes = {p.name for p in sorted(documents_dir.iterdir()) if p.is_dir()}
    logger.info("Found %d hash directories under documents/", len(doc_hashes))

    document_orphans = sorted(doc_hashes - input_hashes)
    print()
    print(f"=== Document orphans ({len(document_orphans)}) ===")
    print("(hash directories whose source PDF is no longer in the input set)")
    if not document_orphans:
        print("  (none)")
    else:
        for h in document_orphans:
            summary_file = documents_dir / h / "summary.json"
            last_known = ""
            if summary_file.exists():
                try:
                    with summary_file.open() as f:
                        s = json.load(f)
                    rels = s.get("relative_paths") or []
                    if rels:
                        last_known = f"  (last known: {rels[0]}"
                        if len(rels) > 1:
                            last_known += f" + {len(rels) - 1} more"
                        last_known += ")"
                except Exception:
                    pass
            print(f"  {h}{last_known}")

    vector_orphans: List[str] = []
    vector_db_path = output_dir / "vector_db" / "lancedb"
    if vector_db_path.is_dir():
        try:
            import lancedb  # type: ignore
            db = lancedb.connect(str(vector_db_path))
            if "document_chunks" in db.list_tables():
                table = db.open_table("document_chunks")
                # Schema: pipeline/embed.py:ChunkMetadata stores the
                # hash as ``metadata.pdf_hash`` (nested), not a flat
                # ``hash`` column. select(["metadata"]) returns the
                # whole nested struct.
                hashes_in_table = {
                    (row.get("metadata") or {}).get("pdf_hash")
                    for row in table.search().select(["metadata"]).limit(10**9).to_list()
                    if (row.get("metadata") or {}).get("pdf_hash")
                }
                vector_orphans = sorted(hashes_in_table - doc_hashes)
        except ImportError:
            logger.info("lancedb not importable; skipping vector-index audit")
        except Exception as e:
            logger.warning("Could not audit LanceDB: %s", e)

    print()
    print(f"=== Vector-index orphans ({len(vector_orphans)}) ===")
    print("(LanceDB hashes with no documents/<HASH>/ directory)")
    if not vector_orphans:
        print("  (none)")
    else:
        for h in vector_orphans:
            print(f"  {h}")

    total = len(document_orphans) + len(vector_orphans)
    print()
    print(f"Total orphans: {total}. Audit is read-only — nothing was deleted.")
    return total


def create_output_structure(output_dir: Path):
    """Create the output directory structure."""
    documents_dir = output_dir / "documents"
    vector_db_dir = output_dir / "vector_db"
    
    documents_dir.mkdir(parents=True, exist_ok=True)
    vector_db_dir.mkdir(parents=True, exist_ok=True)
    
    return documents_dir, vector_db_dir


def get_relative_paths(pdf_paths: List[Path], input_dir: Path) -> List[str]:
    """Get relative paths of PDFs from input directory."""
    return [str(path.relative_to(input_dir)) for path in pdf_paths]

def create_summary_json(
    pdf_hash_full: str,
    pdf_paths: List[Path],
    input_dir: Path,
    hash_dir: Path,
    processing_summary: Dict,
):
    """Write ``summary.json`` for one document.

    Records both the short directory prefix and the full SHA-256 so that
    ``--resume`` can verify that a re-encountered prefix refers to the same
    PDF (not a hash-prefix collision).
    """
    relative_paths = get_relative_paths(pdf_paths, input_dir)

    summary = {
        "pdf_hash": short_hash(pdf_hash_full),
        "pdf_hash_full": pdf_hash_full,
        "hash_algorithm": "sha256",
        "input_dir": str(input_dir),
        "relative_paths": relative_paths,
        "total_copies_found": len(pdf_paths),
        "processing_summary": processing_summary,
        "output_directory": str(hash_dir),
    }

    # Atomic write via per-writer tmp + rename. Per-writer tmp name
    # (``.tmp.<pid>.<ns>``) keeps two concurrent writers in the same
    # hash_dir from corrupting each other's payload — see #55. Without
    # the rename, a reader hitting the partial write also saw garbage.
    summary_file = hash_dir / "summary.json"
    tmp = summary_file.with_suffix(
        f"{summary_file.suffix}.tmp.{os.getpid()}.{time.monotonic_ns()}"
    )
    with open(tmp, "w") as f:
        json.dump(stamp_artifact(summary), f, indent=2)
    tmp.replace(summary_file)

    return summary_file


def refresh_summary_source_paths(pdf_hash_full: str, pdf_paths: List[Path],
                                 input_dir: Path, hash_dir: Path) -> bool:
    """Refresh cheap source inventory even when every extraction stage skips.

    Directory moves and added/removed identical copies do not change metadata's
    basename input, but they do change the paths stored in the vector index.
    Keep extraction receipts/timings intact and avoid rewriting a true no-op.
    """
    summary = json.loads((hash_dir / "summary.json").read_text(encoding="utf-8"))
    if (summary.get("relative_paths") == get_relative_paths(pdf_paths, input_dir)
            and summary.get("input_dir") == str(input_dir)
            and summary.get("total_copies_found") == len(pdf_paths)):
        return False
    create_summary_json(pdf_hash_full, pdf_paths, input_dir, hash_dir,
                        summary.get("processing_summary") or {})
    return True


def _verify_or_raise_collision(hash_dir: Path, pdf_hash_full: str) -> Optional[bool]:
    """If ``hash_dir/summary.json`` exists, verify its recorded full hash
    matches ``pdf_hash_full``. Returns True if it matches (resume-safe), False
    if no summary is present (fresh dir), and raises ``RuntimeError`` on a
    real hash-prefix collision.
    """
    summary_file = hash_dir / "summary.json"
    if not summary_file.exists():
        return False
    try:
        with open(summary_file, "r") as f:
            existing = json.load(f)
    except Exception as e:
        logger.warning("Could not read %s (%s); treating as incomplete", summary_file, e)
        return False
    existing_full = existing.get("pdf_hash_full")
    if existing_full is None:
        # Legacy summary from before we recorded full hashes. Trust it but warn.
        logger.warning(
            "Existing summary at %s has no pdf_hash_full; cannot verify "
            "against prefix collision. Treating as a match.",
            summary_file,
        )
        return True
    if existing_full != pdf_hash_full:
        raise RuntimeError(
            f"Hash-prefix collision detected at {hash_dir}: "
            f"existing summary records full hash {existing_full!r} but this "
            f"PDF hashes to {pdf_hash_full!r}. Increase HASH_PREFIX_LEN or "
            f"investigate duplicate inputs."
        )
    return True
