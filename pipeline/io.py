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


def find_all_pdfs(input_dir: Path) -> Dict[str, List[Path]]:
    """Recursively find all PDFs under ``input_dir`` and group by full SHA-256.

    Returns a dict mapping full hex digest → list of paths that share it
    (duplicates). The caller derives the short directory name via
    :func:`short_hash`.
    """
    pdf_map: Dict[str, List[Path]] = {}
    for pdf_path in input_dir.rglob("*.pdf"):
        if not pdf_path.is_file():
            continue
        try:
            full_hash = calculate_pdf_hash(pdf_path)
            pdf_map.setdefault(full_hash, []).append(pdf_path)
        except Exception as e:
            logger.warning("Could not hash %s: %s", pdf_path, e)
    return pdf_map


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
    input_pdf_map = find_all_pdfs(input_dir)
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
    vector_db_path = output_dir / "vector_db"
    if vector_db_path.is_dir():
        try:
            import lancedb  # type: ignore
            db = lancedb.connect(str(vector_db_path))
            if "document_chunks" in db.table_names():
                table = db.open_table("document_chunks")
                hashes_in_table = {
                    row["hash"]
                    for row in table.search().select(["hash"]).limit(10**9).to_list()
                    if row.get("hash")
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

    summary_file = hash_dir / "summary.json"
    with open(summary_file, "w") as f:
        json.dump(stamp_artifact(summary), f, indent=2)

    return summary_file


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

