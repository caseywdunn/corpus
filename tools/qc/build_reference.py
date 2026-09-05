"""Snapshot and compare build evidence without accepting a new baseline (#187).

This is an operator-side regression reference, not a fidelity score: unchanged
output can still be wrong. JSON fingerprints ignore only top-level producer
versions and relocate paths under the build root. References are also retained
field by field so equal counts cannot hide corrupted titles or authors.
Binary hashes are diagnostic (PDF timestamps and TEI IDs can change on rebuild).
Use the independent gold scorers beside this report.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
from pathlib import Path

from tools.qc.index_reference import figure_pixels, index_snapshot

SCHEMA_VERSION = 2
ARTIFACTS = ("text.json", "metadata.json", "references.json", "figures.json",
             "chunks.json", "taxa.json", "intext_citations.json")
BINARY_ARTIFACTS = ("processed.pdf", "grobid.tei.xml")
MANIFEST_FACTS = ("paper_count", "figure_count", "chunk_count", "embedding_model",
                  "embedding_dimension")


def _read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                    separators=(",", ":")).encode()).hexdigest()


def _file_digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _relocate(value, root):
    if isinstance(value, dict):
        return {key: _relocate(val, root) for key, val in value.items()}
    if isinstance(value, list):
        return [_relocate(val, root) for val in value]
    if isinstance(value, str) and value.startswith(str(root) + "/"):
        return "<build>/" + value[len(str(root)) + 1:]
    return value


def _logical(payload, root):
    return _relocate({k: v for k, v in payload.items()
                      if k not in ("pipeline_version", "schema_version")}, root)


def snapshot(build: Path) -> dict:
    """Read evidence only; missing/invalid primary files are explicit problems."""
    build = build.resolve()
    documents = build / "documents"
    if not documents.is_dir():
        raise ValueError(f"No documents directory: {documents}")
    papers, problems = {}, []
    for hd in sorted(documents.iterdir()):
        if not hd.is_dir():
            continue
        entry = {"artifacts": {}, "binary_sha256": {}, "references": None,
                 "quality_flags": [], "stage_failures": []}
        for name in ARTIFACTS:
            path = hd / name
            if not path.exists():
                entry["artifacts"][name] = None
                if name in ("text.json", "metadata.json", "references.json",
                            "figures.json", "chunks.json"):
                    problems.append(f"{hd.name}/{name}: missing")
                continue
            data = _read(path)
            logical = _logical(data, build)
            entry["artifacts"][name] = {
                "sha256": _digest(logical), "schema_version": data.get("schema_version"),
                "pipeline_version": data.get("pipeline_version"),
            }
            if name == "references.json":
                entry["references"] = logical.get("references", [])
        for name in BINARY_ARTIFACTS:
            path = hd / name
            entry["binary_sha256"][name] = _file_digest(path) if path.exists() else None
        summary = hd / "summary.json"
        if summary.exists():
            ps = _read(summary).get("processing_summary") or {}
            for key in ("quality_flags", "stage_failures"):
                entry[key] = _relocate(ps.get(key) or [], build)
        else:
            problems.append(f"{hd.name}/summary.json: missing")
        papers[hd.name] = entry
        entry["figure_pixels"] = figure_pixels(hd / "figures")
    if not papers:
        problems.append("No document artifacts")
    manifest_path = build / "_serve/bundle_manifest.json"
    manifest = _read(manifest_path) if manifest_path.exists() else {}
    return {
        "schema_version": SCHEMA_VERSION, "build_root": str(build),
        "producer": {key: manifest.get(key) for key in
                     ("pipeline_version", "pipeline_git_sha")},
        "manifest": {key: manifest.get(key) for key in MANIFEST_FACTS},
        "problems": problems, "documents": papers, "indexes": index_snapshot(build),
    }


def _reference_diff(before, after):
    """Ordinal alignment is deliberately conservative: reordering needs review."""
    changes = []
    for index in range(max(len(before), len(after))):
        old = before[index] if index < len(before) else None
        new = after[index] if index < len(after) else None
        if old == new:
            continue
        fields = {key: {"before": (old or {}).get(key), "after": (new or {}).get(key)}
                  for key in sorted(set(old or {}) | set(new or {}))
                  if (old or {}).get(key) != (new or {}).get(key)}
        changes.append({"ordinal": index, "before_present": old is not None,
                        "after_present": new is not None, "fields": fields})
    return changes


def compare(before: dict, after: dict) -> dict:
    for report in (before, after):
        if report.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("Incompatible build-reference schema; resnapshot both builds")
    old, new = before["documents"], after["documents"]
    changes, binary = {}, {}
    for sha in sorted(old.keys() & new.keys()):
        a, b = old[sha], new[sha]
        delta = {}
        changed = []
        for name in ARTIFACTS:
            x, y = a["artifacts"].get(name), b["artifacts"].get(name)
            if (x is None) != (y is None) or (
                x is not None and y is not None and
                (x["sha256"], x["schema_version"]) != (y["sha256"], y["schema_version"])
            ):
                changed.append(name)
        if changed:
            delta["artifacts"] = changed
        refs = _reference_diff(a.get("references") or [], b.get("references") or [])
        if refs:
            delta["references"] = refs
        if a["figure_pixels"] != b["figure_pixels"]:
            delta["figure_pixels"] = {"before": a["figure_pixels"], "after": b["figure_pixels"]}
        for key in ("quality_flags", "stage_failures"):
            if a[key] != b[key]:
                delta[key] = {"before": a[key], "after": b[key]}
        if delta:
            changes[sha] = delta
        bd = [name for name in BINARY_ARTIFACTS
              if a["binary_sha256"].get(name) != b["binary_sha256"].get(name)]
        if bd:
            binary[sha] = bd
    manifest = {key: {"before": before["manifest"].get(key),
                      "after": after["manifest"].get(key)}
                for key in MANIFEST_FACTS
                if before["manifest"].get(key) != after["manifest"].get(key)}
    hard_failures = {sha: {"stage_failures": paper["stage_failures"],
                           "quality_flags": [flag for flag in paper["quality_flags"]
                                             if flag.get("severity") == "error"]}
                     for sha, paper in new.items()
                     if paper["stage_failures"] or any(
                         flag.get("severity") == "error" for flag in paper["quality_flags"])}
    indexes = {key: {"before": before["indexes"].get(key), "after": after["indexes"].get(key)}
               for key in sorted(set(before["indexes"]) | set(after["indexes"]))
               if before["indexes"].get(key) != after["indexes"].get(key)}
    result = {"schema_version": SCHEMA_VERSION, "added": sorted(new.keys() - old.keys()),
              "removed": sorted(old.keys() - new.keys()), "changed": changes,
              "index_changes": indexes,
              "manifest_changes": manifest, "binary_changes": binary,
              "baseline_problems": before["problems"], "problems": after["problems"],
              "hard_failures": hard_failures,
              "unchanged_warning_documents": [sha for sha in sorted(new.keys() & old.keys())
                  if new[sha]["quality_flags"] and new[sha]["quality_flags"] == old[sha]["quality_flags"]]}
    result["requires_review"] = bool(changes or manifest or indexes or result["added"] or result["removed"]
                                      or before["problems"] or after["problems"] or hard_failures)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    take = sub.add_parser("snapshot")
    take.add_argument("build", type=Path)
    take.add_argument("--out", required=True, type=Path)
    diff = sub.add_parser("compare")
    diff.add_argument("before", type=Path)
    diff.add_argument("after", type=Path)
    diff.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    try:
        result = snapshot(args.build) if args.command == "snapshot" else compare(_read(args.before), _read(args.after))
        if args.out:
            # Never silently replace a reviewed baseline.
            with args.out.open("x", encoding="utf-8") as stream:
                json.dump(result, stream, ensure_ascii=False, indent=2)
                stream.write("\n")
        else:
            print(json.dumps(result, ensure_ascii=False, indent=2))
        return int(bool(result.get("requires_review") or result.get("problems")))
    except (OSError, ValueError, KeyError, TypeError, sqlite3.DatabaseError) as exc:
        parser.exit(2, f"{exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
