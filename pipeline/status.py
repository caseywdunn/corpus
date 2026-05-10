#!/usr/bin/env python3
"""Corpus status report — single command for build state (#40).

Walks ``<output_dir>/documents/*/summary.json`` and rolls up:

  * Stage completion: success/total per stage (from ``stage_timings[]``)
  * Failures: counts grouped by ``reason_code`` × ``stage``
  * Quality flags: counts by gate (from ``quality_flags[]``)
  * Legacy errors: papers with the old ``errors[]`` field populated
    but no ``stage_failures[]`` (papers processed before #34 landed)

The data shape comes from #34 (structured ``stage_failures[]``,
``stage_timings[]``) and #36 (``quality_flags[]``). Older summary.json
files are still readable; their fields just don't contribute to the
new sections.

Usage:

    python corpus_status.py /path/to/output_dir
    python corpus_status.py /path/to/output_dir --json
    python corpus_status.py /path/to/output_dir --filter-reason timeout
    python corpus_status.py /path/to/output_dir --filter-gate gibberish_after_ocr --list-hashes

Pairs with #28 (granular per-stage resume — once stages are first-class,
this command's "stage completion" rows become the source of truth for
what to re-run) and the per-CLI ``--dry-run`` flags (#41 — show what a
re-run would do).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

logger = logging.getLogger("corpus_status")


# ---------------------------------------------------------------------------
# JSON helper (used by aggregation + stale-fingerprint detection)
# ---------------------------------------------------------------------------


def _safe_load_json(path: Path) -> Any:
    try:
        if path.exists() and path.stat().st_size > 0:
            with path.open(encoding="utf-8") as f:
                return json.load(f)
    except Exception as e:
        logger.warning("Could not read %s: %s", path, e)
    return None


# ---------------------------------------------------------------------------
# Stale-fingerprint detection (#29)
# ---------------------------------------------------------------------------


def _file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda: f.read(64 * 1024), b""):
            h.update(buf)
    return h.hexdigest()


# Stale-fingerprint detection lived here in v0.1; v0.2 moves the
# signal into ``pipeline_state.json`` (per-stage completion records
# stamped with PIPELINE_VERSION + input_fingerprint), so a plain
# ``--resume`` already re-runs whichever stages disagree with the
# current input. No extra detection step is required.


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


def _iter_summaries(documents_dir: Path) -> Iterable[Tuple[str, Dict]]:
    """Yield ``(hash, summary_dict)`` for every documents/<HASH>/summary.json
    that loads successfully.
    """
    for hash_dir in sorted(documents_dir.iterdir()):
        if not hash_dir.is_dir():
            continue
        summary = _safe_load_json(hash_dir / "summary.json")
        if isinstance(summary, dict):
            yield hash_dir.name, summary


def aggregate(documents_dir: Path) -> Dict[str, Any]:
    """Build the rollup. Pure function over summary.json contents."""
    rollup: Dict[str, Any] = {
        "documents_dir": str(documents_dir),
        "total_documents": 0,
        "documents_with_summary": 0,
        "stages": {},               # stage → {ok: N, fail: N, total: N}
        "failures_by_reason_stage": Counter(),  # (reason_code, stage) → N
        "papers_by_reason_stage": defaultdict(set),  # (reason, stage) → {hashes}
        "quality_flags": Counter(),  # gate → N
        "papers_by_gate": defaultdict(set),  # gate → {hashes}
        "papers_with_legacy_errors_only": [],
        "papers_with_failures": set(),
        "papers_with_quality_flags": set(),
    }

    for h, summary in _iter_summaries(documents_dir):
        rollup["total_documents"] += 1
        rollup["documents_with_summary"] += 1
        ps = summary.get("processing_summary") or {}

        timings = ps.get("stage_timings") or []
        failures = ps.get("stage_failures") or []
        quality_flags = ps.get("quality_flags") or []
        legacy_errors = ps.get("errors") or []

        for t in timings:
            stage = t.get("stage", "?")
            row = rollup["stages"].setdefault(stage, {"ok": 0, "fail": 0, "total": 0})
            row["total"] += 1
            if t.get("ok"):
                row["ok"] += 1
            else:
                row["fail"] += 1

        for f in failures:
            key = (f.get("reason_code", "?"), f.get("stage", "?"))
            rollup["failures_by_reason_stage"][key] += 1
            rollup["papers_by_reason_stage"][key].add(h)
            rollup["papers_with_failures"].add(h)

        for q in quality_flags:
            gate = q.get("gate", "?")
            rollup["quality_flags"][gate] += 1
            rollup["papers_by_gate"][gate].add(h)
            rollup["papers_with_quality_flags"].add(h)

        if legacy_errors and not failures and not timings:
            rollup["papers_with_legacy_errors_only"].append(h)

    return rollup


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _bar(n: int, total: int, width: int = 30) -> str:
    if total == 0:
        return ""
    filled = int(width * n / total)
    return "█" * filled + "░" * (width - filled)


def render_artifacts(output_dir: Path) -> str:
    """Cross-paper artifact presence section for `--report` (#57).

    Lists the four corpus-level outputs and ✓/✗ for each. The MCP
    server can run with any subset present; missing ones surface here.
    """
    out: List[str] = ["Cross-paper artifacts:"]
    for rel in ("biblio_authority.sqlite", "taxon_mentions.sqlite",
                "taxonomy.sqlite", "vector_db/lancedb"):
        p = output_dir / rel
        mark = "✓" if p.exists() else "✗"
        out.append(f"  {mark} {rel}")
    return "\n".join(out) + "\n"


def render_text(rollup: Dict[str, Any]) -> str:
    out: List[str] = []
    n_total = rollup["total_documents"]
    out.append(f"Corpus status — {rollup['documents_dir']} ({n_total} documents)")
    out.append("")

    # Stage completion
    stages = rollup["stages"]
    if stages:
        out.append("Stages (success / total):")
        # Order: by total desc, then alpha
        for stage, row in sorted(stages.items(), key=lambda kv: (-kv[1]["total"], kv[0])):
            ok, total = row["ok"], row["total"]
            pct = (100.0 * ok / total) if total else 0.0
            out.append(f"  {stage:<28s} {ok:>5d} / {total:<5d}  ({pct:5.1f}%)  {_bar(ok, total)}")
    else:
        out.append("Stages: no stage_timings recorded yet (older summary.json files; #34 was added in v0.2).")
    out.append("")

    # Failures
    failures = rollup["failures_by_reason_stage"]
    n_fail_papers = len(rollup["papers_with_failures"])
    if failures:
        n_fail_total = sum(failures.values())
        out.append(f"Failures ({n_fail_papers} papers, {n_fail_total} stage_failures):")
        for (reason, stage), count in failures.most_common():
            out.append(f"  {reason} × {stage:<28s} {count:>5d}")
    else:
        out.append("Failures: none recorded.")
    out.append("")

    # Quality flags
    qf = rollup["quality_flags"]
    n_qf_papers = len(rollup["papers_with_quality_flags"])
    if qf:
        n_qf_total = sum(qf.values())
        out.append(f"Quality flags ({n_qf_papers} papers, {n_qf_total} flags):")
        for gate, count in qf.most_common():
            out.append(f"  {gate:<32s} {count:>5d}")
    else:
        out.append("Quality flags: none recorded.")
    out.append("")

    # Legacy errors
    legacy = rollup["papers_with_legacy_errors_only"]
    if legacy:
        out.append(
            f"Note: {len(legacy)} paper(s) have legacy errors[] only "
            "(processed before #34). Re-run those to populate the structured "
            "stage_failures[] / stage_timings[]."
        )

    return "\n".join(out)


def render_json(rollup: Dict[str, Any]) -> str:
    """JSON output. Sets are converted to sorted lists; Counter / defaultdict
    to plain dicts so the output is stable.
    """
    j = dict(rollup)
    # Counter(tuple → int) → list of {reason_code, stage, count}
    j["failures_by_reason_stage"] = [
        {"reason_code": r, "stage": s, "count": c}
        for (r, s), c in rollup["failures_by_reason_stage"].most_common()
    ]
    j["papers_by_reason_stage"] = [
        {"reason_code": r, "stage": s, "hashes": sorted(hs)}
        for (r, s), hs in rollup["papers_by_reason_stage"].items()
    ]
    j["quality_flags"] = [
        {"gate": g, "count": c} for g, c in rollup["quality_flags"].most_common()
    ]
    j["papers_by_gate"] = {g: sorted(hs) for g, hs in rollup["papers_by_gate"].items()}
    j["papers_with_failures"] = sorted(rollup["papers_with_failures"])
    j["papers_with_quality_flags"] = sorted(rollup["papers_with_quality_flags"])
    return json.dumps(j, indent=2)


# ---------------------------------------------------------------------------
# Filters → hash list
# ---------------------------------------------------------------------------


def filtered_hashes(
    rollup: Dict[str, Any],
    filter_stage: Optional[str] = None,
    filter_reason: Optional[str] = None,
    filter_gate: Optional[str] = None,
) -> List[str]:
    """Return the sorted list of hashes that match the active filter(s).

    With multiple filters, hashes must match all of them (intersection).
    With no filter, returns the union of papers-with-failures-or-flags.
    """
    selected: Optional[set] = None

    def _intersect(s: set):
        nonlocal selected
        selected = s if selected is None else (selected & s)

    if filter_stage or filter_reason:
        matched: set = set()
        for (r, s), hs in rollup["papers_by_reason_stage"].items():
            if filter_reason and r != filter_reason:
                continue
            if filter_stage and s != filter_stage:
                continue
            matched |= hs
        _intersect(matched)

    if filter_gate:
        _intersect(rollup["papers_by_gate"].get(filter_gate, set()))

    if selected is None:
        selected = (
            rollup["papers_with_failures"] | rollup["papers_with_quality_flags"]
        )
    return sorted(selected)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "output_dir", type=Path,
        help="Corpus output directory (contains documents/<HASH>/ subdirs)",
    )
    parser.add_argument(
        "--json", action="store_true",
        help="Emit the rollup as JSON instead of the text report.",
    )
    parser.add_argument(
        "--report", action="store_true",
        help="Print the full report (default rollup + cross-paper artifact "
             "presence) — also written to <output_dir>/run.log on `corpus run` "
             "completion. (#57)",
    )
    parser.add_argument(
        "--list-hashes", action="store_true",
        help="Print one hash per line for papers matching --filter-* "
             "(suitable for `xargs`). Combine filters to narrow.",
    )
    parser.add_argument(
        "--filter-stage", default=None,
        help="Only papers with a failure at this stage "
             "(e.g. metadata_extraction).",
    )
    parser.add_argument(
        "--filter-reason", default=None,
        help="Only papers with this reason_code (e.g. timeout, crash, "
             "external_unavailable, unsupported_format, corrupted, "
             "quality_gate).",
    )
    parser.add_argument(
        "--filter-gate", default=None,
        help="Only papers with this quality_flag gate "
             "(e.g. gibberish_after_ocr, empty_text).",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    documents_dir = args.output_dir / "documents"
    if not documents_dir.is_dir():
        logger.error("Not a corpus output dir: %s (no documents/)", args.output_dir)
        return 1

    rollup = aggregate(documents_dir)

    if args.list_hashes:
        for h in filtered_hashes(
            rollup,
            filter_stage=args.filter_stage,
            filter_reason=args.filter_reason,
            filter_gate=args.filter_gate,
        ):
            print(h)
        return 0

    if args.json:
        print(render_json(rollup))
    else:
        print(render_text(rollup))
        if args.report:
            print()
            print(render_artifacts(args.output_dir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
