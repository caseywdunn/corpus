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
# #54 — worst-first triage workflow
# ---------------------------------------------------------------------------


def quality_flag_counts_per_paper(rollup: Dict[str, Any]) -> Counter:
    """Return ``{paper_hash: n_quality_flags}`` from the rollup.

    Backs ``--sort-by quality_flag_count --tail N`` (#54).
    """
    counts: Counter = Counter()
    for gate, hs in rollup["papers_by_gate"].items():
        for h in hs:
            counts[h] += 1
    return counts


def stage_failure_counts_per_paper(rollup: Dict[str, Any]) -> Counter:
    """Return ``{paper_hash: n_stage_failures}`` from the rollup."""
    counts: Counter = Counter()
    for (reason, stage), hs in rollup["papers_by_reason_stage"].items():
        for h in hs:
            counts[h] += 1
    return counts


def render_worst_first(
    rollup: Dict[str, Any],
    metric: str,
    tail: int,
) -> str:
    """Worst-N papers by ``metric``. Used by ``corpus status --sort-by …
    --tail N`` (#54)."""
    if metric == "quality_flag_count":
        per_paper = quality_flag_counts_per_paper(rollup)
        label = "quality flags"
    elif metric == "stage_failure_count":
        per_paper = stage_failure_counts_per_paper(rollup)
        label = "stage failures"
    else:
        return (
            f"unknown --sort-by metric: {metric!r}. Supported: "
            "quality_flag_count, stage_failure_count. Per-paper QC metrics "
            "(citation_resolution_rate, dictionary_hit_rate, …) deferred "
            "to follow-up; the gate count is the available proxy."
        )
    if not per_paper:
        return f"no papers carry {label}; nothing to triage."
    out: List[str] = [
        f"Worst {tail} paper(s) by {label}:",
    ]
    for h, n in per_paper.most_common(tail):
        out.append(f"  {h}  {n} {label}")
    return "\n".join(out)


def render_propose_skips(
    rollup: Dict[str, Any],
    output_dir: Path,
    min_flags: int = 2,
) -> str:
    """Papers flagged by ≥ ``min_flags`` quality gates that aren't yet
    ``works.serve = 0``. BibTeX-paste-ready output (#54)."""
    counts = quality_flag_counts_per_paper(rollup)
    biblio = output_dir / "biblio_authority.sqlite"
    skipped = _load_skipped_hashes(biblio)
    candidates = [
        (h, n) for h, n in counts.most_common() if n >= min_flags and h not in skipped
    ]
    if not candidates:
        return (
            f"no propose-skip candidates "
            f"(threshold: ≥{min_flags} quality flags, not already serve=0)."
        )
    out = [
        f"Propose-skip candidates (≥{min_flags} quality flags, not yet serve=0):",
        "",
        "% paste these `serve = {false}` lines into the matching BibTeX",
        "% entries; `corpus bib import` will write them to works.serve.",
        "",
    ]
    for h, n in candidates:
        gates_for_h = [g for g, hs in rollup["papers_by_gate"].items() if h in hs]
        out.append(f"  % {h}  {n} flags: {', '.join(sorted(gates_for_h))}")
        out.append(f"  serve = {{false}},")
        out.append("")
    return "\n".join(out)


def render_skipped(output_dir: Path) -> str:
    """Currently excluded papers + their reasons, grouped by reason (#54)."""
    biblio = output_dir / "biblio_authority.sqlite"
    if not biblio.is_file():
        return f"no biblio_authority.sqlite at {biblio}; can't list skipped papers."
    import sqlite3 as _sql
    try:
        conn = _sql.connect(f"file:{biblio}?mode=ro", uri=True)
        try:
            cols = {r[1] for r in conn.execute("PRAGMA table_info(works)")}
            if "serve" not in cols:
                return ("biblio_authority.sqlite predates v0.3 — works.serve "
                        "column missing. Re-run `corpus run` to migrate.")
            rows = conn.execute(
                "SELECT corpus_hash, COALESCE(serve_reason, '(no reason)'), title "
                "FROM works WHERE serve = 0 AND corpus_hash IS NOT NULL "
                "ORDER BY serve_reason, corpus_hash"
            ).fetchall()
        finally:
            conn.close()
    except _sql.Error as e:
        return f"could not read biblio_authority.sqlite: {e}"
    if not rows:
        return "no papers currently flagged works.serve = 0."
    grouped: Dict[str, List[tuple]] = defaultdict(list)
    for h, reason, title in rows:
        grouped[reason].append((h, title or "(no title)"))
    out = [f"Skipped papers ({len(rows)} total, grouped by reason):"]
    for reason in sorted(grouped):
        out.append(f"\n  [{reason}] ({len(grouped[reason])})")
        for h, title in grouped[reason]:
            out.append(f"    {h}  {title[:80]}")
    return "\n".join(out)


def _load_skipped_hashes(biblio_path: Path) -> set:
    """Mirror of mcpsrv.bundle._load_skipped_hashes; kept here so the
    status report doesn't depend on the mcpsrv package."""
    if not biblio_path.is_file():
        return set()
    import sqlite3 as _sql
    try:
        conn = _sql.connect(f"file:{biblio_path}?mode=ro", uri=True)
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
    except _sql.Error:
        return set()


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
    # #54 — worst-first triage workflow.
    parser.add_argument(
        "--sort-by", default=None,
        choices=["quality_flag_count", "stage_failure_count"],
        help="Worst-N papers by metric. Pair with --tail. Per-paper QC "
             "metrics (citation_resolution_rate, dictionary_hit_rate, …) "
             "land in a #54 follow-up; gate counts are the proxy for now.",
    )
    parser.add_argument(
        "--tail", type=int, default=20,
        help="With --sort-by, return the N worst papers (default 20).",
    )
    parser.add_argument(
        "--propose-skips", action="store_true",
        help="Papers flagged by ≥ N quality gates that aren't yet "
             "works.serve = 0 (--min-flags, default 2). BibTeX-paste-ready. "
             "(#54)",
    )
    parser.add_argument(
        "--min-flags", type=int, default=2,
        help="Threshold for --propose-skips (default 2).",
    )
    parser.add_argument(
        "--skipped", action="store_true",
        help="Currently excluded papers (works.serve = 0), grouped by "
             "serve_reason. Useful for comparing new candidates against "
             "precedent. (#54)",
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

    # #54 — triage modes short-circuit the default rollup output.
    if args.sort_by:
        print(render_worst_first(rollup, args.sort_by, args.tail))
        return 0
    if args.propose_skips:
        print(render_propose_skips(rollup, args.output_dir, args.min_flags))
        return 0
    if args.skipped:
        print(render_skipped(args.output_dir))
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
