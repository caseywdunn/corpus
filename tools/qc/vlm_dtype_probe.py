#!/usr/bin/env python3
"""Compare local-VLM panel detection across weight dtypes (#258).

`figures.vision_dtype` defaults to `auto`, which is float32 everywhere
except CUDA — ~30.4 GB of weights for Qwen2.5-VL-7B against ~15.2 GB at
half precision, so on Apple Silicon it does not fit a 32 GB machine at
all. Apple Silicon supports bf16 and fp16. What nobody has established is
whether Qwen2.5-VL is numerically *sound* in either on Metal, and which
of the two this model prefers.

The default does not move on a mocked test. It moves on this script's
output.

**ROI equivalence is the criterion, not memory.** Half precision that
quietly finds fewer or sloppier panels is worse than float32 that does
not fit: the corpuscle still builds, every gate passes, and nothing says
the panels got worse. So the report leads with per-figure ROI agreement
against the float32 baseline and treats memory and wall clock as
secondary.

Usage
-----
    python tools/qc/vlm_dtype_probe.py <corpuscle_output_dir>

    # fewer figures / different dtypes / save the report
    python tools/qc/vlm_dtype_probe.py <dir> --figures 6 \\
        --dtypes float32,bfloat16,float16 --out vlm_dtype_report.md

Each dtype loads the model once and runs every selected figure through
it, so cost is one model load per dtype plus a few seconds per figure. A
dtype that cannot load is recorded as a finding rather than a crash —
"float32 does not fit" is one of the answers this is looking for.

Paste the report into
https://github.com/caseywdunn/corpus/issues/258.
"""
from __future__ import annotations

import argparse
import json
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

DEFAULT_DTYPES = ("float32", "bfloat16", "float16")


# ── environment ────────────────────────────────────────────────────────

def machine_report() -> Dict[str, str]:
    """What a reader of the report needs to interpret the numbers."""
    info = {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": platform.python_version(),
    }
    try:
        info["cpu"] = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"],
            capture_output=True, text=True, timeout=10,
        ).stdout.strip() or "unknown"
        total = subprocess.run(
            ["sysctl", "-n", "hw.memsize"],
            capture_output=True, text=True, timeout=10,
        ).stdout.strip()
        info["ram_gb"] = f"{int(total) / 1e9:.0f}" if total.isdigit() else "unknown"
    except Exception:
        info.setdefault("cpu", "unknown")
        info.setdefault("ram_gb", "unknown")
    try:
        import torch
        info["torch"] = torch.__version__
        info["mps_available"] = str(
            hasattr(torch.backends, "mps") and torch.backends.mps.is_available()
        )
        info["cuda_available"] = str(torch.cuda.is_available())
    except Exception as exc:
        info["torch"] = f"unavailable: {exc}"
    try:
        import transformers
        info["transformers"] = transformers.__version__
    except Exception as exc:
        info["transformers"] = f"unavailable: {exc}"
    return info


def peak_rss_gb() -> float:
    """Peak resident set size of this process, in GB.

    `ru_maxrss` is bytes on macOS and kilobytes on Linux — the one place
    this script has to know which platform it is on.
    """
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(raw / 1e9 if sys.platform == "darwin" else raw / 1e6, 2)


# ── figure selection ───────────────────────────────────────────────────

def select_figures(output_dir: Path, limit: int) -> List[Dict]:
    """Figures worth asking a panel detector about, deterministically.

    Prefers figures whose caption already names panels: those have an
    expected label set, so a dtype that loses a panel is visible as a
    miss rather than as a judgement call. Sorted by (hash, figure_id) so
    every dtype sees the same figures in the same order.
    """
    candidates, fallback = [], []
    documents = output_dir / "documents"
    if not documents.is_dir():
        raise SystemExit(
            f"No documents/ under {output_dir}.\n"
            f"Build a corpuscle first, e.g.:\n"
            f"    cd demo && corpus run --only extract"
        )
    for hd in sorted(documents.iterdir()):
        figures_json = hd / "figures.json"
        if not figures_json.is_file():
            continue
        try:
            figures = json.loads(figures_json.read_text()).get("figures", [])
        except (OSError, ValueError):
            continue
        for fig in figures:
            name = fig.get("filename") or Path(fig.get("file_path") or "").name
            if not name:
                continue
            image = hd / "figures" / name
            if not image.is_file():
                continue
            entry = {
                "hash": hd.name,
                "figure_id": fig.get("figure_id") or name,
                "image": image,
                "caption": fig.get("caption_text") or "",
                "expected": list(fig.get("panels_from_caption") or []),
                # Whatever produced this corpuscle's ROIs — in the shipped
                # fixture, production Pass 3b on an H200 at bfloat16. A
                # same-machine float32 baseline says "half precision agrees
                # with float32 here"; this says "it agrees with the answer
                # the cluster actually shipped", which is the stronger claim.
                "reference_rois": [
                    {"label": r.get("label"), "roi_px": r.get("roi_px")}
                    for r in (fig.get("rois") or [])
                ],
                "reference_source": next(
                    (r.get("source") for r in (fig.get("rois") or [])
                     if r.get("source")), None),
            }
            (candidates if entry["expected"] else fallback).append(entry)
    chosen = sorted(candidates, key=lambda e: (e["hash"], e["figure_id"]))[:limit]
    if len(chosen) < limit:
        chosen += sorted(fallback, key=lambda e: (e["hash"], e["figure_id"]))[
            : limit - len(chosen)]
    if not chosen:
        raise SystemExit(f"No figure images found under {documents}")
    return chosen


# ── one dtype ──────────────────────────────────────────────────────────

def probe_dtype(dtype: str, figures: List[Dict], device: Optional[str]) -> Dict:
    """Load the model at ``dtype`` and detect panels on every figure."""
    from pipeline.vision import LocalVLMBackend

    result: Dict = {"dtype": dtype, "loaded": False, "error": None,
                    "load_seconds": None, "detect_seconds": None,
                    "peak_rss_gb": None, "figures": {}}
    t0 = time.monotonic()
    try:
        backend = LocalVLMBackend(dtype=dtype, device=device)
    except BaseException as exc:      # OOM here is a finding, not a crash
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["peak_rss_gb"] = peak_rss_gb()
        return result
    result["loaded"] = True
    result["load_seconds"] = round(time.monotonic() - t0, 1)
    result["device"] = backend._device

    t1 = time.monotonic()
    for entry in figures:
        key = f"{entry['hash']}/{entry['figure_id']}"
        try:
            rois = backend.detect_figure_panels(
                entry["image"], entry["caption"], entry["expected"])
            result["figures"][key] = {
                "rois": [{"label": r.get("label"), "roi_px": r.get("roi_px")}
                         for r in rois],
                "error": None,
            }
        except BaseException as exc:
            result["figures"][key] = {"rois": None,
                                      "error": f"{type(exc).__name__}: {exc}"}
    result["detect_seconds"] = round(time.monotonic() - t1, 1)
    result["peak_rss_gb"] = peak_rss_gb()
    return result


# ── comparison ─────────────────────────────────────────────────────────

def _iou(a, b) -> float:
    """Intersection over union of two ``[x0, y0, x1, y1]`` boxes.

    Anything unusable scores 0.0 rather than raising: a dtype that
    returns a malformed ROI has to show up as *disagreement*, not as a
    traceback that loses the whole run — the run is expensive and the
    machine may not be available again.
    """
    try:
        ax0, ay0, ax1, ay1 = (float(v) for v in list(a)[:4])
        bx0, by0, bx1, by1 = (float(v) for v in list(b)[:4])
    except (TypeError, ValueError):
        return 0.0
    ix0, iy0 = max(ax0, bx0), max(ay0, by0)
    ix1, iy1 = min(ax1, bx1), min(ay1, by1)
    if ix1 <= ix0 or iy1 <= iy0:
        return 0.0
    inter = (ix1 - ix0) * (iy1 - iy0)
    union = ((ax1 - ax0) * (ay1 - ay0) + (bx1 - bx0) * (by1 - by0) - inter)
    return round(inter / union, 3) if union > 0 else 0.0


def compare(baseline: Dict, other: Dict) -> Dict:
    """Per-figure agreement between two dtype runs.

    Boxes are matched by panel label where both runs produced one,
    because that is the identity a reader cares about — "panel B moved"
    is a different statement from "there is one fewer box".
    """
    rows = []
    for key, base in (baseline.get("figures") or {}).items():
        theirs = (other.get("figures") or {}).get(key) or {}
        b_rois = base.get("rois") or []
        o_rois = theirs.get("rois") or []
        by_label = {r.get("label"): r.get("roi_px") for r in o_rois}
        ious = [
            _iou(r.get("roi_px"), by_label.get(r.get("label")))
            for r in b_rois if r.get("label") in by_label
        ]
        rows.append({
            "figure": key,
            "baseline_rois": len(b_rois),
            "other_rois": len(o_rois),
            "matched_labels": len(ious),
            "min_iou": min(ious) if ious else None,
            "mean_iou": round(sum(ious) / len(ious), 3) if ious else None,
            "error": theirs.get("error") or base.get("error"),
        })
    same_count = sum(1 for r in rows if r["baseline_rois"] == r["other_rois"])
    tight = [r for r in rows if r["min_iou"] is not None and r["min_iou"] >= 0.95]
    return {"rows": rows, "figures": len(rows),
            "same_roi_count": same_count, "tight_boxes": len(tight)}


# ── report ─────────────────────────────────────────────────────────────

def render(machine: Dict, runs: List[Dict], comparisons: Dict,
           against_reference: Optional[Dict] = None,
           reference_source: Optional[str] = None) -> str:
    out = ["# Local VLM dtype probe (#258)", "",
           "## Machine", ""]
    out += [f"- **{k}**: {v}" for k, v in machine.items()]
    out += ["", "## Did it load, and what did it cost", "",
            "| dtype | loaded | device | load s | detect s | peak RSS GB | error |",
            "|---|---|---|---|---|---|---|"]
    for r in runs:
        out.append(
            f"| {r['dtype']} | {'yes' if r['loaded'] else '**no**'} | "
            f"{r.get('device', '-')} | {r['load_seconds'] or '-'} | "
            f"{r['detect_seconds'] or '-'} | {r['peak_rss_gb'] or '-'} | "
            f"{(r['error'] or '')[:80]} |")
    out += ["", "## ROI agreement against float32 — the criterion", ""]
    if not comparisons:
        out.append("_No baseline to compare against: float32 did not load._")
    for dtype, cmp in comparisons.items():
        out += [f"### {dtype} vs float32", "",
                f"- figures compared: **{cmp['figures']}**",
                f"- same ROI count: **{cmp['same_roi_count']}/{cmp['figures']}**",
                f"- all boxes within IoU 0.95: **{cmp['tight_boxes']}/{cmp['figures']}**",
                "",
                "| figure | float32 ROIs | this ROIs | matched | min IoU | mean IoU | error |",
                "|---|---|---|---|---|---|---|"]
        for row in cmp["rows"]:
            out.append(
                f"| {row['figure']} | {row['baseline_rois']} | "
                f"{row['other_rois']} | {row['matched_labels']} | "
                f"{row['min_iou'] if row['min_iou'] is not None else '-'} | "
                f"{row['mean_iou'] if row['mean_iou'] is not None else '-'} | "
                f"{(row['error'] or '')[:60]} |")
        out.append("")
    if against_reference:
        out += ["## ROI agreement against the recorded reference", "",
                f"Reference ROIs were produced by `{reference_source or 'unknown'}`.",
                "This comparison does not depend on float32 loading here, and",
                "matching it is the stronger claim: it says half precision on",
                "this machine agrees with the answer production shipped.", "",
                "| dtype | figures | same ROI count | boxes within IoU 0.95 |",
                "|---|---|---|---|"]
        for dtype, cmp in against_reference.items():
            out.append(f"| {dtype} | {cmp['figures']} | "
                       f"{cmp['same_roi_count']}/{cmp['figures']} | "
                       f"{cmp['tight_boxes']}/{cmp['figures']} |")
        out.append("")
    out += ["## Verdict", "",
            "Fill this in. The default moves to half precision on MPS only if a",
            "half-precision dtype matches float32's ROI counts and keeps boxes",
            "tight (IoU >= 0.95). If it does not, say which dtype and which",
            "figures diverged — that is the comment #258 asks for, and it should",
            "name the failing op if the logs reveal one, so nobody re-litigates",
            "this.", "",
            "If float32 could not load at all on this machine, that is itself a",
            "finding: it means the current default makes the local VLM backend",
            "unusable here, and the comparison has to be run against whichever",
            "dtype did load.", ""]
    return "\n".join(out)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output_dir", type=Path,
                        help="A built corpuscle output directory (has documents/)")
    parser.add_argument("--dtypes", default=",".join(DEFAULT_DTYPES),
                        help=f"Comma-separated (default: {','.join(DEFAULT_DTYPES)})")
    parser.add_argument("--figures", type=int, default=8,
                        help="How many figures to probe (default 8)")
    parser.add_argument("--device", default=None,
                        help="Force a torch device; default autodetects")
    parser.add_argument("--out", type=Path, default=None,
                        help="Write the markdown report here as well as stdout")
    args = parser.parse_args(argv)

    dtypes = [d.strip() for d in args.dtypes.split(",") if d.strip()]
    figures = select_figures(args.output_dir.resolve(), args.figures)
    print(f"probing {len(dtypes)} dtype(s) over {len(figures)} figure(s)",
          file=sys.stderr)
    for entry in figures:
        print(f"  {entry['hash']}/{entry['figure_id']} "
              f"expected={entry['expected'] or '-'}", file=sys.stderr)

    reference = {
        "dtype": "recorded (reference)",
        "loaded": True,
        "figures": {
            f"{e['hash']}/{e['figure_id']}": {"rois": e["reference_rois"],
                                              "error": None}
            for e in figures if e["reference_rois"]
        },
        "source": next((e["reference_source"] for e in figures
                        if e.get("reference_source")), None),
    }

    runs = []
    for dtype in dtypes:
        print(f"\n=== {dtype} ===", file=sys.stderr)
        run = probe_dtype(dtype, figures, args.device)
        if run["error"]:
            print(f"  did not load: {run['error']}", file=sys.stderr)
        else:
            print(f"  loaded in {run['load_seconds']}s, detected in "
                  f"{run['detect_seconds']}s, peak {run['peak_rss_gb']} GB",
                  file=sys.stderr)
        runs.append(run)

    baseline = next((r for r in runs
                     if r["dtype"] == "float32" and r["loaded"]), None)
    if baseline is None:
        baseline = next((r for r in runs if r["loaded"]), None)
    comparisons = {}
    if baseline is not None:
        for run in runs:
            if run["loaded"] and run["dtype"] != baseline["dtype"]:
                comparisons[run["dtype"]] = compare(baseline, run)
    # And every dtype against the recorded reference, which is the
    # comparison that does not depend on float32 having loaded at all.
    against_reference = {}
    if reference["figures"]:
        for run in runs:
            if run["loaded"]:
                against_reference[run["dtype"]] = compare(reference, run)

    report = render(machine_report(), runs, comparisons,
                    against_reference, reference.get("source"))
    print(report)
    if args.out:
        args.out.write_text(report, encoding="utf-8")
        print(f"\nwrote {args.out}", file=sys.stderr)
    # Raw data alongside the report, for anything the tables flatten.
    raw = (args.out.with_suffix(".json") if args.out
           else Path("vlm_dtype_probe.json"))
    raw.write_text(json.dumps(
        {"machine": machine_report(), "runs": runs,
         "comparisons": comparisons,
         "against_reference": against_reference,
         "reference_source": reference.get("source")},
        indent=2, default=str), encoding="utf-8")
    print(f"wrote {raw}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
