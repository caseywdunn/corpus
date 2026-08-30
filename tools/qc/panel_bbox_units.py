#!/usr/bin/env python3
"""Are Pass 3b's panel bboxes normalized 0–1 floats, or pixels? (#253)

The prompt demands "each coordinate is a float in 0.0 .. 1.0"
(pipeline/vision.py), and `_norm_to_px` multiplies by the image dimensions
on that assumption. #253 observed response heads carrying `[50, 30, 250,
500]` and `[17, 808, 150, 1365]` — absolute pixels. When that happens
`_norm_to_px` clamps `x1` with `min(1.0, ...) * w` while `x0` becomes
`50 * w`, so `x1 <= x0` and the guard drops the panel silently, with no
warning on that path at all.

That matters because it lands in the *same* counters as the truncation bug
(#253's main defect): both show up as missing ROIs and `no_labels_found`.
Fixing truncation alone and re-measuring would look like the fix
underperforming, when the residue is a second, unrelated cause.

This script does not fix anything and does not need a GPU. It reads the
`figures.json` a build already produced and reports the `pass3_status`
breakdown bucketed by caption-derived panel count.

**On the document count:** `figure_passes.py` writes `pass3_status` only on
figures whose caption yields more than one panel, so a document with no
multi-panel figure carries no status at all. The "documents with >=1
multi-panel figure" line is therefore the *eligible* population, not a
partial sample — Pass 3b ran over everything and had nothing to do on the
rest. Reported that way because the earlier wording ("reached Pass 3b")
read as though 84% of a finished build had been skipped.

Read-only. Safe against a build that is still running — it skips documents
that have not reached Pass 3b yet.

Usage:
    tools/qc/panel_bbox_units.py --corpuscle ~/corpora/<name>_<stamp>/output
    tools/qc/panel_bbox_units.py --corpuscle <dir> --json out.json
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path


def _read_json(path: Path):
    try:
        if path.exists() and path.stat().st_size > 0:
            with path.open(encoding="utf-8") as f:
                return json.load(f)
    except Exception:
        pass
    return None


def _classify(bbox, w: int, h: int) -> str:
    """What units does this bbox look like?

    A 4-list of floats all within 0..1 is normalized. Anything with a
    component above 1 cannot be normalized; if it is also within the image
    dimensions it is very likely pixels. Above both is neither.
    """
    if not (isinstance(bbox, list) and len(bbox) == 4):
        return "malformed"
    try:
        vals = [float(v) for v in bbox]
    except (TypeError, ValueError):
        return "malformed"
    if all(0.0 <= v <= 1.0 for v in vals):
        return "normalized"
    if w and h and vals[0] <= w and vals[2] <= w and vals[1] <= h and vals[3] <= h:
        return "pixels"
    return "out_of_range"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpuscle", required=True, type=Path,
                    help="build output dir (the one containing documents/)")
    ap.add_argument("--json", type=Path, default=None,
                    help="also write the full result as JSON")
    ap.add_argument("--examples", type=int, default=5,
                    help="how many example bboxes to print per class")
    args = ap.parse_args()

    docs = args.corpuscle / "documents"
    if not docs.is_dir():
        print(f"error: no documents/ under {args.corpuscle}", file=sys.stderr)
        return 2

    status = Counter()
    units = Counter()
    by_panels = {}          # panel_count bucket -> Counter of statuses
    examples = {}
    n_docs = n_pass3b = 0
    surviving_rois = 0

    for hd in sorted(docs.iterdir()):
        if not hd.is_dir():
            continue
        n_docs += 1
        figs = _read_json(hd / "figures.json")
        if not isinstance(figs, dict):
            continue
        seen_pass3b = False
        for fig in figs.get("figures") or []:
            if not isinstance(fig, dict):
                continue
            st = fig.get("pass3_status")
            if st is None:
                continue          # not reached by Pass 3b yet
            seen_pass3b = True
            status[st] += 1

            n_panels = int(fig.get("panel_count_from_caption") or 0)
            bucket = ("2-3" if n_panels <= 3 else "4-5" if n_panels <= 5
                      else "6-7" if n_panels <= 7 else "8-9" if n_panels <= 9
                      else "10+")
            by_panels.setdefault(bucket, Counter())[st] += 1

            rois = fig.get("rois") or []
            surviving_rois += len(rois)
            size = fig.get("image_size_px") or [0, 0]
            w, h = (int(size[0]), int(size[1])) if len(size) == 2 else (0, 0)

            # The surviving ROIs carry pixel bboxes already converted by
            # _norm_to_px, so they cannot tell us about the input units.
            # What can: any raw *_norm field the backend recorded.
            for roi in rois:
                if not isinstance(roi, dict):
                    continue
                for key in ("panel_bbox_norm", "label_bbox_norm", "bbox_norm"):
                    if key in roi:
                        cls = _classify(roi[key], w, h)
                        units[cls] += 1
                        examples.setdefault(cls, [])
                        if len(examples[cls]) < args.examples:
                            examples[cls].append(
                                {"doc": hd.name, "figure": fig.get("figure_id"),
                                 key: roi[key], "image_size_px": [w, h]})
        if seen_pass3b:
            n_pass3b += 1

    print(f"corpuscle : {args.corpuscle}")
    print(f"documents : {n_docs} total, {n_pass3b} with >=1 multi-panel figure")
    print("            (Pass 3b writes pass3_status only on figures whose "
          "caption\n"
          "             yields >1 panel, so the rest are not a missing "
          "sample --\n"
          "             Pass 3b ran and had nothing to do on them.)")
    print(f"rois kept : {surviving_rois}")
    print()
    print("pass3_status:")
    for k, v in status.most_common():
        print(f"  {v:6d}  {k}")

    print("\nby caption-derived panel count:")
    hdr = ["2-3", "4-5", "6-7", "8-9", "10+"]
    keys = [k for k, _ in status.most_common()]
    print(f"  {'panels':<8}" + "".join(f"{k[:18]:>20}" for k in keys))
    for b in hdr:
        c = by_panels.get(b)
        if not c:
            continue
        print(f"  {b:<8}" + "".join(f"{c.get(k, 0):>20}" for k in keys))

    print("\nraw *_norm bbox units, where the backend recorded them:")
    if not units:
        print("  (none recorded — see the note below)")
    for k, v in units.most_common():
        print(f"  {v:6d}  {k}")
    for cls, ex in examples.items():
        print(f"\n  examples [{cls}]:")
        for e in ex:
            print(f"    {e}")

    if not units:
        print(
            "\nNOTE: figures.json stores only the *converted* pixel bboxes, so\n"
            "raw normalized values are not recoverable from artifacts alone.\n"
            "The units question has to be answered from the Pass 3b logs\n"
            "instead — see the grep in the issue thread. The pass3_status\n"
            "table above is still the useful half: it is the coverage\n"
            "gradient #253 measured, recomputed on this build."
        )

    if args.json:
        args.json.write_text(json.dumps({
            "corpuscle": str(args.corpuscle),
            "documents_total": n_docs,
            "documents_with_pass3b": n_pass3b,
            "rois_kept": surviving_rois,
            "pass3_status": dict(status),
            "by_panel_bucket": {k: dict(v) for k, v in by_panels.items()},
            "bbox_units": dict(units),
            "examples": examples,
        }, indent=2))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
