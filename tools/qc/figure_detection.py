#!/usr/bin/env python3
"""Score figure *detection* against the gold transcription set (#194).

Two questions, and they are not the same one:

  * **Recall** — is every figure on the page found?
  * **Precision** — is publisher furniture being called a figure? A journal
    logo, an ORCID icon, a rule.

Caption association is a separate scorer (#195), and OCR of text *inside* a
figure is explicitly not a target — ``fidelity.py`` reports that as
``figure_coverage`` beside prose coverage, and it stays reported-only.

WHY COUNTING, AND WHY PER PAGE
==============================
The gold set marks ``[FIGURE]`` and ``[PLATE]`` blocks but records no bounding
boxes, so an individual gold block cannot be matched to an individual
``figures.json`` entry. Per-page counts can be, and they are enough to answer
both questions above.

Counting raw entries corpus-wide is *not* enough. Logical plate children and
split panels make that total a category error, while document-level misses and
surpluses can cancel even among physical detections. **The per-page
disagreement is the measurement; the total hides it.**

So per page: ``matched = min(gold, found)``, ``missed = gold - matched``,
``surplus = found - matched``. Recall and precision are the corpus sums of
those, which is a deliberately conservative bipartite approximation — it
cannot tell a page that found the right two figures from one that found two
wrong ones, and no page-level count can.

THE FURNITURE QUESTION IS A FILTER QUESTION
===========================================
``figures.json`` already carries ``figure_type``. If ``graphical_element``
were a reliable furniture predicate, dropping those entries would raise
precision without costing recall — and that is a thing this script measures
rather than assumes, by scoring the same corpuscle under each candidate
filter and reporting what each one buys and costs.

Usage::

    tools/qc/figure_detection.py --gold <transcriptions/> --corpuscle <corpuscle/>
        [--out figures.json] [--quiet]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import sys
from collections import Counter, defaultdict

# `tools/` is not an importable package, so the sibling harness is loaded by
# path. Sharing it rather than reimplementing matters: the gold-markup parser
# is subtle (a page prints brackets, and notes discuss brackets), and a second
# copy would drift from the one the tests cover.
_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("qc_fidelity", _HERE / "fidelity.py")
fid = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fid)

# The default MCP surface is part of the measurement. Import its semantic
# type set from the producer rather than maintaining a lookalike filter here;
# otherwise a server policy change can silently leave the reported "served"
# precision describing a surface that no longer exists.
sys.path.insert(0, str(_HERE.parent.parent))
from pipeline.figures import EVIDENCE_FIGURE_TYPES  # noqa: E402

# Only these are figures. `[TABLE]` is body content — `fidelity.py` scores it
# as prose for the same reason — and docling reports tables separately from
# pictures, so counting them here would compare two different populations.
GOLD_FIGURE_TAGS = ("FIGURE", "PLATE")


def gold_blocks_by_page(gold_dir):
    """``{page: Counter({"FIGURE": n, "PLATE": m})}`` for one document.

    Top-level blocks only. A ``[FIGURE]`` quoted inside a ``[NOTE:]`` is
    commentary about a figure, not a second figure, and the shared span
    scanner already declines to yield it.
    """
    out = defaultdict(Counter)
    for gf in sorted(gold_dir.glob("page_*.txt")):
        page = int(gf.stem.split("_")[1])
        text = gf.read_text(encoding="utf-8")
        depth = 0
        for a, b in fid.spans(text):
            inner = text[a + 1:b - 1]
            if inner in fid.STRUCTURAL_OPEN:
                if depth == 0 and inner in GOLD_FIGURE_TAGS:
                    out[page][inner] += 1
                depth += 1
            elif inner in fid.STRUCTURAL_CLOSE:
                depth = max(0, depth - 1)
        out.setdefault(page, Counter())
    return out


# --- candidate furniture filters ---------------------------------------------
#
# Each answers "which figures.json entries should count as a figure at all?".
# Scored side by side so the choice is made on evidence rather than on which
# field name sounds most like furniture.

def _captioned(f):
    return bool((f.get("caption_text") or "").strip())


FILTERS = {
    "all entries": lambda f: True,
    "drop graphical_element": lambda f: f.get("figure_type") != "graphical_element",
    "drop uncaptioned graphical_element":
        lambda f: f.get("figure_type") != "graphical_element" or _captioned(f),
    "served MCP types": lambda f: f.get("figure_type") in EVIDENCE_FIGURE_TYPES,
    "captioned only": _captioned,
}


def _is_image_sharing_record(figure):
    """Whether this is a logical child of one physical image host.

    Caption-derived plate children predate compound/vision materialization and
    store the host's Docling index in ``shares_image_with``. Newer Pass 3c and
    vision-discovery children store a figure ID in ``image_shared_with``.
    They encode the same fact for this scorer: neither is another physical
    detection.
    """
    return (
        figure.get("shares_image_with") is not None
        or figure.get("image_shared_with") is not None
    )


def collapse_panels(entries):
    """Return physical figures after collapsing derivative records.

    A figure with panels produces one image *per panel* — `fig_99_a.png` and
    `fig_99_b.png` — while the gold records one `[FIGURE]` block carrying both,
    with the panels enumerated inside its caption. Comparing raw entry counts
    against gold block counts therefore penalises the pipeline for splitting
    panels correctly (#211), and it flips the sign of the error: `Totton1965a`
    reads as over-counting by 3 entries when it is in fact under-counting by 4
    figures.

    Image-sharing plate/compound children are logical retrieval records, not
    additional physical detections. They carry ``shares_image_with`` or
    ``image_shared_with`` and are excluded here: their one host remains
    countable. This distinction became load-bearing when #195 materialized
    hundreds of caption-bound children from historical plate legends.

    Panel siblings are identified by the same typed
    ``(page, namespace, figure_number)`` key used by ``dedupe_figures``.
    ``plate:10`` and ``figure:10`` therefore remain separate physical objects,
    while ordinary figures and their ``subpanel`` entries share the ``figure``
    namespace. Non-evidence types are never merged merely because a proximity
    heuristic happened to attach the same number.
    """
    by_key, loose = defaultdict(list), []
    for f in entries:
        if _is_image_sharing_record(f):
            continue
        num = f.get("figure_number")
        figure_type = f.get("figure_type")
        namespace = (
            "plate" if figure_type == "plate"
            else "figure" if figure_type in {"figure", "subpanel"}
            else None
        )
        if f.get("page") and num is not None and str(num).strip() \
                and namespace is not None:
            by_key[(int(f["page"]), namespace, str(num))].append(f)
        else:
            # No typed number to group on. Counted as its own physical record,
            # which is conservative: it cannot silently merge two real figures
            # or hide a numbered item from the furniture review population.
            loose.append(f)
    figures = [v[0] for v in by_key.values()] + loose
    panel_groups = {k: len(v) for k, v in by_key.items() if len(v) > 1}
    return figures, panel_groups


def score(gold_pages, entries, keep):
    """Page-level recall/precision counts for one document under one filter.

    Counts *figures*, not entries — panel siblings are collapsed first, so a
    correctly split multi-panel figure counts once, as the gold counts it.
    """
    kept_entries = [f for f in entries if keep(f)]
    collapsed, _ = collapse_panels(kept_entries)
    found = Counter()
    for f in collapsed:
        if f.get("page"):
            found[int(f["page"])] += 1
    tally = Counter()
    for page in set(gold_pages) | set(found):
        g = sum(gold_pages.get(page, Counter()).values())
        j = found.get(page, 0)
        m = min(g, j)
        tally["gold"] += g
        tally["found"] += j
        tally["matched"] += m
        tally["missed"] += g - m
        tally["surplus"] += j - m
        tally["pages_exact"] += 1 if g == j else 0
        tally["pages"] += 1
    return tally


def _rates(t):
    recall = t["matched"] / t["gold"] if t["gold"] else None
    precision = t["matched"] / t["found"] if t["found"] else None
    f1 = (2 * recall * precision / (recall + precision)
          if recall and precision else None)
    return recall, precision, f1


def build_report(gold_root, corpuscle_root):
    bound, unmatched = fid.bind_documents(gold_root, corpuscle_root)
    docs, per_filter = {}, {name: Counter() for name in FILTERS}
    segments = {axis: defaultdict(lambda: Counter()) for axis in ("era", "file_type")}
    kind_tally = Counter()
    surplus_profile = Counter()

    for stem, sha, gold_dir, corpus_dir in bound:
        gold_pages = gold_blocks_by_page(gold_dir)
        figs = (fid._read_json(corpus_dir / "figures.json") or {}).get("figures") or []
        meta = fid._read_json(corpus_dir / "metadata.json") or {}
        scan = fid._read_json(corpus_dir / "scan_detection.json") or {}
        # #188 — figure `page` is a position in the subset when a selection
        # is active; the gold was transcribed over the whole file. Restate
        # the gold in subset coordinates and drop the pages the operator
        # excluded, so a deliberately-dropped plate atlas title page is not
        # counted as a figure corpus failed to find.
        gold_pages = fid.rebase_gold_pages(gold_pages, scan)

        for page, kinds in gold_pages.items():
            kind_tally.update(kinds)

        base = score(gold_pages, figs, FILTERS["all entries"])
        _, panel_groups = collapse_panels(figs)
        for name, keep in FILTERS.items():
            t = score(gold_pages, figs, keep)
            per_filter[name].update(t)
            if name == "all entries":
                era = fid.era_bucket(meta.get("year"))
                segments["era"][era].update(t)
                segments["file_type"][str(scan.get("file_type"))].update(t)

        # What sits on pages carrying more physical figures than the gold has
        # blocks — the population furniture would be drawn from. Raw logical
        # children and split-panel images must not distort this diagnostic any
        # more than they distort the headline score.
        by_page = defaultdict(list)
        physical_figures, _ = collapse_panels(figs)
        for f in physical_figures:
            if f.get("page"):
                by_page[int(f["page"])].append(f)
        for page, on_page in by_page.items():
            g = sum(gold_pages.get(page, Counter()).values())
            if len(on_page) > g:
                for f in on_page:
                    surplus_profile[(f.get("figure_type"), _captioned(f))] += 1

        r, p, f1 = _rates(base)
        docs[stem] = {
            "year": meta.get("year"), "era": fid.era_bucket(meta.get("year")),
            "file_type": scan.get("file_type"),
            "gold_figures": base["gold"], "found": base["found"],
            "missed": base["missed"], "surplus": base["surplus"],
            "recall": round(r, 4) if r is not None else None,
            "precision": round(p, 4) if p is not None else None,
            "pages_exact": base["pages_exact"], "pages": base["pages"],
            # Reported beside detection rather than folded into it: a split
            # panel is a *correct* extra image, and #195 scores whether the
            # split matches the panels the gold caption enumerates.
            "entries": len(figs),
            "shared_image_records": sum(
                _is_image_sharing_record(f) for f in figs
            ),
            "physical_figures": len(physical_figures),
            "panel_groups": len(panel_groups),
            "panel_images": sum(panel_groups.values()),
        }

    def pack(t):
        r, p, f1 = _rates(t)
        return {
            **{k: t[k] for k in ("gold", "found", "matched", "missed", "surplus",
                                 "pages", "pages_exact")},
            "recall": round(r, 4) if r is not None else None,
            "precision": round(p, 4) if p is not None else None,
            "f1": round(f1, 4) if f1 is not None else None,
        }

    return {
        "schema_version": 2,
        "question": "are all figures found (recall), and is furniture being "
                    "called a figure (precision)",
        "method": "per-page counts of gold [FIGURE]/[PLATE] blocks against "
                  "physical figures.json detections; image-sharing logical "
                  "plate children are excluded, typed panel siblings collapse "
                  "to one figure, and matched = min(gold, found) per page",
        "caveat": "A page-level count cannot tell a page that found the right "
                  "two figures from one that found two wrong ones. The gold set "
                  "records no bounding boxes, so nothing finer is available.",
        "gold_root": str(gold_root), "corpuscle_root": str(corpuscle_root),
        "documents_bound": len(bound), "documents_unmatched": unmatched,
        "gold_block_kinds": dict(kind_tally),
        "filters": {name: pack(t) for name, t in per_filter.items()},
        "segments": {axis: {k: pack(v) for k, v in sorted(buckets.items())}
                     for axis, buckets in segments.items()},
        "surplus_page_profile": {f"{k[0]}/{'captioned' if k[1] else 'uncaptioned'}": v
                                 for k, v in surplus_profile.most_common()},
        "documents": docs,
    }


def _r(value, width=8):
    """Rate or a dash, right-aligned — rates are legitimately None when a
    document has no gold figures, and a printed 0.000 would read as failure."""
    return f"{value:>{width}.3f}" if isinstance(value, float) else f"{'-':>{width}}"


def print_summary(report, stream=None):
    # Resolved here, not as a default argument: a default binds
    # `sys.stdout` at definition time, so the function would ignore
    # any later redirection — including a test capturing it.
    w = (stream or sys.stdout).write
    w(f"\n{report['documents_bound']} documents bound on sha256; gold blocks "
      f"{report['gold_block_kinds']}\n")
    pg = sum(d["panel_groups"] for d in report["documents"].values())
    pi = sum(d["panel_images"] for d in report["documents"].values())
    en = sum(d["entries"] for d in report["documents"].values())
    sr = sum(d["shared_image_records"] for d in report["documents"].values())
    pf = sum(d["physical_figures"] for d in report["documents"].values())
    w(f"{en} figures.json entries contain {sr} image-sharing logical records; "
      f"the remainder collapse to {pf} physical figures ({pg} multi-panel "
      f"groups holding {pi} images) — counted as the gold counts them\n")

    w("\n-- does a furniture filter help? " + "-" * 38 + "\n")
    w(f"{'filter':38}{'found':>7}{'recall':>8}{'precis':>8}{'F1':>8}"
      f"{'missed':>8}{'surplus':>8}\n")
    for name, t in report["filters"].items():
        w(f"{name[:38]:38}{t['found']:>7}{_r(t['recall'])}{_r(t['precision'])}"
          f"{_r(t['f1'])}{t['missed']:>8}{t['surplus']:>8}\n")

    w("\n-- what sits on pages with more physical figures than the gold has "
      + "-" * 3 + "\n")
    for k, v in report["surplus_page_profile"].items():
        w(f"  {k:44} {v}\n")

    for axis, buckets in report["segments"].items():
        w(f"\n-- by {axis} " + "-" * (52 - len(axis)) + "\n")
        w(f"{'':20}{'gold':>7}{'found':>7}{'recall':>8}{'precis':>8}\n")
        for name, t in buckets.items():
            w(f"{name[:20]:20}{t['gold']:>7}{t['found']:>7}"
              f"{_r(t['recall'])}{_r(t['precision'])}\n")

    w("\n-- documents furthest from the gold count " + "-" * 28 + "\n")
    w(f"{'':26}{'gold':>6}{'found':>7}{'missed':>8}{'surplus':>9}{'recall':>8}{'precis':>8}\n")
    ranked = sorted(report["documents"].items(),
                    key=lambda kv: -(kv[1]["missed"] + kv[1]["surplus"]))
    for stem, d in ranked[:12]:
        w(f"{stem[:26]:26}{d['gold_figures']:>6}{d['found']:>7}{d['missed']:>8}"
          f"{d['surplus']:>9}{_r(d['recall'])}{_r(d['precision'])}\n")
    w("\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gold", required=True)
    ap.add_argument("--corpuscle", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    gold_root = pathlib.Path(a.gold).expanduser().resolve()
    corpuscle_root = pathlib.Path(a.corpuscle).expanduser().resolve()
    if not (corpuscle_root / "documents").is_dir() and \
            (corpuscle_root / "output" / "documents").is_dir():
        corpuscle_root = corpuscle_root / "output"

    report = build_report(gold_root, corpuscle_root)
    if a.out:
        out = pathlib.Path(a.out).expanduser()
        out.write_text(json.dumps(report, indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"wrote {out}")
    if not a.quiet:
        print_summary(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
