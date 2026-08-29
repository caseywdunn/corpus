#!/usr/bin/env python3
"""Score caption association against the gold transcription set (#195).

Detection asks whether the figure *objects* are right (`figure_detection.py`).
This asks the next question: is the caption attached to the figure it belongs
to? `dev_docs/OVERVIEW.md` calls caption association "the highest-value
annotation per figure and the hardest in historical layouts"; this turns that
from an assertion into a number.

WHY THIS BINDS ON THE FIGURE NUMBER
===================================
The obvious measure — match each gold caption to a `caption_text` by string
similarity — was built, run, and thrown away. It reports 44% bound and is
mostly artifact:

  * **Bilingual captions.** One paper prints every caption twice, once in
    Chinese and once in English. The extraction captured the Chinese; a
    matcher handed the English scores 0 of 10 while every figure is in fact
    bound correctly.
  * **A label is not a caption.** Plate pages carry `FIG. 1` and nothing else
    — two tokens, which no similarity threshold can bind — while the prose
    caption sits on another page entirely.
  * **A document can be its own translation**, so `Fig. 1` legitimately
    appears twice with different text.

A figure *number* is language-independent and is what a client actually asks
with. So the primary measure is whether the number the pipeline reports for a
page matches a number printed on a figure on that page.

WHAT THE DENOMINATOR HAS TO BE, AND WHY IT IS EASY TO GET WRONG
===============================================================
Gold pages are full of figure numbers that are *references* — "see Fig. 18",
"figured by Bigelow (op. cit., fig. 34)". Counting those as figures the
pipeline failed to find gives 939 numbers across the set and a recall of
0.296, which is meaningless: most of them are not on that page at all.

Restricted to numbers printed *inside* a `[FIGURE]`/`[PLATE]` block — the
figures actually on the leaf — the denominator is 482 and the measure means
what it says.

Usage::

    tools/qc/caption_binding.py --gold <transcriptions/> --corpuscle <corpuscle/>
        [--out captions.json] [--quiet]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import sys
from collections import Counter, defaultdict
from difflib import SequenceMatcher

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("qc_fidelity", _HERE / "fidelity.py")
fid = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fid)

sys.path.insert(0, str(_HERE.parent.parent))
from pipeline.figures import parse_figure_number  # noqa: E402

GOLD_FIGURE_TAGS = ("FIGURE", "PLATE")

# A block with at least this many tokens after its label is carrying a prose
# caption; below it, the label is all that is printed. The distinction matters
# because only the first is a caption-binding test — scoring a bare `FIG. 1`
# for caption text would report a failure where the page has nothing to bind.
_PROSE_CAPTION_MIN_TOKENS = 8


def gold_blocks(gold_dir, scan=None):
    """Yield ``(page, kind, body)`` for each top-level figure/plate block.

    ``page`` is restated in subset coordinates when a `keeppages` selection
    is active (#188), and blocks on excluded pages are not yielded at all —
    the operator said those pages are not the paper, so a caption printed on
    one is not a caption corpus failed to bind.
    """
    mapping = fid.keeppages_map(scan)
    for gf in sorted(gold_dir.glob("page_*.txt")):
        page = int(gf.stem.split("_")[1])
        if mapping:
            if page not in mapping:
                continue
            page = mapping[page]
        text = gf.read_text(encoding="utf-8")
        depth, kind, start = 0, None, None
        for a, b in fid.spans(text):
            inner = text[a + 1:b - 1]
            if inner in fid.STRUCTURAL_OPEN:
                if depth == 0:
                    kind, start = inner, b
                depth += 1
            elif inner in fid.STRUCTURAL_CLOSE:
                depth -= 1
                if depth == 0 and start is not None:
                    if kind in GOLD_FIGURE_TAGS:
                        yield page, kind, text[start:a]
                    start = None


def classify_block(body):
    """What kind of evidence this gold block offers.

    Only ``prose_caption`` blocks can test caption *text*; the rest are
    counted so no rate is computed over blocks it does not apply to.
    """
    stripped = fid.strip_markup(body).text
    numbers, first_line = set(), None
    for line in stripped.splitlines():
        line = line.strip()
        num = parse_figure_number(line)
        if num:
            numbers.add(num)
            if first_line is None:
                first_line = line
    tokens = fid.tokens(stripped)
    if not tokens:
        return "nothing_printed", numbers, ""
    if not numbers:
        return "lettering_only", numbers, ""
    label_tokens = fid.tokens(first_line or "")
    if len(label_tokens) >= _PROSE_CAPTION_MIN_TOKENS:
        return "prose_caption", numbers, first_line
    return "bare_label", numbers, first_line


def score_document(gold_dir, corpus_dir):
    figs = (fid._read_json(corpus_dir / "figures.json") or {}).get("figures") or []
    scan = fid._read_json(corpus_dir / "scan_detection.json") or {}
    by_page = defaultdict(list)
    for f in figs:
        if f.get("page"):
            by_page[int(f["page"])].append(f)

    kinds = Counter()
    gold_nums, found_nums, matched = defaultdict(set), defaultdict(set), 0
    captions = []                      # (gold_caption, extracted_caption)

    for page, _kind, body in gold_blocks(gold_dir, scan):
        kind, numbers, first_line = classify_block(body)
        kinds[kind] += 1
        gold_nums[page] |= numbers
        if kind == "prose_caption" and numbers:
            for f in by_page.get(page, []):
                if str(f.get("figure_number")) in numbers:
                    captions.append((first_line, f.get("caption_text") or ""))
                    break

    for page, on_page in by_page.items():
        found_nums[page] |= {str(f["figure_number"]) for f in on_page
                             if f.get("figure_number")}

    g = sum(len(v) for v in gold_nums.values())
    j = sum(len(v) for v in found_nums.values())
    matched = sum(len(gold_nums[p] & found_nums.get(p, set())) for p in gold_nums)

    sims = []
    for gold_cap, got_cap in captions:
        if got_cap:
            sims.append(SequenceMatcher(None, fid.tokens(gold_cap),
                                        fid.tokens(got_cap),
                                        autojunk=False).ratio())
    return {
        "block_kinds": dict(kinds),
        "gold_numbers": g, "found_numbers": j, "matched_numbers": matched,
        "number_recall": round(matched / g, 4) if g else None,
        "number_precision": round(matched / j, 4) if j else None,
        "captions_compared": len(sims),
        "median_caption_similarity": fid._median(sims),
    }


def build_report(gold_root, corpuscle_root):
    bound, unmatched = fid.bind_documents(gold_root, corpuscle_root)
    docs, totals, kinds = {}, Counter(), Counter()
    segments = defaultdict(lambda: Counter())
    for stem, sha, gold_dir, corpus_dir in bound:
        d = score_document(gold_dir, corpus_dir)
        meta = fid._read_json(corpus_dir / "metadata.json") or {}
        d["era"] = fid.era_bucket(meta.get("year"))
        docs[stem] = d
        kinds.update(d["block_kinds"])
        for k in ("gold_numbers", "found_numbers", "matched_numbers",
                  "captions_compared"):
            totals[k] += d[k]
            segments[d["era"]][k] += d[k]
    return {
        "schema_version": 1,
        "question": "is the caption bound to the figure it belongs to",
        "method": "figure numbers printed INSIDE a gold [FIGURE]/[PLATE] block, "
                  "compared per page against figures.json; caption text scored "
                  "only within number-matched pairs",
        "caveat": "Naive caption-text similarity reports 44% and is mostly "
                  "artifact — bilingual captions, bare labels, and documents "
                  "that are their own translation. Numbers are the "
                  "language-independent signal.",
        "documents_bound": len(bound), "documents_unmatched": unmatched,
        "gold_block_kinds": dict(kinds),
        "totals": {
            **dict(totals),
            "number_recall": round(totals["matched_numbers"] / totals["gold_numbers"], 4)
            if totals["gold_numbers"] else None,
            "number_precision": round(totals["matched_numbers"] / totals["found_numbers"], 4)
            if totals["found_numbers"] else None,
        },
        "segments": {k: dict(v) for k, v in sorted(segments.items())},
        "documents": docs,
    }


def print_summary(report, stream=None):
    w = (stream or sys.stdout).write
    t = report["totals"]
    w(f"\n{report['documents_bound']} documents bound on sha256\n")
    w(f"gold figure blocks by what they print: {report['gold_block_kinds']}\n")
    w(f"\nfigure numbers printed on a figure : {t['gold_numbers']}\n")
    w(f"numbers the pipeline reports       : {t['found_numbers']}\n")
    w(f"  binding recall                   : {t['number_recall']}\n")
    w(f"  binding precision                : {t['number_precision']}\n")

    w("\n-- worst documents by binding recall " + "-" * 30 + "\n")
    w(f"{'':26}{'gold':>6}{'found':>7}{'recall':>8}{'precis':>8}{'capsim':>8}\n")
    ranked = sorted(report["documents"].items(),
                    key=lambda kv: (kv[1]["number_recall"] is None,
                                    kv[1]["number_recall"] or 0))
    for stem, d in ranked[:12]:
        def f(v):
            return f"{v:>8.3f}" if isinstance(v, float) else f"{'-':>8}"
        w(f"{stem[:26]:26}{d['gold_numbers']:>6}{d['found_numbers']:>7}"
          f"{f(d['number_recall'])}{f(d['number_precision'])}"
          f"{f(d['median_caption_similarity'])}\n")
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
        pathlib.Path(a.out).expanduser().write_text(
            json.dumps(report, indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"wrote {a.out}")
    if not a.quiet:
        print_summary(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
