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
pipeline failed to find inflates the denominator with objects that are not on
that page. Restricted to numbers printed *inside* a `[FIGURE]`/`[PLATE]`
block — the figures actually on the leaf — the measure means what it says.

Usage::

    tools/qc/caption_binding.py --gold <transcriptions/> --corpuscle <corpuscle/>
        [--out captions.json] [--quiet]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import re
import sys
from collections import Counter, defaultdict
from difflib import SequenceMatcher

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("qc_fidelity", _HERE / "fidelity.py")
fid = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fid)

sys.path.insert(0, str(_HERE.parent.parent))
from pipeline.figures import (  # noqa: E402
    EVIDENCE_FIGURE_TYPES,
    caption_evidence_summary,
)

GOLD_FIGURE_TAGS = ("FIGURE", "PLATE")

# A block with at least this many tokens after its label is carrying a prose
# caption; below it, the label is all that is printed. The distinction matters
# because only the first is a caption-binding test — scoring a bare `FIG. 1`
# for caption text would report a failure where the page has nothing to bind.
_PROSE_CAPTION_MIN_TOKENS = 8


# Gold labels are parsed independently of the production extractor. Importing
# parse_figure_number/caption_figure_entries here made the yardstick move when
# those functions changed: the Totton source-citation fix changed the gold
# denominator from 496 to 486 even though no transcription changed. The gold
# is human transcription, so it needs only the printed, non-OCR-damage opener
# vocabulary plus lists/ranges. Keep this deliberately narrower and fixture-
# backed; production's fuzzy OCR spellings belong only on the measured side.
_GOLD_OPENER = (
    r"fig(?:ure|ur|s)?|figg|figuren|abb(?:ildung)?|pl(?:ate)?|plate|рис(?:унок)?"
    r"|taf(?:el)?|tab(?:ula)?|lám(?:ina)?|tav(?:ola)?|bild|image|illustration"
    r"|text[\-\s]?fig(?:ure)?"
)
_GOLD_NUMBER = r"(?:\d+(?:\.\d+)?[a-z]?|[IVXLCDMivxlcdm]+(?![A-Za-z]))"
_GOLD_CONNECTOR = r"(?:,\s*(?:&|und|and|et|u\.?)?|&|und|and|et|u\.?|[-\u2013\u2014])"
_GOLD_ENTRY_RE = re.compile(
    r"(?<![A-Za-z])(?:" + _GOLD_OPENER + r")\s*\.?\s*"
    r"(?P<first>" + _GOLD_NUMBER + r")"
    r"(?P<tail>(?:\s*" + _GOLD_CONNECTOR + r"\s*" + _GOLD_NUMBER + r")*)",
    re.IGNORECASE,
)
_GOLD_ENTRY_START_RE = re.compile(r"^[\s.\u00b7\u2022]*" + _GOLD_ENTRY_RE.pattern,
                                  re.IGNORECASE)
_GOLD_TAIL_RE = re.compile(
    r"(?P<connector>" + _GOLD_CONNECTOR + r")\s*"
    r"(?P<number>" + _GOLD_NUMBER + r")",
    re.IGNORECASE,
)
_GOLD_ROMAN_RE = re.compile(
    r"^m{0,4}(cm|cd|d?c{0,3})(xc|xl|l?x{0,3})(ix|iv|v?i{0,3})$",
    re.IGNORECASE,
)
_GOLD_ROMAN_VALUES = {"i": 1, "v": 5, "x": 10, "l": 50,
                      "c": 100, "d": 500, "m": 1000}


def _gold_number(token):
    """Canonical number from a hand-transcribed printed label."""
    value = (token or "").strip()
    if not value or not _GOLD_ROMAN_RE.match(value):
        return value.lower()
    total, previous = 0, 0
    for char in reversed(value.lower()):
        current = _GOLD_ROMAN_VALUES[char]
        total += -current if current < previous else current
        previous = current
    return str(total)


def gold_caption_entries(line):
    """Figure numbers asserted by one line of hand-transcribed gold.

    A line must open with a printed label. Later entries require a sentence or
    semicolon boundary; comma-delimited terms without another opener remain a
    supported list, while source citations such as ``..., fig. 36`` do not
    become additional figures.
    """
    text = line or ""
    if not _GOLD_ENTRY_START_RE.match(text):
        return []
    matches = []
    for match in _GOLD_ENTRY_RE.finditer(text):
        if not matches:
            matches.append(match)
            continue
        before = text[:match.start()]
        if not before[before.rfind("\n") + 1:].strip() \
                or before.rstrip()[-1:] in ".;":
            matches.append(match)

    numbers, seen = [], set()
    for match in matches:
        values = [_gold_number(match.group("first"))]
        previous = values[0]
        for term in _GOLD_TAIL_RE.finditer(match.group("tail") or ""):
            connector = term.group("connector").lower().rstrip(".")
            number = _gold_number(term.group("number"))
            if connector in {"-", "\u2013", "\u2014"} \
                    and previous.isdigit() and number.isdigit():
                start, end = int(previous), int(number)
                if 0 < end - start <= 100:
                    values.extend(str(n) for n in range(start + 1, end + 1))
                elif number not in values:
                    values.append(number)
            elif number not in values:
                values.append(number)
            previous = number
        for number in values:
            if number not in seen:
                seen.add(number)
                numbers.append(number)
    return numbers


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
        # A gold block can contain descriptive prose and references. Only a
        # line that opens as a caption is allowed to contribute; once it does,
        # split any additional entries/lists on that same printed line.
        entries = gold_caption_entries(line)
        numbers.update(entries)
        if first_line is None and entries:
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


def _layout_bucket(structural_kinds):
    """Document-level gold layout used for a precision-capable segment."""
    present = {k for k, count in structural_kinds.items() if count}
    if present == {"FIGURE"}:
        return "figure_blocks"
    if present == {"PLATE"}:
        return "plate_blocks"
    if present:
        return "mixed_figure_plate"
    return "no_figure_blocks"


def _number_key(value):
    value = str(value)
    return (0, int(value)) if value.isdigit() else (1, value)


def _pack_counts(counts):
    """Add rates to one aggregate without hiding undefined precision."""
    gold = counts["gold_numbers"]
    found = counts["found_numbers"]
    matched = counts["matched_numbers"]
    recall = matched / gold if gold else None
    precision = matched / found if found else None
    f1 = None
    if recall is not None and precision is not None:
        f1 = (
            2 * recall * precision / (recall + precision)
            if recall + precision else 0.0
        )
    return {
        **{k: counts[k] for k in (
            "gold_numbers", "found_numbers", "matched_numbers",
            "captions_compared",
        )},
        "number_recall": round(recall, 4) if recall is not None else None,
        "number_precision": round(precision, 4) if precision is not None else None,
        "f1": round(f1, 4) if f1 is not None else None,
    }


def score_document(gold_dir, corpus_dir, figure_types=None):
    figs = (fid._read_json(corpus_dir / "figures.json") or {}).get("figures") or []
    if figure_types is not None:
        figs = [f for f in figs if f.get("figure_type") in figure_types]
    scan = fid._read_json(corpus_dir / "scan_detection.json") or {}
    by_page = defaultdict(list)
    for f in figs:
        if f.get("page"):
            by_page[int(f["page"])].append(f)

    kinds, structural_kinds = Counter(), Counter()
    gold_nums, found_nums, matched = defaultdict(set), defaultdict(set), 0
    gold_blocks_by_page = defaultdict(list)
    captions = []                      # (gold_caption, extracted_caption)

    for page, structural_kind, body in gold_blocks(gold_dir, scan):
        kind, numbers, first_line = classify_block(body)
        kinds[kind] += 1
        structural_kinds[structural_kind] += 1
        gold_nums[page] |= numbers
        gold_blocks_by_page[page].append({
            "structural_kind": structural_kind,
            "caption_kind": kind,
            "numbers": sorted(numbers, key=_number_key),
        })
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
    page_diagnostics = {}
    for page in sorted(set(gold_nums) | set(found_nums)):
        gold_on_page = gold_nums.get(page, set())
        found_on_page = found_nums.get(page, set())
        figures = []
        for figure in by_page.get(page, []):
            summary = caption_evidence_summary(figure)
            figures.append({
                "figure_id": figure.get("figure_id"),
                "figure_type": figure.get("figure_type"),
                "figure_number": figure.get("figure_number"),
                "caption_preview": (
                    figure.get("caption_text") or figure.get("caption") or ""
                )[:240],
                "caption_source": figure.get("caption_source"),
                **summary,
            })
        page_diagnostics[str(page)] = {
            "gold_blocks": gold_blocks_by_page.get(page, []),
            "gold_numbers": sorted(gold_on_page, key=_number_key),
            "found_numbers": sorted(found_on_page, key=_number_key),
            "matched_numbers": sorted(
                gold_on_page & found_on_page, key=_number_key,
            ),
            "missing_numbers": sorted(
                gold_on_page - found_on_page, key=_number_key,
            ),
            "surplus_numbers": sorted(
                found_on_page - gold_on_page, key=_number_key,
            ),
            "reported_figures": figures,
        }

    evidence = Counter()
    for figure in figs:
        summary = caption_evidence_summary(figure)
        evidence[f"status/{summary['caption_status']}"] += 1
        evidence[f"confidence/{summary['caption_confidence']}"] += 1
        evidence[f"kind/{summary['caption_kind']}"] += 1

    return {
        "block_kinds": dict(kinds),
        "structural_kinds": dict(structural_kinds),
        "layout": _layout_bucket(structural_kinds),
        "gold_numbers": g, "found_numbers": j, "matched_numbers": matched,
        "number_recall": round(matched / g, 4) if g else None,
        "number_precision": round(matched / j, 4) if j else None,
        "captions_compared": len(sims),
        "median_caption_similarity": fid._median(sims),
        "caption_evidence_counts": dict(sorted(evidence.items())),
        "pages": page_diagnostics,
    }


def build_report(gold_root, corpuscle_root):
    bound, unmatched = fid.bind_documents(gold_root, corpuscle_root)
    docs, totals, kinds = {}, Counter(), Counter()
    segments = {
        axis: defaultdict(Counter)
        for axis in ("era", "file_type", "layout")
    }
    served_totals = Counter()
    served_segments = {
        axis: defaultdict(Counter)
        for axis in ("era", "file_type", "layout")
    }
    for stem, sha, gold_dir, corpus_dir in bound:
        d = score_document(gold_dir, corpus_dir)
        meta = fid._read_json(corpus_dir / "metadata.json") or {}
        scan = fid._read_json(corpus_dir / "scan_detection.json") or {}
        d["era"] = fid.era_bucket(meta.get("year"))
        d["file_type"] = scan.get("file_type")
        docs[stem] = d
        kinds.update(d["block_kinds"])
        for k in ("gold_numbers", "found_numbers", "matched_numbers",
                  "captions_compared"):
            totals[k] += d[k]
            for axis in segments:
                segments[axis][str(d[axis])][k] += d[k]
        served = score_document(
            gold_dir, corpus_dir, figure_types=EVIDENCE_FIGURE_TYPES,
        )
        for axis in ("era", "file_type"):
            served[axis] = d[axis]
        for k in ("gold_numbers", "found_numbers", "matched_numbers",
                  "captions_compared"):
            served_totals[k] += served[k]
            for axis in served_segments:
                served_segments[axis][str(served[axis])][k] += served[k]

    def pack_segments(source):
        return {
            axis: {
                key: _pack_counts(counts)
                for key, counts in sorted(buckets.items())
            }
            for axis, buckets in source.items()
        }

    return {
        "schema_version": 3,
        "question": "is the caption bound to the figure it belongs to",
        "method": "figure numbers parsed independently from labels printed "
                  "INSIDE a gold [FIGURE]/[PLATE] block, compared per page "
                  "against figures.json; caption text scored only within "
                  "number-matched pairs",
        "caveat": "Naive caption-text similarity reports 44% and is mostly "
                  "artifact — bilingual captions, bare labels, and documents "
                  "that are their own translation. Numbers are the "
                  "language-independent signal.",
        "documents_bound": len(bound), "documents_unmatched": unmatched,
        "gold_block_kinds": dict(kinds),
        "totals": _pack_counts(totals),
        "segments": pack_segments(segments),
        "surfaces": {
            "all entries": {
                "totals": _pack_counts(totals),
                "segments": pack_segments(segments),
            },
            "default MCP types": {
                "totals": _pack_counts(served_totals),
                "segments": pack_segments(served_segments),
            },
        },
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

    w("\n-- by retrieval surface " + "-" * 30 + "\n")
    w(f"{'surface':24}{'gold':>7}{'found':>7}{'recall':>9}{'precis':>9}\n")
    for name, surface in report["surfaces"].items():
        values = surface["totals"]
        def surface_rate(value):
            return f"{value:>9.3f}" if isinstance(value, float) else f"{'-':>9}"
        w(f"{name[:24]:24}{values['gold_numbers']:>7}"
          f"{values['found_numbers']:>7}{surface_rate(values['number_recall'])}"
          f"{surface_rate(values['number_precision'])}\n")

    for axis in ("era", "layout"):
        w(f"\n-- all entries by {axis} "
          + "-" * max(1, 36 - len(axis)) + "\n")
        w(f"{'bucket':24}{'gold':>7}{'found':>7}{'recall':>9}{'precis':>9}\n")
        for bucket, values in report["segments"][axis].items():
            def rate(value):
                return f"{value:>9.3f}" if isinstance(value, float) else f"{'-':>9}"
            w(f"{bucket[:24]:24}{values['gold_numbers']:>7}"
              f"{values['found_numbers']:>7}{rate(values['number_recall'])}"
              f"{rate(values['number_precision'])}\n")

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
