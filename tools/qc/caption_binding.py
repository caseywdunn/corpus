#!/usr/bin/env python3
"""Score caption association against the gold transcription set (#195).

Detection asks whether the figure *objects* are right (`figure_detection.py`).
This asks the next questions: is the caption attached to the figure it belongs
to, and does its declared panel set match? `dev_docs/OVERVIEW.md` calls caption
association "the highest-value annotation per figure and the hardest in
historical layouts"; this turns that from an assertion into numbers.

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
that page. Restricted to numbers printed in a `[FIGURE]`/`[PLATE]` block plus
a plate heading immediately adjacent to that block — the figures actually on
the leaf — the measure means what it says. Plate blocks may express engraved
numbers as explicit inventories, standalone digits, or `F. N` labels; those
are gold-format conventions, not production OCR heuristics.

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

# Panel labels in the gold transcription are parsed independently of
# pipeline.figures.parse_panels_from_caption for the same reason figure
# numbers are: importing the implementation under test would make the
# yardstick move with it.  This intentionally covers only the production
# parser's stated A--L, lettered-panel contract. Numeric parts of a grouped
# plate are figure targets, not panels (#203).
_GOLD_PANEL_RANGE_RE = re.compile(
    r"(?<![A-Za-z])([A-L])\s*[-\u2013\u2014]\s*([A-L])(?=[.,:;\s)]|$)"
)
_GOLD_PANEL_PAREN_RE = re.compile(
    r"(?<![A-Za-z])\(\s*([A-L])\s*\)(?![A-Za-z])"
)
_GOLD_PANEL_PERIOD_RE = re.compile(r"(?<![A-Za-z.])([A-L])\.(?=\s+\S)")
_GOLD_PANEL_COMMA_RE = re.compile(r"(?<![A-Za-z])([A-L])\s*,(?=\s)")
_GOLD_INITIAL_SURNAME_RE = re.compile(r"^\s+[A-Z][a-z\u00c0-\u00ff]+")
_GOLD_CREDIT_CONTEXT_RE = re.compile(
    r"(?:credit|courtesy|photo(?:graph)?s?\s+by|drawn\s+by|"
    r"\bby(?:\s+the\s+late)?|after|from)"
    r"[A-Za-z\u00c0-\u00ff.\s]{0,32}$",
    re.IGNORECASE,
)
_GOLD_AUTHORITY_TAIL_RE = re.compile(
    r"^\s+[A-Z][a-z\u00c0-\u00ff]+\s*(?:\)|,\s*\d{4})"
)
# ``Pl.M.`` is the conventional mouth-plate key inside siphonophore drawings,
# not an abbreviation of "Plate M". The generic opener grammar cannot tell
# those strings apart after converting Roman M to 1000, so reject the compact
# anatomical spelling before applying it. A printed plate label uses ``Plate``
# or a spaced abbreviation such as ``Pl. IV`` in the gold set.
_GOLD_MOUTH_PLATE_LABEL_RE = re.compile(r"^\s*Pl\.M\.?(?:\s|$)", re.IGNORECASE)

_NUMBER_COUNT_KEYS = (
    "gold_numbers", "found_numbers", "matched_numbers",
    "reported_pair_capacity", "captions_compared",
)
_PANEL_COUNT_KEYS = (
    "gold_panel_figures", "reported_panel_figures", "exact_panel_figures",
    "gold_panel_labels", "reported_panel_labels", "matched_panel_labels",
    "gold_unpanelled_figures", "number_matched_unpanelled_figures",
    "unexpected_panel_figures",
)
_COUNT_KEYS = _NUMBER_COUNT_KEYS + _PANEL_COUNT_KEYS


# Gold labels are parsed independently of the production extractor. Importing
# parse_figure_number/caption_figure_entries here made the yardstick move when
# those functions changed: the Totton source-citation fix changed the gold
# denominator from 496 to 480 even though no transcription changed. The gold
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
_GOLD_PLATE_NUMBER_LIST_RE = re.compile(
    r"^\s*Figure numbers on the plate(?:,[^:]*)?\s*:\s*(?P<numbers>.+?)\.?\s*$",
    re.IGNORECASE,
)
_GOLD_BARE_PLATE_NUMBER_RE = re.compile(
    r"^\s*(?:F\.\s*)?(\d+[a-z]?)\.?\s*$",
    re.IGNORECASE,
)
_GOLD_PLATE_HEADING_RE = re.compile(
    r"^\s*(?:plate|pl\.?)\s+" + _GOLD_NUMBER + r"(?:\s|$)",
    re.IGNORECASE,
)


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
    if _GOLD_MOUTH_PLATE_LABEL_RE.match(text):
        return []
    if not _GOLD_ENTRY_START_RE.match(text):
        return []
    matches = []
    for match in _GOLD_ENTRY_RE.finditer(text):
        if _GOLD_MOUTH_PLATE_LABEL_RE.match(text[match.start():]):
            continue
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


def gold_plate_number_list(line):
    """Numbers in the gold set's explicit bare-plate inventory sentence.

    Historical plate transcriptions record engraved labels using a controlled
    sentence inside ``[PLATE]``. It is annotation about source-visible labels,
    not caption prose, so it contributes to the number denominator without
    turning a bare-label block into a prose-caption test.
    """
    match = _GOLD_PLATE_NUMBER_LIST_RE.match(line or "")
    if not match:
        return []
    numbers = []
    for raw in match.group("numbers").split(","):
        token = raw.strip()
        if not re.fullmatch(_GOLD_NUMBER, token, re.IGNORECASE):
            return []
        number = _gold_number(token)
        if number not in numbers:
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
        depth, kind, start, opening = 0, None, None, None
        for a, b in fid.spans(text):
            inner = text[a + 1:b - 1]
            if inner in fid.STRUCTURAL_OPEN:
                if depth == 0:
                    kind, start, opening = inner, b, a
                depth += 1
            elif inner in fid.STRUCTURAL_CLOSE:
                depth -= 1
                if depth == 0 and start is not None:
                    if kind in GOLD_FIGURE_TAGS:
                        body = text[start:a]
                        # The transcription set has both faithful markup
                        # shapes: ``PLATE XVII`` can be inside [PLATE] or on
                        # the immediately preceding nonblank line. Only that
                        # adjacent line is eligible; a plate reference
                        # elsewhere on the page is not host identity.
                        preceding = text[:opening].rstrip().splitlines()
                        plate_heading = (
                            preceding[-1].strip() if kind == "PLATE"
                            and preceding
                            and _GOLD_PLATE_HEADING_RE.match(preceding[-1])
                            else ""
                        )
                        if plate_heading and plate_heading not in {
                                    line.strip() for line in body.splitlines()
                                }:
                            body = f"{plate_heading}\n{body}"
                        yield page, kind, body
                    start = opening = None


def gold_block_number_sets(body, structural_kind=None):
    """Return figure and plate numbers without collapsing their namespaces."""
    stripped = fid.strip_markup(body).text
    figure_numbers, plate_numbers = set(), set()
    first_figure_line = first_plate_line = None
    for line in stripped.splitlines():
        line = line.strip()
        if structural_kind == "PLATE" and _GOLD_PLATE_HEADING_RE.match(line):
            entries = gold_caption_entries(line)
            plate_numbers.update(entries)
            if first_plate_line is None and entries:
                first_plate_line = line
            continue
        listed_plate_figures = gold_plate_number_list(line)
        if listed_plate_figures:
            figure_numbers.update(listed_plate_figures)
            # This controlled transcription sentence describes the labels
            # printed on the image; it is not source caption text.
            continue
        if structural_kind == "PLATE":
            bare_number = _GOLD_BARE_PLATE_NUMBER_RE.fullmatch(line)
            if bare_number:
                figure_numbers.add(_gold_number(bare_number.group(1)))
                # A standalone engraved label is number evidence but not
                # source caption prose. Restrict this to [PLATE]: standalone
                # numerals in [FIGURE] commonly belong to chart axes/scales.
                continue
        # A gold block can contain descriptive prose and references. Only a
        # line that opens as a caption is allowed to contribute; once it does,
        # split any additional entries/lists on that same printed line.
        entries = gold_caption_entries(line)
        figure_numbers.update(entries)
        if first_figure_line is None and entries:
            first_figure_line = line
    return figure_numbers, plate_numbers, first_figure_line, first_plate_line


def classify_block(body, structural_kind=None):
    """What kind of evidence this gold block offers.

    Only ``prose_caption`` blocks can test caption *text*; the rest are
    counted so no rate is computed over blocks it does not apply to.
    """
    stripped = fid.strip_markup(body).text
    figure_numbers, plate_numbers, first_figure_line, first_plate_line = \
        gold_block_number_sets(body, structural_kind)
    numbers = figure_numbers | plate_numbers
    first_line = first_figure_line or first_plate_line
    tokens = fid.tokens(stripped)
    if not tokens:
        return "nothing_printed", numbers, ""
    if not numbers:
        return "lettering_only", numbers, ""
    # The controlled inventory sentence is transcription metadata for labels
    # visible on an otherwise uncaptioned plate. Later ``Fig. N`` lines record
    # lettering inside those drawings, even when a bracketed uncertainty makes
    # such a line long enough to resemble prose. The block still offers no
    # source caption text to compare.
    if any(gold_plate_number_list(line.strip()) for line in stripped.splitlines()):
        return "bare_label", numbers, first_line
    label_tokens = fid.tokens(first_line or "")
    if len(label_tokens) >= _PROSE_CAPTION_MIN_TOKENS:
        return "prose_caption", numbers, first_line
    return "bare_label", numbers, first_line


def gold_panel_labels(body):
    """Return independently parsed lettered panels from a gold caption.

    Figure blocks also contain lettering printed *inside* the image.  Only
    text beginning at the first caption opener is eligible, otherwise a row
    of engraved ``A`` / ``B`` labels would be mistaken for a caption that
    enumerates panels. A result must begin at A and contain at least two labels:
    this scorer measures a panel *split*, not an isolated letter that might be
    an abbreviation. Printed gaps are retained; Totton Figure 74 deliberately
    runs A--H, K, L rather than A--L.
    """
    stripped = fid.strip_markup(body).text
    lines = stripped.splitlines()
    start = next(
        (i for i, line in enumerate(lines)
         if not _GOLD_MOUTH_PLATE_LABEL_RE.match(line)
         and _GOLD_ENTRY_START_RE.match(line.strip())),
        None,
    )
    if start is None:
        return []
    caption = "\n".join(lines[start:])

    labels = set()
    range_spans = []
    for match in _GOLD_PANEL_RANGE_RE.finditer(caption):
        first, last = match.group(1), match.group(2)
        if ord(last) <= ord(first):
            continue
        labels.update(chr(code) for code in range(ord(first), ord(last) + 1))
        range_spans.append((match.start(), match.end()))

    def in_range(position):
        return any(first <= position < last for first, last in range_spans)

    for pattern in (
        _GOLD_PANEL_PAREN_RE,
        _GOLD_PANEL_PERIOD_RE,
        _GOLD_PANEL_COMMA_RE,
    ):
        for match in pattern.finditer(caption):
            if in_range(match.start()):
                continue
            # An author/photo credit initial is not a panel.  Keep this
            # independent from production's helper and limit it to contexts
            # that are unambiguous in a hand transcription.
            if pattern is _GOLD_PANEL_PERIOD_RE \
                    and (
                        (
                            _GOLD_CREDIT_CONTEXT_RE.search(
                                caption[:match.start()]
                            )
                            and not re.search(
                                r"[a-z\u00e0-\u00ff]\.\s*$",
                                caption[:match.start()],
                            )
                        )
                        or re.search(r"(?:\(|&)\s*$", caption[:match.start()])
                        or re.match(
                            r"\s*&\s*(?:[a-z]\.\s*)?[A-Z]\.",
                            caption[match.end():],
                        )
                        or (
                            _GOLD_INITIAL_SURNAME_RE.match(caption[match.end():])
                            and _GOLD_AUTHORITY_TAIL_RE.match(
                                caption[match.end():]
                            )
                        )
                    ):
                continue
            labels.add(match.group(1))

    ordered = sorted(labels)
    if len(ordered) < 2 or ordered[0] != "A":
        return []
    return ordered


def score_panel_case(page, figure_number, expected_labels, by_page):
    """Score one gold panel enumeration against same-page, same-number data."""
    candidates = [
        figure for figure in by_page.get(page, [])
        if str(figure.get("figure_number") or "") == str(figure_number)
    ]
    reported = {
        str(panel.get("label") or "").upper()
        for figure in candidates
        for panel in (figure.get("panels_from_caption") or [])
        if re.fullmatch(r"[A-L]", str(panel.get("label") or "").upper())
    }
    expected = set(expected_labels)
    return {
        "figure_number": str(figure_number),
        "number_bound": bool(candidates),
        "expected_labels": sorted(expected),
        "reported_labels": sorted(reported),
        "matched_labels": sorted(expected & reported),
        "missing_labels": sorted(expected - reported),
        "surplus_labels": sorted(reported - expected),
        "reported": len(reported) >= 2,
        "exact": reported == expected,
    }


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


def _identity_text(identity):
    namespace, number = identity
    return f"{namespace}:{number}"


def _identity_key(identity):
    namespace, number = identity
    return (namespace, *_number_key(number))


def _pack_counts(counts):
    """Add rates to one aggregate without hiding undefined precision."""
    gold = counts["gold_numbers"]
    found = counts["found_numbers"]
    matched = counts["matched_numbers"]
    recall = matched / gold if gold else None
    precision = matched / found if found else None
    capacity = counts.get("reported_pair_capacity", 0)
    capacity_rate = capacity / gold if gold else None
    f1 = None
    if recall is not None and precision is not None:
        f1 = (
            2 * recall * precision / (recall + precision)
            if recall + precision else 0.0
        )
    panel_figures = counts.get("gold_panel_figures", 0)
    panel_labels = counts.get("gold_panel_labels", 0)
    reported_labels = counts.get("reported_panel_labels", 0)
    reported_panel_figures = counts.get("reported_panel_figures", 0)
    unexpected_panel_figures = counts.get("unexpected_panel_figures", 0)
    return {
        **{k: counts.get(k, 0) for k in _COUNT_KEYS},
        "number_recall": round(recall, 4) if recall is not None else None,
        "number_precision": round(precision, 4) if precision is not None else None,
        "reported_pair_capacity_rate": (
            round(capacity_rate, 4) if capacity_rate is not None else None
        ),
        "f1": round(f1, 4) if f1 is not None else None,
        "panel_figure_recall": (
            round(reported_panel_figures / panel_figures, 4)
            if panel_figures else None
        ),
        "panel_figure_precision": (
            round(
                reported_panel_figures
                / (reported_panel_figures + unexpected_panel_figures),
                4,
            )
            if reported_panel_figures + unexpected_panel_figures else None
        ),
        "panel_exact_rate": (
            round(counts.get("exact_panel_figures", 0) / panel_figures, 4)
            if panel_figures else None
        ),
        "panel_label_recall": (
            round(counts.get("matched_panel_labels", 0) / panel_labels, 4)
            if panel_labels else None
        ),
        "panel_label_precision": (
            round(counts.get("matched_panel_labels", 0) / reported_labels, 4)
            if reported_labels else None
        ),
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
    gold_nums, found_nums = defaultdict(set), defaultdict(set)
    gold_ids, found_ids = defaultdict(set), defaultdict(set)
    gold_blocks_by_page = defaultdict(list)
    gold_panels_by_page = defaultdict(list)
    unexpected_panels_by_page = defaultdict(list)
    panel_counts = Counter()
    captions = []                      # (gold_caption, extracted_caption)

    for page, structural_kind, body in gold_blocks(gold_dir, scan):
        kind, numbers, first_line = classify_block(body, structural_kind)
        figure_numbers, plate_numbers, _first_figure, _first_plate = \
            gold_block_number_sets(body, structural_kind)
        identities = (
            {("figure", number) for number in figure_numbers}
            | {("plate", number) for number in plate_numbers}
        )
        kinds[kind] += 1
        structural_kinds[structural_kind] += 1
        gold_nums[page] |= numbers
        gold_ids[page] |= identities
        gold_blocks_by_page[page].append({
            "structural_kind": structural_kind,
            "caption_kind": kind,
            "numbers": sorted(numbers, key=_number_key),
            "identities": [
                _identity_text(identity)
                for identity in sorted(identities, key=_identity_key)
            ],
        })
        expected_panels = gold_panel_labels(body)
        # A grouped legend naming several numbered figures is deliberately
        # outside the panel denominator.  Its numeric children are scored by
        # the binding measure above (#203).
        if len(numbers) == 1:
            case = score_panel_case(
                page, next(iter(numbers)), expected_panels, by_page,
            )
            if len(expected_panels) >= 2:
                gold_panels_by_page[page].append(case)
                panel_counts["gold_panel_figures"] += 1
                panel_counts["reported_panel_figures"] += int(case["reported"])
                panel_counts["exact_panel_figures"] += int(case["exact"])
                panel_counts["gold_panel_labels"] += len(case["expected_labels"])
                panel_counts["reported_panel_labels"] += len(case["reported_labels"])
                panel_counts["matched_panel_labels"] += len(case["matched_labels"])
            else:
                panel_counts["gold_unpanelled_figures"] += 1
                panel_counts["number_matched_unpanelled_figures"] += int(
                    case["number_bound"]
                )
                panel_counts["unexpected_panel_figures"] += int(case["reported"])
                if case["reported"]:
                    unexpected_panels_by_page[page].append(case)
        if kind == "prose_caption" and figure_numbers:
            for f in by_page.get(page, []):
                if f.get("figure_type") != "plate" \
                        and str(f.get("figure_number")) in figure_numbers:
                    captions.append((first_line, f.get("caption_text") or ""))
                    break

    for page, on_page in by_page.items():
        found_nums[page] |= {str(f["figure_number"]) for f in on_page
                             if f.get("figure_number")}
        found_ids[page] |= {
            (
                "plate" if f.get("figure_type") == "plate" else "figure",
                str(f["figure_number"]),
            )
            for f in on_page if f.get("figure_number")
        }

    g = sum(len(v) for v in gold_ids.values())
    j = sum(len(v) for v in found_ids.values())
    matched = sum(len(gold_ids[p] & found_ids.get(p, set())) for p in gold_ids)
    # This is the best recall possible if the existing count of distinct typed
    # identities on each page is held fixed and every identity may be
    # relabelled. It diagnoses missing upstream evidence; it is not a
    # theoretical ceiling on a rebuild that discovers additional labels.
    reported_pair_capacity = sum(
        min(len(identities), len(found_ids.get(page, set())))
        for page, identities in gold_ids.items()
    )

    sims = []
    for gold_cap, got_cap in captions:
        if got_cap:
            sims.append(SequenceMatcher(None, fid.tokens(gold_cap),
                                        fid.tokens(got_cap),
                                        autojunk=False).ratio())
    page_diagnostics = {}
    for page in sorted(set(gold_ids) | set(found_ids)):
        gold_on_page = gold_nums.get(page, set())
        found_on_page = found_nums.get(page, set())
        gold_identities = gold_ids.get(page, set())
        found_identities = found_ids.get(page, set())
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
            "gold_identities": [
                _identity_text(identity)
                for identity in sorted(gold_identities, key=_identity_key)
            ],
            "found_identities": [
                _identity_text(identity)
                for identity in sorted(found_identities, key=_identity_key)
            ],
            "matched_identities": [
                _identity_text(identity)
                for identity in sorted(
                    gold_identities & found_identities, key=_identity_key,
                )
            ],
            "missing_identities": [
                _identity_text(identity)
                for identity in sorted(
                    gold_identities - found_identities, key=_identity_key,
                )
            ],
            "surplus_identities": [
                _identity_text(identity)
                for identity in sorted(
                    found_identities - gold_identities, key=_identity_key,
                )
            ],
            "reported_figures": figures,
            "gold_panel_cases": gold_panels_by_page.get(page, []),
            "unexpected_panel_cases": unexpected_panels_by_page.get(page, []),
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
        "reported_pair_capacity": reported_pair_capacity,
        "number_recall": round(matched / g, 4) if g else None,
        "number_precision": round(matched / j, 4) if j else None,
        "reported_pair_capacity_rate": (
            round(reported_pair_capacity / g, 4) if g else None
        ),
        "captions_compared": len(sims),
        **{key: panel_counts[key] for key in _PANEL_COUNT_KEYS},
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
        for k in _COUNT_KEYS:
            totals[k] += d[k]
            for axis in segments:
                segments[axis][str(d[axis])][k] += d[k]
        served = score_document(
            gold_dir, corpus_dir, figure_types=EVIDENCE_FIGURE_TYPES,
        )
        for axis in ("era", "file_type"):
            served[axis] = d[axis]
        for k in _COUNT_KEYS:
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
        "schema_version": 7,
        "question": "is the caption bound to the figure it belongs to, and "
                    "does its declared panel set match",
        "method": "figure numbers parsed independently from labels printed "
                  "inside a gold [FIGURE]/[PLATE] block, including controlled "
                  "plate inventories and standalone engraved-number lines; "
                  "an immediately adjacent plate heading is included as "
                  "plate identity evidence. Compared per page against "
                  "typed figures.json identities (plate:N is distinct from "
                  "figure:N); caption text is scored only within "
                  "number-matched pairs",
        "caveat": "Naive caption-text similarity reports 44% and is mostly "
                  "artifact — bilingual captions, bare labels, and documents "
                  "that are their own translation. Numbers are the "
                  "language-independent signal. Reported-pair capacity is "
                  "the maximum same-page typed-identity recall possible if "
                  "the current number of distinct reported identities on "
                  "each page is held fixed; it diagnoses missing evidence, "
                  "not a ceiling on a rebuild that discovers more labels.",
        "panel_method": "lettered panels parsed independently from the "
                        "gold caption text after its figure opener; only "
                        "explicit multi-label sets beginning at A are scored, "
                        "printed gaps are retained, and grouped numeric "
                        "plates are excluded; declaration precision uses "
                        "single-number gold blocks with a same-page, "
                        "same-number extracted record as its negative set",
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
    w(f"\ntyped identities in gold           : {t['gold_numbers']}\n")
    w(f"typed identities pipeline reports : {t['found_numbers']}\n")
    w(f"  binding recall                   : {t['number_recall']}\n")
    w(f"  binding precision                : {t['number_precision']}\n")
    w(f"  fixed reported-pair capacity     : "
      f"{t['reported_pair_capacity_rate']}\n")

    w("\n-- by retrieval surface " + "-" * 30 + "\n")
    w(f"{'surface':24}{'gold':>7}{'found':>7}{'recall':>9}{'precis':>9}\n")
    for name, surface in report["surfaces"].items():
        values = surface["totals"]
        def surface_rate(value):
            return f"{value:>9.3f}" if isinstance(value, float) else f"{'-':>9}"
        w(f"{name[:24]:24}{values['gold_numbers']:>7}"
          f"{values['found_numbers']:>7}{surface_rate(values['number_recall'])}"
          f"{surface_rate(values['number_precision'])}\n")

    panel = report["surfaces"]["default MCP types"]["totals"]
    w("\n-- caption-enumerated panel splits " + "-" * 24 + "\n")
    w(f"gold panelled figures             : {panel['gold_panel_figures']}\n")
    w(f"panel declarations reported       : {panel['reported_panel_figures']}\n")
    w(f"exact label sets                  : {panel['exact_panel_figures']}\n")
    w(f"  panel figure recall             : {panel['panel_figure_recall']}\n")
    w(f"  panel figure precision          : {panel['panel_figure_precision']}\n")
    w(f"  exact-set rate                  : {panel['panel_exact_rate']}\n")
    w(f"  label recall / precision        : "
      f"{panel['panel_label_recall']} / {panel['panel_label_precision']}\n")
    w(f"unexpected declarations on "
      f"{panel['number_matched_unpanelled_figures']} same-page, same-number "
      f"non-panel figure blocks: "
      f"{panel['unexpected_panel_figures']}\n")

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
