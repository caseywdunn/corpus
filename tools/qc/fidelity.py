#!/usr/bin/env python3
"""Score a built corpuscle's extraction against an independent gold transcription set (#193).

Every other quality signal this pipeline has measures it against *itself*.
The soft consistency rates in ``tests/test_corpus_wide.py`` are
corpus-internal agreement, the per-document quality gates are plausibility
checks, and fingerprint diffing compares one build to the next. None of them
can say whether the text is *right*.

The gold set can. It is a set of documents transcribed page by page from a
rendered page image only, with the transcriber forbidden from opening any
software extraction of the page they are working on. That independence
guarantee is the whole value: an extractor scored against this set is not
being compared to a cleaned-up version of itself.

This script lives here rather than in the library repo that holds the
transcriptions so that the evaluator is versioned with the extractor it
grades.

WHICH SIDE IS ON TRIAL, BECAUSE IT IS THE EASIEST THING TO GET BACKWARDS
========================================================================
The library repo ships ``scripts/crosscheck.py``, which compared gold against
a poppler text layer **in order to validate the gold**, and a companion
``CROSSCHECK_REPORT.md`` documenting five ways that signal misleads. Read that
report before trusting any number here — but read it knowing that the thing on
trial has changed. There, poppler was the yardstick and the transcription was
the subject; here the gold is the yardstick and **the corpuscle is the
subject**.

Note what does *not* change: the arithmetic. ``coverage`` is
``|gold ∩ other| / |gold|`` in both scripts and ``recall`` is
``|gold ∩ other| / |other|`` in both — the measures are symmetric and gold is
the reference term on both sides. Nothing is mirrored, transposed or
recomputed. What changes is two judgement calls that the arithmetic cannot
make for you:

  1. **Which measure you read first.** ``recall`` leads there, because "poppler
     saw words the gold lacks" is the shape of a transcription omission.
     ``coverage`` leads here, because "the page holds words the pipeline
     missed" is the shape of an extraction failure.
  2. **What an unscorable page counts as** — see the status handling in
     ``score_page`` below, and the note on it further down.

Those two calls are why several of the report's false positives are this
script's true positives:

  * A garbage text layer (roman OCR over Fraktur, ~0.05-0.15) tells the report
    nothing about the gold page. Here it is a real extraction failure — OCR is
    exactly what the pipeline is supposed to bring to that page — so it is
    scored, not excluded.
  * The gold correctly holding *more* than the comparison layer penalised the
    gold there. Here it is a real miss: whole-page plates carry engraved
    lettering that exists only as image, and the pipeline is meant to OCR it.
  * Text hallucinated out of image texture was noise in the comparison layer
    there. Here it is a pipeline defect, and ``excess_wordlike`` finds it.

Consequently ``coverage`` is the primary measure in this script, where the
report correctly names ``recall`` primary for its own purpose. Coverage asks
how much of the true page the pipeline recovered; recall asks how much of what
the pipeline emitted is actually on the page.

The second judgement call is the one with teeth, and it is worth stating what
it costs to get wrong. The report excludes a page whose comparison layer
carries no usable signal, correctly: an empty poppler layer says nothing about
the gold. Here an empty extraction is not an absence of evidence, it *is* the
evidence, so it is scored 0.0 and stays in every median. Measured on the
35-document set, adopting the report's exclusion policy instead would drop 57
of 761 pages and move the median coverage from 0.891 to 0.908 — hiding, with
precision, exactly the pages that need work.

Two of the report's caveats do carry over unchanged, and no amount of
arithmetic here escapes them:

  * **This cannot detect invention.** Fluent invented prose does not diff as
    an error against a page that simply lacks it.
  * **A low volume ratio is not evidence of loss.** It is as often the other
    side carrying more noise.

Alignment needs no pipeline change. ``text.json`` carries a flat ``text``
string and a page *count* only, but ``docling_doc.json`` retains
``prov[].page_no`` on every text item, so per-page extracted text is
reconstructible from the persisted artifact and directly comparable to the
gold's ``page_NNN.txt``. Reconstruction follows ``body.children`` reading
order, which is the order ``export_to_markdown()`` uses to build the flat
string that ends up in ``text.json`` — so the reading order being scored is
the one a consumer actually sees.

Documents are bound **on the source PDF's checksum, never on its filename**.
Corpuscle document directories are named for the first 12 hex characters of
that checksum, and the gold set's ``sources.json`` records the full digest.
Filenames are not unique across a library: two editions of Lery's *Histoire
d'un Voyage* were both once named ``Lery1594.pdf`` and set the same *Physalia*
passage on different folios, with only one of them transcribed.

Usage::

    tools/qc/fidelity.py --gold <transcriptions/> --corpuscle <corpuscle/> \
        [--out fidelity.json] [--quiet]

The full run is manual and release-time. ``tests/test_fidelity_harness.py``
covers the arithmetic against a committed fixture, with no corpus dependency.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys
import unicodedata
from collections import Counter, defaultdict
from difflib import SequenceMatcher, get_close_matches

# Below this many extracted characters on a page that the gold says carries
# text, there is nothing to compare and the extraction has simply failed.
# Note the difference from the library-side cross-check, which treats the same
# condition as "no signal" and excludes the page: there, an empty comparison
# layer says nothing about the gold; here, it is the finding.
MIN_EXTRACTED_CHARS = 50

# --- gold markup ------------------------------------------------------------
#
# The stripping rules below are ported from the library repo's crosscheck.py so
# that the two sets of numbers stay comparable. Divergence here would make a
# score from this script silently incommensurable with the ones in
# CROSSCHECK_REPORT.md.

# Markers whose bracketed payload is transcriber commentary, not printed text.
# These must be dropped before diffing or they inflate the disagreement.
DROP_WHOLE = ("NOTE:", "PAGE ", "BLANK PAGE", "illegible", "no lettering")

# Markers that label printed text: drop the tag, keep what follows.
KEEP_AFTER_TAG = ("RUNNING HEAD", "FOOTNOTE", "MARGIN", "STAMP", "HANDWRITTEN")

# Paired structural tags. The tag itself is dropped and the block contents are
# kept, but we track the nesting so a page's figure/plate/table share of the
# text can be reported separately — see `_StrippedGold.structural`.
STRUCTURAL_OPEN = ("FIGURE", "PLATE", "TABLE")

# Which of those count as "inside a figure" when prose is separated from figure
# text. Tables are body content — text the pipeline is expected to get right —
# so they stay on the prose side. Measured either way the corpus-wide prose
# coverage moves by 0.002, so this is a naming choice rather than a lever.
FIGURE_TAGS = ("FIGURE", "PLATE")
STRUCTURAL_CLOSE = tuple("/" + t for t in STRUCTURAL_OPEN)

# `[PAGE 8]` records the printed folio; `[PAGE 12: unnumbered]` records that the
# leaf carries none. Kept per page because page-number provenance on served
# figures is an open question (#188) and this set is the only place where the
# printed folio and the PDF page index are both known.
_PAGE_MARKER = re.compile(r"^PAGE\s+(\d+)(?::\s*(.*))?$")


# Every keyword that can legally open a marker. A ``[`` followed by anything
# else is not markup — it is a bracket printed on the page, and this set is what
# tells the two apart.
#
# Counting nesting over *every* bracket instead, as a naive scanner does, breaks
# on two things this corpus actually contains. Both were found by checking the
# 761 gold pages for brackets that never close:
#
#   * A note that quotes a bracket while explaining the convention —
#     `[NOTE: ... Everywhere else on this page "[" is my marker.]` — leaves the
#     scanner one level deep, so the note never closes and its entire text
#     leaks into the compared page as though it had been printed there. 13 of
#     one document's 17 pages, and the reason that document posted 0.767
#     coverage against 0.998 recall: the gold appeared to hold text no
#     extractor could ever find, because the page does not contain it.
#   * A transcription typo — `[continued opposite` with no closing bracket —
#     swallows the `[/FIGURE]` that follows it, so the figure block never ends.
MARKER_PREFIXES = (
    "NOTE", "PAGE", "BLANK PAGE", "illegible", "no lettering",
    "RUNNING HEAD", "FOOTNOTE", "MARGIN", "/MARGIN", "STAMP", "HANDWRITTEN",
    "FIGURE", "/FIGURE", "PLATE", "/PLATE", "TABLE", "/TABLE",
    "?",
)


def _marker_at(text, i):
    """True if a marker keyword opens at ``text[i]``."""
    return text.startswith("[", i) and any(
        text.startswith(m, i + 1) for m in MARKER_PREFIXES)


# A note's prose can contain a bracket that is not markup, and it comes in two
# forms that need opposite handling. Both were found by parsing all 761 pages
# and asking which notes failed to close:
#
#   * A bracketed expression the note quotes from the page — `[sic]` in a
#     translator's insertion, `[21]` in a citation the note reproduces. It has
#     a closing bracket, and that bracket must not be mistaken for the note's
#     own: doing so ends the note early, and every marker after it is then read
#     at the wrong nesting level.
#   * A mention of the bracket *character*, always written quoted:
#     `Everywhere else on this page "[" is my marker.` There is no partner to
#     absorb, so it must simply be ignored.
#
# A quote immediately after the `[` is what tells them apart, and it is a
# property of the convention rather than a guess. The width cap is a second
# guard: a bracketed expression a note quotes is short, and refusing to skip a
# long one keeps a single stray bracket from swallowing a whole note.
_LITERAL_MAX_CHARS = 40
_QUOTES = "\"'“”‘’"


def _literal_pair_end(text, j):
    """Index past the ``]`` of a non-marker bracket pair at ``j``, else ``None``."""
    if j + 1 < len(text) and text[j + 1] in _QUOTES:
        return None                      # `["` — a mention of the character
    stop = min(len(text), j + 1 + _LITERAL_MAX_CHARS)
    k = j + 1
    while k < stop:
        if text[k] == "]":
            return k + 1
        if text[k] == "[" and _marker_at(text, k):
            return None                  # a real marker opens first
        k += 1
    return None


def _span_end(text, i):
    """Index just past the ``]`` closing the bracketed span opening at ``i``.

    ``None`` when the span never closes, which is the case the caller must not
    paper over: consuming an unterminated ``[`` swallows every marker after it.
    """
    if _marker_at(text, i):
        # Markers nest inside markers — notes routinely quote a `[FIGURE]` —
        # but only *markers* nest, and a non-marker bracket in the note's prose
        # is handled by `_literal_pair_end` rather than counted as a level.
        depth, j = 1, i + 1
        while j < len(text):
            ch = text[j]
            if ch == "[":
                if _marker_at(text, j):
                    depth += 1
                    j += 1
                    continue
                k = _literal_pair_end(text, j)
                j = k if k is not None else j + 1
                continue
            if ch == "]":
                depth -= 1
                if depth == 0:
                    return j + 1
            j += 1
        return None
    # A non-marker bracket: printed content, or ad-hoc commentary. It does not
    # nest, and if a real marker opens before its closing bracket then it never
    # had one — it is a stray `[` and must be left as literal text.
    j = i + 1
    while j < len(text):
        if text[j] == "]":
            return j + 1
        if text[j] == "[" and _marker_at(text, j):
            return None
        j += 1
    return None


def spans(text):
    """Yield ``(start, end)`` of bracketed spans, in order and non-overlapping.

    A ``[`` that opens nothing — an unterminated one — is skipped rather than
    consumed, so the markers after it are still seen.
    """
    i, n = 0, len(text)
    while i < n:
        if text[i] != "[":
            i += 1
            continue
        end = _span_end(text, i)
        if end is None:
            i += 1
            continue
        yield i, end
        i = end


def unterminated_brackets(text):
    """Count ``[`` characters that open no span — a gold-integrity signal.

    Reported per page rather than silently tolerated. On this corpus the count
    is a mix of the harmless (a note quoting a bracket, now handled) and the
    genuine (a transcription typo), and only inspection tells them apart.
    """
    n, i = 0, 0
    while i < len(text):
        if text[i] == "[":
            end = _span_end(text, i)
            if end is None:
                n += 1
                i += 1
                continue
            i = end
        else:
            i += 1
    return n


class _StrippedGold:
    """A gold page reduced to what was printed on it.

    ``text`` is everything printed; ``structural`` is the subset that sat
    inside a ``[FIGURE]``/``[PLATE]``/``[TABLE]`` block. The split matters
    because a page whose text is mostly engraved plate lettering fails in a
    completely different way from a page of prose, and a single coverage number
    over both is not actionable.
    """

    __slots__ = ("text", "structural", "prose", "printed_page", "unnumbered")

    def __init__(self, text, structural, prose, printed_page, unnumbered):
        self.text = text
        self.structural = structural
        self.prose = prose
        self.printed_page = printed_page
        self.unnumbered = unnumbered


def strip_markup(text):
    """Reduce a gold page to just the text printed on the page."""
    out, structural_out, prose_out, pos, depth = [], [], [], 0, 0
    fig_depth = 0
    printed_page, unnumbered = None, False

    def emit(chunk):
        out.append(chunk)
        # `structural` is everything inside any block, kept for the
        # figure/prose share; `prose` excludes only figure and plate blocks.
        if depth:
            structural_out.append(chunk)
        if not fig_depth:
            prose_out.append(chunk)

    for a, b in spans(text):
        emit(text[pos:a])
        pos = b
        inner = text[a + 1:b - 1]
        if inner.startswith("?"):
            # `[?reading]` is an uncertain reading. Keep the transcriber's best
            # guess: it is what they actually saw on the page.
            emit(inner[1:].split(":", 1)[-1] if inner.startswith("?reading:") else inner[1:])
        elif inner in STRUCTURAL_OPEN:
            depth += 1
            if inner in FIGURE_TAGS:
                fig_depth += 1
        elif inner in STRUCTURAL_CLOSE:
            depth = max(0, depth - 1)
            if inner[1:] in FIGURE_TAGS:
                fig_depth = max(0, fig_depth - 1)
        elif any(inner.startswith(d) for d in DROP_WHOLE):
            m = _PAGE_MARKER.match(inner)
            if m and printed_page is None:
                printed_page = int(m.group(1))
                unnumbered = bool(m.group(2) and "unnumbered" in m.group(2))
                if unnumbered:
                    printed_page = None
        elif any(inner.startswith(k) for k in KEEP_AFTER_TAG):
            pass  # the tag is its own span; the text it labels is outside it
    emit(text[pos:])
    return _StrippedGold("".join(out), "".join(structural_out), "".join(prose_out),
                         printed_page, unnumbered)


# --- tokenisation -----------------------------------------------------------

# CJK ranges plus kana. These scripts are not written with spaces, so
# whitespace tokenisation would yield a handful of enormous tokens and the
# similarity would be meaningless; they are compared character by character.
CJK = re.compile(r"[぀-ヿ㐀-䶿一-鿿豈-﫿ｦ-ﾟ]")
CYRILLIC = re.compile(r"[Ѐ-ӿ]")
GREEK = re.compile(r"[Ͱ-Ͽ]")
NONLATIN = re.compile(
    r"[぀-ヿ㐀-䶿一-鿿豈-﫿ｦ-ﾟЀ-ӿͰ-Ͽ]"
)
# Vowels across the Latin, Cyrillic and Greek alphabets in this set. Used to
# tell a real word from OCR salad; an ASCII-only vowel class would classify
# every Russian and Greek word as noise.
_VOWELS = re.compile(r"[aeiouyаеиоуяыэюі"
                     r"αεηιουω]")


def tokens(text):
    """Normalise aggressively — we are measuring agreement on words, not on
    typography. OCR mangles case, spacing and punctuation constantly and those
    differences are not what we want to measure.

    Must be Unicode-aware: this set spans 13 languages, and an ASCII-only
    filter silently deletes the Cyrillic, Greek and CJK documents' actual
    content, leaving only stray Latin species names to compare. That produces a
    confident-looking score computed from a few percent of the page.
    """
    text = unicodedata.normalize("NFKD", text.lower())
    text = "".join(c for c in text if not unicodedata.combining(c))
    text = re.sub(r"[^\w\s]+", " ", text, flags=re.UNICODE)
    text = text.replace("_", " ")
    out = []
    for tok in text.split():
        if CJK.search(tok):
            out.extend(tok)
        else:
            out.append(tok)
    return out


def dominant_script(text):
    """Name the writing system a gold page is mostly set in.

    Measured from the page rather than declared per document, because the
    property is genuinely per page: at least one document in the set is an
    English translation followed by the vertical Japanese original, and a
    document-level label would mis-segment both halves.
    """
    counts = {
        "cjk": len(CJK.findall(text)),
        "cyrillic": len(CYRILLIC.findall(text)),
        "greek": len(GREEK.findall(text)),
    }
    latin = sum(1 for c in text if c.isalpha()) - sum(counts.values())
    counts["latin"] = max(0, latin)
    total = sum(counts.values())
    if not total:
        return "none"
    script, n = max(counts.items(), key=lambda kv: kv[1])
    # A page is only called non-Latin when that script genuinely carries it;
    # a Latin page quoting one Greek word is a Latin page.
    return script if n >= 0.15 * total else "latin"


def script_missing(gold_text, extracted_text):
    """True when the gold page is substantially non-Latin and the extraction
    contains almost none of that script.

    On the library side this meant "no signal, say nothing about the gold". On
    this side it is a finding: the page is in a script the extraction never
    represented, which is an OCR language-pack question (#176) rather than a
    text-quality one, and it is reported under its own status so it does not
    get averaged into pages that were merely extracted badly.
    """
    g = NONLATIN.findall(gold_text)
    if len(g) < 0.15 * max(len(gold_text.split()), 1):
        return False           # gold is not meaningfully non-Latin
    return len(NONLATIN.findall(extracted_text)) < 0.02 * len(g)


# --- corpuscle side ---------------------------------------------------------


def _item_text(item, kind):
    """The text of one docling item, whatever shape that item stores it in.

    Tables are the case that matters. A table item has no ``text`` field at
    all — its content lives in ``data.table_cells[].text`` — so a walk that
    reads only ``text`` silently discards every table on the page. Measured on
    this gold set before the fix, that cost 26% of one document's tokens and
    would have been scored as the pipeline losing text it had in fact
    extracted.
    """
    if kind == "tables":
        cells = ((item.get("data") or {}).get("table_cells") or [])
        # Row-major, matching how the table reads and how markdown export
        # linearises it. Cell order in the JSON is not guaranteed to be either.
        cells = sorted(cells, key=lambda c: (c.get("start_row_offset_idx") or 0,
                                             c.get("start_col_offset_idx") or 0))
        return " ".join((c.get("text") or "").strip() for c in cells).strip()
    return (item.get("text") or "").strip()


def page_texts_from_docling(doc):
    """Reconstruct per-page extracted text from a persisted ``docling_doc.json``.

    Walks ``body.children`` so the text is assembled in the same reading order
    ``export_to_markdown()`` uses for ``text.json``. Items without a ``prov``
    entry carry no page attribution and are dropped rather than guessed at.
    """
    resolvable = {
        "texts": doc.get("texts") or [],
        "tables": doc.get("tables") or [],
        "pictures": doc.get("pictures") or [],
    }
    groups = doc.get("groups") or []
    pages = defaultdict(list)
    seen = set()

    def resolve(ref):
        # Refs look like `#/texts/12`, `#/groups/0`, `#/body`.
        parts = ref.lstrip("#/").split("/")
        if len(parts) != 2:
            return None, None
        kind, idx = parts[0], int(parts[1])
        if kind == "groups":
            return "groups", (groups[idx] if idx < len(groups) else None)
        items = resolvable.get(kind)
        return kind, (items[idx] if items and idx < len(items) else None)

    def walk(node, depth=0):
        if node is None or depth > 64:
            # Depth guard: a malformed parent/child cycle would otherwise
            # recurse until the interpreter gives up.
            return
        for child in node.get("children") or []:
            ref = child.get("$ref") if isinstance(child, dict) else None
            if not ref or ref in seen:
                continue
            seen.add(ref)
            kind, item = resolve(ref)
            if item is None:
                continue
            if kind == "groups":
                walk(item, depth + 1)
                continue
            text = _item_text(item, kind)
            prov = item.get("prov") or []
            if text and prov and isinstance(prov[0], dict) and prov[0].get("page_no"):
                pages[int(prov[0]["page_no"])].append(text)
            walk(item, depth + 1)

    walk(doc.get("body") or {})

    # Fallback for a document whose body traversal reached nothing — an older
    # or truncated artifact. Flat list order is a worse reading order than
    # body order, but a worse order still scores better than no text at all,
    # and the alternative is reporting a total extraction failure that isn't.
    if not pages:
        for kind in ("texts", "tables"):
            for item in resolvable[kind]:
                text = _item_text(item, kind)
                prov = item.get("prov") or []
                if text and prov and isinstance(prov[0], dict) and prov[0].get("page_no"):
                    pages[int(prov[0]["page_no"])].append(text)

    return {n: "\n".join(chunks) for n, chunks in pages.items()}


def _read_json(path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def bind_documents(gold_root, corpuscle_root):
    """Match gold documents to corpuscle document directories on sha256.

    Returns ``(bound, unmatched)`` where ``bound`` is a list of
    ``(stem, sha256, gold_dir, corpus_dir)`` and ``unmatched`` names the gold
    documents absent from this corpuscle.
    """
    sources = _read_json(gold_root / "sources.json")
    if sources is None:
        raise SystemExit(f"no readable sources.json under {gold_root}")
    docs_dir = corpuscle_root / "documents"
    if not docs_dir.is_dir():
        raise SystemExit(f"no documents/ directory under {corpuscle_root}")
    by_prefix = {p.name: p for p in docs_dir.iterdir() if p.is_dir()}

    bound, unmatched = [], []
    for stem, entry in sorted(sources.items()):
        sha = (entry or {}).get("sha256") or ""
        gold_dir = gold_root / stem
        if not gold_dir.is_dir() or not any(gold_dir.glob("page_*.txt")):
            continue
        # Bind on the checksum prefix the corpuscle uses for its directory
        # names. Never fall back to the stem: two documents in this library
        # have shared a filename.
        corpus_dir = by_prefix.get(sha[:12])
        if corpus_dir is None:
            unmatched.append(stem)
            continue
        bound.append((stem, sha, gold_dir, corpus_dir))
    return bound, unmatched


def era_bucket(year):
    """Coarse era bands. Typography, not chronology, is what these separate:
    long-s and Fraktur before 1800, hot metal through the nineteenth century,
    photo-offset scans after, born-digital at the end.
    """
    if year is None:
        return "unknown"
    if year < 1800:
        return "pre-1800"
    if year < 1900:
        return "1800-1899"
    if year < 1950:
        return "1900-1949"
    if year < 2000:
        return "1950-1999"
    return "2000-"


# --- scoring ----------------------------------------------------------------


# A taxonomic name has to be at least this long to count. Shorter strings are
# abbreviations and rank words — `sp`, `var`, `f` — which collide with ordinary
# vocabulary and would inflate the rate with tokens nobody is retrieving on.
_MIN_TAXON_TOKEN = 4

# Below this many taxon tokens on a page, the rate is noise and is not
# reported. The gold set spans all of nature while a corpuscle's taxonomy is
# one clade: the 801-taxon siphonophore snapshot labels 58 tokens in the whole
# of *Systema Naturae*. A recall over three tokens is a coin flip printed to
# four decimals.
_MIN_TAXON_TOKENS_FOR_RATE = 10

# A page whose taxon coverage trails its prose coverage by more than this is
# worth looking at: the headline says the page came out well and the names it
# is retrieved by did not.
_TAXON_TRAILS_BY = 0.20


def taxon_vocabulary(corpuscle_root):
    """Every word of every name in the corpuscle's taxonomy, normalised.

    Word-level rather than name-level because the gold is compared token by
    token: a page printing *Agalmopsis elegans* contributes two tokens, and
    partial recovery of a binomial is a real, gradable outcome.

    Returns an empty set when there is no taxonomy — the rate is then omitted
    rather than reported as zero.
    """
    db = pathlib.Path(corpuscle_root) / "taxonomy.sqlite"
    if not db.exists():
        return frozenset()
    try:
        import sqlite3
        conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
        try:
            rows = conn.execute("SELECT name_lowercase FROM names").fetchall()
        finally:
            conn.close()
    except Exception as e:                       # pragma: no cover
        print(f"warning: could not read {db}: {e}", file=sys.stderr)
        return frozenset()
    vocab = set()
    for (name,) in rows:
        for word in tokens(name or ""):
            if len(word) >= _MIN_TAXON_TOKEN:
                vocab.add(word)
    return frozenset(vocab)


def score_page(gold_raw, extracted, taxon_vocab=frozenset()):
    """Score one page of extraction against one gold page.

    ``coverage`` is the headline: the fraction of the gold page's words that
    appear anywhere in the extraction. ``similarity`` is order-sensitive and
    will be depressed by reading-order differences alone — a multi-column page
    can score low on similarity with the text entirely correct, so the two are
    only meaningful read together:

        low similarity + high coverage -> reading-order difference only
        low similarity + low  coverage -> genuine content loss

    ``taxon_coverage`` answers a different question than ``coverage`` and is
    reported beside it rather than folded in. Coverage weights every token
    equally, and for this literature that undervalues exactly what retrieval
    keys on: replacing `por` with `eng` on a Portuguese paper costs 0.010 on
    English-wordlist tokens and 0.129 on everything else, which is where the
    binomials live. A change can look neutral in the headline while damaging
    the corpus's whole purpose. See dev_docs/OCR_LANGUAGES.md.
    """
    gold = strip_markup(gold_raw)
    gt = tokens(gold.text)
    et = tokens(extracted)
    st = tokens(gold.structural)
    pt = tokens(gold.prose)

    rec = {
        "printed_page": gold.printed_page,
        "printed_page_unnumbered": gold.unnumbered,
        "gold_words": len(gt),
        "extracted_words": len(et),
        "gold_structural_share": round(len(st) / len(gt), 4) if gt else None,
        "gold_prose_words": len(pt),
        "gold_figure_words": len(st),
        # The two coverages that matter separately. Body text is the pipeline's
        # job; lettering engraved inside a plate is a different problem with a
        # different value, and averaging them produces a number that answers
        # neither question. Both use the gold side's own split — the extraction
        # is a single blob per page, but that is fine for coverage, which only
        # asks whether each true word turned up somewhere.
        "prose_coverage": None,
        "figure_coverage": None,
        # Reported with its denominator, never as a bare rate: a corpuscle's
        # taxonomy covers one clade and the gold spans all of nature, so most
        # pages carry too few taxon tokens for a rate to mean anything.
        "taxon_tokens": 0,
        "taxon_coverage": None,
        "script": dominant_script(gold.text),
        "similarity": None,
        "coverage": None,
        "recall": None,
        "excess_wordlike": None,
        "excess_novel": None,
        "volume_ratio": round(len(gt) / len(et), 3) if et else None,
        "uncertain": len(re.findall(r"\[\?", gold_raw)),
        "illegible": len(re.findall(r"\[illegible", gold_raw)),
        # A gold-integrity signal, reported rather than absorbed. A bracket that
        # opens nothing is either a transcription typo or a note mentioning the
        # character; only inspection tells them apart, so the count is surfaced
        # instead of being silently tolerated.
        "unterminated_brackets": unterminated_brackets(gold_raw),
        "status": "ok",
    }

    if not gt:
        # A genuinely blank leaf, or a plate with no lettering at all. There is
        # nothing to recover, so scoring it would reward the extractor for the
        # page being empty. Counted, not scored.
        rec["status"] = "gold_blank"
        return rec

    if not et or len("".join(et)) < MIN_EXTRACTED_CHARS:
        # The extraction is empty on a page the gold says carries text. This is
        # a total failure and it is scored as one — coverage 0.0, not None.
        # Leaving it unscored would quietly delete the pipeline's worst pages
        # from the average.
        rec["status"] = "extraction_empty"
        rec["coverage"] = 0.0
        rec["similarity"] = 0.0
        if pt:
            rec["prose_coverage"] = 0.0
        if st:
            rec["figure_coverage"] = 0.0
        if taxon_vocab:
            tx = [w for w in gt if w in taxon_vocab]
            rec["taxon_tokens"] = len(tx)
            if len(tx) >= _MIN_TAXON_TOKENS_FOR_RATE:
                # A page the pipeline lost entirely lost its taxa with it.
                rec["taxon_coverage"] = 0.0
        return rec

    if script_missing(gold.text, extracted):
        rec["status"] = "script_missing"

    rec["similarity"] = round(SequenceMatcher(None, gt, et, autojunk=False).ratio(), 4)

    gc, ec = Counter(gt), Counter(et)
    inter = sum((gc & ec).values())
    # THE HEADLINE. How much of what is actually printed on the page did the
    # pipeline recover? Order-insensitive on purpose, so it is not contaminated
    # by column interleaving the way `similarity` is.
    rec["coverage"] = round(inter / len(gt), 4)
    if pt:
        rec["prose_coverage"] = round(sum((Counter(pt) & ec).values()) / len(pt), 4)
    if taxon_vocab:
        tx = [w for w in gt if w in taxon_vocab]
        rec["taxon_tokens"] = len(tx)
        if len(tx) >= _MIN_TAXON_TOKENS_FOR_RATE:
            rec["taxon_coverage"] = round(
                sum((Counter(tx) & ec).values()) / len(tx), 4)
    if st:
        rec["figure_coverage"] = round(sum((Counter(st) & ec).values()) / len(st), 4)
    # The opposite question: how much of what the pipeline emitted is really on
    # the page? Low recall means the extraction carries material the page does
    # not — either text bled in from elsewhere, or invented from image texture.
    rec["recall"] = round(inter / len(et), 4)

    # Low recall has two utterly different causes and conflating them wastes
    # the finding. Word-likeness separates them cheaply: real words are
    # alphabetic, of reasonable length, and contain a vowel. Character salad
    # extracted from a marbled endpaper is not.
    excess = list((ec - gc).elements())
    real = [w for w in excess if len(w) >= 3 and w.isalpha() and _VOWELS.search(w)]
    rec["excess_wordlike"] = round(len(real) / len(excess), 4) if excess else 0.0
    # Final discriminator. A word-like excess token is still not evidence that
    # the pipeline read something new if it is merely a MANGLED FORM of a word
    # the gold already has. Long-s typography guarantees this on eighteenth-
    # century pages, where an extractor reads "ſ" as "f" and nearly every word
    # mismatches while the content is right. `excess_novel` counts only excess
    # words with no near-match in the gold page: genuinely spurious text looks
    # like whole unfamiliar words, OCR damage looks like near-misses.
    goldset = set(gt)
    novel = [w for w in real if not get_close_matches(w, goldset, n=1, cutoff=0.75)]
    rec["excess_novel"] = round(len(novel) / len(real), 4) if real else 0.0
    return rec


def _median(values):
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return None
    mid = len(vals) // 2
    if len(vals) % 2:
        return round(vals[mid], 4)
    return round((vals[mid - 1] + vals[mid]) / 2, 4)

def _mean(values):
    vals = [v for v in values if v is not None]
    return round(sum(vals) / len(vals), 4) if vals else None


def _percentile(values, q):
    """Nearest-rank percentile, for measures whose median is saturated."""
    vals = sorted(v for v in values if v is not None)
    if not vals:
        return None
    return round(vals[min(len(vals) - 1, int(q * len(vals)))], 4)



def _aggregate(pages):
    """Summarise a set of page records.

    Median rather than mean: this set was deliberately built around its hardest
    cases, so a handful of unreadable plates would drag a mean somewhere no
    individual page is. The status counts beside it are what keep the median
    honest — a good median over few scored pages is not a good result.
    """
    scored = [p for p in pages if p["coverage"] is not None]
    statuses = Counter(p["status"] for p in pages)
    return {
        "pages": len(pages),
        "scored": len(scored),
        "median_coverage": _median([p["coverage"] for p in scored]),
        "median_prose_coverage": _median([p.get("prose_coverage") for p in scored]),
        # Mean and a low percentile, deliberately NOT the median used
        # everywhere else. 53% of qualifying pages recover every taxon token,
        # so the median sits at exactly 1.000 and cannot move — a metric
        # pinned at its ceiling cannot detect the regressions it exists for.
        # The tail is where the signal is: p10 = 0.643 on the same set.
        "mean_taxon_coverage": _mean([p.get("taxon_coverage") for p in scored]),
        "p10_taxon_coverage": _percentile(
            [p.get("taxon_coverage") for p in scored], 0.10),
        "taxon_scored_pages": sum(1 for p in scored
                                  if p.get("taxon_coverage") is not None),
        # The actionable count: pages the headline calls fine while the tokens
        # retrieval depends on are being lost. This is the shape #244 was
        # filed about, and it is real — Hosiaetal2024 p6 scores prose 0.976
        # and taxon 0.140.
        "taxon_trails_prose": sum(
            1 for p in scored
            if p.get("taxon_coverage") is not None
            and p.get("prose_coverage") is not None
            and p["prose_coverage"] - p["taxon_coverage"] > _TAXON_TRAILS_BY),
        "median_figure_coverage": _median([p.get("figure_coverage") for p in scored]),
        "prose_pages": sum(1 for p in scored if p.get("prose_coverage") is not None),
        "figure_pages": sum(1 for p in scored if p.get("figure_coverage") is not None),
        "median_recall": _median([p["recall"] for p in scored]),
        "median_similarity": _median([p["similarity"] for p in scored]),
        "pages_below_0.5_coverage": sum(1 for p in scored if p["coverage"] < 0.5),
        "status_counts": dict(sorted(statuses.items())),
        "gold_words": sum(p["gold_words"] for p in pages),
        "taxon_tokens": sum(p.get("taxon_tokens", 0) for p in pages),
        "extracted_words": sum(p["extracted_words"] for p in pages),
        "uncertain": sum(p["uncertain"] for p in pages),
        "illegible": sum(p["illegible"] for p in pages),
        "unterminated_brackets": sum(p.get("unterminated_brackets", 0) for p in pages),
    }


def keeppages_map(scan):
    """``{source_page: subset_page}`` for a document, or ``{}`` if unselected.

    The gold set was transcribed over the whole file — `Beklemishev1969`'s
    gold page 1 is an ownership endpaper, `Kawamura1911a`'s is twelve pages
    of English typescript — while a `keeppages` document's extracted page
    numbers are positions in the subset (#188). Every scorer here binds gold
    to extraction by page number, so without this the two coordinate systems
    silently disagree and a selected document scores near zero: a numbering
    artifact that looks exactly like a catastrophic regression.

    `scan_detection.json`'s `keeppages_selected` is the forward map (subset
    page i is source page ``selected[i-1]``); this inverts it.
    """
    selected = (scan or {}).get("keeppages_selected") or []
    return {src: i + 1 for i, src in enumerate(selected)}


def rebase_gold_pages(gold_by_page, scan):
    """Restate a ``{source_page: X}`` gold mapping in subset coordinates.

    Pages the operator excluded are dropped rather than scored as misses:
    corpus was told they are not the paper, so counting a deliberately
    removed library title page as a figure it failed to find would penalise
    the selection for working. Returns the input unchanged when no selection
    is active, which is the ordinary case.
    """
    mapping = keeppages_map(scan)
    if not mapping:
        return gold_by_page
    return {mapping[p]: v for p, v in gold_by_page.items() if p in mapping}


def score_document(gold_dir, corpus_dir, taxon_vocab=frozenset()):
    """Score every gold page of one document against the corpuscle's extraction."""
    doc_json = _read_json(corpus_dir / "docling_doc.json")
    extracted = page_texts_from_docling(doc_json) if doc_json else {}
    meta = _read_json(corpus_dir / "metadata.json") or {}
    scan = _read_json(corpus_dir / "scan_detection.json") or {}

    # #188 — a `keeppages` selection makes the extracted page numbers
    # positions in a *subset*, while the gold set was transcribed over the
    # whole file: `Beklemishev1969`'s gold page 1 is an ownership endpaper
    # the operator has since declared not part of the paper. Binding by
    # position would shift every page by one and score the document at
    # roughly zero — a numbering artifact indistinguishable from a
    # catastrophic regression.
    #
    # `keeppages_selected` is the map: extracted page i is source page
    # selected[i-1]. Invert it to look the other way.
    selected = scan.get("keeppages_selected") or []
    gold_to_extracted = keeppages_map(scan)

    records, excluded = [], []
    for gf in sorted(gold_dir.glob("page_*.txt")):
        n = int(gf.stem.split("_")[1])
        if selected:
            if n not in gold_to_extracted:
                # The operator said this page is not the paper. Scoring it as
                # a miss would penalise the selection for doing its job — but
                # it is counted and reported, never silently dropped, or the
                # harness looks like it scored pages it never opened.
                excluded.append(n)
                continue
            extracted_page = gold_to_extracted[n]
        else:
            extracted_page = n
        rec = score_page(gf.read_text(encoding="utf-8"),
                         extracted.get(extracted_page, ""), taxon_vocab)
        rec["page"] = n
        if selected:
            rec["subset_page"] = extracted_page
        records.append(rec)

    year = meta.get("year")
    summary = _aggregate(records)
    summary.update({
        "year": year,
        "era": era_bucket(year),
        "file_type": scan.get("file_type"),
        "needs_ocr": scan.get("needs_ocr"),
        "tesseract_packs": scan.get("tesseract_packs"),
        "detected_language": scan.get("detected_language"),
        # Documents mix scripts page by page, so the document-level label is a
        # tally rather than a single name.
        "scripts": dict(sorted(Counter(r["script"] for r in records).items())),
        "docling_pages_seen": len(extracted),
    })
    if selected:
        summary["keeppages"] = scan.get("keeppages")
        summary["pages_excluded_by_selection"] = excluded
        summary["gold_pages_scored"] = len(records)
    return summary, records


def segment(all_pages):
    """Group page records by the axes the gold set was built to span.

    Reported separately because a corpus-wide mean over 13 languages and five
    centuries is not actionable: it cannot distinguish a pipeline that handles
    born-digital PDFs well and Fraktur not at all from one that is mediocre
    everywhere, and those two call for completely different work.
    """
    out = {}
    for axis in ("script", "era", "file_type"):
        buckets = defaultdict(list)
        for p in all_pages:
            buckets[str(p.get(axis))].append(p)
        out[axis] = {k: _aggregate(v) for k, v in sorted(buckets.items())}
    return out


def _pdf_name(summary, stem):
    """The PDF basename for a scored document, for bib lookup."""
    return (summary.get("filename") or f"{stem}.pdf").lower()


def _bib_keeppages(bib_path):
    """``{pdf_basename_lower: keeppages}`` from a .bib, or ``{}``.

    Read to *cross-check*, never to shift page numbers. The shift uses
    `keeppages_selected` from `scan_detection.json`, which is the directive
    already resolved against the document's real page count — clamped,
    deduplicated, ordered. Re-parsing the raw string here would mean a second
    implementation of that resolution, and the two would disagree the first
    time either changed. Same argument as #215's mapping table.
    """
    if not bib_path:
        return {}
    import re
    try:
        text = pathlib.Path(bib_path).expanduser().read_text(
            encoding="utf-8", errors="replace")
    except OSError as e:
        print(f"warning: could not read {bib_path}: {e}", file=sys.stderr)
        return {}
    out = {}
    for entry in re.split(r"\n(?=@)", text):
        f = re.search(r"\bfile\s*=\s*\{([^}]*)\}", entry)
        k = re.search(r"\bkeeppages\s*=\s*\{([^}]*)\}", entry)
        if f and k:
            out[pathlib.Path(f.group(1).strip()).name.lower()] = k.group(1).strip()
    return out


def build_report(gold_root, corpuscle_root, bib_path=None):
    bound, unmatched = bind_documents(gold_root, corpuscle_root)
    bib_keeppages = _bib_keeppages(bib_path)
    taxon_vocab = taxon_vocabulary(corpuscle_root)
    documents, all_pages, detail, stale = {}, [], [], []
    for stem, sha, gold_dir, corpus_dir in bound:
        summary, records = score_document(gold_dir, corpus_dir, taxon_vocab)
        # Cross-check the directive against what the build actually applied.
        # The dangerous case is scoring a corpuscle built *before* the
        # directive was written: nothing in its artifacts says a selection
        # was intended, the page numbers line up 1:1 with the gold, and the
        # scores look perfectly reasonable — they are simply answering a
        # different question than the operator now thinks they are.
        declared = bib_keeppages.get(_pdf_name(summary, stem))
        if declared and not summary.get("keeppages"):
            stale.append({"document": stem, "bib_keeppages": declared})
            summary["keeppages_not_applied"] = declared
        documents[stem] = summary
        for r in records:
            # Carry the document's axes onto each page so segmentation is a
            # flat group-by rather than a nested walk.
            r = dict(r, era=summary["era"], file_type=summary["file_type"], document=stem)
            all_pages.append(r)
        detail.append({"document": stem, "sha256": sha, "pages": records})

    return {
        "schema_version": 1,
        "method": "token-level SequenceMatcher plus order-insensitive coverage and "
                  "recall; gold markup stripped, case/diacritic/punctuation "
                  "normalised, CJK compared per character",
        "subject": "the corpuscle's extraction, measured against gold as the "
                   "yardstick; coverage (how much of the true page was recovered) "
                   "is the primary measure, and a page the pipeline failed to "
                   "extract at all is scored 0.0 rather than excluded",
        "caveat": "Cannot detect invention, and a low volume ratio is not evidence "
                  "of loss. Read the gold set's CROSSCHECK_REPORT.md before acting "
                  "on any number, noting that it puts the transcription on trial "
                  "rather than the extractor: same arithmetic, different subject, "
                  "so several of its false positives are findings here.",
        "gold_root": str(gold_root),
        "corpuscle_root": str(corpuscle_root),
        # Non-empty means the corpuscle predates a keeppages directive now in
        # the bib. The scores are still internally consistent, they just
        # describe the untrimmed documents.
        "keeppages_declared_but_not_applied": stale,
        "documents_bound": len(bound),
        "documents_unmatched": unmatched,
        "corpus_wide": _aggregate(all_pages),
        "segments": segment(all_pages),
        "documents": documents,
        "detail": detail,
    }


# --- reporting --------------------------------------------------------------


def _fmt(value, width=7):
    return f"{value:>{width}.3f}" if isinstance(value, float) else f"{'-':>{width}}"


def print_summary(report, stream=None):
    # Resolved here, not as a default argument: a default binds
    # `sys.stdout` at definition time, so the function would ignore
    # any later redirection — including a test capturing it.
    w = (stream or sys.stdout).write
    w(f"\ngold {report['gold_root']}\n")
    w(f"corpuscle {report['corpuscle_root']}\n")
    w(f"{report['documents_bound']} documents bound on sha256")
    if report["documents_unmatched"]:
        w(f"; {len(report['documents_unmatched'])} not in this corpuscle: "
          f"{', '.join(report['documents_unmatched'])}")
    w("\n")

    stale = report.get("keeppages_declared_but_not_applied") or []
    if stale:
        w(f"\n!! {len(stale)} document(s) carry a keeppages directive this "
          f"build did not apply.\n")
        w("   The scores below describe the untrimmed documents — the gold "
          "pages and the\n   extraction still line up 1:1, so nothing is "
          "misaligned, but the front matter\n   the operator excluded is "
          "still being scored. Rebuild to measure what was asked for.\n")
        for e in stale[:10]:
            w(f"     {e['document']:28s} keeppages={e['bib_keeppages']}\n")
        if len(stale) > 10:
            w(f"     ... and {len(stale) - 10} more\n")

    selected_docs = [d for d, v in (report.get("documents") or {}).items()
                     if v.get("keeppages")]
    if selected_docs:
        dropped = sum(len(report["documents"][d].get("pages_excluded_by_selection") or [])
                      for d in selected_docs)
        w(f"\nkeeppages active on {len(selected_docs)} document(s); "
          f"{dropped} gold page(s) excluded from scoring as not-the-paper.\n"
          "  Scores are therefore over a different, smaller page set than an "
          "untrimmed build's.\n")

    cw = report["corpus_wide"]
    w(f"\n{cw['pages']} gold pages, {cw['scored']} scored, "
      f"{cw['pages'] - cw['scored']} not comparable\n")
    w(f"median coverage: PROSE {_fmt(cw['median_prose_coverage'])} over "
      f"{cw['prose_pages']} pages | figure text {_fmt(cw['median_figure_coverage'])} "
      f"over {cw['figure_pages']} | combined {_fmt(cw['median_coverage'])}\n")
    w("prose is the measure that matters; figure text is reported, not optimised. "
      "Read the segments, not these.\n")
    if cw.get("mean_taxon_coverage") is not None:
        w(f"TAXON coverage: mean {_fmt(cw['mean_taxon_coverage'])} "
          f"p10 {_fmt(cw['p10_taxon_coverage'])} over "
          f"{cw['taxon_scored_pages']} page(s) with >= "
          f"{_MIN_TAXON_TOKENS_FOR_RATE} taxon tokens "
          f"({cw['taxon_tokens']} tokens corpus-wide)\n")
        w("  mean and p10, not median: half these pages recover every taxon "
          "token, so a median\n  sits at 1.000 and cannot move. A separate "
          "question from coverage, not a component.\n")
        if cw.get("taxon_trails_prose"):
            w(f"  !! {cw['taxon_trails_prose']} page(s) score >"
              f"{_TAXON_TRAILS_BY} better on prose than on taxa — the "
              f"headline calls them fine\n     while the names they are "
              f"retrieved by are being lost.\n")
    elif cw.get("taxon_tokens"):
        w(f"taxon tokens seen: {cw['taxon_tokens']}, but no page carries "
          f"{_MIN_TAXON_TOKENS_FOR_RATE} — no rate reported\n")
    w(f"status: {cw['status_counts']}\n")
    if cw.get("unterminated_brackets"):
        offenders = sorted(
            (stem for stem, agg in report["documents"].items()
             if agg.get("unterminated_brackets")))
        w(f"gold integrity: {cw['unterminated_brackets']} bracket(s) opening no "
          f"marker, in {', '.join(offenders)} — inspect before trusting those "
          f"pages\n")

    for axis, buckets in report["segments"].items():
        w(f"\n-- by {axis} " + "-" * (58 - len(axis)) + "\n")
        w(f"{'':22}{'pages':>7}{'PROSE':>8}{'figure':>8}{'recall':>8}{'all':>8}\n")
        for name, agg in buckets.items():
            w(f"{name[:22]:22}{agg['pages']:>7}"
              f"{_fmt(agg['median_prose_coverage'], 8)}"
              f"{_fmt(agg['median_figure_coverage'], 8)}"
              f"{_fmt(agg['median_recall'], 8)}{_fmt(agg['median_coverage'], 8)}\n")

    w("\n-- worst documents by PROSE coverage " + "-" * 35 + "\n")
    ranked = sorted(
        report["documents"].items(),
        key=lambda kv: (kv[1]["median_prose_coverage"] is None,
                        kv[1]["median_prose_coverage"] or 0),
    )
    w(f"{'':26}{'year':>6}{'pages':>7}{'PROSE':>8}{'figure':>8}{'recall':>8}\n")
    for stem, agg in ranked[:15]:
        w(f"{stem[:26]:26}{str(agg['year'] or '-'):>6}{agg['pages']:>7}"
          f"{_fmt(agg['median_prose_coverage'], 8)}"
          f"{_fmt(agg['median_figure_coverage'], 8)}"
          f"{_fmt(agg['median_recall'], 8)}\n")
    w("\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gold", required=True,
                    help="the gold transcriptions directory (holds sources.json "
                         "and one directory of page_NNN.txt per document)")
    ap.add_argument("--corpuscle", required=True,
                    help="a built corpuscle; either its root or the output/ "
                         "directory that holds documents/")
    ap.add_argument("--bib", default=None,
                    help="optional .bib; cross-checks that a keeppages "
                         "directive was actually applied by the build being "
                         "scored, and warns when it was not")
    ap.add_argument("--out", default=None, help="write the full JSON report here")
    ap.add_argument("--quiet", action="store_true", help="suppress the stdout summary")
    a = ap.parse_args(argv)

    gold_root = pathlib.Path(a.gold).expanduser().resolve()
    corpuscle_root = pathlib.Path(a.corpuscle).expanduser().resolve()
    # Accept either the corpuscle root or its output/ directory, because both
    # are things an operator reasonably has on the clipboard.
    if not (corpuscle_root / "documents").is_dir() and (corpuscle_root / "output" / "documents").is_dir():
        corpuscle_root = corpuscle_root / "output"

    report = build_report(gold_root, corpuscle_root, bib_path=a.bib)
    if a.out:
        out = pathlib.Path(a.out).expanduser()
        out.write_text(json.dumps(report, indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"wrote {out}")
    if not a.quiet:
        print_summary(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
