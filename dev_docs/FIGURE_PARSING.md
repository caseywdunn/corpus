# Figure parsing: what the gold set says

Companion to [OCR_LANGUAGES.md](OCR_LANGUAGES.md). Everything here is measured
against the 35-document gold transcription set — 761 pages, 1594–2026, whose
pages were transcribed from rendered images alone, so the figures counted are
the ones actually printed rather than the ones an extractor found.

Three separate questions, deliberately not averaged into one:

1. **Are all the figures found, and is furniture being called a figure?**
   `tools/qc/figure_detection.py`
2. **Is each caption bound to the figure it belongs to?**
   `tools/qc/caption_binding.py`
3. **Is the text inside a figure recovered?** Reported by
   `tools/qc/fidelity.py` as `figure_coverage`, and explicitly *not optimised* —
   lettering engraved on a plate is a different problem with a different value.

## Guidance

**The default served surface favors precision; the raw record is noisier than
it looks.** Every MCP tool filters on the shared evidence-type set (`figure`,
`plate`, `subpanel`), so neither `graphical_element` nor the `unclassified`
review bucket reaches a client unless it explicitly asks for `include_all`.
Read the row for the surface you mean: excluding only `graphical_element` is a
useful classifier diagnostic, but it is not the default MCP surface.

**If you are building a corpus of modern papers, expect furniture.** Precision
on born-digital documents is 0.586 raw against 0.946 for scans, because
publisher logos, ORCID icons and journal rules are figures as far as a layout
model is concerned. The `graphical_element` filter handles it; historical scans
barely need it.

**Do not read a per-document rate over a handful of figures.** Half the
documents here carry fewer than ten gold figures, and two carry one. The
corpus-wide numbers are meaningful; `Chun1882c` at precision 0.25 is three
surplus records on a one-figure document.

**Front matter costs figures twice.** In the classifier diagnostic that drops
only `graphical_element`, trimming it with `keeppages` improved both recall
and precision, because a plate atlas's duplicate captures and a bound volume's
title page are surplus figures that no classifier can recognise as
not-the-paper. Page selection is therefore part of the measurement, and the
scorer rebases the gold into selected-page coordinates before counting.

## Where it stands

Clean release candidate, 35 documents, 376 gold figure blocks (316
`[FIGURE]`, 60 `[PLATE]`):

| filter | recall | precision | F1 |
| --- | --- | --- | --- |
| all physical detections | 0.883 | 0.865 | 0.874 |
| drop `graphical_element` | 0.867 | 0.985 | 0.922 |
| drop uncaptioned `graphical_element` | 0.867 | 0.973 | 0.917 |
| **shared evidence types** (default MCP surface) | **0.827** | **1.000** | **0.905** |
| captioned only | 0.827 | 0.987 | 0.900 |

Scorer v2 makes a distinction that #195 made load-bearing: a `figures.json`
entry is not necessarily another physical detection. The clean build contains
653 entries, of which 261 are logical records sharing a plate/compound image;
after excluding those and collapsing eight split-panel groups, 384 physical
figures remain. The preceding clean artifact had 422 entries and only 30
image-sharing logical records, but the same 384 physical figures and exactly
the same score in every filter row. Caption expansion therefore did not
regress physical detection. The earlier 0.936 / 0.876 headline counted logical
children as newly detected images and is invalid.

Caption binding is scored on typed page/identity pairs over the default MCP
evidence types: `plate:10` and `figure:10` are distinct even when printed on
the same page. The v7 scorer corrected material gold-denominator bugs: 25
plate blocks explicitly inventory their engraved numbers, and many more list
those numbers as standalone `1`, `F. 1.`, etc. The earlier parser counted only
lines beginning `Fig.` and collapsed same-numbered plate and child identities,
therefore reporting **480 gold pairs where the transcription actually contains
839**. It also treated a correct `PLATE N` host as false when the transcriber
placed that heading immediately outside `[PLATE]`. The earlier 0.596 / 0.911
headline is not comparable and must not be used as a release baseline.

The complete clean source-PDF build reports **538 correct of 545 identities
against 839 gold: recall 0.641, precision 0.987**, with fixed-population
capacity 544/839 (0.648). Totton, whose complete legends are printed on the
leaf *before* each plate, moves from the pre-repair 181/184 correct reported
identities to **406/411** (472 gold). Exact `PLATE N`-to-next-page-`PLATE N`
matching creates logical children with their own source captions while
preserving Plate X and Figure 10 as distinct identities. The clean build also
persists the unconditioned-vision probe's nine correct labels and no false
labels on the four plates left eligible after deterministic expansion.

The seven surplus identities were inspected individually. Three already carry
explicit low-confidence `uncertain` status. Two are faithful bindings to OCR
that damaged the printed number (`1` → `7`, `34` → `3`) and cannot be safely
rewritten without independent number evidence. One is a gold-markup omission:
the transcription note and printed caption both identify Figure 45, but only
Figures 46 and 47 have structural blocks. The final case exposed a real typed
identity defect: a following-page Figure 16 legend could overwrite a
same-page Plate XVI link. The selector now forbids that cross-namespace
replacement; its post-fix clean measurement remains the final #195 gate.

Panel splitting is a separate acceptance measure rather than an inference
from figure-number binding. The independent gold parser finds **98 captions
that explicitly enumerate lettered panels**. The clean build reports a panel
declaration for **92**, an exact label set for **89**, and label recall /
precision of **0.946 / 0.997**. It declares no panels for the 175 number-matched
members of the 202 gold identities classified as non-panelled.

The scorer parses hand-transcribed gold labels independently of the production
parser and reports the raw artifact and default MCP surface separately. That
separation is load-bearing: an earlier version reused production's fuzzy OCR
parser, so changing the extractor silently changed the supposed gold
denominator. Scorer v7 additionally understands the transcription's explicit
plate inventories, standalone engraved numbers, `F. N` labels, and adjacent
plate headings, and keeps plate and child-figure identities distinct. Those
are gold-format rules, not production OCR heuristics.
The compact anatomical key `Pl.M.` remains explicitly excluded: inside these
figures it means mouth-plate, not Roman-numeral Plate 1000.

**Panel labels have their own independent parser and denominator.** Only text
after the gold block's first figure-caption opener is eligible, so an `A` or
`B` engraved inside the picture cannot manufacture a caption declaration.
The yardstick accepts period, parenthesized, comma and range styles, retains
printed gaps such as Totton Figure 74's A–H, K, L, and excludes numeric grouped
plates. It reports declaration recall/precision, exact label sets, and
label-level recall/precision. Production supports those same letter styles,
but uses geometry and a growing label set to join Docling-split continuation
cells; it stops period-marker scanning at abbreviation glossaries so
`C. rad.lat` cannot become a fictitious panel C.

## How it is measured, and why in that shape

**Detection is counted per page, not matched figure by figure.** The gold
records what is printed; it carries no bounding boxes, so there is no
correspondence to match against. Per page, `matched = min(gold, found)`,
`missed = gold - matched`, `surplus = found - matched`. That is the strongest
claim the data supports, and it is why a document can score 1.0 while having
the right *number* of wrong figures.

**Captions bind on the figure number, never on caption text.** A naive
text-similarity match reports 44% bound and is mostly artifact:
`Chenetal2015` scores 0 of 10 because the document prints every caption twice,
Chinese and English, and the extraction took the Chinese while the matcher was
handed the English — every figure is in fact bound correctly.
`Carre1969_Nanomia_tr` scores 1 of 18 because its plate pages carry only
`FIG. 1` as printed matter, two tokens no threshold can bind, with the prose
caption on a different page. Numbers are language-independent, and #16's
Roman↔Arabic normalisation makes them comparable across centuries.

**Gold blocks are classified before scoring, so no rate is computed over
blocks it does not apply to:**

| kind | count | can it test caption text? |
| --- | --- | --- |
| `bare_label` | 254 | no — "Fig. 3" and nothing else |
| `prose_caption` | 105 | yes |
| `lettering_only` | 12 | no — non-number lettering engraved on the plate |
| `nothing_printed` | 5 | no |

**Most "captions" in this literature are not captions.** 254 of 376 blocks are
a bare label. Caption *text* similarity is computable for only 92
number-matched pairs in the current clean candidate, which is why it is
reported and not headlined.

## Where corpus does well

**Mid-century scanned monographs, which is most of this material.**
1950–1999: recall 0.975, precision 0.960. `Totton1965a` — 226 pages, 195 gold
figures — scores recall 0.974 and precision 0.995.

**Scans generally**: 0.870 / 0.946 raw, against born-digital at 0.962 / 0.586.
That inversion is worth internalising. The modern papers are easier to *read*
and harder to *count*.

**Figure numbers, after #205.** `parse_figure_number` fired for only 32% of
figures; the cause was pleasing once found — the opener word is OCR-damaged,
`Fic.` and `Frc.` for `Fig.`, and across 320 captions the damaged spellings
were **more common than the correct one**. Coverage is now 67.6% at 98.2%
precision.

## Where it has trouble

**Publisher furniture on modern papers.** `Ahuja_etal2026`: 6 gold figures,
31 physical records. Logos, icons and rules. The `graphical_element` filter
removes them from the served surface, but the raw record still carries them,
and `figures_report.html` shows them.

**Historical plates carrying several engravings under one legend.** On the
pre-repair served bundle, `Vanhoeffen1906` scores 0.781 recall / 0.943
precision. The clean v1.3 candidate scores **0.891 / 0.966**: enumerating
caption blocks are split into per-number
entries, duplicate picture assignments are reconciled when the counts form a
complete bijection, and collected prose on the following page enriches
already-found bare labels instead of being cloned onto the wrong plate.
Pass 2.5 now also admits the shared plate to ROI detection with numeric figure
targets kept separate from panel letters. Pass 3 runs once on the plate and
puts each detected region on its logical figure record; the default OCR path
uses an exact caption-derived number allow-list, while vision consumes the
same targets. The caption-binding rate measures ownership and numbering, not
crop geometry; ROI accuracy requires a separate geometric yardstick.

The other two conspicuous acceptance failures move in the same replay:
`Hosiaetal2024` goes from 0.692 / 0.900 to **0.923 / 1.000**, and
`Ahuja_etal2026` from 1.000 / 0.500 to **1.000 / 1.000**; both gains are
confirmed in the clean build.

**Caption binding in 1800–1899: recall 0.030** (5 of 169 numbers) in the clean
candidate. The former 5-of-33 denominator omitted
standalone engraved labels. The numbers are printed on the plates as
lettering, not as extracted text, so there is often nothing for a text-only
parser to find. Four documents score zero:
`Bernstein1934`, `Chenetal2015`, `Eschscholtz1825`, and `Tilesius1814`.
Pass 3b now has a deliberately narrow image-evidence path for this layout:
only a confidently bound bare plate is admitted, and it materializes nothing
unless the VLM returns at least two distinct Arabic number+region candidates
at confidence ≥0.80. The host retains every accepted/rejected decision;
derived figure records share the plate image and remain caption-unbound. The
Totton probe validates the mechanism, but the pre-1900 plates still require a
clean run and inspection; do not generalize one monograph's result across
those layouts.

## Furniture: `graphical_element` is an actionable predicate

Established by scoring the corpuscle under each candidate filter rather than
by assuming which field name sounded right. Dropping `graphical_element` takes
precision from 0.865 to 0.985 for 0.016 of recall. The default MCP surface also
hides `unclassified`: precision reaches 1.000, but recall falls to 0.827. That
is an intentional review-state policy, not evidence that those records are all
furniture; the raw artifact and `include_all` retain them for inspection.

Position recurrence is the other signal, and it separates cleanly:
`furniture_positions()` treats a bbox recurring on ≥5 pages as running
furniture. **No real figure in the 35 documents repeats a position even
twice**; the journal logo in one born-digital paper repeats on 24 of 25 pages.
That evidence is what lets `classify_figure` trust a caption on a small item
without also trusting one the caption-proximity heuristic wrongly attached to
a masthead.

## A lesson from #203 and #231

#203 gave each plate-legend entry its own record and shipped on a measured
gain: corpus-wide figure recall 0.849 → 0.917. Scoring the rebuild showed it
also cost precision, 0.970 → 0.892 on the served surface, and **all of the
loss was one document**. The legend scan took a figure number from anywhere in
a text item, so in a 226-page monograph of running prose it read
cross-references — a species heading `Plate XX, figures 1, 2` — as legend
entries and cloned the page's real figure under the referenced number. 33 of
that document's 232 records shared an image file with another record under a
different number: asking for figure 20 returned the picture of figure 53.

Two things worth carrying forward. The harness reports recall *and* precision,
and the summary quoted recall. And the defect was found by opening the page
image after three rounds of JSON analysis had pointed the wrong way.

## Known gaps

- **ROI detection is not extraction-time physical splitting.** A grouped
  plate's logical records point at one source image, with per-record ROI crops
  available only when Pass 3 locates their printed numbers. A multi-panel
  modern figure likewise remains one record with optional panel ROIs. #207
  found the whole-figure/subpanel branch of `dedupe_figures` unreachable — its
  overlap measure is symmetric and the looser threshold always fires first —
  which is corroborated by no figure in the reference corpuscle carrying
  `figure_type: subpanel`.
- **Counting figures measures the definition of "figure" as much as the
  extraction.** `Ahuja_etal2026` has 6 gold blocks against 31 physical records
  because Docling counts logos; `Vanhoeffen1906` has 67 against 33 because the
  gold counts engravings and Docling counts plates. Caption binding is the
  better posed question, and it is the one with the lower score.
- **The metrics weight every figure equally**, which — as with prose coverage
  and taxonomic tokens (#244) — is not how they are valued downstream.
