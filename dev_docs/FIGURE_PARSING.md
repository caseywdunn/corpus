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
on born-digital documents is 0.622 raw against 0.941 for scans, because
publisher logos, ORCID icons and journal rules are figures as far as a layout
model is concerned. The `graphical_element` filter handles it; historical scans
barely need it.

**Do not read a per-document rate over a handful of figures.** Half the
documents here carry fewer than ten gold figures, and two carry one. The
corpus-wide numbers are meaningful; `Chun1882c` at precision 0.25 is three
surplus records on a one-figure document.

**Front matter costs figures twice.** In the classifier diagnostic that drops
only `graphical_element`, trimming it with `keeppages` improved both recall
*and* precision (0.894 → 0.923, 0.962 → 0.967), because a plate atlas's
duplicate captures and a bound volume's title page are surplus figures that no
classifier can recognise as not-the-paper. Those are the release-candidate
comparison values, not the stricter default-MCP row below.

## Where it stands

Release candidate, 35 documents, 376 gold figure blocks (316 `[FIGURE]`,
60 `[PLATE]`):

| filter | recall | precision | F1 |
| --- | --- | --- | --- |
| all entries | 0.936 | 0.876 | 0.905 |
| drop `graphical_element` | 0.920 | 0.972 | 0.945 |
| drop uncaptioned `graphical_element` | 0.926 | 0.967 | 0.946 |
| **shared evidence types** (default MCP surface) | **0.875** | **0.985** | **0.927** |
| captioned only | 0.880 | 0.979 | 0.927 |

Caption binding on the full 35-document regression replay, scored on figure
*numbers* over the default MCP evidence types: **recall 0.596, precision
0.911** (286 matches, 314 reported, 480 gold). Before the v1.3 association
repairs, replaying the same persisted Docling documents produced **0.558 /
0.890**; the immediately preceding parser produced **0.588 / 0.910** (282 /
310). No extraction or gold input changed in this comparison. The last four
matches come from a grouped `Figg. 2-5` legend and one structurally linked
caption whose damaged label follows a species heading. Direct `Fic.` captions
are deliberately not treated as plate legends: a looser trial recovered real
numbers but also cloned same-page running-text cross-references.

The clean source-PDF candidate before that final parser change measured
**0.588 / 0.916** (282 / 308). Replaying the six locally retained clean
Docling artifacts that contain the conspicuous failures, while leaving the
other 29 records unchanged, measures **0.608 / 0.918** (292 / 318): ten added
correct page/number pairs, no removals and no false additions. This targeted
clean-artifact result includes `Fic, 11`, `Fic. ro` → 10, `Fic. ror` → 101,
leading dash/underscore noise, the embedded label, and the `Figg.` range. It
is evidence for the bounded repair, not a substitute for the complete clean
source-PDF rebuild required at the release gate.

The **0.596 is end-to-end page/number coverage, not the accuracy of the final
candidate selector in isolation**. The current replay reports only 314
page/number pairs for 480 gold pairs, and their per-page distribution can cover
at most 308 gold pairs (0.642) without recovering another number-bearing
record. Figure-only documents already match 102/124 gold pairs, with 102/104
reported pairs correct; mixed figure/plate documents match 179/323, with
179/204 reported pairs correct; and plate-only documents match 5/33, with 5/6
reported pairs correct. The selector is therefore near its observed input
ceiling on ordinary layouts. Raising the overall rate requires upstream label
discovery and logical plate expansion, not just another caption-ranking
heuristic.

Panel splitting is now a separate acceptance measure rather than an inference
from figure-number binding. The independent gold parser finds **98 captions
that explicitly enumerate lettered panels**. Replaying current caption
selection and panel parsing over the persisted Docling artifacts reports a
panel declaration for **87**, an exact label set for **84**, and label recall /
precision of **0.895 / 0.997**. It declares no panels for the 176 same-page,
same-number gold figure blocks that the gold parser classifies as non-panelled.
This metadata replay deliberately retains the older entry/classification
population, so its figure-number totals are not a replacement for the 0.588 /
0.910 full regression result above. The clean release rebuild must confirm
both measures together.

The scorer now parses hand-transcribed gold labels independently of the
production parser and reports the raw artifact and default MCP surface
separately. That separation is load-bearing: an earlier version reused
production's fuzzy OCR parser, so changing the extractor silently changed the
supposedly fixed gold denominator from 496 to 480. The compact anatomical key
`Pl.M.` is explicitly excluded: inside these figures it means mouth-plate, not
Roman-numeral Plate 1000.

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
| `bare_label` | 229 | no — "Fig. 3" and nothing else |
| `prose_caption` | 113 | yes |
| `lettering_only` | 29 | no — labels engraved on the plate |
| `nothing_printed` | 5 | no |

**Most "captions" in this literature are not captions.** 229 of 376 blocks are
a bare label. Caption *text* similarity is computable for only 85
number-matched pairs in the served-bundle replay, which is why it is reported
and not headlined.

## Where corpus does well

**Mid-century scanned monographs, which is most of this material.**
1950–1999: recall 0.975, precision 0.964. `Totton1965a` — 226 pages, 195 gold
figures — scores recall 0.974 and precision 0.995.

**Scans generally**: 0.932 / 0.941 raw, against born-digital at 0.962 / 0.622.
That inversion is worth internalising. The modern papers are easier to *read*
and harder to *count*.

**Figure numbers, after #205.** `parse_figure_number` fired for only 32% of
figures; the cause was pleasing once found — the opener word is OCR-damaged,
`Fic.` and `Frc.` for `Fig.`, and across 320 captions the damaged spellings
were **more common than the correct one**. Coverage is now 67.6% at 98.2%
precision.

## Where it has trouble

**Publisher furniture on modern papers.** `Ahuja_etal2026`: 6 gold figures,
27 records. Logos, icons and rules. The `graphical_element` filter removes
them from the served surface, but the raw record still carries them, and
`figures_report.html` shows them.

**Historical plates carrying several engravings under one legend.** On the
pre-repair served bundle, `Vanhoeffen1906` scores 0.781 recall / 0.943
precision. The v1.3 regression replay over its stored Docling artifact scores
**0.922 / 0.967**: enumerating caption blocks are split into per-number
entries, duplicate picture assignments are reconciled when the counts form a
complete bijection, and collected prose on the following page enriches
already-found bare labels instead of being cloned onto the wrong plate.
Pass 2.5 now also admits the shared plate to ROI detection with numeric figure
targets kept separate from panel letters. Pass 3 runs once on the plate and
puts each detected region on its logical figure record; the default OCR path
uses an exact caption-derived number allow-list, while vision consumes the
same targets. ROI accuracy still needs a corpus run with the selected backend;
the caption-binding rate measures ownership and numbering, not crop geometry.

The other two conspicuous acceptance failures move in the same replay:
`Hosiaetal2024` goes from 0.692 / 0.900 to **0.923 / 1.000**, and
`Ahuja_etal2026` from 1.000 / 0.500 to **1.000 / 1.000**. These are targeted
artifact replays, not a substitute for the full-corpus rebuild required before
release.

**Caption binding before 1900: recall 0.152** (5 of 33 numbers). The numbers
are printed on the plates as engraved lettering, not as text, so there is
nothing for a text-based parser to find. Four documents score zero:
`Bernstein1934`, `Chenetal2015`, `Eschscholtz1825`, and `Tilesius1814`.

**Small-denominator eras look worse than they are.** pre-1800 precision 0.400
is four gold figures against ten records; 1800–1899 precision 0.438 is nine
against sixteen. Real, but not a rate.

## Furniture: `graphical_element` is an actionable predicate

Established by scoring the corpuscle under each candidate filter rather than
by assuming which field name sounded right. Dropping `graphical_element` takes
precision from 0.876 to 0.972 for 0.016 of recall. The default MCP surface also
hides `unclassified`: precision reaches 0.985, but recall falls to 0.875. That
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
  extraction.** `Ahuja_etal2026` has 6 gold blocks against 27 records because
  docling counts logos; `Vanhoeffen1906` has 67 against 56 because the gold
  counts engravings and docling counts plates. Caption binding is the better
  posed question, and it is the one with the lower score.
- **The metrics weight every figure equally**, which — as with prose coverage
  and taxonomic tokens (#244) — is not how they are valued downstream.
