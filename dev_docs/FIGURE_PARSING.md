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

**The served surface is in good shape; the raw record is noisier than it
looks.** Every MCP tool filters on `_REAL_FIGURE_TYPES` (`figure`, `plate`,
`subpanel`), so `graphical_element` never reaches a client. Read the served
numbers, not the raw ones — the gap between them is 0.10 of precision.

**If you are building a corpus of modern papers, expect furniture.** Precision
on born-digital documents is 0.607 raw against 0.919 for scans, because
publisher logos, ORCID icons and journal rules are figures as far as a layout
model is concerned. The `graphical_element` filter handles it; historical scans
barely need it.

**Do not read a per-document rate over a handful of figures.** Half the
documents here carry fewer than ten gold figures, and two carry one. The
corpus-wide numbers are meaningful; `Chun1882c` at precision 0.25 is three
surplus records on a one-figure document.

**Front matter costs figures twice.** Trimming it with `keeppages` improved
both recall *and* precision (0.894 → 0.923, 0.962 → 0.967), because a plate
atlas's duplicate captures and a bound volume's title page are surplus figures
that no classifier can recognise as not-the-paper.

## Where it stands

Release candidate, 35 documents, 376 gold figure blocks (316 `[FIGURE]`,
60 `[PLATE]`):

| filter | recall | precision | F1 |
| --- | --- | --- | --- |
| all entries | 0.939 | 0.872 | 0.904 |
| **drop `graphical_element`** (the served surface) | **0.923** | **0.967** | **0.944** |
| drop uncaptioned `graphical_element` | 0.928 | 0.961 | 0.945 |
| drop `graphical_element` + `unclassified` | 0.875 | 0.979 | 0.924 |
| captioned only | 0.880 | 0.974 | 0.925 |

Caption binding, scored on figure *numbers*: **recall 0.574, precision 0.878**.

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
a bare label. Caption *text* similarity is computable for 86 pairs out of 465
gold numbers, which is why it is reported and not headlined.

## Where corpus does well

**Mid-century scanned monographs, which is most of this material.**
1950–1999: recall 0.975, precision 0.964. `Totton1965a` — 226 pages, 195 gold
figures — scores recall 0.974 and precision 0.995.

**Scans generally**: 0.935 / 0.919, against born-digital at 0.974 / 0.607.
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

**Historical plates carrying several engravings under one legend.**
`Vanhoeffen1906` still scores recall 0.806 — its plates print six
separately-numbered figures under a single legend and docling extracts each
plate as *one* picture. #203 recovered 23 of those by giving each legend entry
its own record sharing the plate image; the engravings are still not cropped
apart, because locating them needs OCR of the lettering on the plate itself.

**Caption binding before 1900: recall 0.091** (3 of 33 numbers). The numbers
are printed on the plates as engraved lettering, not as text, so there is
nothing for a text-based parser to find. Six documents score zero:
`Bernstein1934`, `Carre1968_Hippopodius_tr`, `Chenetal2015`,
`Eschscholtz1825`, `Tilesius1814`, `Totton1965b`.

**Small-denominator eras look worse than they are.** pre-1800 precision 0.400
is four gold figures against ten records; 1800–1899 precision 0.438 is nine
against sixteen. Real, but not a rate.

## Furniture: `graphical_element` is an actionable predicate

Established by scoring the corpuscle under each candidate filter rather than
by assuming which field name sounded right. Dropping `graphical_element` takes
precision from 0.872 to 0.967 for 0.016 of recall. Taking `unclassified` as
well is over-reach — it holds real figures and costs five points of recall for
one of precision.

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

- **Panels are not split.** A plate's engravings share one image; a
  multi-panel modern figure is one record. #207 found the
  whole-figure/subpanel branch of `dedupe_figures` unreachable — its overlap
  measure is symmetric and the looser threshold always fires first — which is
  corroborated by no figure in the reference corpuscle carrying
  `figure_type: subpanel`.
- **Counting figures measures the definition of "figure" as much as the
  extraction.** `Ahuja_etal2026` has 6 gold blocks against 27 records because
  docling counts logos; `Vanhoeffen1906` has 67 against 56 because the gold
  counts engravings and docling counts plates. Caption binding is the better
  posed question, and it is the one with the lower score.
- **The metrics weight every figure equally**, which — as with prose coverage
  and taxonomic tokens (#244) — is not how they are valued downstream.
