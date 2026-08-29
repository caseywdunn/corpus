# PLAN.md — Corpus pipeline (v1.2)

**Prior cycles are recorded elsewhere, not here.** A minor or major
release entry in [CHANGELOG.md](../CHANGELOG.md) opens with the cycle's
organizing theme — required going forward by
[CONTRIBUTING.md](../CONTRIBUTING.md)'s release ritual, step 2, and
present on every entry back to 0.4.0 — and each cycle's punch list is
preserved in its own tag's copy of this file
([v0.6.0](https://github.com/caseywdunn/corpus/blob/v0.6.0/dev_docs/PLAN.md),
[v1.0.0](https://github.com/caseywdunn/corpus/blob/v1.0.0/dev_docs/PLAN.md)).
In one line each: v0.1 shipped the extraction → annotation → indexing →
MCP-serving stack; v0.2 hardened internals; v0.3 collapsed the user
surface into one CLI plus a per-corpuscle `config.yaml`; v0.4 closed
silent-failure modes and gated the next cycle on tiered CI; v0.5 was the
served-bundle quality cycle; **v0.6** froze the MCP surface at **38
tools**, which every cycle since has held; **v1.0** made a green CI badge
mean a fresh install works, validated by 1,769 of 1,769 documents on
Bouchet with 0 stage failures; **v1.1.1** made `CITATION.cff` schema-valid,
which is what had been silently blocking Zenodo archival since v0.3.0.

**v1.1 (2026-08-07) is the one v1.2 grows directly out of** — its theme
paragraph is in the CHANGELOG; what matters here is the inheritance. Its
699-paper viburnum validation build opened
[#186](https://github.com/caseywdunn/corpus/issues/186),
[#187](https://github.com/caseywdunn/corpus/issues/187) and
[#188](https://github.com/caseywdunn/corpus/issues/188), three of the
four items v1.2 picks up. And `ocrlang`
([#176](https://github.com/caseywdunn/corpus/issues/176)) established the
mechanism the feature half of v1.2 extends: a per-document bib directive,
carried by per-stage fingerprints descended from `processed.pdf`.

**v1.2 is the extraction-fidelity cycle.** Every quality signal the
pipeline has ever had measures it against *itself*: the #185 soft rates
are corpus-internal consistency, the quality gates are per-document
plausibility, and #187's fingerprints compare one build to the next.
None of them can tell whether the text is *right*. The
[`dunnlab/siphonophores`](https://github.com/dunnlab/siphonophores)
library repo closed that gap on 2026-08-24 with `transcriptions/`:
**35 documents, 761 pages, 1594–2026, 13 languages**, every page
transcribed from a rendered page image only. The protocol forbids the
transcriber from opening any software extraction of the page it is
working on, and that independence guarantee is the whole value — an
extractor scored against this set is not being compared to a cleaned-up
version of itself. The set was chosen to span the collection's hardest
axes (Fraktur, Cyrillic, CJK, vertical Japanese, plate-only atlases,
documents with no text layer at all), and it carries 348 `[FIGURE]`
blocks, 76 `[PLATE]` and 65 `[TABLE]` with their printed labels and
captions verbatim, so it is ground truth for figure and caption
extraction as much as for text.

So the cycle is: **measure OCR, docling, Grobid and figure/caption
extraction against truth, then fix what the measurement finds.** The
feature half is the same shape — per-document bib directives that let an
operator correct a document the pipeline gets confidently wrong, which is
where `ocrlang` left off. Scope is deliberately not fixed beyond the six
goals below; as in v1.1, the material is expected to redirect it.

Doc map unchanged: architectural background in
[OVERVIEW.md](OVERVIEW.md); per-feature history in
[CHANGELOG.md](../CHANGELOG.md); the API contract in
[API_STABILITY.md](API_STABILITY.md); HPC operations in
[BOUCHET.md](BOUCHET.md); deployment in [DEPLOY.md](../DEPLOY.md);
platform-portability criteria in
[PLATFORM_SMOKE.md](PLATFORM_SMOKE.md). Open work is tracked in
[GitHub issues](https://github.com/caseywdunn/corpus/issues).

## Standing gates

The gate v1.0 established is permanent now, and not a checklist item:

> **A clean-room install from `environment.yaml` must be verified by CI,
> not by hand, before a release is tagged.**

v1.0.0 was the first release held to it — T3 ran on the release PR and
was green on the release commit before the tag existed.

Now that [`clean-room.yml`](../.github/workflows/clean-room.yml) (**T3**)
lives on the default branch, all three of its triggers work: the weekly
`schedule:` is the standing drift detector, `workflow_dispatch` is
available (it returned `HTTP 404` while the lane sat on a feature
branch), and **a pull request targeting `main` — which *is* the release
proposal — runs the lane automatically**. That last one is the path that
satisfies the gate before the merge rather than depending on someone
remembering to dispatch it.

[`dev_docs/ec2_smoke.sh`](ec2_smoke.sh) (**T3-bare**) stays manual and
pre-release: it covers the one thing T3 cannot, the bare-host bootstrap
(apt, miniforge install) on a real Ubuntu EC2 instance, against
[PLATFORM_SMOKE.md](PLATFORM_SMOKE.md)'s criteria.

The per-push tiers (T0, T1/T2, T1-compose) and the full tier table live
in [CONTRIBUTING.md](../CONTRIBUTING.md).

**"Re-measure the baselines" was the wrong ask, and is retired.** The
recorded numbers (T0 at 623 passed / 2 skipped, `corpus_required` at 163
/ 14 / 4, taken at `1dbb69a`) were treated as a regression reference.
They cannot be one. T0 went 755 → 799 in a single afternoon of v1.1
because tests were *added*, and a number that moves whenever someone
writes a test is a changelog rather than a detector. The property that
matters for T0 — zero failures — is already enforced by CI on every push,
so the count adds nothing on top of it.

What a real reference records is pipeline *output*: the #185 soft rates,
quality-gate counts, and bundle manifest counts for a fixed corpus,
diffed against a rebuild of that same corpus.
[#187](https://github.com/caseywdunn/corpus/issues/187) specifies it and
targets the **gold corpuscle** (v1.2 §5) — big enough for a rate to mean
something, small enough to rebuild per release. That comparison was run
by hand against a viburnum rebuild during v1.1 and is what proved the
3.4 GB → 2.3 GB drop was #184's re-encoding rather than lost content.

The `--deselect` question that paragraph raised is answered: **none**.
[#167](https://github.com/caseywdunn/corpus/issues/167) removed all three
flags, T1 now runs bare `-m corpus_required`, and T2 additionally ignores
`test_reference_extraction.py` because Grobid is disabled there. Both
workflow files state the reasoning inline.

## v1.2 — extraction fidelity, measured against ground truth — **released**

Kept here through the release commit because the sections below record *why*
each decision went the way it did, and several were reversed by measurement
along the way. The CHANGELOG says what shipped; prune this section at the
head of v1.3.

**The release under work.** Six goals, numbered below and referred to as
§1–§6 elsewhere in this file. §1 is the instrument; §2 and §3 point it at
figures and references; §4 is the operator's recourse when the instrument
finds something no code change can fix; §5 rebases the drift reference on
the same corpus; §6 is cleanup this cycle owes. No wave ordering — §1 is
a precondition for *reading* §2 and §3's numbers, not for writing their
code.

### 1. The fidelity harness — [#193](https://github.com/caseywdunn/corpus/issues/193), **landed**

`tools/qc/fidelity.py` scores a built corpuscle against the gold tree.
Ships in this repo, not the library repo, so the evaluator is versioned
with the extractor it grades. Manual pre-release tier **T5**;
`tests/test_fidelity_harness.py` covers its arithmetic in T0 against a
committed fixture.

- Resolve gold ↔ corpuscle documents through `transcriptions/sources.json`
  (stem → `library/<LETTER>/<stem>.pdf` + sha256) and the corpuscle's
  content-addressed hash dirs. **Bind on the checksum, not the stem.**
  The library holds two editions of Léry's *Histoire d'un Voyage*, both
  once named `Lery1594.pdf`, setting the same *Physalia* passage on
  different folios; only the 1594 is transcribed.
- **Per-page alignment needs no pipeline change.** `text.json` carries a
  flat `text` string and a `pages` *count* only, but `docling_doc.json`
  keeps `prov[].page_no` on every text item, so per-page extracted text
  is reconstructible from the persisted artifact and comparable to
  `page_NNN.txt`.
- Reuse the measures `siphonophores/scripts/crosscheck.py` established:
  token-level `SequenceMatcher` similarity **plus** the order-insensitive
  `coverage` and `recall`, gold markup stripped, case/diacritic/
  punctuation normalised.
- **Read `transcriptions/CROSSCHECK_REPORT.md` before trusting a
  number.** It documents five ways this signal misleads, each found by
  inspecting pages: a garbage text layer scores ~0.05–0.15 and says
  nothing; scrambled reading order is punished as hard as wrong content
  (Hosiaetal2024: 0.59 sequence, 0.99 coverage); the gold correctly holds
  *more* than any extractor can (Ahuja p8 — 1007 gold words against 54,
  because the rest is engraved lettering); extractors hallucinate text
  from image texture; and long-s typography guarantees mismatch on
  18th-century pages. Its phase-1 funnel went 83 flagged pages → **0**
  confirmed gold omissions. A harness that ignores this will spend the
  cycle chasing pages that are already correct.
- Report per document and per page, **segmented by the axes the set was
  built to span** — script, era, scanned vs born-digital. A corpus-wide
  mean over 13 languages is not actionable.
- Commit a two- or three-page fixture under `tests/fixtures/` so T0
  covers the scorer's own arithmetic with no corpus dependency. The full
  run is manual and release-time, alongside T3-bare/T4 — add the row to
  CONTRIBUTING.md's tier table when it lands.

**Which side is on trial is the part that had to be got right — and it
is not a matter of arithmetic.** `crosscheck.py` used a poppler text
layer as the yardstick to validate the *transcription*; this harness uses
the gold as the yardstick to validate the *extractor*. The measures
themselves are symmetric and identical in both (`coverage` is
`|gold ∩ other| / |gold|` either way), so nothing is mirrored. What
changes is two judgement calls the arithmetic cannot make: which measure
leads — `recall` there, `coverage` here — and what an unscorable page
counts as. The report excludes a page whose comparison layer carries no
signal, correctly for its purpose; here an empty extraction *is* the
finding and is scored 0.0. Adopting the report's exclusion policy instead
would drop 57 of 761 pages and lift the median from 0.891 to 0.908,
hiding precisely the pages that need work. Those two calls are why three
of the report's five false positives are this harness's true positives —
a garbage text layer, a plate whose lettering exists only as image, and
text hallucinated from image texture. The harness docstring tabulates the
mapping.

**First run, against the 35-document gold corpuscle** (761 pages, ~7 s).
Corpus-wide median coverage 0.891, but the segments are the result:

| axis | pages | median coverage | pages < 0.5 |
| --- | --- | --- | --- |
| CJK | 19 | **0.344** | 11 |
| Cyrillic | 32 | 0.924 | 3 |
| Latin | 710 | 0.897 | 98 |
| pre-1800 | 26 | **0.645** | 4 |
| 1800–1899 | 80 | 0.807 | 24 |
| 1900–1949 | 99 | 0.881 | 19 |
| 1950–1999 | 463 | 0.920 | 57 |
| 2000– | 93 | 0.901 | 8 |
| born-digital | 75 | 0.882 | 7 |
| scanned | 686 | 0.892 | 105 |

Three things this says immediately, and none of them was knowable before:

- **CJK is the worst axis by a wide margin**, and vertical Japanese is
  most of it. `Kawamura1911a` is bilingual — an English translation
  followed by the 1911 vertical original — and scores 0.344 overall with
  one page flagged `script_missing`. This is where #186's OCR-*mode*
  override lands, and it is the strongest argument yet for scheduling it.
- **The pipeline is dramatically better than poppler on Fraktur, and now
  there is a number.** The library-side cross-check recorded roman OCR
  over Fraktur at ~0.05–0.15 with nothing triageable;
  `Eschscholtz1825` scores 0.830, `Keferstein_Ehlers1860` 0.801,
  `DeHaan1827` 0.824. The OCR path is doing exactly what it was built
  for. That comparison is only legible because the measures were kept
  identical to `crosscheck.py`'s.
- **54 pages extract to nothing at all**, concentrated in plate-heavy
  documents (`Quoy_Gaimard1834Plates` at 0.216 median coverage is the
  floor). That is §2's material, sized: it is figure and plate lettering,
  not prose, and `gold_structural_share` per page separates the two.

**Prose is the measure; figure text is reported, not optimised.**
Scoring both together answered neither question. Text inside a
`[FIGURE]`/`[PLATE]` block is engraved plate lettering and panel labels —
12.3% of the gold set's words but about 70% of its worst pages — so a
combined number is dragged down by the material the pipeline is least
expected to recover, while body text, which is its actual job, vanishes
into the average. Split: **0.946 prose**, 0.731 figure text, 0.898
combined. `[TABLE]` counts as prose; measured the other way the corpus
figure moves 0.002.

The split changes what the run says to work on. The 1800–1899 band reads
0.812 prose against 0.114 figure text, so its apparent weakness was
almost entirely plate lettering, and `Chenetal2015` moves 0.667 → 0.812
the same way. Two findings survive, and one was being hidden:

| segment / document | prose | figure text |
| --- | --- | --- |
| CJK | **0.351** | 0.406 |
| pre-1800 | **0.645** | 0.616 |
| `Linnaeus1735` | **0.628** | 0.634 |

**What the two survivors turned out to be** — investigated, and none of
it was the long-s typography the cross-check report had primed us for:

- **Vertical CJK has no model selected**
  ([#196](https://github.com/caseywdunn/corpus/issues/196)).
  `jpn_vert`, `chi_sim_vert`, `chi_tra_vert` and `kor_vert` are all
  installed; `scan.py` can reach none of them. `Kawamura1911a` isolates
  it exactly — same PDF, same scan, same OCR call — with its English
  translation at ≈0.99 prose coverage over pages 1–13 and its vertical
  Japanese at ≈0.25 over pages 14–23. The fix has a precedent in the same
  file: `deu` → `deu_latf` on Fraktur, which is why Fraktur scores
  0.80–0.85 here against poppler's recorded 0.05–0.15.
- **`tesseract_packs` is recorded empty on three documents**
  ([#197](https://github.com/caseywdunn/corpus/issues/197)) whose
  language detection was untrusted. OCR itself was fine — `Linnaeus1735`
  ran with `eng+deu+deu_latf+fra+lat+spa+por` — so the operator-facing
  record disagrees with the invocation, which `scan.py` asserts cannot
  happen. That record is what `corpus status` and #176's `ocrlang`
  workflow tell an operator to consult.
- **The untrusted fallback union is not the pre-1800 problem — tested
  and refuted.** The hypothesis was that seven competing packs cost more
  than they buy, since `_compose_ocr_langs`'s docstring says accuracy
  degrades as packs multiply. `Linnaeus1735` was rebuilt with
  `ocrlang = {lat}` pinned (#176), same config, same CPU accelerator, and
  scored against the gold set: prose coverage **0.628 → 0.549**, recall
  0.630 → 0.518, similarity 0.562 → 0.443, pages under 0.5 coverage 2 → 4.
  Worse on 12 of 13 pages, with `excess_novel` up on 12 of 13. Pinning the
  single correct language made it worse, so for antiqua-set early-modern
  Latin the extra packs are evidently supplying character coverage `lat`
  alone lacks. The fallback is doing its job, and whatever depresses
  pre-1800 prose — 0.607, 0.628, 0.644 across the three worst — is not
  pack over-selection. First real use of the #176 directive against ground
  truth, and it cost one afternoon to kill a plausible hypothesis.

  Unmeasured but suggestive: OCR ran 39 s on one pack against roughly 11
  minutes on seven. Not comparable as they stand — the seven-pack run was
  one of seven documents extracting in parallel — but if pack count
  dominates, the union buys +0.08 coverage for a large wall-clock
  multiple. Worth measuring properly before widening it.

  Its score also depends on the table decision more than any other
  document's:
  0.628 with `[TABLE]` counted as prose against 0.400 with tables on the
  figure side. *Systema Naturae* is largely tabular — 18,358 characters
  of table-cell text — so it is the one document where that classification
  is a lever rather than a rounding difference, and any number quoted for
  it has to say which convention it used. Under the shipped convention the
  worst pre-1800 document is `Hjortberg1769` at 0.607, and the worst prose
  anywhere in the set is `Kawamura1911a` at 0.351.

- **docling picks a GPU the pinned torch cannot use**
  ([#198](https://github.com/caseywdunn/corpus/issues/198)), found while
  running the experiment above. Same machine, same pinned set, 20 days
  apart: `Accelerator device: 'cpu'` on 2026-08-06, `'cuda:0'` and 20 ×
  `CUDA error: no kernel image is available` on 2026-08-26. The card is a
  GTX 1080 at compute capability 6.1 and the pinned torch ships `sm_75`
  and up. Nothing in the project changed — the GPU became visible to
  torch, and `torch.cuda.is_available()` alone is what docling decides on.
  A driver appearing turned a working install into a broken one, which is
  precisely the reproducibility the #98 pins exist to provide. Nothing in
  a corpuscle records which accelerator produced it.

**Two defects the harness found on its way in**, and the reason to
validate a measuring instrument against a second source before trusting
it.

The first was in reading the gold. The transcription convention uses `[`
for markers, but pages print brackets too and notes discuss them, so a
scanner that counts every `[` as a nesting level mis-parses three
constructs the set actually contains — a note quoting the bracket
character, a `[sic]` or `[21]` quoted inside a note, and one unterminated
`[continued opposite`. The first alone leaked whole notes into 13 of one
document's 17 pages, which is why it posted 0.767 coverage against 0.998
recall: gold apparently holding text no extractor could find. Gating
marker recognition on a known keyword vocabulary brings the structural
tag counts to exactly the **348 `[FIGURE]` / 76 `[PLATE]` / 65 `[TABLE]`**
the gold set documents — against 341 / 76 / 60 before — with no page left
unbalanced. That agreement is the check: the parser and the set now count
the same thing. The affected document moves to 0.927 and no other moves
by more than 0.002. **This was a precondition for §2**, which scores
against those very block boundaries.

The second was in reading the corpuscle: a
docling table item has no `text` field at all — its content lives only in
`data.table_cells[].text` — so a per-page walk reading `text` discarded
every table on the page. Checking reconstruction against `text.json`
rather than assuming it caught 8 of 35 documents short, the worst by 26%
of its tokens. Unfixed, those pages would have been reported as the
pipeline losing prose it had in fact extracted. All 35 now reconstruct
≥95.8% of `text.json`'s tokens.

### 2. Figure fidelity — [#194](https://github.com/caseywdunn/corpus/issues/194) detection **landed**; captions [#195](https://github.com/caseywdunn/corpus/issues/195) **blocked on [#205](https://github.com/caseywdunn/corpus/issues/205)**

Score `figures.json` — `caption_text`, `page`, `figure_id`,
`panels_from_caption` — against the gold `[FIGURE]` / `[PLATE]` blocks:
figures found vs missed, captions bound to the right figure, panel splits
right, and plate lettering, which on an engraved plate is often the only
text that exists. The code to work with rather than around is
`pipeline/figures.py` — `extract_caption_info` (structural docling link,
then the proximity heuristic with its cross-page penalty) and
`parse_panels_from_caption`, including the archaic `Plate IV.` / `Taf.
III.` handling.

This is the first real measurement of the caption heuristic. OVERVIEW.md
already calls caption association "the highest-value annotation per
figure and the hardest in historical layouts"; the gold set is what turns
that from an assertion into a number.

**§1's run put the damage here**, which is what makes this the next
piece of work rather than §3. Figure and plate text is 12.3% of gold
tokens (36,470 of 295,745) and roughly 70% of the failures: 47 of the 54
pages that extracted to nothing are >80% inside a structural block, as
are 78 of the 112 pages under 0.5 coverage. `Totton1965a` contributes 40
empty pages, every one a plate.

**Detection is measured: recall 0.833, precision 0.841.** 71 gold
figures have no entry; 67 entries have no gold block. And
`figure_type == "graphical_element"` turns out to be an actionable
furniture predicate — established by scoring the corpuscle under each
candidate filter rather than assuming which field name sounds right:

| filter | recall | precision | F1 |
| --- | --- | --- | --- |
| all entries | 0.833 | 0.841 | 0.837 |
| **drop `graphical_element`** | **0.811** | **0.964** | **0.881** |
| drop uncaptioned `graphical_element` | 0.818 | 0.940 | 0.875 |
| drop `graphical_element` + `unclassified` | 0.708 | 0.987 | 0.824 |
| captioned only | 0.724 | 0.959 | 0.825 |

Surplus falls 67 → 13 for nine real figures lost. Taking `unclassified`
as well is over-reach: it holds 49 real figures and recall drops ten
points. The two halves are split out because they have nothing in common:
[#204](https://github.com/caseywdunn/corpus/issues/204) is the furniture
decision, [#203](https://github.com/caseywdunn/corpus/issues/203) the
recall hole — where `Vanhoeffen1906` alone is 33 of the 71 misses, its
pages of six separately-numbered engravings each collapsing to one
2430×4002 region captioned with one of the six labels and no panels
split. Segments show why no single number would do — born-digital
precision 0.562 against 0.890 for scans, because modern papers carry
publisher furniture historical scans do not; and 1900–1949 recall 0.477,
driven by `Vanhoeffen1906` (34 of 67 found) and `Kawamura1911a` (6 of
16). **What to do with the filter is a decision, not a measurement** —
excluding those entries from the served bundle, or fixing classification
upstream — and it is left open.

**Counting figures measures the definition of "figure", not the
extraction.** 420 entries against 424 gold blocks corpus-wide — but the
totals agreeing is a coincidence. `Ahuja_etal2026` has 6 gold blocks
against 31 entries because docling counts logos and icons;
`Vanhoeffen1906` has 67 against 34. Caption binding is the measure that
is well posed.

**Caption binding is blocked, and the block is worth stating.** §2's
plan was to bind on the figure *number* rather than the caption text,
because a number is language-independent. It only exists for a third of
the figures: `figure_number` is populated on **135 of 420 entries
(32.1%)**, with **15 of 35 documents carrying none at all**. Where it does
fire it is right — 131 of 135 (**97.0%**) are genuinely printed on that
page, allowing #16's Roman↔Arabic normalisation, which works. The gold
says at least 88% of figure blocks carry a printed number, so the gap is
coverage, not correctness.

Binding on a signal present for a third of the population would give the
other two thirds a silent pass, which is precisely the error shape this
cycle keeps catching. So
[#205](https://github.com/caseywdunn/corpus/issues/205) came first, and
has **landed**: coverage is now **67.6%** with precision *up* at 98.2%.
The largest cause was pleasing — the opener word is OCR-damaged, `Fic.`
and `Frc.` for `Fig.`, and across 320 captions the damaged spellings are
**more common than the correct one**.

It also surfaced a latent defect by amplifying it. `dedupe_figures`
grouped same-numbered figures on number alone, and its bbox tests carry no
page, so two figures at similar coordinates on different leaves read as
redundant crops — routine in a document that is its own translation.
Numbering more figures fed more of them into that grouping, and
`Carre1969_Nanomia_tr` fell from 18 extracted figures to 13. Keying on
`(page, number)` takes it to **22 — every figure docling detects, and
exactly what the gold records** — with 16 numbered against 7. That is a
direct dent in §2's recall hole
([#203](https://github.com/caseywdunn/corpus/issues/203)), found from the
number side rather than the detection side.

Tracing it also turned up
[#207](https://github.com/caseywdunn/corpus/issues/207): the
whole-figure/subpanel branch of `dedupe_figures` is unreachable, because
its overlap measure is symmetric and step 1's looser threshold always
fires first. No figure in the reference corpuscle carries
`figure_type: subpanel`, which corroborates it. Its largest single cause is pleasing: the opener word is
OCR-damaged, `Fic.` and `Frc.` for `Fig.`, consistently and document-wide
(53 of 185 captioned-but-unnumbered entries), on documents whose captions
are otherwise extracted perfectly. Anchoring the label to the caption
start costs another 59, and `Figg.` and `第九圖` are unhandled spellings.

**A naive caption-text match was tried and it lies — 44% is not the
answer.** Token similarity at a 0.6 threshold reported 44% bound, 7%
wrong page, 49% unmatched. Reading the pages showed most of that is
artifact, in the same shape CROSSCHECK_REPORT.md warns about:
`Chenetal2015` scored 0 of 10 because the document prints every caption
twice and the extraction took the Chinese while the matcher was handed
the English — every figure is in fact bound correctly;
`Carre1969_Nanomia_tr` scored 1 of 18 because its plate pages carry only
`FIG. 1` as printed matter, two tokens no threshold can bind, with the
prose caption on another page. #194 therefore binds on the **figure
number** first, which is language-independent, and classifies each gold
block — prose caption, bare label, plate lettering, nothing printed —
before scoring, so no rate is computed over blocks it does not apply to.
One genuine finding did fall out of the exercise: `Chenetal2015`'s figure
numerals OCR as `图 ;` and `图 <` for 3 and 4, which belongs to the OCR
axis.

**The cycle's fixes were rebuilt and rescored, and one of them was wrong.**
A full 35-document rebuild on `dev` including
[#203](https://github.com/caseywdunn/corpus/issues/203) confirms the
recall it was written for — corpus-wide figure recall **0.849 → 0.917**,
`Vanhoeffen1906` alone recovering 23 of its 33 missing figures. It also
cost precision, 0.861 → 0.807, and **the whole of that loss is one
document**: `Totton1965a` gained 31 spurious records against 3 real ones.

Every one of them is a cross-reference read as a legend entry.
`plate_legend_entries` took a figure number from anywhere in a text item,
and a 226-page monograph of running prose is full of them — a species
heading reading `Plate XX, figures 1, 2`, a parenthetical `Text-figure
106 (see p. 170)`, a citation of someone else's plate. Two on a page
cleared the threshold and the page's one real text-figure was cloned
under the referenced number. **33 of that document's 232 records shared
an image file with another record under a different number**, so asking
for figure 20 returned the picture of figure 53 — a confident wrong
answer, which is worse than a miss.

The evidence that settled it was the page image, not the JSON. Three
rounds of counting said the surplus was concentrated in one document;
rendering p. 105 said why, in one look.
[#231](https://github.com/caseywdunn/corpus/issues/231) anchors a legend
entry to a line that *opens* with a figure label, with the number
required to follow it. On the served surface (`drop graphical_element`):

| build | recall | precision | F1 |
| --- | --- | --- | --- |
| pre-#203 | 0.833 | 0.970 | 0.896 |
| #203 as built | 0.901 | 0.892 | 0.897 |
| **#203 + anchor** | **0.894** | **0.962** | **0.927** |

Caption binding moves the same way — recall 0.587 → 0.581, precision
**0.795 → 0.870**, which is exactly the precision the corpus had before
#203, with binding recall up from 0.527. Prose fidelity is untouched at
0.947 median coverage, as it should be: none of this cycle's figure work
goes anywhere near the text.

Two things this is worth recording for. #203 shipped on a real measured
gain and still regressed the corpus, because the measurement it shipped
on was recall and the cost landed in precision — the harness reports both
and the summary quoted one. And requiring *description* after the number
as well was tried as a second guard and rejected on measurement: it costs
8 of `Vanhoeffen1906`'s 24 real entries for no further precision.

### 3. Grobid output and reference consolidation — two of three answered

Three questions that sound like one, and the code gives them different
answers:

- **`consolidateCitations` — answered, and now reachable.** It had never
  run: `grobid_client.py` defaulted it to 0, `metadata.py` overrode
  nothing, and config exposed only `url` and `disable`. Measured off
  against on, same PDFs: `Ahuja_etal2026` 86.1% → 88.9% DOIs at 3.1s →
  6.4s; `Stepanjants2014` 0% → 6.9%; `Bernstein1934` 0% → 3.6%;
  `Benasso_Stroiazzo1976` 0% → 0%. **Six DOIs across 194 references for
  1.4–2× the Grobid time**, so the default stays off — but weighed rather
  than unexamined. The split is by era, not by the flag: corpus-wide 719
  references are 32.4% DOI-bearing, 69–86% in papers from 2020 on and
  **0% in every document before 1980**, because the historical works are
  not in CrossRef at all. Both flags are now settable from `grobid:`,
  since that arithmetic belongs to the library.

  Original framing: **`consolidateCitations` has never run.**
  `pipeline/grobid_client.py` defaults `consolidate_header=1` and
  `consolidate_citations=0`, `pipeline/metadata.py` calls
  `process_fulltext()` overriding neither, and the `grobid:` block in
  `pipeline/config.template.yaml` exposes only `url` and `disable`. So
  reference consolidation against CrossRef is off and always has been.
  Measure what enabling it costs (a network round-trip per reference) and
  buys, then expose both flags in config with the measured default.
- **In-corpus references — answered: the looseness is justified.** The
  gold carries the reference sections as printed, so reference *parsing*
  can be measured directly, separately from corpus cohesion. Of 659
  references Grobid extracted, **94.8% have a title appearing in that
  document's gold text** — it is not inventing references. Modern papers
  score 100%, `Totton1965a` 95.5%, `Benasso1976` 83.9%,
  `Eschscholtz1825` 50% and `Linnaeus1735` 0%: the same Fraktur and
  pre-1800 axes that limit everything else. So the loose ceiling on
  `references_match_corpus_papers` is right — parsing is good, and what
  that check varies with is how tightly a corpus cites itself, which is
  an assembly property rather than a pipeline one. The evidence is now
  recorded beside the ceiling it justifies.
- **Out-of-corpus references** — this is
  [#155](https://github.com/caseywdunn/corpus/issues/155).
  `get_missing_references` is dominated by unreconciled citation-string
  variants of works already in the corpus; `resolve_reference "Bigelow
  1911 The Siphonophorae"` returns 53 matches, one canonical and ~50
  variants. Full reconciliation is a cycle of its own (DOI normalization,
  block-and-cluster canonicalization, alternate-DOI aliasing, junk
  filtering, probably an auditable LLM adjudication pass) and stays
  unscheduled. **The cheap slice was scoped, measured, and does not
  survive contact with the data**
  ([#225](https://github.com/caseywdunn/corpus/issues/225)). Collapsing
  on normalized (title, year) merges *different papers*: six viburnum
  works titled `phytochemistry` in 1989 are six papers by six author
  teams, and the rule loses five of them. Even with the author
  component it needs order-insensitive author-set agreement and
  DOI-aware canonical selection — which is the full reconciliation
  cycle #155 puts out of scope.

  The reported symptom is not reproducible on either corpus here: the
  gold 35 show 4 duplicate groups in 865 works, and viburnum's high-row
  cases (`Wang 2020` → 70 rows) are seventy different people sharing a
  surname, not one work. A fix wants measuring against the
  1,769-document production corpus that exhibits it.

  Upstream of it, and probably the larger contributor:
  [#226](https://github.com/caseywdunn/corpus/issues/226) — journal
  names recorded as work titles, `phytochemistry` ×62, `nature` ×7,
  `kb` ×19 across 229 works — which manufactures the duplicate groups a
  de-duplication rule would then have to clean up.

### 4. Per-document bib directives

`ocrlang` (#176) established the mechanism: a flat BibTeX field on the
entry whose `file = {…}` matches the PDF, parsed in `bib/parser.py`,
carried through `bib/importer.py` and `bib/export.py`, documented in the
README's directive table. Two more fields follow the same path.

- [x] **Per-document page-range selection**
  ([#188](https://github.com/caseywdunn/corpus/issues/188)) — **landed.**
  `keeppages` applies before scan detection by rebinding the temp PDF the
  later stages already read, so `scan.py` needed no change. The page-number
  question was settled as the issue proposed: `page` stays subset-relative
  because it indexes the artifacts on disk, `source_page` rides beside it,
  and the resolved selection recorded in `scan_detection.json` *is* the map.
  #174 turned out not to be a precondition — `keeppages` rides the same
  per-document fingerprint path `ocrlang` established, across every
  OCR-dependent stage. Original framing: scanner
  front matter, prepended translations and blank runs skew scan detection
  before they waste OCR. The open design question is page-number
  provenance on served figures, and the gold set is unusually well placed
  to settle it: `[PAGE n]` opens every gold page and records the
  **printed** folio where one exists, `[PAGE n: unnumbered]` where it
  does not.
- **Precondition worth stating:** a bib directive that changes extraction
  is a silent no-op unless editing the bib invalidates the fingerprint.
  That is [#174](https://github.com/caseywdunn/corpus/issues/174)
  (fingerprints over bib entries, filenames and config keys), listed
  below as unscheduled. #188 either pulls it in or ships with a
  documented "re-run with `--force`" caveat.
- **Decision pending: [#186](https://github.com/caseywdunn/corpus/issues/186)**,
  the per-document OCR-*mode* override. `redo_ocr` preserves a corrupt
  digital text layer — legacy symbolic Greek fonts read as spaced
  punctuation — so OCR output is discarded and no `ocrlang` value can
  reach it. Not scoped into v1.2, but flagged here because A and B will
  land on precisely those documents, and the gold set contains several.

### 5. The gold corpuscle, and #187 rebased on it

- [x] **The gold 35 become the smoke-test corpuscle** — **landed.**
  `siphonophores_sample` has no references left in the tree; BOUCHET.md and
  both SLURM scripts point at `siphonophore_gold_YYYYMMDD`, and the cycle
  built it a dozen times. Original framing:
  `siphonophores_sample/library` (30 PDFs, **zero overlap** with the gold
  set) is retired. The selection criteria already coincide —
  [BOUCHET.md](BOUCHET.md) picked its 30 to span "born-digital modern +
  historical scans + German Fraktur to exercise the OCR / language /
  figure paths"; the gold 35 were picked to span the same axes, more
  deliberately, and come with truth attached. One corpuscle then serves
  three jobs: pre-production rehearsal, drift reference, and accuracy
  scoring ([#192](https://github.com/caseywdunn/corpus/issues/192), which
  supersedes the corpus-selection half of #187). It is also a *cheaper*
  build than the sample: 761 pages / ~644 needing OCR, against 1,290 /
  ~916. Totton1965a is half that OCR load on its own and is kept
  deliberately — #192 records the trigger for splitting the rehearsal off
  if it dominates wall-clock.
- [ ] **Fingerprint-based regression reference**
  ([#187](https://github.com/caseywdunn/corpus/issues/187)) — **slipped to
  v1.3**, and worth saying so rather than quietly dropping the box. The
  accuracy half of this section landed and was used all cycle; the drift half
  was never started, and `tools/` still has no fingerprint diff. **A measures
  accuracy against truth, #187 measures drift between builds** — the second is
  what tells you a rebuild changed something you did not intend, which this
  cycle answered by hand each time.
- **Recalibration this forces.** `_SOFT_RATE_CEILINGS` in
  `tests/test_corpus_wide.py` is calibrated across two deliberately
  unalike corpora. A hardest-cases set will legitimately post worse rates
  than either, and the right response is the existing
  `CORPUS_SOFT_RATE_CEILINGS` env override — built for exactly this third
  corpus — not looser shipped defaults.

### 6. Housekeeping, in-cycle

- [x] **16 source comments cited PLAN.md sections that no longer exist** —
  `§10` (×8), `§9` (×4), `§7` (×1), `§3` (×3) across `mcpsrv/`,
  `pipeline/`, `tests/` and BOUCHET.md, pointing at a v0.1-era numbering
  this file has been pruned past several times. Repointed at the stable
  docs (OVERVIEW.md, DEPLOY.md) or dropped. A roadmap that gets renumbered
  every cycle is the wrong anchor for a code comment; **don't add new
  ones.**

## v1.3 — skills and usage

Scoped now so v1.2 can stay narrow. Where v1.2 asks whether the corpus is
*right*, v1.3 asks whether anyone can *use* it — the packaging and
recipes that turn 38 frozen tools into an answer to a question someone
actually has. Additive by construction: skills live beside the server,
not inside the API contract.

- [ ] **A `skills/` plugin directory**
  ([#178](https://github.com/caseywdunn/corpus/issues/178)) — where a
  corpuscle's task recipes live and how a client discovers them.
- [ ] **A clade-monograph skill**
  ([#179](https://github.com/caseywdunn/corpus/issues/179)) — the first
  real one, and the synthesis recipe Q2 has been missing since v0.5.
- [ ] **A README quick start**
  ([#180](https://github.com/caseywdunn/corpus/issues/180)) — the
  shortest path from clone to a served answer.
- [ ] **Fingerprint-based regression reference**
  ([#187](https://github.com/caseywdunn/corpus/issues/187)) — carried from
  v1.2 §5, where the accuracy half landed and this did not.
- [ ] **`corpus bib inspect-pages`**
  ([#217](https://github.com/caseywdunn/corpus/issues/217)) — read-only
  per-page evidence for curating `keeppages` and `doclang`. Deferred in v1.2
  and now has a working prototype with a caller: the siphonophore library's
  `scripts/inspect_pages.py`, written free of clade-specific knowledge so
  upstreaming is a move rather than a rewrite. It produced the 515 `keeppages`
  and 1,773 `pagemap` values v1.2's #188 consumes.
- [ ] **The reconciliation rework**
  ([#240](https://github.com/caseywdunn/corpus/issues/240)) — **a release
  theme of its own, not an item here.** v1.2 measured what the current
  cascade does and where it fails: 58% of duplicate ghosts with *identical*
  titles escape the `(first author, year)` block, nothing ever re-reconciles,
  and a corrupted DOI short-circuits the whole cascade
  ([#239](https://github.com/caseywdunn/corpus/issues/239)). The fix is a
  data-model change — observations separated from canonical works — plus a
  model dependency, touching the served MCP surface. #226 and #239 should land
  before it; both corrupt the evidence any clustering pass would consume.
- **Also in scope for the cycle, if it earns its way in:**
  **`export_to_disk` + `suggest_command`**
  ([#88](https://github.com/caseywdunn/corpus/issues/88) Part 2) and the
  **`corpus export` CLI**
  ([#93](https://github.com/caseywdunn/corpus/issues/93)) — the bulk
  export surface for "the LLM isn't the consumer" workflows, which is a
  usage question wearing a tooling hat.

## Unscheduled

Not claimed by a cycle, and none blocks another. Carried so nothing is
silently dropped. Split because the two halves get picked up for
different reasons: a known defect is picked up when it bites someone, an
unbuilt feature when something makes it worth building.

### Open defects

Issue-backed, in dependency-free groups.

**Served-surface correctness**

- [ ] **Taxonomic authority linking assumes zoological authorship**
  ([#175](https://github.com/caseywdunn/corpus/issues/175)), so
  `get_original_description` is structurally dead for any botanical
  corpus — 889 of 913 viburnum taxa had authorship, 0 with a year.
  Overlaps PR #144 below, from the opposite end.
- [ ] **Lexicon translations match only uninflected forms**
  ([#165](https://github.com/caseywdunn/corpus/issues/165)), zeroing
  anatomy coverage on German papers: Eschscholtz prints `Luftblasen`,
  the lexicon has `Luftblase`.
- [ ] **A hub work's depth-1 `get_citation_graph` payload can exceed MCP
  transport limits** ([#166](https://github.com/caseywdunn/corpus/issues/166))
  while `truncated: false` stays accurate.

**Extraction quality**

- [ ] **Abbreviated genus binomials**
  ([#164](https://github.com/caseywdunn/corpus/issues/164)) —
  `Ph. pelagica` resolves to nothing.
- [ ] **`visual_script` is never populated**
  ([#172](https://github.com/caseywdunn/corpus/issues/172)), leaving the
  text-layer/visual cross-check inert.
- [ ] **Move the docling pin forward**
  ([#98](https://github.com/caseywdunn/corpus/issues/98) follow-up).
  Still `docling==2.94.0`. Reproduce on an arm64 Mac, determine whether
  2.95/2.96 broke MPS extraction via an API change or an upstream bug,
  then advance deliberately. Needs Apple-Silicon hardware. **v1.2 §1 gives
  this a criterion it never had** — "better or worse" against the gold
  set rather than against impressions.

**Operator surface**

- [ ] **Progress heartbeat during long per-document stages**
  ([#170](https://github.com/caseywdunn/corpus/issues/170)).
- [ ] **`--filter-gate` is silently ignored without `--list-hashes`**
  ([#169](https://github.com/caseywdunn/corpus/issues/169)) — the hint
  was fixed in 1.0, the underlying flag was not.
- [ ] **Surface the naive-chunker fallback in `corpus status`**
  ([#168](https://github.com/caseywdunn/corpus/issues/168)).
- [ ] **Extend per-stage fingerprints to bib entries, filenames and
  config keys** ([#174](https://github.com/caseywdunn/corpus/issues/174)),
  so editing the input bib stops being a silent no-op. See v1.2 §4 — #188
  depends on this or on a documented caveat.

**Housekeeping**

- [ ] **`pipeline/CITATION.cff` is a hand-maintained copy**
  ([#173](https://github.com/caseywdunn/corpus/issues/173)) with nothing
  enforcing sync — currently byte-identical, so this is a guard rather
  than a fix. `tests/test_citation_cff.py` (v1.1.1) now validates both
  files independently, which is adjacent but not the same check.
- [ ] **README is ambiguous about where `instructions.md` lives**
  ([#171](https://github.com/caseywdunn/corpus/issues/171)).
- [x] **`tools/` is in the pyflakes gate.** Done alongside
  [#193](https://github.com/caseywdunn/corpus/issues/193), which added
  another script there. `tests/test_no_undefined_names.py` had linted
  `pipeline/`, `mcpsrv/` and `bib/` only, so operator scripts never got
  the NameError check [#75](https://github.com/caseywdunn/corpus/issues/75)
  built it for — and those are run by hand at release time, where a
  NameError costs a whole manual run rather than a fast test failure.

**External contribution**

- [ ] **[PR #144](https://github.com/caseywdunn/corpus/pull/144)** from
  @beroe — original-description linking against in-corpus works
  (423 → 505 of 598 ctenophore taxa). Mergeable and substantive, but it
  is a 4-commit external change that has never run CI, and it touches
  `bib/authority.py`, which #154 rewrote during the 1.0 cycle. It wants
  a review and a CI run, not a fast merge. Note it overlaps #175: both
  are `parse_authority`, from opposite ends.

### Features awaiting motivation

Net-new, safe to add after 1.0 without breaking the frozen surface —
held because nothing has yet made them worth the cost.

- **Container distribution image.** The only channel that can ship the
  full native toolchain in one artifact: Docker is already a
  prerequisite for Grobid, `docker-compose.yml` could bring up grobid +
  corpus together, and #153's HPC user already runs Apptainer, which
  pulls straight from a Docker registry. Costs: bind-mounting the PDF
  directory, GPU passthrough for the local VLM, and a large image with
  torch in it. 1.0's install path is verified continuously now, which was
  the precondition — worth its own issue.
- **`verify_claim`** ([#123](https://github.com/caseywdunn/corpus/issues/123)).
  Per-claim ledger anchoring as a thin similarity-only wrapper over
  `get_chunks_for_topic`. New tool — post-freeze by construction.
- **Drift detection** ([#80](https://github.com/caseywdunn/corpus/issues/80)).
  Pre-run diff of resolved config + per-input content SHAs. Overlaps
  #174 and #187; pick it up with whichever of those moves first.
- **Column-store shape for `lexicon_matrix`**
  ([#83](https://github.com/caseywdunn/corpus/issues/83)). Token saving
  at large-matrix scale; held pending a prompt-suite analysis showing it
  matters.
- **Figure-number extraction: the non-caption cases.**
  [#16](https://github.com/caseywdunn/corpus/issues/16) is **closed** —
  it landed the parsing half (`Taf. III.`, `Tab. XII.`, `Plate IV.`,
  Roman→Arabic normalization, 49 tests) of the ~538-of-1,787-paper gap.
  What remains is papers with no caption at all, or a caption not near
  its image, which needs vision OCR or a body-text-mention fallback.
  **Untracked** — file an issue if picked up. v1.2 §2 is what would size it.
- **Vision pass corpus-scale validation.**
  [#11](https://github.com/caseywdunn/corpus/issues/11) is **closed** as
  carried-out-in-code: `corpus run` invokes the vision pass whenever
  `figures.panel_detection` selects a vision backend and the host
  capability check passes. The corpus-scale run happened in v1.0 — 934 of
  934 eligible figures through the vision pass on Bouchet. What remains
  is the **figure-coverage audit**: count figures with `pass3c_status`
  set, sum `missing_figures[]` lengths, and establish what "eligible"
  excluded. v1.2 §2 subsumes the accuracy half of this on 35 documents; the
  corpus-scale count is still untracked release-validation work.
- **Evaluate Cloud Run vs the EC2+ALB stack**
  ([#89](https://github.com/caseywdunn/corpus/issues/89)). Deployment
  decision, not part of the MCP API contract.

## Out of scope (longer horizon)

- [#124](https://github.com/caseywdunn/corpus/issues/124) — whether to
  expand the server-side LLM-call surface at all. `translate_chunk` was
  removed in v0.6, leaving the MCP server a purely deterministic
  retrieval layer. Re-opening that is a direction question, not a task.
- [#14](https://github.com/caseywdunn/corpus/issues/14) — Trait
  extraction + identification keys (Q3). Substantial enough to warrant
  its own plan section when picked up.
- [#13](https://github.com/caseywdunn/corpus/issues/13) — Geographic
  mention layer. Deferred to v2.0+; the mention-layer surface is likely
  to be reworked at the major-version boundary.
- [#38](https://github.com/caseywdunn/corpus/issues/38) — Embedding
  model migration path. Design-only; implement when a model swap is
  actually needed.
- [#39](https://github.com/caseywdunn/corpus/issues/39) — MCP server
  lazy index loading. Premature at 1.8K papers; documented as a known
  scaling cliff; revisit when a corpuscle pushes 10K+.
- [#5](https://github.com/caseywdunn/corpus/issues/5) — Streamable HTTP
  transport with OAuth. Deferred indefinitely. SSE + bearer-token works
  for the ~20-collaborator deploy target.
- Multi-region failover, autoscaling, or blue/green deploys for the AWS
  served bundle. Single-instance per corpuscle is fine until it isn't.
- Authentication beyond bearer tokens / OAuth (Cognito, institutional
  SSO).
- Mirror to Cloudflare R2 / Backblaze B2 for cost. Defer until S3 egress
  shows up on a bill.
- A thin HTML/web UI on top of the MCP server. Out of scope until the
  MCP-only experience has actual non-Claude-Desktop users.

## Reference: target queries

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`.

| # | Pattern | Status entering v1.2 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer ([#13](https://github.com/caseywdunn/corpus/issues/13), deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; citation-trust gap closed by [#79](https://github.com/caseywdunn/corpus/issues/79) in v0.5, provenance preservation by [#142](https://github.com/caseywdunn/corpus/issues/142) / [#152](https://github.com/caseywdunn/corpus/issues/152) in v1.0; synthesis recipe scoped for v1.3 as the clade-monograph skill ([#179](https://github.com/caseywdunn/corpus/issues/179)) |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred ([#14](https://github.com/caseywdunn/corpus/issues/14)) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Indices in place; the corpus-scale vision run landed with v1.0's Bouchet production run (934/934 eligible figures); figure coverage becomes measurable against truth in v1.2 §2 |
| Q5 | "Summarize `<author X>`'s comments about `<author Y>`." | Indices in place |
| Q6 | "Summarize `<topic>` across the corpus." | Indices in place; cache cost addressed by dossier tools [#76](https://github.com/caseywdunn/corpus/issues/76) in v0.5 |
| Q7 | "Plot species described per decade." | Indices in place |
| Q8 | "Summarize what is known about `<anatomy>`." | Indices in place; figure-retrieval synonym blindness fixed by [#143](https://github.com/caseywdunn/corpus/issues/143) in v1.0, though translations still miss inflected forms ([#165](https://github.com/caseywdunn/corpus/issues/165), unscheduled) |

## Reference: versioning + release ritual

`__version__` in [pipeline/version.py](../pipeline/version.py) is the
single source of truth and is stamped into every persistent artifact
(bundle manifest, MCP `bundle_info`).
[CONTRIBUTING.md](../CONTRIBUTING.md) covers the branching model and
release ritual. Note **step 8** — pruning this document at release time.
It was skipped at v0.6.0, which is why this file described a finished
cycle in the present tense for two months; carried out at v1.0.0, skipped
again across v1.1.0 and v1.1.1, and carried out at the head of the v1.2
cycle instead. Twice skipped, twice caught late: do it in the release
commit, not afterward.
