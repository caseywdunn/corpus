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

## v1.2 — extraction fidelity, measured against ground truth

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

**A defect the harness found on its way in**, and the reason to validate
a measuring instrument against a second source before trusting it: a
docling table item has no `text` field at all — its content lives only in
`data.table_cells[].text` — so a per-page walk reading `text` discarded
every table on the page. Checking reconstruction against `text.json`
rather than assuming it caught 8 of 35 documents short, the worst by 26%
of its tokens. Unfixed, those pages would have been reported as the
pipeline losing prose it had in fact extracted. All 35 now reconstruct
≥95.8% of `text.json`'s tokens.

### 2. Figure and caption fidelity, **to file**

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

### 3. Grobid output and reference consolidation, **to file**

Three questions that sound like one, and the code gives them different
answers:

- **`consolidateCitations` has never run.**
  `pipeline/grobid_client.py` defaults `consolidate_header=1` and
  `consolidate_citations=0`, `pipeline/metadata.py` calls
  `process_fulltext()` overriding neither, and the `grobid:` block in
  `pipeline/config.template.yaml` exposes only `url` and `disable`. So
  reference consolidation against CrossRef is off and always has been.
  Measure what enabling it costs (a network round-trip per reference) and
  buys, then expose both flags in config with the measured default.
- **In-corpus references.** Do references to papers that are themselves
  in the corpuscle resolve to the canonical work? `references_match_corpus_papers`
  is flagged in `tests/test_corpus_wide.py` as "the outlier, deliberately
  set" — that looseness is a measurement of something, and the gold
  reference sections can now say whether it is justified.
- **Out-of-corpus references** — this is
  [#155](https://github.com/caseywdunn/corpus/issues/155).
  `get_missing_references` is dominated by unreconciled citation-string
  variants of works already in the corpus; `resolve_reference "Bigelow
  1911 The Siphonophorae"` returns 53 matches, one canonical and ~50
  variants. Full reconciliation is a cycle of its own (DOI normalization,
  block-and-cluster canonicalization, alternate-DOI aliasing, junk
  filtering, probably an auditable LLM adjudication pass) and stays
  unscheduled. **v1.2 takes the cheap slice**: collapse candidates by
  normalized (author, year, title) before ranking, so one work stops
  appearing as N rows. Makes the tool honest without pretending to fix
  reconciliation, and it is additive under the freeze.

### 4. Per-document bib directives

`ocrlang` (#176) established the mechanism: a flat BibTeX field on the
entry whose `file = {…}` matches the PDF, parsed in `bib/parser.py`,
carried through `bib/importer.py` and `bib/export.py`, documented in the
README's directive table. Two more fields follow the same path.

- [ ] **Per-document page-range selection**
  ([#188](https://github.com/caseywdunn/corpus/issues/188)) — scanner
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

- [ ] **The gold 35 become the smoke-test corpuscle.**
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
  ([#187](https://github.com/caseywdunn/corpus/issues/187)), retargeted
  from `siphonophore_sample_YYYYMMDD` to `siphonophore_gold_YYYYMMDD`.
  Complements A rather than duplicating it: **A measures accuracy against
  truth, #187 measures drift between builds.**
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
