# PLAN.md — Corpus pipeline (v1.4)

**Prior cycles are recorded elsewhere, not here.** A minor or major
release entry in [CHANGELOG.md](../CHANGELOG.md) opens with the cycle's
organizing theme — required going forward by
[CONTRIBUTING.md](../CONTRIBUTING.md)'s release ritual, step 2, and
present on every entry back to 0.4.0 — and each cycle's punch list is
preserved in its own tag's copy of this file
([v0.6.0](https://github.com/caseywdunn/corpus/blob/v0.6.0/dev_docs/PLAN.md),
[v1.0.0](https://github.com/caseywdunn/corpus/blob/v1.0.0/dev_docs/PLAN.md),
[v1.2.0](https://github.com/caseywdunn/corpus/blob/v1.2.0/dev_docs/PLAN.md)).
In one line each: v0.1 shipped the extraction → annotation → indexing →
MCP-serving stack; v0.2 hardened internals; v0.3 collapsed the user
surface into one CLI plus a per-corpuscle `config.yaml`; v0.4 closed
silent-failure modes and gated the next cycle on tiered CI; v0.5 was the
served-bundle quality cycle; **v0.6** froze the MCP surface at **38
tools**, which every cycle since has held; **v1.0** made a green CI badge
mean a fresh install works, validated on the full Bouchet corpus with 0
stage failures; **v1.1.1** made `CITATION.cff` schema-valid,
which is what had been silently blocking Zenodo archival since v0.3.0; and
**v1.2** built the gold transcription set and scored the extractor against
it — the first signal this project has had that measures the pipeline
against something other than itself; and **v1.3** made evidence auditable —
captions carry provenance, references are separated from canonical works, a
re-run cannot leave stale evidence, and an operator can inspect a parse
against the original page ([v1.3.0](https://github.com/caseywdunn/corpus/blob/v1.3.0/dev_docs/PLAN.md),
Zenodo [10.5281/zenodo.22647843](https://doi.org/10.5281/zenodo.22647843)).

**v1.3 (2026-09-07) is the one v1.4 grows directly out of.** Its theme
paragraph is in the CHANGELOG; what matters here is the inheritance. Evidence
is now auditable — captions carry provenance, references are separated from
canonical works, and a re-run cannot leave stale evidence behind. Two
instruments came with it and should keep being used: the gold set from v1.2,
against which any extraction change can be scored before and after, and the
build-comparison harness (`tools/qc/build_reference.py`) that made the update
contract checkable at all.

What v1.3 did *not* do is fix the defects that are wrong without saying so,
which is why v1.4 exists. It also demonstrated the lesson v1.2 left, and which
has now recurred in every cycle that looked for it: **a measurement is not a
result until you have looked at what it is measuring.** The
plate-legend fix shipped on a recall gain and quietly cost precision;
#254's OCR page-blanking was invisible to every quality gate the pipeline
had, because each one measured the document rather than the machine that
built it; and a stray file sat in the repo root through four green CI runs
and two release PRs because nothing looks there.

**v1.4 is the cleanup cycle.** v1.3 proved the pipeline can show where
evidence came from. What survives behind that is the subject here: results that
are wrong without saying so, build hazards that cost hours per run, and a
tracker that no longer describes the product. It is deliberately broad rather
than deep — the aim is that the open tracker afterwards holds only new
capability and direction questions.

The ordering against v1.5 is deliberate and was argued from v1.3's own
evidence: a skill that generates a monograph from a corpuscle with silent gaps
produces a document that reads authoritative and is wrong, which is the exact
failure class v1.3 spent itself on.

**Standing gates** — the rules that outlive any one cycle, including what
"an update is correct" means and how a release is validated — are at the
end of this document, with the other durable reference material.

Doc map unchanged: architectural background in
[OVERVIEW.md](OVERVIEW.md); per-feature history in
[CHANGELOG.md](../CHANGELOG.md); the API contract in
[API_STABILITY.md](API_STABILITY.md); HPC operations in
[BOUCHET.md](BOUCHET.md); deployment in [DEPLOY.md](../DEPLOY.md);
platform-portability criteria in
[PLATFORM_SMOKE.md](PLATFORM_SMOKE.md). Open work is tracked in
[GitHub issues](https://github.com/caseywdunn/corpus/issues).

## v1.4 — silent wrongs, operational hazards, and a clean tracker

v1.3 proved the pipeline can show where evidence came from. This cycle is
about the defects that survive *behind* that: results that are wrong without
saying so, build hazards that cost hours per run, and a tracker that no longer
describes the product. It is deliberately broad rather than deep — the goal is
to close as much of the backlog as is reasonable, so v1.5's client layer is
built on evidence it can trust instead of workarounds for known gaps.

The ordering is not arbitrary. A skill that generates a monograph from a
corpuscle with silent gaps produces a document that reads authoritative and is
wrong, which is the exact failure class v1.3 spent itself on. Shipping a
generator over it would undo the cycle's point.

### 1. Silent wrongs — results that are wrong without saying so

Highest priority, because nothing surfaces them. An operator cannot act on a
gap they cannot see.

- [x] **Route mojibake non-Latin scans to OCR without a pin**
  ([#266](https://github.com/caseywdunn/corpus/issues/266)), and populate
  `visual_script` on every detection path
  ([#172](https://github.com/caseywdunn/corpus/issues/172)). The unpinned
  Chinese paper is now OCR'd under `chi_sim` and recovers 451 Han characters
  where the old pack choice recovered none. Two findings reshaped the fix.
  First, the raster check had already stopped the `born_digital`
  misclassification the issue describes — what survived was *pack* selection,
  read off the mojibake layer, so the document was OCR'd with `eng`. Second,
  a whole population the issue does not name: a text layer of **unmappable
  glyph indices**, invisible to the gibberish score because glyph indices are
  not letters. Hunt et al. 2001 held 896 usable letters in 59,056 characters
  and was classified `clean_text_layer`; re-OCR recovers 42,398. Six such
  documents. New `ocr.unmappable_char_max`.
- [x] **Make lexicon translations match inflected forms**
  ([#165](https://github.com/caseywdunn/corpus/issues/165)). +14,111 anatomy
  mentions (+9.7%) corpus-wide, 515 documents gaining, none losing, 59
  rescued from zero. Of the two documents the issue names, though, one gains
  a single mention and the other correctly stays at zero — its text is about
  electric organs of fish. The value is on documents with partial coverage.
- [x] **Expand abbreviated genus binomials**
  ([#164](https://github.com/caseywdunn/corpus/issues/164)). +31,041 mentions
  (+15.6%), +2,237 unique taxa, 838 documents gaining, none losing, and 449
  ambiguities reported rather than guessed.
- [x] **Decide the botanical authority policy**
  ([#175](https://github.com/caseywdunn/corpus/issues/175)). Decided:
  the capability is reported unsupported, ICN authorship parses, and the
  author-only matching path is **declined on measurement** — it yields 3
  candidates from 889 Viburnum taxa and 2 of the 3 are wrong. Reasoning lives
  in `_record_authority_convention`, not here.
- [x] **Stop the CLI lying by omission**
  ([#169](https://github.com/caseywdunn/corpus/issues/169),
  [#168](https://github.com/caseywdunn/corpus/issues/168)). A filter now
  implies the listing the report's own hint promises, and a filter that
  cannot apply is named in a warning instead of dropped. The naive-chunker
  fallback is a `naive_chunker_fallback` quality gate.

**Acceptance:** for each, a document that previously produced a silent gap now
either produces the right answer or reports the gap at error or warning
severity. No fix here is complete while its failure mode is still quiet.

**Met, and the section left a method behind.** Every item here was sized
against the reference library before being built, and in four of the five the
measurement contradicted the issue's own framing — a wrong threshold, a
document whose zero was correct, a heuristic that would have written wrong
protologues, and a corroboration rule of my own that turned out circular
(OCRing a Latin page under `rus` transcribes its letters as Cyrillic
lookalikes, so the characters meant to confirm the verdict were manufactured
by the check). Bare Tesseract OSD calls **424 of 1,580** Latin-text-layer
documents non-Latin; only CJK verdicts may now override a text layer, and
only corroborated. Measure the population before writing the fix.

### 2. Operational hazards — each has already cost hours

- [x] **Make a GPU allocation fail instead of degrading to CPU**
  ([#270](https://github.com/caseywdunn/corpus/issues/270)).
  `compute.accelerator: require` and `corpus run --require-gpu`, resolved
  before any step starts. The SLURM GPU scripts pass it; the card-type pin
  stays, because this makes the failure loud rather than making an
  unsupported card work. Tested against real unusable hardware — this
  workstation's GTX 1080 is a card the pinned torch ships no kernels for.
- [x] **Let concurrent builds share a cluster**
  ([#279](https://github.com/caseywdunn/corpus/issues/279)). Per-job port
  pair derived from the job ID by one shared function, so client and server
  cannot drift. **Both Dropwizard connectors have to move** — there is an
  admin connector at 8071 as well, and overriding only the application port
  still dies, which is the trap the issue's own suggested fix would have hit.
  Verified against the real Grobid image: three instances side by side,
  byte-identical TEI from the alternate port.
- [x] **Bound build memory**
  ([#182](https://github.com/caseywdunn/corpus/issues/182)). A `docling`
  block plus `compute.num_threads` bound extraction, which is where the
  memory actually goes; `embeddings.batch_size` is reachable at last; and
  INSTALL.md documents the cgroup cap, which is the outer bound. Turned up a
  separate silent wrong on the way: `pipeline.embed` took no `--config` at
  all, so `compute.accelerator` was honoured by Stage 1 and ignored by
  Stage 2.
- [ ] **Load the local VLM in half precision on MPS**
  ([#258](https://github.com/caseywdunn/corpus/issues/258)) — float32 needs
  ~30 GB for a 7B model, which shuts Apple Silicon out entirely. **Blocked on
  hardware, deliberately not attempted.** This document already records that
  a mocked dtype-selection test does not establish that half precision is
  numerically and operationally sound on MPS, and there is no Apple Silicon
  in reach. Needs a real 7B run on a real M-series machine.
- [x] **Surface why a re-run is doing more work than expected**
  ([#80](https://github.com/caseywdunn/corpus/issues/80)). Reassessed as this
  document asked, and v1.3 had already built the computation — what remained
  was reading it. The renderer printed one line per affected document, so 699
  lines on the Viburnum corpuscle. Rolled up by reason, the case that matters
  is one line: `docling_extraction: pipeline_version — all of 699 documents`,
  which is the sentence #281 needed. Also attached to `corpus run --dry-run`.

**Section status:** four of five done. #258 is hardware-blocked, not
deferred — see its note. Two of the four fixes turned up a defect the issue
did not name, which is now the expected outcome rather than a surprise.

### 3. Cheap, and better done at a version boundary

- [x] **Rename `_serve/`** ([#273](https://github.com/caseywdunn/corpus/issues/273)).
  Now `corpus_bundle/`. An existing `_serve/` is read and updated in place, so
  no corpuscle needs rebuilding and no client breaks; the migration is one
  `mv`, on the operator's schedule. The redundant basename check is removed
  rather than renamed — the manifest was always the robust signal.
- [ ] **Stop warning about a vision downgrade on phases that never run vision**
  ([#263](https://github.com/caseywdunn/corpus/issues/263)).
- [ ] **Fix served-bundle absolute-path audit false positives**
  ([#183](https://github.com/caseywdunn/corpus/issues/183)) — noise inside a
  release gate teaches operators to ignore the gate.
- [ ] **Bound the depth-1 `get_citation_graph` payload**
  ([#166](https://github.com/caseywdunn/corpus/issues/166)) so `truncated`
  actually covers it.
- [ ] **Emit a progress heartbeat during long per-document stages**
  ([#170](https://github.com/caseywdunn/corpus/issues/170)).
- [x] **Retire `siphonophores_sample` for the 35 transcribed documents**
  ([#192](https://github.com/caseywdunn/corpus/issues/192)). Already done in
  this repo, and verified rather than assumed: `siphonophores_sample` has no
  references left, BOUCHET.md documents the gold corpuscle as the smoke test,
  and #187 already targets it (see Standing gates). BOUCHET.md's
  checksum-verified materialisation recipe was run — 35 PDFs, all sha256s
  matched — and every figure it and the issue cite re-measures exactly: 761
  pages, ~644 needing OCR, 127 MB, one document over 100pp, Totton1965a at
  314 pages carrying 49% of the OCR load.
- [x] **Column-store shape for `lexicon_matrix`**
  ([#83](https://github.com/caseywdunn/corpus/issues/83)). **Declined on
  measurement**, and the grid bounded instead. The saving is real (16.4-20.0%)
  but applies only to the opt-in `detail=True` grid, which runs 382 kB over
  1,775 rows — 19% off an undeliverable payload is not a fix. The default view
  callers use is 469-1,606 bytes. What the measurement did find is that #88
  made the grid opt-in without ever bounding it, so it now has a ceiling and
  honest counts.
- [x] **Cap figure resolution** ([#184](https://github.com/caseywdunn/corpus/issues/184)).
  Decided and shipped: a flat `figures.max_pixels_long_side`, default 3000,
  applied at build time on both write paths and available as a backfill flag.
  Re-measured on the current tree (21,521 figures / 11.88 GiB): the cap
  recovers 2.74 GiB (23%) and costs the median panel-detected figure **0%**,
  because 95% of them are already under it. Two things the measurement
  corrected. `max_dpi` cannot substitute — the byte mass is full plate pages
  at an ordinary 400 dpi, so a density cap leaves a 16,000 px figure at
  12,000. And the selective rule's advantage evaporates: capping only
  no-panel figures recovers 0.03 GiB more, while `rois == 0` turns out to mean
  "detection never ran" for 95% of figures, not "no panels" — weaker than the
  gate here assumed.
  Analyses (A), (C) and (D) from the issue are done: pngquant was present, so
  that lever is spent; the OCR path is 58% of the corpus, not a minority.
  (B) is still worth running; (E) is stale archaeology and should be dropped;
  (F) colour depth stacks on top of a cap.

**Section complete.** Two of the eight turned out to be decisions rather than
builds, and both went against the issue's proposal once measured: #83's column
store saves 19% of a payload that is undeliverable either way, and #184's
selective rule is *worse* than the flat cap it was meant to improve on. Two
others were already done and needed verifying rather than doing (#192) or
argued for deleting code rather than adding it (#273's basename check). The
recurring shape: the cheap items were cheap because the previous cycle had
already done the hard part — what was left was reading it.

### 4. Decisions to record rather than defer again

Each of these is a judgment that keeps being re-derived. Write the answer down
and close the issue, or scope the work — either is progress; leaving them open
is not.

- [x] **OCR reproducibility** ([#280](https://github.com/caseywdunn/corpus/issues/280)).
  Decided: normalize CJK whitespace **in the comparison**, keep the criterion
  exclusion, and do *not* pin `--jobs 1` — that hypothesis was tested and
  fails. See Standing gates for the measurements.
- [x] **`get_missing_references` scope**
  ([#155](https://github.com/caseywdunn/corpus/issues/155)). Both, since they
  are complementary. The cheap slice: rows with neither title nor year are
  withheld — 477 of 6,953 at the default threshold, and they outranked real
  gaps (`corpus:|unknown|`, empty-titled with 30 citations, sat 11th, above
  the genuinely-missing Bigelow 1906). And the docstring — the MCP tool
  description a client actually reads — now says the tool is best-effort, names
  the residual 96 title/year-only leads, and points at `resolve_reference` and
  the QC tool. The remaining cases need a per-block LLM pass or a similarity
  threshold loose enough to merge distinct works; neither is a cheap slice.

**Cycle acceptance:** every issue above is closed or has a recorded decision,
and the open tracker contains only new capability and direction questions —
nothing that describes the product being wrong. Issues close when their fix
lands on `dev` (CONTRIBUTING.md, "Closing issues"), so the tracker should
shrink continuously through the cycle rather than in a bulk close at release.

**Met.** All four sections are closed: 22 issues, 34 open → 13, with the
tracker holding only v1.5 skills work (#178, #179, #180, #217), direction
questions (#88, #89, #93, #123), deferred layers (#13, #14, #38, #39), and one
item blocked on hardware rather than on a decision (#258 — Apple Silicon).

Two things this cycle established, worth carrying into v1.5:

**Measure the population before writing the fix.** Nine of the twenty-two
items had their shape changed by measurement, and in six the measurement
contradicted the issue's own proposal — #83's column store saves 19% of a
payload undeliverable either way; #184's selective rule is *worse* than the
flat cap it was meant to improve; #280's `--jobs 1` remedy addresses a
mechanism that does not reproduce; #175's author-only matcher would write
wrong protologues at a 67% error rate; #266's headline symptom was already
fixed while a worse one went unnamed; and #165's two named documents are not
where its value is. A proposal in an issue is a hypothesis, including a
convincing one.

**Reassess before building.** #80, #192 and half of #155 turned out to be
already done, or already built by the previous cycle with only the reading of
it left. #273 argued for deleting a check rather than renaming it. Cheap items
were cheap because someone had already done the hard part.

Explicitly **not** in v1.4: the skills and client layer (v1.5), new extraction
layers (#13, #14), bulk export (#93), `verify_claim` (#123), embedding-model
migration (#38), MCP scaling (#39), and the direction questions (#88, #89,
#124).

## v1.5 — skills and usage

Deferred from v1.4 so it can consume a clean tracker and trustworthy
evidence rather than design around known gaps. v1.5 is the
small client/workflow layer that turns the frozen retrieval surface into a
repeatable answer. Its scope is deliberately limited to one library-building
workflow, one corpus-consuming workflow, and the shortest public path through
them:

- [ ] **A `skills/` plugin directory and library-assembly skill**
  ([#178](https://github.com/caseywdunn/corpus/issues/178)). Skills may import
  public functions from `pipeline/`; the product never imports a skill.
- [ ] **`corpus bib inspect-pages`**
  ([#217](https://github.com/caseywdunn/corpus/issues/217)), the read-only
  pre-build evidence used by the library-assembly workflow to curate
  `keeppages`, `doclang` and related judgments. It inspects the library; it does
  not duplicate the post-build page report from #274.
- [ ] **A clade-monograph skill**
  ([#179](https://github.com/caseywdunn/corpus/issues/179)) that consumes the
  caption/reference provenance shipped by v1.3 and writes deliverables on the
  client, never on the MCP host.
- [ ] **A README quick start**
  ([#180](https://github.com/caseywdunn/corpus/issues/180)) covering that path
  from a library to a served answer.

Bulk export ([#88](https://github.com/caseywdunn/corpus/issues/88) Part 2 and
[#93](https://github.com/caseywdunn/corpus/issues/93)), reconciliation changes,
new MCP tools, unrelated housekeeping and new scientific extraction layers are
not part of v1.5. If the clade-monograph acceptance run exposes another
evidence-integrity defect, fix it as a defect; do not expand the skills cycle
into another pipeline redesign.

## Unscheduled

This is a selected orientation list, not a second issue tracker; GitHub issues
are authoritative. Dependencies that matter are stated inline. Split because
the two halves get picked up for different reasons: a known defect is picked up
when it bites someone, an unbuilt feature when something makes it worth
building.

**v1.4 now owns these, and the section above is where their scheduling lives:**
#80, #83, #155, #164, #165, #166, #168, #169, #170, #172, #175, #182, #183,
#184, #192, #258, #263, #266, #270, #273, #279, #280. Their notes are kept
below because the rationale is still worth reading — but v1.4's list is the
one to work from, not this one. Anything here *not* in that list is genuinely
unscheduled: new extraction layers, bulk export, `verify_claim`,
embedding-model migration, MCP scaling, and the direction questions.

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
- [ ] **Move the docling pin forward**
  ([#98](https://github.com/caseywdunn/corpus/issues/98) follow-up).
  Still `docling==2.94.0`. Reproduce on an arm64 Mac, determine whether
  2.95/2.96 broke MPS extraction via an API change or an upstream bug,
  then advance deliberately. Needs Apple-Silicon hardware. **v1.2's fidelity harness
  ([#193](https://github.com/caseywdunn/corpus/issues/193)) gives this a
  criterion it never had** — "better or worse" against the gold
  set rather than against impressions.

**Operator surface**

- [ ] **Progress heartbeat during long per-document stages**
  ([#170](https://github.com/caseywdunn/corpus/issues/170)).
- [ ] **`--filter-gate` is silently ignored without `--list-hashes`**
  ([#169](https://github.com/caseywdunn/corpus/issues/169)) — the hint
  was fixed in 1.0, the underlying flag was not.
- [ ] **Surface the naive-chunker fallback in `corpus status`**
  ([#168](https://github.com/caseywdunn/corpus/issues/168)).
- [ ] **Do not warn about a downgraded vision pass on phases that never run
  it** ([#263](https://github.com/caseywdunn/corpus/issues/263)).
- [ ] **Bound build memory explicitly**
  ([#182](https://github.com/caseywdunn/corpus/issues/182)); embedding batch
  size and docling resource controls should be reachable from configuration.
- [ ] **Make accelerator requirements explicit**
  ([#270](https://github.com/caseywdunn/corpus/issues/270)). A GPU phase should
  fail before expensive setup when no usable accelerator is present; keep
  scheduler-specific constraints in deployment configuration.

**Housekeeping**

- [ ] **The local VLM loads in float32 on MPS**
  ([#258](https://github.com/caseywdunn/corpus/issues/258)). This is not part of
  the bounded v1.3 hardening tranche: validate the actual 7B model on suitable
  Apple Silicon, together with any move beyond the current docling pin. A
  mocked dtype-selection test alone does not establish that half precision is
  numerically and operationally sound on MPS.
- [x] **CI now looks at the repo root.** Shipped in v1.2.1 as
  `tests/test_repo_root_is_clean.py`, an allowlist over `git ls-files`.
  A stray 9-byte `%PDF-1.4` fragment named `6` sat next to `README.md`
  through four green CI runs and a full release PR, caught by eye at the
  v1.2.0 tag boundary — one merge from a citable Zenodo archive. T0 lints
  `pipeline/`, `mcpsrv/`, `bib/` and `tools/` for undefined names; nothing
  had an opinion about the top level. Note the file was *not* a shell typo
  as first assumed: the test suite regenerated it on every run
  ([#257](https://github.com/caseywdunn/corpus/issues/257)), which is why
  deleting it once did not hold. Same shape as the `tools/` pyflakes gap
  below — the check that would have caught it did not exist because nobody
  had been bitten yet.
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
  Pre-run explanation of why a run will invalidate each stage. v1.3 implements
  #174 and #187; reassess afterward whether #80 still names distinct work.
- **Bulk export outside the MCP response channel**
  ([#88](https://github.com/caseywdunn/corpus/issues/88) Part 2 and
  [#93](https://github.com/caseywdunn/corpus/issues/93)). A local `corpus
  export` or client-side download workflow fits the execution-plane contract;
  an MCP tool that writes the server's filesystem does not.
- **Rename `_serve/`**
  ([#273](https://github.com/caseywdunn/corpus/issues/273)). Naming cleanup,
  not evidence integrity; provide a deprecation/compatibility path because the
  directory appears in operator scripts and documentation.
- **Column-store shape for `lexicon_matrix`**
  ([#83](https://github.com/caseywdunn/corpus/issues/83)). Token saving
  at large-matrix scale; held pending a prompt-suite analysis showing it
  matters.
- **Figure-number extraction: the non-caption cases.**
  [#16](https://github.com/caseywdunn/corpus/issues/16) is **closed** —
  it landed the parsing half (`Taf. III.`, `Tab. XII.`, `Plate IV.`,
  Roman→Arabic normalization and fixture-backed tests) of the corpus-wide gap.
  What remains is papers with no caption at all, or a caption not near
  its image, which needs vision OCR or a body-text-mention fallback.
  **Untracked** — file an issue if picked up. v1.2's figure-fidelity
  scoring ([#194](https://github.com/caseywdunn/corpus/issues/194)) is what
  would size it.
- **Vision pass corpus-scale validation.**
  [#11](https://github.com/caseywdunn/corpus/issues/11) is **closed** as
  carried-out-in-code: `corpus run` invokes the vision pass whenever
  `figures.panel_detection` selects a vision backend and the host
  capability check passes. A corpus-scale run happened in v1.0 with every
  eligible figure reaching the vision pass. What remains
  is the **figure-coverage audit**: count figures with `pass3c_status`
  set, sum `missing_figures[]` lengths, and establish what "eligible"
  excluded. v1.2's figure-detection scoring
  ([#194](https://github.com/caseywdunn/corpus/issues/194)) subsumes the
  accuracy half of this on 35 documents; the
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
  lazy index loading. Premature for the current ~2,000-paper reference corpus;
  documented as a known
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

## Standing gates

### What "an update is correct" means

Established in v1.3 and permanent. Do not re-derive it, and do not restate it
as byte equality — that was the old form, and it is not a property this
pipeline has or should chase.

> For every supported change class, an incremental run **re-processes exactly
> the documents whose fingerprinted inputs changed and no others**; documents
> it does not touch are left unchanged; and the resulting current document set,
> cross-paper mappings and vector rows **match a clean rebuild**.

Equality here is **semantic, not byte-for-byte**. Provenance that is not
reproducible by construction is excluded from the comparison, and the exclusion
list is explicit rather than assumed: the taxonomy snapshot's file hash (#278,
since fixed), absolute build paths, local-VLM ROI coordinates, and CJK OCR
whitespace segmentation (#280). Shortening that list is real work; ignoring it
quietly is not.

**On #280 specifically, the decision is made and the remedy is not `--jobs 1`.**
The harness normalizes whitespace *between two CJK characters* before digesting
(`tools/qc/build_reference.normalize_cjk_spacing`), so segmentation noise stops
burying real differences — and it does that in the comparison only, never in the
artifact or the stage fingerprint, because a fingerprint decides what re-runs.
The `jobs=12` mechanism the issue proposed was tested and does not hold: both
affected documents are byte-identical across repeated OCR runs at `--jobs` 1, 4
and 12, and across `OMP_THREAD_LIMIT` 1, 4, 12 and unset, with docling
deterministic on fixed input. Pinning `--jobs 1` would cost real build time for
nothing. The underlying nondeterminism is real but its mechanism is still
unidentified and does not reproduce on a workstation, so the exclusion stands
rather than being traded for a guess.

Two method notes that cost real time to learn, and will again:

- **Incrementals must run in place.** A build tree embeds absolute paths in
  `figures.json`, `summary.json` and `pipeline_state.json`; only `_serve/` is
  scrubbed (#70). A copy relocated to another path reports every document as
  changed, and the entire diff is path noise.
- **One build per machine, or per node.** Two concurrent local builds blew the
  per-page OCR timeout on 3-4 documents each and silently degraded the clean
  side of the comparison. The quality gates did catch it, loudly.

### Validate a release against the full corpus, not only the gold set

v1.3's release validation rebuilt all 1775 siphonophore documents, and that run
found four defects no unit test could reach: #278 (the taxonomy fingerprint
hashed a file containing timestamps), #279 (concurrent Grobid jobs collide on a
fixed port), #280 (CJK OCR whitespace is not reproducible), and #281 (the GPU
vision phase re-extracted every document, turning 1h27m into a projected 35h
that could not finish in one allocation). The gold set had passed cleanly
beforehand. Scale is what surfaced them, so budget a full rebuild before a
release rather than treating the 35-document set as sufficient.

**T3-bare** stays waivable under its own "where platform behavior changed"
clause — v1.3 recorded the first waiver — but re-run it at the next release
that touches apt packages, the miniforge bootstrap, or a runtime dependency
pin. None of those moved over `v1.2.1..HEAD`, which is why the waiver held.

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
targets the **gold corpuscle** (built in v1.2) — big enough for a rate to mean
something, small enough to rebuild per release. That comparison was run
by hand against a viburnum rebuild during v1.1 and is what proved the
3.4 GB → 2.3 GB drop was #184's re-encoding rather than lost content.

The `--deselect` question that paragraph raised is answered: **none**.
[#167](https://github.com/caseywdunn/corpus/issues/167) removed all three
flags, T1 now runs bare `-m corpus_required`, and T2 additionally ignores
`test_reference_extraction.py` because Grobid is disabled there. Both
workflow files state the reasoning inline.

## Reference: target queries

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`.

| # | Pattern | Status entering v1.3 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer ([#13](https://github.com/caseywdunn/corpus/issues/13), deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; v1.3 repairs caption/reference evidence and makes uncertainty inspectable before v1.4 adds the clade-monograph skill ([#179](https://github.com/caseywdunn/corpus/issues/179)) |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred ([#14](https://github.com/caseywdunn/corpus/issues/14)) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Indices in place; a corpus-scale vision run landed in v1.0 and figure detection became measurable against truth in v1.2. On the clean gold corpuscle, physical detection is 0.883 recall / 0.865 precision raw and 0.827 / 1.000 on the default MCP type surface. Caption identity binding is 0.641 / 0.987 against the corrected 839-identity yardstick; the remaining recall gap is predominantly absent upstream number evidence rather than selector error |
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
commit, not afterward. v1.2.0 did exactly that and it worked.

**What v1.2.0 got wrong instead was the far end of the ritual.** Its
`release:` commit (`e90ad09`) did steps 2 and 8 — version string, dated
CHANGELOG entry with its theme, INSTALL.md pin, this file pruned — and then
steps 3–5 never ran. `main` sat on v1.1.1 and no tag existed for a full day,
while every file in the tree said 1.2.0 and the INSTALL.md pin named a tag
that did not resolve. A commit titled `release: vX.Y.Z` reads like a release
and is not one. **The only evidence a version shipped is
`git ls-remote --tags origin` and `gh release list`** — check those, not the
log, before assuming a version is out. Two things were then caught in the
gap: #254, which belonged in 1.2.0's section rather than `## [Unreleased]`,
and a stray file in the repo root that would otherwise have been archived to
Zenodo permanently.
