# PLAN.md — v1.5 correctness; v1.6 skills

## Release decision and preserved work

The September interface audit exposed errors in bibliographic identity,
scientific text, figures and bounded retrieval. Repairing that evidence takes
priority over the skills that consume it. **v1.5.0 is the correctness release;
the skills release moves to v1.6.0.** Additional response metadata and optional
pagination require a minor release under [API_STABILITY.md](API_STABILITY.md).
Existing fields, defaults and successful calls remain compatible; fixes must
not silently redesign the frozen MCP surface.

The enhancement work through `44977fea1657985f075f7e3bafc7833d2ab914cd` is
preserved on **`enhancements/v1.6`**, and `dev` remains at that history.
**`release/v1.5-correctness` starts from the published `v1.4.0` tag.** Work in
isolated issue branches/worktrees, integrate reviewed fixes into the correctness
branch, and carry them forward to the enhancement branch before v1.6 resumes.
No reset or force-push of `dev` is part of this plan. The release proposal will
be a PR from the correctness branch to `main`, subject to the usual gates.

The release branch is the temporary integration target for this cycle. Apply
CONTRIBUTING's issue-branch, evidence and closure rules to that target; distinguish
implemented fixes from remaining source/corpus acceptance in both the tracker
and the checklist. Do not close an issue on a partial fix. GitHub is the
issue-status authority; this file records sequence and acceptance, not a claim
that every checked code change has shipped.

### Execution planes

Each item names the plane whose data it writes: **build** for extraction,
annotation, reconciliation and bundled evidence; **serve** for bounded,
read-only response formatting/filtering; **library-curation** for input edits;
**client** for synthesis. Infrastructure and documentation have no plane.
Serve corrections below must not introduce OCR, reconciliation, general LLM
calls, external enrichment or writes to the immutable bundle. Curator changes
belong in source inputs and flow through a rebuild, never artifact patches.

## v1.5 — ordered correctness work

The ordering is by dependency and consequence, not issue number. Independent
tracks can run concurrently. Read each issue's full acceptance criteria before
implementation. The audited bundle was built from `1.4.0.dev0`, so a historical
symptom is not by itself a fresh reproduction on this release baseline.
Separate source-verified failures, current-code reproductions and unproven
mechanisms. Detector hits are review populations, not error counts.

### 1. Establish the release baseline and regression cases

- [x] Preserve the enhancement branch and create the correctness branch from
  v1.4.0. Bring forward the existing root-replacement warning and honest help
  text for **#298** — `[plane:build]`; this does not add multi-root ingestion.
- [x] **#285** — remove the unused Chrome repository from the affected CI
  runners before apt update; verify the actual runner path. Both workflow
  steps now identify/remove existing Chrome `.list` or `.sources` files and
  retry transient downloads. Hosted integration run `35375772899` passed and
  confirms removal of `/etc/apt/sources.list.d/google-chrome.sources`;
  clean-room run `35375772962` and T0 run `35375772896` also passed. *No plane.*
- [ ] Inventory the audit cases against the existing gold source manifest and
  add explicit expectations for the failure mechanisms. *No plane; validation.*
  Reuse the existing 35-document siphonophore gold corpuscle and independent
  transcriptions. Chun1898b and Hosia2024 are already present. A document's
  presence does not establish that its citation identity or clearance is tested.
  The hash-to-manifest inventory and remaining source gates are recorded in
  [V1_5_SOURCE_ACCEPTANCE.md](V1_5_SOURCE_ACCEPTANCE.md); the remaining extraction
  families still need their issue-specific expectations.

### 2. Restore bibliographic identity and authoritative fields

- [ ] **#296** — preserve build-time BibTeX origin, ordered authors and canonical
  field precedence through materialization, reconciliation and import/export.
  Implementation and source-derived bibliography regressions landed; the
  named production merge-history review remains pending.
  `[plane:build]`
- [ ] **#299** — reject incompatible reconciliation of a curated document onto
  another same-author/year work; repair identity and graph membership, not just
  the displayed title. Depends on retaining #296's evidence. Implemented with
  source-derived Chun controls; rebuilt-corpus mappings remain to verify.
  `[plane:build]`
- [ ] **#300** — distinguish publication parts from shared book DOIs and short
  key collisions while preserving genuine duplicate scans. Define deterministic
  identity and rebuild/migration behavior before changing matching. `[plane:build]`
  Implemented with all seven source Delle Chiaje entries, Moore records,
  duplicate controls and clean/incremental checks; full-corpus replay remains.
- [ ] **#301** — preserve volume, issue, pages and article locators through the
  complete parser → authority → bundle → formatter round trip. Source-derived
  round trips pass; regenerated bundle acceptance remains. `[plane:build]`
- [x] **#310** — share citation-query author parsing so additional correct
  author text cannot turn a match into false absence. Title-only search is
  outside this fix. Named examples, ambiguity/absence controls and twelve
  additional real-work queries pass against the retained v1.2.1 authority;
  the rebuilt candidate's audit replay remains a release gate. `[plane:serve]`
- [ ] **#311** — parse comma-separated taxonomic authors and expose supported
  original-description candidates without asserting an ambiguous match.
  Implemented with the source-supported Apolemia example and ambiguous
  Physalia/Nectopyramis controls. Candidate evidence stays separate from
  curator-reviewed originals; rebuilt taxonomy-link acceptance remains.
  `[plane:build]` (bounded candidate presentation: `[plane:serve]`)
- [ ] **#313**, then **#314** — preserve raw reference observations while
  distinguishing parse debris and resolving source-supported publication-year
  conflicts before ranking missing works. Depends on trustworthy identities;
  do not use broad title-only merges. Implemented with conservative fragment
  quarantine, visible uncertainty, source-supported per-observation year
  adjudication and shared-edge deduplication. Raw evidence remains unchanged;
  rebuilt-corpus acceptance remains. `[plane:build]`

Acceptance includes reversed ingestion/merge order, curated/extracted conflicts,
Unicode and ordered authors, duplicate-scan controls, distinct parts sharing an
identifier, fresh builds and unchanged refreshes. All bibliography and graph
routes must agree on the requested document's identity. Keep observations and
decision provenance; document how old bundles obtain repaired mappings.

### 3. Restore source text and citation evidence

- [ ] **#306**, **#312**, **#315**, **#316** — reproduce and distinguish Big5-like
  text-layer corruption, other mojibake, surname substitutions and spacing
  diacritics/truncated names. Use source images and correct negative controls;
  no global replacement or guessed correction is sufficient. #306 and #316
  now have same-region byte corroboration and geometric accent recovery,
  source replay, idempotence and conflicting-number controls. Incomplete
  surnames remain reviewable without guessed aliases or reverse-order merges.
  #312 now retains original damaged-layer evidence before OCR and requires
  agreeing regional raster readings plus prepared-word geometry. A real
  extraction/materialization/chunk/query replay recovers the German example;
  eight other suspect tokens remain unresolved. Producer changes reprepare
  dependent artifacts and retire obsolete receipts. #315 uses curated
  author/year candidates and unhinted regional OCR agreement; a source sample
  admits ten correct repairs and leaves four unresolved. The complete current
  catalog reaches extraction and figure resets, with both resume gates and
  audits tracking consumed evidence. Actual Mapstone extraction/chunk/query
  replay recovers both tested names; fresh full-corpus citation mapping and
  broader rebuilt-source acceptance remain pending. `[plane:build]`
- [ ] **#303** — preserve scientific signs, units and exponents through stored
  text, chunks and embedding input. Word coverage alone cannot validate these
  semantics because its normalization removes punctuation. Source-glyph and
  raster-supported repairs are implemented, with the falsely encoded Kidwai
  dash retained as a negative control; integrated release replay remains.
  `[plane:build]`
- [ ] **#304**, then **#319** — preserve multi-column reading order and enclosing
  species context, and expose diagnosis passages through a documented route.
  Wrong element order and missing headings need separate checks. Source order,
  treatment propagation and context-preserving chunk boundaries are integrated.
  Ordinary dated headings cannot create species treatments; global Results
  headings clear the prior context. Served provenance has explicit per-row and
  shared optional-evidence limits. Complete rebuilt acceptance remains.
  `[plane:build]`
- [ ] **#307**, **#308**, **#334** — preserve table/key cells and branch destinations,
  remove artificial repeated cell text, and recover source-supported word
  boundaries without splitting legitimate compounds. Logical-cell/key
  serialization and source-supported Hissmann/DuClos spacing are implemented.
  The source-printed Mapstone p200 typography is preserved; Daniel's disputed
  spelling remains explicit. Rebuilt retrieval acceptance remains. `[plane:build]`
- [x] **#309**, **#317** — validate complete citation spans and preserve the
  source paragraph across grouped citations; trace original TEI as well as
  extracted artifacts. Link validation and text preservation are separate
  acceptance boundaries. Fresh source TEI reproduces the named defects;
  coordinate-backed repair and complete-span regression checks pass while
  ambiguous/unmatched targets retain explicit status and raw observations.
  Production/gold replay remains a release gate. `[plane:build]`

Corrected text must invalidate downstream chunks, annotations, references and
embeddings as appropriate. Source-grounded assertions must cover relationships
and scientific meaning, not only token presence. Run targeted current builds
before claiming an old OCR-routing defect persists.

### 4. Restore figure content, captions and clearance

- [ ] **#305** — record and transform ROI coordinate frames correctly. Fix the
  deterministic Claude resize bug, investigate the deployed Qwen failures
  separately, and validate actual panel content after rebuild. Bounds checks
  and the Claude fix alone do not close the issue. Both backend frame fixes
  and provenance are implemented; cached Qwen processor dimensions corroborate
  the double-resize mechanism. The source-pilot capture tool now preserves actual
  processor frames, raw responses and production crops for review. Fresh
  source-pilot inference and scientific panel review remain.
  `[plane:build]`
- [x] **#324**, **#322**, **#329** — parse panels beyond L with specific
  descriptions, preserve caption fragments and incomplete-binding evidence,
  and retain the edge species label in the source-verified figure. Source
  geometry replay passes for Hosia, Sutherland and the clipped Erenna figure;
  integrated rebuilt-corpus validation remains the release gate. `[plane:build]`
- [ ] **#302** — materialize figure-specific rights exclusions and provenance;
  apply them before inherited publication clearance at every strict delivery
  boundary. Whole figures, panels, fallbacks, URLs and HTTP must agree.
  Implemented; #322 now recovers the actual Hosia exclusion from source
  caption fragments. Rebuilt gold/bundle replay remains pending.
  `[plane:build]` (enforcement: `[plane:serve]`)
- [x] **#321** — caption matches precede paper-only mentions in both taxon
  figure routes, with deterministic ties and preserved legacy scores.
  Boundary/tie regressions pass; integrated corpus replay remains in the
  release acceptance gate. `[plane:serve]`
- [ ] **#323** — match unambiguous caption abbreviations without inventing
  associations. Post-build caption links and evidence are now materialized from
  final figures, chunks and taxonomy with content receipts. Both figure routes
  use those links; legacy bundles expose unavailable provenance. Source Hosia
  checks pass; rebuilt bundle acceptance remains. `[plane:build]`
- [x] **#327** — expose consistent structured refusal reasons while preserving
  successful MCP image responses. Actual MCP result conversion is covered by
  the regression suite; integrated audit replay remains a release gate.
  `[plane:serve]`
- [ ] **#332** — recompute record totals after expansion/removal and validate
  them before bundling; consistently read legacy artifacts. Implemented and
  integration-tested; full rebuilt count concordance remains. `[plane:build]`

Acceptance uses source-reviewed complete panels, labels and scale context,
permitted/restricted/unknown rights controls, mixed-panel inheritance, and
append/remove/count consistency. No query-time interpretation of captions as
new licensing evidence and no cache writes into the bundle.

### 5. Make bounded query results interpretable

- [x] **#318** — expose aggregate scope and selected/available paper counts
  without silently changing existing aggregate meanings. `[plane:serve]`
- [x] **#325** — add deterministic excerpt pagination, available/returned counts,
  truthful truncation and response-byte bounds; distinguish markers from unique
  paragraphs and retrieve every eligible row without gaps. `[plane:serve]`
- [x] **#326** — share one total graph-edge budget across both directions and
  return a structured invalid-argument error for an unsupported direction.
  `[plane:serve]`
- [x] **#328** — select default top lexicon terms after paper/year filtering;
  retain caller-specified term order and deterministic ties. `[plane:serve]`
- [x] **#331** — correct unavailable tool names, parameters and stale examples.
  *No plane; documentation.*
- [ ] **#320** — after text/context corrections, add a bounded independent
  diagnostic/key retrieval evaluation and measure repetitive-table crowding.
  Fix supported crowding mechanisms against that evaluation; broader ranking
  experiments require explicit follow-up scope. `[plane:build]` for stored
  retrieval units; bounded ranking, if warranted: `[plane:serve]`.
  The [source-graded evaluator](RETRIEVAL_EVALUATION.md), fixed audit calls,
  controls, deterministic independent sampler and pre-tuning acceptance targets
  are implemented. Prose controls now have independently reviewed source labels.
  The rebuilt independent population
  and actual reference/candidate captures remain; no retrieval improvement is
  claimed from the evaluator's unit tests.

### 6. Release acceptance and return to enhancements

- [ ] Relevant regressions, contract tests and established CI lanes pass on the
  integrated candidate. Additions follow the minor API policy; existing
  successful calls and fields remain compatible.
- [ ] Rebuild and score the existing gold corpuscle; review build-reference
  differences and source-based expectations. Use small synthetic cases for
  logic/boundaries and add source material only for a demonstrated coverage gap.
  **Hydractinia is not a second fixture corpus.**
- [ ] Demonstrate clean/incremental semantic equivalence with the standing
  exclusions, artifact invalidation and unchanged-document checks below.
- [ ] Run the full reference-corpus release validation, then replay the audit
  calls through the served candidate. Gold-only checks cannot prove corpus-scale
  citation degree, concurrency or reconciliation behavior.
- [ ] Record repaired cases and remaining investigations accurately, publish
  rebuild/migration instructions, and prepare the release PR. The
  [candidate upgrade procedure](V1_5_MIGRATION.md) and draft PR #335 are ready
  for review; release acceptance remains pending. Production
  replacement requires the rebuilt bundle; server-only updates cannot repair
  already stored evidence.
- [ ] Bring correctness commits into the preserved enhancement history before
  resuming v1.6. The goal is the accepted fixes and verified gates, not an empty
  tracker or an indefinite extractor redesign.

## v1.6 — preserved skills and usage

The implementation through `44977fe` remains on `enhancements/v1.6`. It contains
assemble-library and its helper scripts, plugin scaffolding, configuration path
expansion and associated documentation/tests. Other workflows below remain
unfinished. Resume them after correctness acceptance; validate plugin installation
as well as workflow output. **#297** (the plugin's accidental development MCP
configuration) stays here, because the plugin is not in the correctness release.
*No plane; packaging.*

- [ ] **A `skills/` plugin directory and library-assembly skill** — *library
  curation* ([#178](https://github.com/caseywdunn/corpus/issues/178)). Skills may
  import public functions from `pipeline/`; the product never imports a skill.
  The plugin scaffolding itself has no plane; `assemble-library` produces a
  library, which is what the tag tracks.

  **Keep the plugin installable on its own.** A plugin install needs neither the
  `corpus` package nor an MCP server, and that independence is what lets #286 and
  #179 reach someone who has access to a served corpuscle and has never installed
  Python. It is easy to break by accident — one shared helper that imports
  `pipeline` and the plugin needs a conda environment. The three runtime profiles
  are tabulated in the issue; document the matrix, not a single install sequence.
- [ ] **`corpus bib inspect-pages`** — *library curation*
  ([#217](https://github.com/caseywdunn/corpus/issues/217)), the read-only
  pre-build evidence used by the library-assembly workflow to curate
  `keeppages`, `doclang` and related judgments. It inspects the library; it does
  not duplicate the post-build page report from #274.
- [ ] **A corpuscle-summary skill** — *client/agent*
  ([#286](https://github.com/caseywdunn/corpus/issues/286)). Markdown by default,
  a LaTeX fragment on request. Broken out of #179, where it was Appendix A of the
  monograph. It is two MCP calls, a plot script and a template — `corpus_summary`
  and `bundle_info` already return every field it needs, so it adds no MCP
  surface.
- [ ] **A build-and-triage skill** — *build/materialization*
  ([#287](https://github.com/caseywdunn/corpus/issues/287)). The judgment layer
  over `corpus check` / `run` / `status`: what a `timeout` versus a `corrupted`
  versus a `quality_gate` failure means for this collection. It must not
  re-implement `corpus status`; if it ends up only printing that output, cut it.
  It also emits the `.mcp.json` entry for the corpuscle it just built — the path
  is known here and nowhere else, the plugin cannot ship one, and it is the join
  that lets #288 chain a build into a summary.
- [ ] **A clade-monograph skill** — *client/agent*
  ([#179](https://github.com/caseywdunn/corpus/issues/179)) that consumes the
  caption/reference provenance shipped by v1.3 and writes deliverables on the
  client, never on the MCP host. Appendix A now comes from #286 rather than being
  written here.
- [ ] **A quickstart orchestrator skill** — *client/agent*
  ([#288](https://github.com/caseywdunn/corpus/issues/288)), sequencing the three
  skills above with a gate per step. Its value over a README prompt is that a
  gate becomes an exit code instead of a sentence a model can rationalize past.
- [ ] **A README quick start** — *no plane; documentation*
  ([#180](https://github.com/caseywdunn/corpus/issues/180)) covering that path
  from a library to a served answer. **It now ends at the corpuscle summary, not
  a monograph** — a quick start should not build a book, and stopping at the
  summary takes LaTeX out of a first-time reader's path entirely.

## Unscheduled — closing the improvement loop

Not scheduled; recorded together because they are one body of work and the
ordering inside it matters. The skills cycle builds the forward path — library to build to
bundle to answer. This is the **return edge**: reading a built corpuscle to
find what its inputs got wrong, and proposing the reviewed edit that makes the
next build better. See [OVERVIEW.md](OVERVIEW.md#execution-planes-and-data-ownership),
which now documents the loop, its iterative nature, and the fact that its ends
usually sit on different machines.

**The finding that shapes this cycle: the return edge is missing its
instruments, not its skills.** Two of the three workflows below are blocked on
measurement that does not exist, and building the measurement is `pipeline/`
work, not skill work, under the AGENTS.md tiering. Once it exists each skill is
thin — read the report, propose the edit, re-run. Only the bibliography lap is
ready today, because `get_missing_references` is the one instrument already
built.

Instruments first:

- [ ] **Record taxon-name candidates the snapshot does not have** — *build*
  ([#289](https://github.com/caseywdunn/corpus/issues/289)). `pipeline/taxa.py`
  drops every unresolved candidate silently, so the corpus cannot answer "which
  names does the literature use that my taxonomy lacks" — it computes the answer
  on every run and discards it. Measure the noise floor before choosing what to
  keep.
- [ ] **Lexicon coverage and content validation** — *build*
  ([#290](https://github.com/caseywdunn/corpus/issues/290)). A declared term
  that matched nothing is invisible: `n_terms_hit` has no denominator, and the
  bundle does not ship `lexicon.yaml` for the server to diff against. Deferred
  once already at [QC.md](QC.md) "What's missing in v0.3"; this is that pass.
  Dead terms only — gap detection is a different problem.

Then the laps:

- [ ] **Find papers to add** — *library curation*
  ([#293](https://github.com/caseywdunn/corpus/issues/293)). The ready one.
  Its trap is that `in_corpus = 0` depends on ghost reconciliation, so an
  unreconciled variant of a work already held appears as a lead — verification
  is the skill, not garnish.
- [ ] **DwC-A audit** — *library curation*
  ([#294](https://github.com/caseywdunn/corpus/issues/294)), blocked on #289.
  Reports discrepancies; never asserts a synonymy. The snapshot cannot even
  justify its own, having dropped `nomenclaturalStatus` at ingest.
- [ ] **Build and revise a lexicon** — *library curation*
  ([#295](https://github.com/caseywdunn/corpus/issues/295)), blocked on #290.
  The cheapest lap available, which makes it the best place to show the loop is
  worth taking.

Adjacent, found while scoping the above:

- [ ] **Vernacular names are exported but never ingested** — *build*
  ([#291](https://github.com/caseywdunn/corpus/issues/291)). The export selects
  a `name_type` no ingest path writes, so the extension is always empty. Fixing
  it properly widens `name_set()` and would make vernaculars match as taxon
  mentions — a behavior change, and #178's "kalina" homonym trap is the reason
  to think hard about it.
- [ ] **Gate lexicon annotation per category** — *build*
  ([#292](https://github.com/caseywdunn/corpus/issues/292)). Editing one
  category re-runs all of them plus `taxa.json`; the docs claimed otherwise and
  were corrected. Per-lap cost is what decides whether users keep going round,
  and #178 generates multi-category lexicons by default.

## Other deferred work

- **#330** — optional lexicon surface-form dossier support. `[plane:serve]`
- **#333** — optional taxonomy miss diagnostics and matched-name presentation;
  preserve truthful snapshot scope. `[plane:serve]`
- **#283** — dependency upgrade, held unless a confirmed defect requires it and
  the platform/fidelity checks justify it. *No plane; dependencies.*
- **#88 / #93**, **#123** — bulk export and new claim-verification capability;
  keep separate from repairing existing behavior. `[plane:client]`
- **PR #144** — review any applicable original-description fixes against #311
  and the new identity controls; do not merge an old unvalidated reconciliation
  policy as a shortcut. `[plane:build]`

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

**The local-VLM exclusion is now measured, not assumed.** v1.4 ran the #258
probe twice on one H200 at production `bfloat16`, same weights, same images,
greedy decoding: ROI *counts* were identical on 4/4 figures, but only 3/4 held
boxes within IoU 0.95, one landing at 0.899. So the geometry is not reproducible
run to run on a fixed accelerator, not merely across accelerators. That is the
mechanism behind the v1.4 rebuild's figure churn — 40 documents gained ROIs and
42 lost them against the v1.3.0 build, a symmetry that is the signature of noise
rather than a regression. The cause is greedy decoding over coordinate *digit*
tokens: a reduction-order difference too small to matter flips a digit and a
coordinate moves. Do not treat small ROI deltas between builds as a defect
without first reproducing them.

**Not covered by that measurement: discovery-mode figure *counts*.** The probe
exercises panel detection against a known label set, where counts held. Bare-plate
discovery asks the model how many figures it can see, and on `AgassizL1862ab`
that went 56 -> 35 between the two builds with byte-identical input images.
Whether that is the same nondeterminism amplified by an open-ended count, or
something else, is unmeasured — do not assume it is noise on the strength of the
panel-detection result.

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

v1.3's release validation rebuilt the full siphonophore corpus, and that run
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
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; v1.3 repairs caption/reference evidence and makes uncertainty inspectable before v1.6 adds the clade-monograph skill ([#179](https://github.com/caseywdunn/corpus/issues/179)) |
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
