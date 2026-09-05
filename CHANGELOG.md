# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Reference evidence is now independent of canonical works (#240,
  deterministic core).** Every bibliography occurrence is retained as a
  content-addressed, append-only observation; replaceable source-set pointers
  identify the current evidence without deleting superseded raw or parsed
  citations. A separate observation-to-work relation records the match method,
  score and rule-producer version. When the current observation set changes,
  the complete derived mapping is rebuilt in stable order, so an incremental
  paper addition converges with a clean build and an unchanged rerun leaves
  the evidence and mapping rows untouched. Corpus-paper reconciliation records
  its own scored decision, and maintenance merges redirect observation
  mappings before removing a canonical duplicate. The existing `citations`
  table remains the frozen MCP compatibility materialization; no tool name,
  input or response shape moves.

- **Build regression references expose equal-count evidence changes (#187).**
  The operator-side snapshot/comparison tool fingerprints primary JSON content,
  reports reference-field differences, separates unchanged warnings from hard
  failures, and records manifest facts and diagnostic PDF/TEI hashes. It never
  overwrites an existing reference or automatically accepts a changed build.

- **Grobid reference observations retain raw citation strings.** Fulltext
  requests now explicitly request them; Stage 1 and TEI-cache receipts force a
  one-time migration of older builds. Structured fields remain best-effort
  interpretations of that evidence, not replacements for it. Source-backed
  reference assertions catch lost title words and author surnames even when
  the bibliography count is unchanged.

- **A read-only reference-reconciliation audit makes the missing-work ranking
  inspectable (#155).** `tools/qc/reference_reconciliation.py` reports current,
  mapped and unmapped observation populations; compatibility-edge collapse;
  mapping methods and producer versions; bounded raw/parsed evidence; and a
  deliberately separate same-title/year review signal with exact author-set
  agreement marked. It never mutates the authority graph or MCP surface.

- **Per-document `ocrmode` makes a wrong OCR routing decision correctable
  (#186).** A BibTeX entry may set `ocrmode = {force}`, `{redo}`, or
  `{skip-text}` alongside `ocrlang`. A valid directive forces OCR to run even
  when detection called the PDF born-digital, preserves the detector's
  original verdict in `scan_detection.json`, round-trips through the
  bibliographic authority database, and fingerprints every descendant stage
  so edits cannot be skipped by resume. Unknown values are recorded and
  ignored with a warning.

- **An on-demand, self-contained page audit joins the evidence needed to
  diagnose extraction and caption failures (#274).**
  `python -m pipeline.page_report <output>/documents/<HASH> --pages 12-13`
  renders the processed page beside selectable Docling text, with toggleable
  PDF-word, figure, ROI, chosen-caption and rejected-caption overlays plus
  page-level counts and provenance. It is re-runnable, excluded from served
  bundles, and bounded to 200 selected pages unless explicitly overridden.
  The old unconditional `visualizations/*.png` pass has been removed, so a
  normal build no longer creates a second raster set for every document.

### Fixed

- **Grobid fallback metadata now recovers when capability returns (#174).**
  Both resume gates distinguish deliberate disablement, incomplete extraction
  and complete evidence. Startup outages and per-paper request/parse failures
  retry on recovery without redoing OCR or figures; valid cached evidence
  survives an outage. Curated BibTeX headers no longer conceal missing reference
  extraction. Service-version and optional declared-producer changes invalidate
  TEI; malformed responses cannot become successful caches. Status reports
  recorded Grobid outcomes. Explicit disablement archives active TEI and removes
  Grobid-derived data until reenabled. The Stage 1 CLI now honors configured
  URL/disablement and applies the configured request timeout. Legacy capability
  receipts require one metadata refresh, with Grobid available to regenerate
  unverified TEI. Custom model contents are not remotely verified.

- **Stage 1 configuration changes invalidate their consumers and descendants
  (#174, configuration tranche).** OCR/probe controls, raster settings,
  fallback chunking, Grobid consolidation, panel mode/explicit model and quality
  thresholds now participate in both resume gates. Stale TEI is archived and
  revalidated against the prepared PDF and consolidation inputs. Panel changes
  rebuild an unsplit figure base; full extraction replaces obsolete images
  and sidecars. Interrupted/failed producers cannot retain old success receipts,
  including a failed standalone vision overlay. Chunk/figure links no longer
  accumulate obsolete IDs. Status reports configuration differences read-only.
  Legacy builds without configuration receipts require one Stage 1 refresh;
  this supersedes the metadata-only migration cost described below. External
  service/model provenance and whole-build update acceptance remain open.

- **Input BibTeX edits and PDF renames invalidate metadata (#174, first
  tranche).** Both resume gates compare the canonical resolved entry for each
  paper, including entry absence, and the filename used for provenance and
  fallback title/year. Unaffected papers retain their receipts; affected papers
  reuse OCR, Docling, chunks and materialized figure/vision results instead of
  overwriting vision ROIs with the CPU floor. Added/removed identical source copies refresh
  the path inventory even when extraction skips, allowing embedding to update
  its stored paths. Legacy builds need one metadata refresh, not full extraction.

- **Embedding updates replace each document atomically (#271).** Resume now
  verifies a content/metadata fingerprint and the committed row generation,
  count, model and dimension. Changed text cannot be skipped, shortened or
  empty documents cannot leave stale rows, and retries cannot append duplicates.
  Stage 1 no longer writes fake embedding receipts. Legacy receipts trigger a
  one-time re-embedding; status, dry-run and bundling verify actual completion,
  and bundling checks the entire index instead of sampling a marker. Switching
  models requires a whole-index rebuild even when dimensions match.

- **SLURM builds now load the selected checkout through the entire job
  chain.** An alternate `REPO_DIR` previously changed the working directory
  while the installed `corpus` entry point and phase subprocesses could still
  import an older editable installation. The shared setup now exports the
  selected checkout on `PYTHONPATH`. The launcher and every build phase verify
  package paths; phases also reject a commit change after submission. These
  checks prevent a successful bundle stamped by new code from masking old
  extraction code.

- **`bib.authority --rebuild` no longer retains artifact stamps that empty the
  rebuilt database.** The command previously dropped `works` while preserving
  `paper_artifacts_processed`, causing unchanged `metadata.json` files to be
  skipped when the fresh schema was seeded. Rebuild now drops all derived
  tracking and observation tables along with the authority graph while still
  retaining the rate-limited BHL lookup cache.

- **Authority tests now exercise the production database schema (#237).**
  Five fixtures had copied mutually inconsistent `works`, `work_authors`,
  `citations`, and `work_aliases` declarations, so adding a production column
  broke only whichever tests happened to select it. They now call
  `bib.authority.create_schema` and insert small datasets into the real schema.
  A repository test rejects any future hand-written `works` declaration under
  `tests/`, making schema additions visible everywhere immediately.

- **Reference ingestion no longer manufactures title evidence (#226, #239).**
  A Grobid `monogr/title[@level='j']` is now retained only as the journal;
  only an analytic title or a non-journal monograph title can become the cited
  work's title. The authority builder also clears exact legacy
  title-equals-journal duplicates before matching. When an exact DOI lookup
  misses, narrowly shaped OCR variants (hyphen insertion/loss or a long
  alphabetic suffix glued onto the DOI) may resolve to an existing work only
  when its title independently passes the established token-set and straight
  similarity thresholds. The original DOI is not silently rewritten, and
  neither DOI resemblance nor title similarity can merge a work alone.

- **High-confidence citation variants now reach the in-corpus work across a
  damaged first-author block (#155, #225).** DOI normalization strips
  `info:doi/` and decodes percent escapes; a dangling-parenthesis DOI is only a
  candidate when independent title evidence passes. Cross-block matching
  requires the same year, the complete order-insensitive author-surname set,
  substantive titles, one unique in-corpus candidate, and the established
  exact or dual fuzzy-title thresholds. A conflicting DOI is retained in the
  raw observation and named in the mapping method. On the 2026-09-01 reference
  bundle, the Mapstone false node fell from 54 citation edges to 3 while the
  canonical node held 54 after replay; canonical *Siphonophore biology*
  reached 109, and
  encoded/prefixed DOI ghosts ceased contributing. No resolver-safe
  title/year/author-set candidate remained in the missing list, while 96
  title/year-only review leads remained; #155 stays open for those damaged or
  genuinely ambiguous cases.

- **Unmappable observations no longer force every unchanged authority run to
  rematerialize (#240).** The current reference bundle exposed a legacy schema
  limitation: 17 citing documents that share canonical work identities with
  another corpus document contribute 2,509 observations but cannot be
  represented by scalar `works.corpus_hash`. The producer-validity check now
  compares against the mappable active population, taking the unchanged pass
  from a repeated 338-second rebuild to a `(0, 0)` no-op in 9.22 seconds. The
  QC report exposes the unmapped population; a first-class one-work/many-
  document relation remains required rather than being hidden by this guard.

- **The figure-detection scorer now measures physical detections and the actual
  default MCP type filter (#194).** The documentation called “drop
  `graphical_element`” the served surface, but default retrieval also excludes
  the `unclassified` review bucket. Both paths now consume one shared evidence
  type set. Scorer v2 also excludes caption/vision children that deliberately
  share a plate or compound image and collapses typed panel siblings. On the
  clean 35-document gold corpuscle, 653 entries contain 261 image-sharing
  logical records and collapse to 384 physical detections: 0.883 recall /
  0.865 precision raw, 0.867 / 0.985 after dropping `graphical_element`, and
  0.827 / 1.000 on the default MCP type surface. The preceding clean artifact
  has the same 384 physical figures and identical scores despite carrying only
  30 logical children. The former 0.936 / 0.876 headline counted logical
  retrieval records as newly detected images and is invalid. Raw and
  `include_all` records remain available for review.

- **Whole-document OCR failures now fail visibly instead of satisfying clean
  success paths (#264, #266–#268).** Documents with no content layer or only a
  vendor wrapper use forced OCR rather than the self-defeating `--skip-text`.
  A curated non-Latin `ocrlang` pin against a Latin-only text layer triggers
  the rendered-page script check even below the configurable gibberish floor.
  Quality gates remove Docling `<!-- image -->` placeholders before measuring
  text, so adding figures cannot make an empty document healthier. Finally,
  OCRmyPDF exiting zero with every output page textless is persisted and
  emitted as the error-level `ocr_no_text_recovered` gate, including when the
  condition is derived from an older scan artifact. Verified on the reported
  Lin/Zhang canary: its original symbolic-font layer exposed zero CJK
  characters; `force` with `chi_sim+eng` recovered 447 and readable Chinese
  prose.

- **A BibTeX export could not be imported unchanged when it matched by
  `work_id`.** The matcher already supported that stable identifier, but its
  result counter did not, so the import raised `KeyError` before applying or
  stamping the entry. The counter now covers the existing match route; this
  also makes `ocrmode` export/edit/import round trips usable.

- **Incomplete vision responses can no longer masquerade as successful empty
  detections (#269).** Claude and local Qwen now share a panel-count-sized
  output budget. A provider stop at the token cap retries that figure once at
  twice the budget; a second stop is a hard failure even when a parseable
  prefix survived. Malformed JSON or either missing required list fails the
  same boundary. Only an explicit complete empty response becomes
  `no_labels_found`, and `figures.json` preserves the failure reason in
  `pass3_error` until a clean retry replaces it.

- **Figure captions now carry auditable ownership rather than only plausible
  text (#195, #203).** Caption extraction preserves provenance spans across
  Docling page merges, joins adjacent labels and prose, rejects a next-page
  candidate when a substantial local figure is the better owner, and records
  status, confidence, kind, page distance, and bounded chosen/rejected
  candidates. Grouped plate captions now split lists/ranges into per-number
  records, reconcile duplicate assignments only when the counts form a full
  bijection, and use exact next-page legend matches to enrich preceding bare
  labels. The figure MCP responses expose the summary fields without moving
  association work into the server. On the 35-document regression replay,
  scorer v7 now recognizes explicit plate inventories, standalone engraved
  numbers (`1`, `F. 1.`), and a plate heading immediately outside `[PLATE]`.
  It also scores typed identities, so `plate:10` cannot satisfy `figure:10` on
  the same page. The prior gold parser omitted or collapsed that evidence and
  reported 480 pairs where the transcription contains 839, so its old recall
  headline is not a valid release baseline. On the corrected yardstick, the
  retained clean candidate before facing-page expansion reports 313/318
  correct (0.373 recall / 0.984 precision). The complete clean source-PDF
  build reaches **538/545 correct (0.641 recall / 0.987 precision)** with
  fixed-population capacity 544/839 (0.648). The scorer remains independent of
  production OCR rules, reports raw and default-MCP surfaces plus capacity,
  and excludes the anatomical key `Pl.M.` (mouth-plate) from Roman-numeral
  plate labels.

  Complete legends printed on the leaf before a full-page plate now bind only
  when an explicit `PLATE N` heading exactly matches one plate on the following
  page. OCR-damaged `Fic.`/`Fics.` openers are accepted only inside that
  context. Plate and child-figure numbers use separate deduplication namespaces,
  so Plate X no longer deletes Figure 10. On the persisted Totton replay this
  creates captioned logical records and moves Totton from 181/184 correct
  reported identities to **406/411 against 472 gold** in the clean build.
  Inspection of the seven corpus-wide surpluses found three already marked
  uncertain, two OCR-damaged printed numbers without independent repair
  evidence, and one omitted gold structural block. The remaining defect was a
  typed identity collision: a following-page Figure 16 legend could overwrite
  a same-page Plate XVI link. That cross-namespace replacement is now rejected.

  Panel correctness is now measured independently as well. Caption parsing
  supports the common `A, ...; B, ...` style, preserves strong printed sets
  with gaps through L, stops at abbreviation glossaries, joins geometrically
  adjacent panel-description cells, and can reject a wrong structural link in
  favor of a materially closer same-page numbered caption. In the clean build,
  98 gold captions enumerate panels: 92 receive declarations, 89 label sets
  are exact, and label recall / precision is 0.946 / 0.997. No panels are
  declared on the 175 number-matched members of 202 non-panelled gold
  identities.

  Shared historical plates now proceed beyond logical record expansion:
  Pass 2.5 records their numeric figure targets separately from lettered
  panels, preserves an independent pre-expansion `missing_figures`
  cross-check, and admits the host image to one Pass 3 ROI invocation. OCR
  accepts bare numbers only through the exact caption-derived allow-list;
  vision uses the same target set. Detected regions are distributed back to
  the individual figure records, while the MCP server only reads and crops
  that build-time evidence. A separately scheduled `--only vision` run now
  also executes Pass 3c and rebuilds chunk/figure cross-references, matching
  the inline full-run artifact contract.

  Bare plates with no deterministic legend may use a tightly gated
  unconditioned vision fallback. It requires at least two high-confidence
  Arabic number/region pairs, persists accepted and rejected candidates, and
  never treats a model description as a caption. The 34-plate local-Qwen probe
  produced 216 regions; 215 agree with corrected gold. Its one false label was
  an invented A-H grid copied onto figures 1-8, so that conflicting-grid shape
  is now rejected regardless of self-reported confidence. Deterministic Pass
  2.5 leaves only four of those plates eligible; their probe result adds nine
  correct labels and no false ones, all nine of which persist in the clean
  build.

- **`corpus taxonomy ingest` doubled the `names` table on every re-run, and
  v1.2.1 made it fire automatically (#262).** `names` shipped with no PRIMARY
  KEY and no UNIQUE constraint, and both writes were plain `INSERT`. The only
  dedup was a `set` built inside `insert_records`, which knows nothing about
  rows already on disk. Re-ingesting therefore appended a complete duplicate
  set: 801 names became 1,602, then 2,403.

  Latent for as long as the code existed, because a re-ingest was an operator
  choice. v1.2.1's #251 fix added an unconditional pre-build to
  `slurm/batch_pipeline.sh` — justified by the claim that the ingest no-ops,
  which was asserted without checking a row count and was false — so it began
  firing on every launch, on the exact workflow that release existed to
  unblock.

  `names` now carries a unique index on `(name_lowercase, taxon_id,
  name_type)` and both writes are `INSERT OR IGNORE`. Because
  `CREATE UNIQUE INDEX` would fail on a database that already holds
  duplicates, `create_schema` deduplicates first and logs how many rows it
  removed — so **a corpuscle built with v1.2.1 is repaired in place on its
  next ingest**, rather than needing a rebuild. `n_names` now counts rows
  actually written instead of attempts, so the log stops reporting "801
  names" against a table holding 1,602.

  Lookups were correct throughout — `name_set()` uses `SELECT DISTINCT` and
  `lookup()` chooses among identical rows — which is why nothing surfaced
  this except the file growing. The three places that claimed idempotence now
  say what changed and when.

## [1.2.1] - 2026-08-30

### Fixed

- **Every SLURM job opened its stderr with two alarming, meaningless lines
  (#252).** All four batch scripts began with `module purge`. Because
  `sbatch --export=ALL` propagates `LOADEDMODULES` / `_LMFILES_`, a batch job
  starts believing miniconda is already loaded; `purge` unloads it and the
  modulefile's hook calls `conda`, which is a shell function that is not
  exported and does not exist in a non-interactive shell. Hence
  `conda: command not found` and a `CondaError` at the head of every task,
  on jobs that then ran correctly.

  Cosmetic, but it cost real diagnosis time: during a failed Stage 1 launch
  the actual cause was the missing `taxonomy.sqlite` of #251, and this noise
  sat above it and drew the first hypothesis. `module purge` also does not do
  what its presence implies — `StdEnv` is sticky and survives it — so the
  line neither achieved a clean environment nor was needed for one.
  `module reset` restores the same sticky default, matches YCRC's documented
  convention, and emits one informational line. Verified equivalent on
  Bouchet: same resolved `python`, same successful `docling` + `torch`
  import, including from a shell with no environment active.

- **Pass 3b dropped panel bboxes that arrived as pixels rather than
  normalized floats, and counted some of them as successes (#253).** The
  prompt demands "each coordinate is a float in 0.0 .. 1.0" and both
  backends multiplied by the image dimensions on that assurance.
  Qwen2.5-VL frequently ignores it: measured across every Pass 3b log on
  one cluster, 130 of 142 observable responses carried absolute pixels, and
  100% of them since 2026-05-30.

  That produced silent loss two different ways, neither logged:

  - `[17, 808, 150, 1365]` → `x0 = 17*w = 8432` while `x1 = min(1.0, 150)*w
    = w`, so `x1 <= x0` and the panel was dropped. Any pixel bbox with a
    non-zero origin was lost this way, 100% of the time.
  - `[0, 0, 487, 606]` → the guard passes and you get a "panel" spanning the
    whole figure. Wrong data rather than missing data, and it lands in the
    `completed` counter where nothing looks at it.

  A coordinate above 1.0 cannot be normalized, so the units are recoverable
  without guessing. Pixel bboxes are now used as pixels, which strictly
  dominates dropping them and removes the whole-figure impostor. Every
  disposition — normalized, pixels, out-of-range, malformed, degenerate — is
  counted and logged per figure, so a build reports what its model actually
  emitted instead of discarding the evidence. That reporting is what the
  units question needed: the old code left no trace on any drop path, so the
  rate could only be inferred from the handful of responses that failed to
  parse, which is a biased sample of exactly the wrong population.

  The converter was duplicated in both backends, so this defect had to be
  found twice; there is now one `_bbox_to_px` and a test that fails if it is
  inlined again.

- **Pass 3b logs did not record the corpus version (#253).** They carry the
  GPU and the config path, so a vision pass could not be attributed to a
  build afterwards — the `siphonophore_20260828` run could not be tied to a
  version at all, because `runs/` held only records written after it and
  `run.log` attested to a different, later run.

- **`corpus check` promised a taxonomy that `run --only` will not build, and
  the promise cost two SLURM chains (#251).** A missing `taxonomy.sqlite` was
  reported as "will be created on first run" — true for a full `corpus run`,
  false for the phase-split `run --only <phase>` path the `slurm/` scripts
  use, where the orchestrator hard-errors before any work starts. Two
  consecutive siphonophore builds passed `check` clean and then lost every
  Stage 1 array task about a minute in; `afterok` took Pass 3b, Embed and
  Finalize down with them.

  Fixed at four sites, because the wording alone only helps an operator who
  reads `check` — and both lost builds had its output on screen:

  - `corpus check`'s dwca/dwc branch now says what the orchestrator requires,
    matching what the WoRMS branch has said since #139. It was simply never
    brought into line.
  - `pipeline/config.template.yaml` carried the same claim, and `corpus init`
    copies it verbatim, so a new corpuscle asserted it before `check` did.
  - `slurm/batch_pipeline.sh` pre-builds the taxonomy before the first
    `sbatch`. `corpus taxonomy ingest` no-ops when the sqlite exists, so this
    costs about a second on every subsequent run.
  - The fatal precondition now prints to stdout as well as the log. SLURM
    routes the logger's stderr to a sibling `.err` file, so from the vantage
    points an operator uses — the `.out` file, `squeue`, a `documents/`
    directory filling up — a dead chain looked like slow first documents.
    26 of 28 tasks were `FAILED` while `squeue` still showed `RUNNING`.

### Changed

- **`corpus taxonomy ingest` reads `taxonomy.source`, `path` and `root_id`
  from `config.yaml` (#251).** It required `--source` explicitly, which made
  it the one verb that could not be driven from the corpuscle's own config —
  so the SLURM pre-build above would have had to parse YAML in bash. Explicit
  flags still override, per the house rule. `corpus taxonomy ingest` with no
  arguments now does the right thing inside a corpuscle.

- **Pass 3b recorded a truncated VLM response as "this figure has no panels"
  (#253).** The local backend capped generation at a fixed 1024 tokens while
  the panel prompt asks for one six-field JSON object per panel at roughly
  60–90 tokens each, so panel-rich figures ran out mid-object. `_extract_json`
  found no balanced `{...}`, the backend returned `[]`, and `[]` maps to
  `no_labels_found` — a clean result, indistinguishable from a figure the
  model genuinely found nothing in. The figures it cost most were the ones
  panel ROIs matter most for: measured over 1,772 documents, ROI coverage
  fell from 47.8% at 2–3 panels to 13.6% at 10+, which is the signature of a
  fixed output budget rather than a vision failure.

  Three changes. The token budget now scales with the panel count Pass 3b
  already knows from the caption, with the configured `max_new_tokens` kept
  as a floor so a raised value is still honoured. A response that stops at
  the cap with nothing parseable now raises, landing in the
  `vision_backend_failed` counter instead of the clean one; a response that
  stops at the cap *after* a complete object logs a warning, which previously
  produced a partial ROI set silently. And `_extract_json` no longer returns
  `None` on the first `json.JSONDecodeError` — it keeps scanning, so one
  malformed leading object stops discarding every well-formed one after it.

  Same defect class as #254: silent loss recorded as success. **The budget's
  adequacy is not covered by tests and cannot be** — that needs a real Pass 3b
  run on a GPU. What is covered is the reasoning around the call: budget
  arithmetic, JSON salvage against the response captured in the issue, and the
  truncated-vs-absent distinction at the status boundary.

- **The test suite wrote a stray file into the working directory on every
  run (#257).** `tests/test_ocrlang_override.py`'s ocrmypdf stub wrote to
  `cmd[-1]` on the assumption it was the output path, but `prepare_pdf`
  appended `--jobs N` *after* the positional output argument — so the stub
  created a file named after the resolved worker count. The name varies by
  machine: `6` here, `3` on the 4-CPU allocation that reported it. The tests
  passed throughout.

  **#254 is what made it constant.** Before it, `_resolve_ocr_jobs()`
  returned `None` whenever RAM was not the binding constraint, so on a
  large-memory host the `if ocr_jobs:` branch usually did not fire and
  `cmd[-1]` really was the output path. Making the resolver always pass the
  number turned an occasional stray into one on every run.

  Two fixes, because either alone would leave the trap set: the stub now
  writes the output path it was given rather than deriving it positionally,
  and `--jobs` moved in front of the positionals where it belongs. This also
  corrects the record — `3bd5cf9` removed one of these files calling it "a
  shell typo", which it was not.

- **Nothing in CI looked at the repo root (#257).** One of those stray files
  was committed in `ec63c92` and survived four green CI runs and two release
  PRs, caught by eye at the v1.2.0 tag boundary — one merge from a permanent
  Zenodo archive, since the release archives the tree under a DOI that
  `CITATION.cff`'s concept DOI resolves to. `tests/test_repo_root_is_clean.py`
  now holds an allowlist of the files that belong at the top level, and a
  second test fails if that allowlist rots. An allowlist rather than a
  pattern because the offending name is machine-dependent.

## [1.2.0] - 2026-08-29

### Theme — v1.2 extraction fidelity, measured against ground truth

Every quality signal this pipeline had measured it against itself. The soft
consistency rates are corpus-internal agreement, the per-document quality
gates are plausibility checks, and fingerprint diffing compares one build to
the next. None of them could say whether the text was *right*.

This cycle built an instrument that can. The gold set is 35 documents, 761
pages, 1594–2026, 13 languages, every page transcribed from a rendered image
with the transcriber forbidden to open any software extraction of the page in
front of them. That independence is the whole value: an extractor scored
against it is not being compared to a cleaned-up version of itself.

Then it acted on what the instrument said, which turned out to be the larger
half of the work. Four defects were found only by building and scoring, not by
any test: docling selecting a CUDA device its pinned torch had no kernels for
and failing every page instead of falling back; vertical Japanese read by a
horizontal model, worth half the words on those pages; figure numbers missed
because the OCR-damaged spellings of "Fig." are *more common* than the correct
one; and plate legends whose separately-numbered engravings existed nowhere in
the output.

The instrument also caught its own repairs. The plate-legend fix shipped on a
measured recall gain and quietly cost precision, putting 33 records into one
monograph that served another figure's image under a wrong number — found by
scoring the rebuild, and then by opening a page image after three rounds of
JSON analysis had pointed the wrong way. That pattern recurs often enough to
be the cycle's real lesson: a measurement is not a result until you have
looked at what it is measuring.

A fifth defect came from a different comparison entirely — not the gold set
but two full cluster builds of the same library, run a month apart from an
identical config. OCR had been silently blanking pages on ~9.5% of documents
in each, a different set every time, because ocrmypdf sized its worker pool
from the host's CPU count rather than the SLURM allocation and the pages
starved past their per-page timeout. Nothing failed; every affected document
recorded `status=success`. It is in `### Fixed` below, and it is the same
lesson from the other side: the pipeline's own account of a run is only as
good as what it was built to notice.

Alongside the fixes, the per-document bib directives grew to cover what a
curator knows and detection cannot infer — `keeppages` for the physical pages
that are the paper, `doclang` and `pagemap` for what the document is, and a
public BCP-47 resolver so a library does not keep its own copy of the mapping.
Two reference documents, `dev_docs/OCR_LANGUAGES.md` and
`dev_docs/FIGURE_PARSING.md`, record what the gold set says about choosing OCR
packs and about figures, including the several places where the obvious
explanation was measured and turned out to be wrong.

Where it ended up, on the 35 documents: prose coverage 0.945, figure detection
0.923 recall at 0.967 precision on the served surface, caption binding 0.574
at 0.878, and reference parsing at 94.8%. Those are the first numbers this
project has had that mean anything outside itself.

### Added

- **`tools/qc/fidelity.py` scores a built corpuscle against an independent
  gold transcription set (#193).** Every quality signal this pipeline had
  measured it against itself — the soft consistency rates are corpus-internal
  agreement, the per-document quality gates are plausibility checks, and
  fingerprint diffing compares one build to the next. None of them could say
  whether the text was *right*. The gold set can: each page was transcribed
  from a rendered page image alone, with the transcriber forbidden to open any
  software extraction of the page, so an extractor scored against it is not
  being compared to a cleaned-up version of itself.

  The harness reports per page and per document, segmented by script, era and
  scanned-vs-born-digital, because a single mean over 13 languages and five
  centuries cannot distinguish a pipeline that handles born-digital PDFs well
  and Fraktur not at all from one that is mediocre everywhere. Documents bind
  on the source PDF's sha256, never its filename — two editions of the same
  1594 travel narrative have shared a filename in this library while setting
  the same passage on different folios.

  Alignment needed no pipeline change: `text.json` carries a flat string and a
  page count only, but `docling_doc.json` retains `prov[].page_no` on every
  item, so per-page text is reconstructible from the persisted artifact.
  Reconstruction follows `body.children`, which is the order
  `export_to_markdown()` uses to build the string that reaches `text.json`, so
  the reading order being scored is the one a consumer sees.

  It runs in ~7 s over 761 pages and is a new manual pre-release tier (T5) in
  CONTRIBUTING.md's tier table. `tests/test_fidelity_harness.py` covers the
  arithmetic in T0 against a committed three-page fixture, with no corpus
  dependency.

  The reconstruction is validated against `text.json` rather than trusted:
  all 35 gold documents recover at least 95.8% of its tokens. That check
  earned its place immediately — the first version of the walk read only each
  item's `text` field, and a docling table has none, keeping its content in
  `data.table_cells[]` instead. Every table on the page was being dropped,
  costing one document 26% of its tokens, and those pages would have been
  reported as the pipeline losing prose it had in fact extracted. Nothing in
  the pipeline was wrong; `text.json` had the tables all along.

  **Which side is on trial is easy to get backwards, and it is not the
  arithmetic.** The gold set ships its own `crosscheck.py`, which used a
  poppler text layer as the yardstick to validate the *transcription*; this
  harness uses the gold as the yardstick to validate the *extractor*. The
  measures are symmetric and identical in both — `coverage` is
  `|gold ∩ other| / |gold|` either way — so nothing is mirrored or
  transposed. What differs is two judgement calls: which measure leads
  (`recall` there, `coverage` here), and whether a page that yielded nothing
  is excluded as no-signal or scored as the failure it is. This harness scores
  it; adopting the other policy would drop 57 of 761 pages and lift the median
  coverage from 0.891 to 0.908, hiding exactly the pages that need work. Those
  two calls are why several documented false positives there are findings
  here — a garbage text layer, a plate whose lettering exists only as image,
  and text hallucinated from image texture.

- **The gold-markup parser distinguishes markers from brackets printed on the
  page (#193).** The transcription convention uses `[` for markers, but the
  pages print brackets too — citation numbers, units like `[°C]`, a
  translator's `[sic]` — and notes talk *about* brackets. A scanner counting
  every `[` as a nesting level gets all three wrong, and all three occur in
  the gold set:

  - `[NOTE: ... "[" is my marker.]` never closed, so the whole note leaked
    into the compared page as though printed there. 13 of one document's 17
    pages; it is why that document posted 0.767 coverage against 0.998 recall,
    the signature of gold holding text no extractor could ever find.
  - `[sic]` and `[21]` quoted inside a note closed it early, leaving a
    `[FIGURE]` two paragraphs later parsed at the wrong nesting level.
  - An unterminated `[continued opposite` swallowed the `[/FIGURE]` after it,
    so the block never ended and the rest of the page scored as plate
    lettering.

  A `[` now opens a marker only when a known keyword follows it; a quote
  immediately after it marks a mention of the character rather than an
  expression. Structural tag counts now come out at exactly the 348 `[FIGURE]`
  / 76 `[PLATE]` / 65 `[TABLE]` the gold set documents, against 341 / 76 / 60
  before, and no page is left with unbalanced figure nesting. The affected
  document moves 0.767 → 0.927 median coverage and every other document moves
  by less than 0.002; corpus-wide 0.891 → 0.898, born-digital 0.882 → 0.919.
  Brackets that open no marker are now counted and reported per page as a
  gold-integrity signal rather than silently absorbed — four remain, each a
  transcription detail worth an upstream look.

- **Prose coverage and figure-text coverage are reported separately, and
  prose is the headline (#193).** Scoring both together answered neither
  question. Text inside a `[FIGURE]`/`[PLATE]` block is engraved plate
  lettering and panel labels — 12.3% of the gold set's words but around 70%
  of its worst pages — so a combined number is dragged down by the material
  the pipeline is least expected to recover, and body text, which is the
  pipeline's actual job, disappears into the average.

  Split, corpus-wide median coverage is **0.946 for prose** against 0.731 for
  figure text (0.898 combined). The split changes what the run says to work
  on, not just the digits: the 1800–1899 band is 0.812 prose against 0.114
  figure text, so its apparent weakness was almost entirely plate lettering;
  `Chenetal2015` moves from 0.667 to 0.812. Two findings survive and are the
  real ones — CJK at 0.351 prose, and pre-1800 at 0.645, where long-s
  typography is the known cause. One gets *worse* and had been hidden:
  `Linnaeus1735` reads 0.628 prose against 0.628 combined with its figure text
  at 0.634, so its body text is genuinely as bad as its plates.

  `[TABLE]` counts as prose, not figure: a table is body content the pipeline
  is expected to get right. Measured the other way the corpus-wide figure
  moves by 0.002, so it is a naming choice rather than a lever.

- **Vertically-set CJK: `ocrlang = {jpn_vert}` documented, and a build-time
  hint pointing at it (#196).** Tesseract ships `jpn_vert`, `chi_sim_vert`,
  `chi_tra_vert` and `kor_vert`; detection selects none of them, so a
  vertically-set document is read by a horizontal model. On a 1911 Japanese
  monograph scored against a hand transcription, pinning `jpn_vert` takes the
  vertical pages from 0.246 to **0.574** median prose coverage.

  The obvious fix — a Fraktur-style companion promotion, adding `jpn_vert`
  beside `jpn` the way `deu_latf` is added beside `deu` — was measured and is
  **wrong**: the union scores 0.186, worse than plain `jpn`, because the two
  models compete for the same glyphs. An unconditional swap is worse still,
  taking horizontal Japanese from 0.746 to 0.207. Since writing direction is a
  per-page property and `ocrmypdf` takes one `-l` per document, no automatic
  rule is available in this cycle; the directive is the answer, and the run
  log now says so beside the `langs=` line, naming the pack and warning
  against the union. `_VERTICAL_COMPANION` carries the measurement in a
  comment so a later tidy-up does not merge it with the Fraktur promotion.
- **`tools/qc/figure_detection.py` measures whether the figure *objects* are
  right (#194).** Separate from text fidelity and from caption association
  (#195): is every figure found, and is publisher furniture being called a
  figure? Against the 35-document gold set, **recall 0.833, precision 0.841**
  — 71 gold figures with no entry, 67 surplus entries.

  It counts per page, because the gold records no bounding boxes and because
  the corpus-wide totals are a trap: 424 gold blocks against 420 entries looks
  like agreement and is a coincidence, hiding one paper with 6 gold blocks
  against 31 entries and another with 67 against 34.

  **`figure_type == "graphical_element"` is an actionable furniture
  predicate**, which the scorer establishes by scoring the same corpuscle
  under each candidate filter rather than assuming it:

  | filter | recall | precision | F1 |
  |---|---|---|---|
  | all entries | 0.833 | 0.841 | 0.837 |
  | **drop `graphical_element`** | **0.811** | **0.964** | **0.881** |
  | drop uncaptioned `graphical_element` | 0.818 | 0.940 | 0.875 |
  | drop `graphical_element` + `unclassified` | 0.708 | 0.987 | 0.824 |
  | captioned only | 0.724 | 0.959 | 0.825 |

  Dropping `graphical_element` cuts surplus from 67 to 13 for nine real
  figures lost. Also dropping `unclassified` is over-reach — it is a mixed
  population holding 49 real figures, and recall falls ten points.

  Segmenting shows why one corpus-wide number would mislead: born-digital
  precision is 0.562 against 0.890 for scans, because modern papers carry
  publisher furniture historical scans do not; and 1900–1949 recall is 0.477,
  driven by two documents where figures are simply not found
  (`Vanhoeffen1906` 34 of 67, `Kawamura1911a` 6 of 16).

- **Figure-number extraction measured: 97% precise, 32% available (#205).**
  `parse_figure_number` fires for 135 of 420 figures, with 15 of 35 documents
  getting none at all — but where it fires, 131 of 135 numbers are genuinely
  printed on that page, allowing #16's Roman↔Arabic normalisation. Coverage,
  not correctness, is the gap. The largest single cause is the opener word
  being OCR-damaged — `Fic.` and `Frc.` for `Fig.`, consistently and
  document-wide — on documents whose captions are otherwise extracted
  perfectly. No behaviour change here; this is the measurement, and it blocks
  the caption-binding scorer (#195), which would otherwise have measured the
  number extractor rather than the caption heuristic.

- **Figure numbers are recovered from OCR-damaged captions: coverage
  32.1% → 67.6% (#205).** `parse_figure_number` fired for 135 of 420 captions
  in the reference corpuscle, with 20 of 35 documents getting no figure number
  at all. It now fires for 284, and precision went *up* — 97.0% → **98.2%** of
  extracted numbers are genuinely printed on that page.

  The cause was not a language gap. `FIG.` set in small caps is misread
  document-wide, and across 320 stored captions the damaged spellings are more
  common than the correct one:

  ```
  Fic 65   Fig 53   PLATE 35   Figur 17   FIGURE 9   FIG 8
  Fi   8   FiG   6  Plate  3   Figg   1   Frc    1   Fie 1   Puc 1
  ```

  A caption-only tolerant opener now accepts `F` plus up to five alphanumerics
  before the number. Digits are allowed *inside* the opener because OCR puts
  them there — `Fi1G.` and `Fi16.` for `FIG.` — and without that the match
  stopped at `Fi` and captured the noise as the figure number, which is worse
  than finding nothing. Also picked up: `Figur` (German, which the old prefix
  could not reach), `Figg. 2-5` (Italian plural with a range, yielding the
  first number), and `Puc` (Cyrillic `Рис` read as Latin lookalikes).

  **The tolerance is confined to captions.** `_FIGURE_REF_RE`, which scans
  running body text where a loose prefix would bind ordinary prose to figure
  numbers, is unchanged; the caption regex stays anchored to the caption
  start, so a figure mentioned mid-sentence is still correctly refused.

  One latent bug fixed on the way: `[IVXLCDM]` overlaps with ordinary words,
  so `"Fig5 caption"` read the `c` of "caption" as Roman 100. The old
  docstring predicted this and the tolerant opener made it reachable; the
  Roman branch now requires a non-letter after it.

  `figure_number` feeds `get_figures_for_*` and `get_figure_dossier_*`, so
  more figures become reachable by number — but `figures.json` is persisted,
  so an existing corpuscle needs a rebuild to pick this up.

- **Figures with the same number on different pages are no longer merged
  (#205).** `dedupe_figures` grouped on figure number alone, and both of its
  tests compare bounding boxes — which carry no page. Two figures at similar
  coordinates on different leaves therefore read as redundant crops of one
  another, and one was dropped.

  That is routine here: a document that is its own translation prints
  "Fig. 4" once in the original and again in the translation.
  `Carre1969_Nanomia_tr` lost nine of its twenty-two figures that way. Keying
  on `(page, number)` costs nothing, because every legitimate grouping the
  function performs — coequal panels, whole figure plus subpanels — is within
  a single page.

  Found because the tolerant opener above *amplified* it: numbering more
  figures fed more of them into the faulty grouping, taking that document
  from 18 extracted figures to 13. With both changes it extracts **22 — every
  figure docling detects, and exactly what the gold transcription records** —
  with 16 of them numbered against 7 before.

  `dedupe_figures` had no direct tests; `tests/test_figure_dedupe.py` adds
  nine, including the cross-page case and the panel-grouping behaviour the
  function exists for.
- **`figure_detection.py` counts figures, not panel images (#211).** It
  compared raw `figures.json` entry counts against gold `[FIGURE]` block
  counts, and those differ when a figure has panels — the pipeline emits one
  image per panel, the gold records one block with the panels enumerated in
  its caption. Correct panel splitting therefore scored as a false positive.

  On `Totton1965a` this flipped the sign of the error: 198 entries against 195
  gold blocks reads as over-counting by 3, while 191 distinct figures against
  195 blocks is under-counting by 4. Corpus-wide on the reference corpuscle it
  understated precision — **0.841 → 0.854** for all entries, 0.964 → 0.969
  with the `graphical_element` filter — so #194's reported figures are
  restated accordingly.

  Panel siblings are now collapsed on `(page, figure_number)` before scoring,
  and panel splitting is reported beside detection rather than folded into it,
  which is where #195 will score it against the panels the gold caption
  enumerates. Note `panel_letter` is set by `dedupe_figures` but not persisted
  to `figures.json` — only the filename carries it.

  Found by opening the two PNGs. Every hand-analysis of the JSON pointed the
  wrong way first, including one that concluded the panel branch was inert
  when the filenames showed it had lettered them all along.
- **Vertically-set CJK is now detected and OCR'd with the right model
  (#196).** Tesseract ships `jpn_vert`, `chi_sim_vert`, `chi_tra_vert` and
  `kor_vert`; detection could never reach them, so a vertically-set document
  was read by a horizontal model. On the 1911 Japanese monograph in the gold
  set, its vertical pages go from **0.250 to 0.572** prose coverage and the
  document from 0.352 to **0.640** — with no operator intervention. Its
  English half is byte-identical, which is the control: `--redo-ocr` preserves
  born-digital text so the pack never touches those pages.

  Selection is **geometric, not confidence-based**. The obvious approach — OCR
  a sample with each pack and keep the more confident result — was tried and
  fails: Tesseract reads a vertical column as a stack of single characters and
  is *confident* about each, preferring plain `jpn` at 61.4 against
  `jpn_vert`'s 58.1 on pages where `jpn_vert` is more than twice as accurate.
  A line *box* is tall or wide whatever the glyphs inside were read as, and
  the separation is about three orders of magnitude: median line width/height
  0.04 on vertically-set pages against 21–53 on horizontal ones.

  The vote counts **only raster pages**, because `ocrmypdf` takes one `-l` per
  document while writing direction is per page. Born-digital pages are
  preserved by `--redo-ocr` and cannot be affected by the choice, so letting
  them vote would decide the question on pages the decision cannot reach —
  counting every page gives an ambiguous 40% on that document, counting only
  the pages OCR rewrites gives an unambiguous majority. Both directions are
  equally costly to get wrong, so it is a plain majority with no tuned
  threshold, and horizontal Japanese is verified unchanged to 0.000.

  The vertical pack goes in **alone**. Anything else competes for the same
  glyphs, `eng` included: `jpn_vert+eng` scores 0.176, worse than doing
  nothing, and that was the first version of this fix. The `ocrlang` override
  from the previous release remains as the operator escape hatch.

- **Small figures are no longer misclassified as publisher furniture
  (#204).** `figure_type` is not cosmetic — `_REAL_FIGURE_TYPES` in the served
  layer excludes `graphical_element` from *every* tool that returns figures,
  so misclassifying a real figure makes it unreachable rather than merely
  mislabelled. `classify_figure` condemned any item under 50 pts in either
  dimension on size alone; Vanhoeffen 1906 Fig. 11 is an engraved nectophore
  **49 pts wide**, captioned and numbered, and was invisible to every figure
  tool because of it.

  A caption carrying a parseable figure number now overrides the size floor —
  publishers do not number their own mastheads. That evidence is guarded by
  **position recurrence**: running furniture sits at the same place page after
  page, and in the 35-document reference corpus no real figure repeats a
  position even twice while one paper's logo repeats on 24 of 25 pages.
  Without the guard the relaxation promotes that logo ten times, because
  caption proximity had attached a real figure's caption to it.

  Net effect on the reference corpus: **3 real figures recovered, 0
  regressions**, each verified by opening the image.

- **The tolerant caption opener no longer reads prose as a figure number
  (#204).** A regression from the previous release, found by inspecting a
  promoted image that turned out to be a handwritten marginal scribble: its
  caption began `"from  the  coasts  of  British  Columbia"`, the opener
  matched `fro`, and the capture read the `m` of `from` as Roman numeral M —
  figure number **1000**.

  Correctly spelled prefixes may still be followed by a Roman numeral with no
  separator (`PLATE XXI`, `Figur 23`); an OCR-damaged opener must now be
  followed by a period. Every damaged spelling observed in the corpus carries
  one, so coverage moves only 67.6% → 66.9% while precision holds at 98.2%,
  and the three lost captions are exactly the false matches.
- **`--config` with a relative path no longer silently discards every tuned
  setting (#210).** Reported by @ejedwards against 1.1.1.
  `corpus --config <relative-path> run` kept `input_pdfs`, `bib`, `lexicon`
  and `taxonomy` — the CLI resolves those itself and passes them as absolute
  arguments — while `ocr`, `chunking`, `stage_timeouts`, `huge_document` and
  `quality_gates` were replaced by built-in defaults. The run looked entirely
  normal; one INFO line was the only trace.

  Two independent causes, both fixed. `_resolve_config_path` returned the flag
  value verbatim, and it is forwarded to each stage subprocess, which the
  orchestrator runs from `REPO_ROOT` rather than the operator's directory — so
  a relative path missed. `CORPUS_CONFIG` had the same hole. Both are now
  resolved to absolute, with `~` expanded.

  Separately, a **named** config that cannot be read is now an error rather
  than a fallback. `load_config` treats a missing file as "use defaults",
  which is correct for the implicit `./config.yaml` and wrong for a file the
  operator asked for by name; a typo would otherwise run to completion on
  defaults.

  This had been masking #209: an `ocr.jobs` setting appeared to have no effect
  because the file carrying it was never read.
- **OCR no longer OOM-kills a workstation (#209).** Reported by @ejedwards
  with a process table: `ocrmypdf` was never passed `--jobs`, so it ran one
  Tesseract worker per CPU at ~1.9 GB each. A 12-core host reached for ~20 GB
  of Tesseract, plus 3.4 GB for the Grobid JVM, and was killed on a 532-page
  scan.

  It could not be worked around from outside. `ocrmypdf` takes its worker
  count from `multiprocessing.cpu_count()`, which **ignores CPU affinity** —
  so neither `taskset` nor a cgroup CPU limit reaches it, and #182's
  `systemd-run` memory limit only contains the blast radius while the build
  still dies.

  New `ocr.jobs` config key. Left unset, the cap is derived from host RAM
  using the components measured in that issue — ocrmypdf's own parent process,
  per-worker resident set, and a reserve for the Grobid JVM, docling and the
  page cache — and it only ever *lowers* the worker count: where RAM is not
  the binding constraint it returns nothing and ocrmypdf's default stands, so
  no large host gets slower. An explicit value is honoured verbatim, including
  above the derived cap. The chosen value appears on the `Running OCR` line.

- **A timed-out OCR no longer orphans its Tesseract workers (#209).**
  `subprocess.run(timeout=)` kills the direct child; ocrmypdf's workers are
  grandchildren, so they were reparented to PID 1 and kept ~20 GB allocated
  after the pipeline had given up on the document. OCR now runs in its own
  process group and the whole tree is killed together, with `KeyboardInterrupt`
  forwarded explicitly since a new session no longer receives terminal signals
  by itself. A `SIGKILL` delivered to the pipeline from outside remains
  unhandleable, but is now the only path that orphans the tree.

- **A plate holding several numbered engravings now yields one figure per
  engraving (#203).** Historical plates carry several separately-numbered
  figures under a single legend. docling extracts the plate as *one* picture
  and emits the legend as separate text items; caption association then bound
  whichever labelled line was vertically nearest — on the page this was built
  for, a bare `Fig. 36.` at the foot — and Figures 31 to 35 existed nowhere in
  the output. Asking for Figure 33 returned nothing.

  A page whose text carries several distinct figure numbers, more than were
  extracted from it, now gets one record per number, each with its own caption
  line. Corpus-wide figure recall **0.849 → 0.917**; `Vanhoeffen1906` alone
  goes 0.463 → 0.806 with 23 figures recovered, and no document regresses.

  The records **share the plate's image**, marked `shares_image_with`, rather
  than copying it — cropping an individual engraving needs OCR of the
  lettering printed on the plate and is a separate problem. Sharing rather
  than duplicating matters: written naively it produced six byte-identical
  copies of a 987 KB plate and took one document's figures directory from 13
  MB to 30 MB.

  Only pages where the legend names *more* figures than were extracted are
  expanded, so a modern paper docling has already separated is untouched.

- **`tools/qc/caption_binding.py` measures whether captions are bound to the
  right figure (#195).** The third of the three questions the gold set was
  brought in for, after text fidelity (#193) and figure detection (#194).
  Against the reference corpuscle: **binding recall 0.527, precision 0.870** —
  when the pipeline reports a figure number it is right about the page 87% of
  the time, but it finds a number for only about half the figures that print
  one.

  It binds on the **figure number**, not the caption text. A caption-similarity
  matcher was built first and reports 44%, which is mostly artifact: one paper
  prints every caption twice, in Chinese and English, and scores 0 of 10 while
  every figure is in fact bound correctly; plate pages carry `FIG. 1` and
  nothing else; and a document that is its own translation prints `Fig. 1`
  twice, legitimately.

  The denominator matters as much as the measure. Gold pages are full of
  figure numbers that are *references* — "see Fig. 18", "figured by Bigelow
  (op. cit., fig. 34)". Counting those gives 939 numbers and a recall of
  0.296, which measures nothing. Restricted to numbers printed *inside* a
  `[FIGURE]`/`[PLATE]` block the denominator is 482.

  Each gold block is classified before scoring — prose caption 123, bare label
  **236**, lettering only 59, nothing printed 6 — so no rate is computed over
  blocks it does not apply to. That bare-label majority is itself the finding:
  most figure blocks in this material print a label and nothing more, which is
  why text similarity was never going to work.

- **`dedupe_figures`' whole-figure/subpanel branch is reachable (#207).** Its
  two stages shared one measure. `_bbox_overlap_fraction` divides by the
  *smaller* box and is therefore symmetric, so a fully contained panel scored
  1.0 — tripping stage 1's 0.5 redundancy threshold and being discarded before
  stage 2's 0.8 containment test could classify it. Anything stage 2 would
  have accepted, stage 1 had already thrown away. `FIGURE_TYPE_SUBPANEL` was
  never assigned by that path, which the data confirms: none of the 420
  figures in the reference corpuscle carried it.

  The stages ask different questions and now use different measures. "Are
  these the same box?" is intersection over *union*, which punishes a size
  difference so a nested panel is not mistaken for a duplicate crop. "Does
  this box contain that one?" keeps the original formula.

  **This changes nothing on the reference corpus, and that is stated rather
  than glossed.** Across its 97 same-page picture pairs, zero merge decisions
  differ and zero pairs are nested-but-not-duplicate — docling does not emit a
  whole-figure crop alongside nested panel crops on this material. The fix
  removes a documented mode that could not execute, and will work if a corpus
  does produce that shape; it is not an improvement to any current number.
- **Grobid consolidation is reachable from config, with a measured default
  (PLAN v1.2 §3).** `consolidateCitations` had never run in this project's
  history: `process_fulltext` defaulted it to 0, `metadata.py` called that
  method overriding nothing, and the `grobid:` block exposed only `url` and
  `disable`. The setting existed and could not be changed.

  It stays off, but now on evidence rather than assumption. Same PDFs, flag
  off then on, against the gold corpuscle:

  | document | era | DOI rate | Grobid time |
  |---|---|---|---|
  | `Ahuja_etal2026` | 2026 | 86.1% → 88.9% | 3.1s → 6.4s |
  | `Stepanjants2014` | 2014 | 0% → 6.9% | 2.0s → 2.7s |
  | `Bernstein1934` | 1934 | 0% → 3.6% | 2.3s → 3.4s |
  | `Benasso_Stroiazzo1976` | 1976 | 0% → 0% | 1.7s → 2.6s |

  Six DOIs recovered across 194 references, for 1.4× to 2× the Grobid time.
  The split is by era rather than by anything the flag controls — modern
  reference lists already carry DOIs and CrossRef holds the rest, while the
  historical works this corpus is mostly made of are not in it. Which is why
  both flags are now exposed rather than hard-coded: the arithmetic is a
  property of the library, not of the pipeline. A corpus of modern papers may
  well find it worth the round trips.

  Baseline for context: 719 references across the corpuscle, 32.4% carrying a
  DOI — 69–86% in papers from 2020 on, and **0%** in every document before
  1980.
- **Vendor wrapper detection covers JSTOR and Google Books (#216).**
  `_VENDOR_BOILERPLATE` held three markers and missed the two most common
  wrappers in a scanned library. Measured over the 1,772-document
  siphonophore library, pages 1–2: JSTOR 20 documents, ResearchGate 6, BHL 6,
  blank notice 5, ProQuest 1, Google Books 1 — **34 with any wrapper**.

  **Wrappers only; publisher imprints are deliberately excluded.** A wrapper
  is a cover sheet or rights notice that is not the paper; an imprint is
  branding on the paper's own pages — a ScienceDirect header, a Springer
  footer. In the same library **373 documents carry an imprint** against 34 a
  wrapper. Here a false match costs a wasteful re-OCR; it would be
  destructive for #188, where a marker is evidence to *drop* a page, and
  dropping a ScienceDirect header page deletes the article's first page.

  Effect is small and that is worth stating: of 27 documents matching a new
  string, **2 actually re-route** — the rest are born-digital papers carrying
  a JSTOR cover, with 9.8K–21.8K characters of real text, correctly held back
  by the existing `total_chars < 5000` gate. Rasterising those would be wrong.

  The list is high-precision, and its recall is better than that first
  reading suggested. An independent page-by-page annotation of the same
  library — a reader working from rendered pages rather than strings — found
  34 documents with a vendor wrapper, which is what the list finds. The
  inference that BHL wrapper pages were sitting there as un-OCR'd images was
  wrong, but so is the simplest replacement for it. 220 of these PDFs carry
  BHL or Internet Archive provenance in their *embedded metadata* and only 8
  still have a cover page — because the covers were stripped by hand when the
  library was assembled, not because BHL never shipped them. 210 of the 212
  without one have a ModDate later than their CreationDate, against 4 of the
  8 that kept theirs. So this recall is partly borrowed from someone else's
  editing: a corpus built from fresh BHL downloads would carry 220 cover
  sheets and would need the list to fire on all of them. It would — the 8
  survivors match — but 34-of-1,773 is a fact about a curated library, not a
  bound on what wrappers cost.

  What the annotation does show is that wrappers are the small part of the
  problem — a title page in 391 documents against a wrapper in 34. Front
  covers, flyleaves, bookplates, a bound volume's own title page: no vendor
  string, because there is no vendor. Those pages are the book, and no
  addition to this list can reach them, which is why #188 needs a structural
  signal rather than a longer table.

- **A public BCP-47 → Tesseract pack resolver (#215).** `_ISO_TO_TESSERACT`
  was the bridge between a detected language and an OCR pack, and it was
  private, so anything outside `scan.py` needing that mapping had to
  reimplement it — a library's annotation pass, resolving a per-document
  language into an `ocrlang` directive, already carried a copy marked
  temporary. `bcp47_to_tesseract` is that copy's exit path.

  BCP-47 subsumes what the table held — bare ISO 639-1 codes are already valid
  tags, so nothing regresses — and reaches two packs no key could name before.
  `grc` existed only inside the Greek fallback union; `deu_latf` was reachable
  only via a visual OSD verdict or a `deu` special case, and langdetect can
  only ever say `de`. There was no way to state "German, set in Fraktur".

  **Deliberately not wired into the build.** The OCR language decision stays
  two tiers — an explicit `ocrlang` pin, else detection. Deriving packs from a
  bib field at run time would make the derivation table an input to
  `processed.pdf` without being part of any fingerprint: improve the table and
  nothing invalidates, so documents keep their old OCR while the log reports
  the new `-l`. Resolving at annotation time leaves `ocrlang` a literal,
  directly-fingerprinted value.

  Vertical CJK stays out. BCP-47 describes language and script; vertical
  setting is typesetting, and `jpn_vert` must not be unioned with `jpn` —
  0.574 and 0.246 alone, 0.186 together. A test asserts no `_vert` pack is
  ever returned, because making the CJK entries symmetric with the Fraktur one
  would look like tidying and would regress #196.

- **`AGENTS.md` says where library-facing code goes, and a test enforces it.**
  Three tiers, sorted by whether the code encodes knowledge about PDFs and OCR
  or knowledge about one collection: generic goes in `pipeline/` as a public
  function, read-only inspection in `tools/`, collection-specific judgment in
  `skills/`. `pipeline/` may never import from `tools/` or `skills/`, which
  was true by convention and is now `tests/test_import_direction.py`.

- **Plate legend expansion read cross-references as legend entries (#231).**
  #203 gives a plate one record per figure its legend names, and shipped on a
  measured recall gain (0.849 → 0.917). Scoring the rebuild showed it also
  cost precision — 0.970 → 0.892 on the served surface — and the whole of the
  loss was one 226-page monograph of running prose.

  The scan took a figure number from anywhere in a text item. A legend line
  qualifies; so does a cross-reference, and a monograph is full of them — a
  species heading reading `Plate XX, figures 1, 2`, a parenthetical
  `Text-figure 106 (see p. 170)`, a citation of someone else's plate. Two on
  a page cleared the threshold, and the page's one real text-figure was then
  cloned under the referenced number. **33 of that document's 232 records
  shared an image file with another record under a different number**: asking
  for figure 20 returned the picture of figure 53, which is worse than a miss.

  A legend line opens with the label of the figure it describes. Anchoring to
  that, with the number required to follow the label and only punctuation
  between (so `figured by` is not `fig. 53`), gives back the precision without
  giving back the recall — 0.894 / 0.962 on the served surface, against
  0.833 / 0.970 before #203 and 0.901 / 0.892 as it shipped. The offending
  document goes from 34 spurious records to none, and the plate volume #203
  was written for keeps all 24 of its real ones.

- **Two inert bib fields for curation: `doclang` and `pagemap` (#214).**
  Curating a scanned library means recording two things the pipeline had no
  way to hold. `ocrlang` (#176) is an *instruction* to Tesseract, meaningful
  only when OCR runs; what was missing is the *fact* — this paper is Russian,
  this one is 19th-c. German set in Fraktur, this one is Ancient Greek. That
  is what a person or an annotation pass determines by looking, it is worth
  recording for born-digital papers too, and it is in a different vocabulary
  from Tesseract's. `doclang` holds a BCP-47 tag, which is the only candidate
  that can express `de-Latf`, `zh-Hant` and `grc` — and `de-Latf` decides
  whether a scanned paper OCRs to text or to whitespace.

  `pagemap` is free text describing the scan's physical structure. It exists
  so a page-range directive is reviewable: a bare `keeppages = {3--20}` tells
  the next reader nothing about whether pages 1–2 were a scanner wrapper, a
  blank verso, or a mistake.

  **Both are read by nothing** — no stage, no fingerprint — and that is the
  feature rather than a limitation. Correcting a language label or fixing a
  typo in a note must not reprocess a document, where `ocrlang` rewrites
  `processed.pdf` and rightly invalidates every OCR-dependent stage. Tests
  assert the absence directly, including that the fingerprint builder takes
  no such argument and that no `entry_doclang()` accessor exists, because the
  instinct when adding a bib field is to copy the `ocrlang` template whole and
  that template fingerprints.

  Deliberately not the standard BibTeX `language` field, for the reason
  `ocrlang` isn't either: reference managers populate it by default, and an
  ordinary Zotero export must not silently start steering anything.
- **`tesseract_packs` no longer records "unknown" as "none" (#197).** Three
  documents in the reference corpuscle once recorded `"tesseract_packs": []`
  while OCR ran with seven. `_compose_ocr_langs` returns early when
  `_available_tesseract_langs()` is empty — before the configured fallback
  union is reached — so a failure to *enumerate* the installed languages was
  written down as a resolution that found nothing.

  `scan_detection.json` is the operator-facing record: `corpus status` reads
  it, and #176's `ocrlang` workflow tells an operator to consult it when
  choosing a pack to pin. A record saying "no packs" when seven were used
  sends that diagnosis backwards, on exactly the documents an operator would
  be investigating.

  **The symptom no longer reproduces** — checked across all 35 documents by
  comparing each `scan_detection.json` against the `Running OCR` line in its
  own log, and record and invocation now agree everywhere. So this ships the
  guard rather than a fix: if the probe comes back empty again, the record
  says `tesseract_langs_unavailable` and a warning points at the log line
  that holds the truth, instead of looking like a resolution that found
  nothing.

- **`keeppages`: which physical pages of a file are the paper (#188).** A PDF
  in a scanned library is frequently not just the paper — a library cover
  sheet, a Russian original bound in front of its English translation, runs of
  blank versos. The costs compound rather than add: `detect_scan_type` samples
  pages to choose the OCR mode *and* the language pack, so front matter
  decides how the body behind it is read; then OCR pays full price for the
  filler, and a calibration target becomes a figure.

  Physical 1-based positions, never printed folios — on the documents this
  targets it is precisely the front matter that has no printed number. An
  entry routinely carries both `pages = {41--118}` (the journal pagination)
  and a `keeppages` that disagrees, and that is correct. BibTeX page syntax
  throughout: `2,4,8--20,22--40,55`, and `40--` for "to the end". Normalised
  to a sorted, deduplicated set, so a selection cannot reorder a document. A
  range past the last page is clamped with a warning, recorded in
  `scan_detection.json`; an unparseable range is an error, because silently
  keeping every page looks exactly like success.

  Applied **before** scan detection, by rebinding the temp PDF the later
  stages already read — so `scan.py` needed no change at all. It runs before
  the huge-document gate too, which makes a selection the supported way to
  bring an oversized bound volume into scope, exactly as that gate's own
  error text asks.

  **Page-number provenance is the part that needed care.** Once pages are
  dropped, `page` in `figures.json` and `text.json` is a position in the
  subset, and that is the number served to a client — a figure reported as
  page 3 that is page 44 of the scan is a citation error nothing downstream
  can detect. Both are carried: `page` stays subset-relative because it is
  what indexes the artifacts on disk, and `source_page` says where it came
  from. The resolved selection recorded as `keeppages_selected` *is* the map,
  so there is no second structure to keep in sync.

  Fingerprinted across every OCR-dependent stage, and required rather than
  defaulted at all four call sites — a page-range edit invalidates strictly
  more than an `ocrlang` edit, since it changes not how the PDF is read but
  which pages it consists of.
- **An `ocrlang` pin that contradicts measured page geometry is now
  recorded and warned about.** The vertical-CJK hint went silent as soon as
  a pin was honored, reasoning that the operator had made the call
  themselves. That holds for a hand-written tag and fails completely for a
  derived one.

  An annotation pass deriving `ocrlang` from `doclang` (#214, #215) *cannot*
  get this right: BCP-47 describes language and script, and vertical setting
  is typesetting — deliberately out of `bcp47_to_tesseract`'s scope. So the
  derivation emits `jpn` for a vertically-set Japanese paper every time, the
  pin overrides #196's geometric verdict, and the one mechanism that would
  have said so is disabled *by the thing that caused it*. Found on
  `Kawamura1911a`, the document #196 was written for, where the pin costs
  more than half the words on those pages — `jpn_vert` 0.574 against `jpn`
  0.246.

  The pin still wins, because `ocrlang` is documented to beat every inferred
  signal and silently overriding an explicit instruction is worse than
  obeying a bad one. But the conflict is now `ocrlang_overrides_vertical_cjk`
  in `scan_detection.json` rather than only a log line that scrolled past
  during a 40-minute build, and the warning says why a derived tag gets it
  wrong so the generator gets fixed rather than one bib entry.


- **`dev_docs/OCR_LANGUAGES.md` — what the gold set says about choosing
  Tesseract packs.** The cycle produced a lot of measurement about pack
  selection, spread across the README, four issue threads and a working
  session. Collected into one stable page that code comments can cite.

  It also records a correction. The obvious explanation for why a mismatched
  pack damages this literature is dictionaries — Tesseract exposes
  `language_model_penalty_non_dict_word` and loads six dawgs per language, and
  89–92% of these tokens fall outside any English word list. Measured, that is
  wrong: re-running with every dawg disabled changed **zero of 2,236
  recognised tokens**, because those parameters govern the legacy word
  permuter and Tesseract 5 runs the LSTM engine.

  What differs is the character repertoire the model was trained on. Of the
  tokens `por` recovers and `eng` loses on the same pages, **83.5% carry a
  diacritic** in the original against 11.4% of those both recover — and
  native-pack advantage tracks diacritic density, from Dutch at 0.4% where
  `eng` ties to Portuguese at 3.4% where it loses 0.08.

  The two failure modes that bracket the problem are recorded with numbers:
  models contesting the same glyphs (`jpn_vert+jpn` 0.186 against 0.574 for
  `jpn_vert` alone; seven packs on monolingual Latin below `lat` alone), and
  substituting rather than adding a pack (13× the damage on out-of-wordlist
  vocabulary, which is the vocabulary retrieval depends on).
- **A pin that narrows what detection resolved is now recorded (#245).** An
  `ocrlang` pin does not merely choose packs, it *discards* the ones detection
  had resolved — and that was invisible. Pinning 31 of 35 gold documents from
  a derived `doclang` tag moved corpus-wide prose coverage 0.9474 → 0.9450
  with every language correct, because each pin replaced a union with a single
  pack: `eng+deu+deu_latf+fra+lat+spa+por` → `lat` cost 0.079,
  `swe+cat+fra+eng` → `swe` cost 0.063.

  Narrowing is not uniformly bad — two documents improved — so the pin still
  wins, and the record is `ocrlang_narrowed_from` in `scan_detection.json`
  plus one INFO line. What was wrong is that a directive costing 0.05 looked
  identical to one gaining it.

  The comparison is against the list OCR would actually have used. An earlier
  version compared against *targeted* resolution instead, reasoning that
  pinning over a fallback union is not narrowing but the case `ocrlang` exists
  for. Run against the reference corpuscle that flagged 4 of the 22 documents
  whose pack list the pin changed, and missed the largest regression in the
  set — `Linnaeus1735`, seven packs down to `lat`, −0.079, with no targeted
  resolution recorded at all. The distinction was real but not the one worth
  acting on; it survives as `ocrlang_narrowed_from_targeted`.
- **`fidelity.py` reports taxon-token coverage beside prose coverage (#244).**
  Coverage weighted every token equally, and for this literature that
  undervalues exactly what retrieval keys on: replacing `por` with `eng` on a
  Portuguese paper costs 0.010 on English-wordlist tokens and 0.129 on
  everything else, and the binomials live in the second bucket. Every
  extraction decision this cycle was made by reading this instrument, so a
  systematic bias in it was in all of them.

  Scored against the corpuscle's own `taxonomy.sqlite`, word by word, with the
  denominator always reported and no rate below 10 taxon tokens on a page — a
  corpuscle's taxonomy covers one clade while the gold spans all of nature, and
  the 801-taxon siphonophore snapshot labels 58 tokens in the whole of
  *Systema Naturae*.

  **Summarised by mean and p10, not the median used everywhere else.** 53% of
  qualifying pages recover every taxon token, so a median sits at exactly 1.000
  and cannot move while the tail does: mean 0.885, p10 0.643 on the same 229
  pages.

  It found what it was built to find on the first run. **17 pages score more
  than 0.2 better on prose than on taxa** — `Hosiaetal2024` p6 is prose 0.976
  and taxon 0.140, a page the headline calls near-perfect while 86% of the
  names it is retrieved by are lost.


- **`dev_docs/FIGURE_PARSING.md` — what the gold set says about figures.**
  Companion to `OCR_LANGUAGES.md`, covering the three questions this cycle
  learned to keep separate: are all the figures found and is furniture being
  called one, is each caption bound to the right figure, and is the text
  inside a figure recovered.

  Records where the numbers are strong — 1950–1999 scans at recall 0.975 /
  precision 0.964, `Totton1965a` at 0.974 / 0.995 over 195 figures — and where
  they are not: born-digital precision 0.607 raw against 0.919 for scans,
  because publisher logos are figures as far as a layout model is concerned;
  caption binding before 1900 at recall 0.091, because the numbers are
  engraved on the plates rather than typeset.

  Also why each measure has the shape it does. Detection is counted per page
  rather than matched figure by figure, because the gold carries no bounding
  boxes and per-page counting is the strongest claim the data supports.
  Captions bind on the figure *number*, never on caption text, because a naive
  text match reports 44% and is mostly artifact. And gold blocks are
  classified before scoring — 229 of 376 are a bare label, so caption text is
  computable for 86 pairs out of 465 and is reported rather than headlined.


### Changed

- **The vertical-CJK section of the README predated #196.** It said detection
  never selects the vertical models and that `ocrlang` is "currently the only
  way to get it right"; neither has been true since #196 landed. Rewritten to
  describe what happens — a sample of the scanned pages is rasterised, line
  orientation measured, and the `_vert` pack swapped in on a majority vote —
  and to recast `ocrlang` as the override for disagreeing with that choice.

- **The pyflakes gate now covers `tools/` (#75, #193).** It linted
  `pipeline/`, `mcpsrv/` and `bib/` only, so the operator scripts under
  `tools/` never got the NameError check it was built for — and those are run
  by hand at release time, where a NameError costs a whole manual run rather
  than a fast test failure.

### Fixed

- **OCR sized its worker pool from the host, not the allocation, and
  silently blanked ~10% of a cluster build's documents (#254).**
  `--tesseract-timeout` does not fail a page. Its documented behaviour is to
  give up on OCR, copy the un-OCR'd image into the output and exit 0 — so the
  page survives visually, carries an empty text layer, and the document is
  recorded `status=success` with an empty `stage_failures[]`. The only trace
  was a `WARNING` in the SLURM log of whichever array task happened to
  process it.

  The cap added in #209 only ever fired when *memory* bound. On a Bouchet
  stage1 node — `CPUTot=64`, 991 GiB — an 8-CPU step could afford 390
  workers, so the cap returned `None`, and `None` hands the decision back to
  ocrmypdf's `multiprocessing.cpu_count()`: 64 Tesseract workers inside an
  8-CPU cgroup, every page on a sliver of a core, until the per-page timeout
  fired and blanked it. Nothing OOMs; the pages just starve.

  `_resolve_ocr_jobs` now reads the allocation — `ocr.jobs`, then
  `$SLURM_CPUS_PER_TASK`, then the CPU affinity mask, then a cgroup `cpu.max`
  quota, then the host — with the memory budget derived the same way
  (`$SLURM_MEM_PER_NODE`, cgroup `memory.max`, physical RAM, smallest wins),
  and **always passes `--jobs` when it knows the number**. Returning `None`
  is now reserved for a host nothing could be determined about. On a
  workstation the value it passes *is* ocrmypdf's own default, so nothing
  gets slower.

  Found by comparing two full builds of the same 1,769-document library.
  Because the condition is load-dependent it selects a different set each
  run: 31 documents lost >80% of their text between builds while 28
  independently gained it back. `Johnson_Widder2001.pdf` went from 83,293
  characters to 671, with `OCR completed successfully` logged both times.

- **Blanked pages are recorded and gated instead of counted and discarded
  (#254).** ocrmypdf names the pages it abandons, one line each, page number
  first; `_log_ocr_warnings` reduced that to `stderr.count(needle)` and
  kept only the total, while a separate end-state check found pages with no
  text and conceded in its own docstring that it could not tell why. Joining
  them removes the ambiguity exactly, with no heuristic: a page that is empty
  **and** was named by Tesseract is lost text; a page that is empty and was
  not named is a plate, a blank verso or a figure-only page. The first is a
  new `error`-severity quality gate, `ocr_pages_blanked`; the second stays
  expected, which matters in a corpus this full of plates. The page list is
  persisted to `scan_detection.json` as `pages_blanked`, so it reaches
  `summary.json` and is queryable after the fact rather than living in one
  array task's log. Selectable like any other gate:
  `corpus status --filter-gate ocr_pages_blanked`,
  `corpus run --re-process-flagged ocr_pages_blanked`.

  The gate fails the *document*, not the run. Aborting would kill a 28-task
  array over a transient condition affecting ~10% of documents while the rest
  of that task's documents are fine.

- **The main OCR invocation now goes through `_run_ocr`, so a timeout takes
  the Tesseract workers with it.** #209 added the process-group kill but
  wired it only into the `--redo-ocr` retry path; the common path still used
  `subprocess.run(timeout=)`, which kills the direct child and leaves its
  grandchildren reparented to PID 1 — burning the cores the *next* document
  needs, which is the same starvation #254 is about.

## [1.1.1] - 2026-08-08

### Fixed

- **`CITATION.cff` is now schema-valid, which unblocks Zenodo archival.**
  The file had failed CFF 1.2.0 validation since it was added, on two
  independent counts: a top-level `year` key, which the schema does not
  define (`date-released` is the key), and `license: "(see repository)"`,
  where an SPDX identifier is required. `license` is now `MIT`, matching
  `LICENSE`.

  This is why Zenodo stopped archiving. The correlation is exact —
  `CITATION.cff` landed one day before v0.3.0, v0.1.0 and v0.2.0 archived
  without it, and every release from v0.3.0 through v1.1.0 did not.
  Zenodo reads the file to build deposition metadata, so an invalid one
  fails the deposit *after* its receiver has already answered `202`,
  which matches the symptom exactly: webhook accepted, no record created.
  The webhook was never at fault.

  Consequence while it was broken: the concept DOI in `CITATION.cff`
  always resolves to the newest archived version, so every citation
  pointed at v0.2.0 from May while seven newer releases existed.

## [1.1.0] - 2026-08-08

### Theme — v1.1 post-1.0 correctness

v1.1 was triaged out of the 1.0 release rather than absorbed into it:
the pre-release UX review, the viburnum production build, and two items
carried across several cycles. Nothing in it blocked an install, which
is the bar 1.0 set, and the v0.6 38-tool freeze permitted all of it —
every item is a bug fix or an additive change under
[API_STABILITY.md](dev_docs/API_STABILITY.md).

**What it turned out to be: silent wrongness in OCR, and test signals
that had stopped meaning anything.** A narrower cycle than its candidate
list, with the shape coming from the material rather than the plan. The
OCR half began as one issue about six viburnum papers (#176) and proved
to be three defects, only one of which that issue had identified — the
third being that there was no per-document override at all, so an
operator watching a Korean flora paper OCR as English had no recourse.
The testing half was #185: corpus-wide soft checks firing on 65% of a
production corpus, which is indistinguishable from off.

Validation was a from-scratch **699-paper viburnum build** on macOS
arm64 — 699/699 documents, 0 stage failures, 67,221 chunks, 5,340
figures, 7h57m. It found two resume bugs that would have shipped a
feature that silently did nothing, confirmed #184's figure shrinking at
production scale (3.4 GB → 2.3 GB, −32.4%, against the −33.2% measured
on the siphonophore sample), and opened #186, #187 and #188.

*(Backfilled 2026-08-25, from this cycle's own PLAN.md write-up rather
than reconstructed from the diff. Older entries without a theme section
were deliberately left alone — see CONTRIBUTING.md's release ritual.)*

### Added

- **A paper's bib entry can now pin which Tesseract packs OCR it, via a
  new `ocrlang` field (#176).** Language detection had no per-document
  override: `ocr.ocr_languages_default` is consulted only when targeted
  resolution returns nothing, which a *confident but wrong* detection
  never reaches, so an operator watching a Korean flora paper OCR as
  English had no recourse short of hand-patching `scan_detection.json`.
  The bib is where it belongs — entries are already keyed to individual
  PDFs, already loaded before the run, and already carry per-paper
  operator directives that aren't bibliographic facts (`license`,
  `serve`). Write it the way ocrmypdf spells it, `ocrlang = {ell+eng}`.

  Deliberately *not* the standard BibTeX `language` field, which means
  "language of the work" and which reference managers populate by
  default — reusing it would let an ordinary Zotero export silently start
  steering OCR.

  The pin beats both langdetect and Tesseract OSD, because being outvoted
  by the signal you are correcting would defeat the purpose; detection
  still runs and its verdict stays in `scan_detection.json` alongside
  `ocrlang_requested` / `ocrlang_honored` / `ocrlang_dropped`, so what
  the pipeline believed and what the operator overruled are visible side
  by side. Uninstalled pack names are dropped with a warning rather than
  passed to ocrmypdf, which would fail the paper outright; if none
  survive, the tag is ignored and detection decides. It selects packs
  only — it does not force OCR, so tagging a born-digital paper is a
  no-op.

  Adding, changing or removing a tag re-runs the paper on `--resume` —
  every stage, not just the OCR ones, because re-OCR rewrites
  `processed.pdf` and docling, Grobid, chunking and taxon extraction all
  descend from those bytes. Without that the feature would be a trap in
  two different ways: the paper would be skipped outright, or it would be
  re-OCR'd correctly and then have its new text discarded by stages that
  still counted as complete. Untagged papers keep skipping as before, so
  no existing corpuscle is re-OCR'd. Validated on a 699-paper corpus:
  tagging two papers re-processed exactly those two and skipped 697, and
  a Korean flora paper that had been OCR'd as English cleared its
  `gibberish_after_ocr` gate with 7,372 Hangul characters recovered.

- **Figure PNGs are now written in their smallest lossless encoding
  (#184).** Figures are ~97% of a served bundle — 16.3 GB of the 19 GB
  siphonophore corpuscle — and on a 2,527-figure sample of that corpus
  47.6% of them are greyscale stored as RGB, because a scanned line
  engraving rasterizes to three bit-identical channels. Dropping the two
  redundant ones and re-encoding measured **-33.2%** over the whole
  figure set, with nothing to trade away.

  Two rules keep it lossless. The original is always a candidate: a
  blanket `optimize=True` re-encode made **52.1%** of sampled figures
  *larger*, since line art compresses worse under PIL's filter choices
  than under the encoder that wrote it, so smallest-of-N is what
  separates -33.2% from -18.1%. And the channel drop is verified rather
  than trusted — the candidate is decoded and compared against the source
  pixels before it may replace anything, with any mismatch keeping the
  original. Detection is exact channel equality, so a figure with even one
  genuinely coloured pixel is left alone.

  Bitonal (≤2 grey levels) → mode `1` is deliberately not attempted: it is
  only bit-exact when the two levels are 0 and 255, covers 3.4% of
  figures, and is a rounding error in the savings — not worth owning a
  lossy path.

  Runs once per paper after every producer has written its final bytes,
  rather than at each save site, because in the default `native` mode the
  #121 resolution pass re-renders and overwrites docling's output. Applies
  to new ingests; existing corpuscles keep their current figures until
  rebuilt. No effect on `processed.pdf`, `docling_doc.json`, chunks, text
  or embeddings.

- **`corpus check` now warns when the Grobid container isn't the image
  `docker-compose.yml` specifies.** `/api/isalive` proves something is
  listening on the port, not what — and because the compose service
  carries `restart: unless-stopped`, a container created from an older
  compose file keeps serving that port indefinitely. A macOS arm64 host
  was found running the full DeLFT image months after the compose default
  moved to `lfoppiano/grobid:0.8.1`; the pre-flight reported `[ok] Grobid:
  reachable` throughout, while the image README §Grobid forbids on Apple
  Silicon was parsing every reference in the corpus. The check is a warn,
  not a failure — DeLFT is a supported opt-in on AVX-capable Linux
  x86_64 — and stays silent for a remote or Apptainer Grobid, a host
  without Docker, and a non-clone install with no compose file.

- **A GPU too old for the pinned torch no longer breaks every build
  (#198).** `torch.cuda.is_available()` answers "is there a visible NVIDIA
  GPU", not "can this torch build run kernels on it", and the two differ on
  exactly the hardware a lab workstation has. On the same machine with an
  unchanged dependency set, 20 days apart:

  ```
  2026-08-06  Accelerator device: 'cpu'      → clean build
  2026-08-26  Accelerator device: 'cuda:0'   → CUDA error: no kernel image
                                               is available for execution
  ```

  A GTX 1080 (compute capability 6.1) became visible to torch, whose pinned
  build ships `sm_75` and up. **Nothing in the project changed** — a driver
  appeared and a working install stopped working, which is precisely the
  reproducibility the #98 pins exist to provide.

  `pipeline/accelerator.py` now checks capability rather than availability,
  allowing both an exact binary kernel and forward PTX JIT so working
  hardware is never pushed onto the CPU, and logs which card it rejected and
  why. It is applied at all four call sites that made the same assumption
  independently: docling in `extract.py` and `prefetch.py` (which previously
  set no `accelerator_options` at all, leaving docling on `auto`), the
  `vision-local` host gate in `cli.py`, and the embedding encoder in
  `embeddings.py`. The last of those was only caught by running a real build
  — the docling fix alone left `corpus run` failing at the embed stage.

  New `compute.accelerator` config key (`auto` | `cpu` | `cuda` | `mps`); a
  pinned value is honoured verbatim, since second-guessing it would make the
  knob useless for the case it exists for. `CORPUS_DEVICE` still wins for
  embeddings. The resolved device is recorded in `text.json`, so two
  corpuscles that differ because one ran on CPU and one on CUDA are no longer
  indistinguishable after the fact.

### Changed

- **The corpus-wide soft consistency checks are now asserted as corpus
  rates rather than per paper (#185).** Seven checks compare two derived
  artifacts and read a disagreement as a pipeline defect — a `.bib` title
  that should appear in the body text, a figure number cited in text that
  should have a figure. That premise holds for some material and not for
  the rest, and which one you get is a property of the *corpus*, not the
  code: a curated title legitimately differs from what an offprint with
  no title page prints, and a 19th-century monograph legitimately cites
  plates bound in another volume. Asserted per paper on a production
  corpuscle they produced 1,690 failures across 1,157 papers — 65% of the
  corpus — which cannot be triaged, so in practice the signal was off.

  Each check is now one aggregate assertion against a ceiling, bucketed
  by file type, and a breach names the rate, the denominator, the ceiling
  and the first ten offending documents, so triage still starts from a
  hash. Ceilings are calibrated across two deliberately unalike corpora —
  1,787 documents of marine invertebrate zoology heavy on plate-based
  monographs, and 699 of botany that is mostly modern journal articles —
  and `CORPUS_SOFT_RATE_CEILINGS` points at a JSON file overriding any
  subset of them for a corpus of a different shape.

  `references_match_corpus_papers` is deliberately set where it can only
  catch catastrophe: its premise, that a paper cites other papers *in
  this corpus*, measures corpus cohesion rather than reference parsing,
  and cohesion is not a pipeline property. A corpus assembled by mining
  references outward runs 52.7% / 84.4% as its normal case.

  The hard per-paper checks are unchanged — `has_text`, `has_title`,
  `has_authors`, `has_chunks`, `text_min_length`, `chars_per_page`. They
  agree with the quality gates `run.log` already reports and are few
  enough to act on individually. The whole rate group stands down below
  25 documents, so the demo is unaffected.

- **`corpus status` no longer blames missing language packs for every
  `gibberish_after_ocr`.** The hint said the cause is "usually a missing
  Tesseract language pack ... install the pack and rerun"; on a 699-paper
  build with all 126 packs installed that was wrong for every paper it
  fired on. It now walks three causes cheapest-first — pack absent, pack
  present but a different one chosen (which `ocrlang` overrides), or the
  OCR result discarded because `redo_ocr` / `skip_text` preserved a
  corrupt digital text layer that no pack choice can reach — and names
  the artifact that answers each. It also notes that a table-dense paper
  can score high without being broken, since the score is a text
  heuristic and numeric tables read as noise to it.


### Fixed

- **A non-Latin Tesseract OSD verdict no longer discards a
  high-confidence language detection (#176).** Pack resolution branched
  `if visual_script … elif detected_iso`, so one OSD call could throw
  away langdetect's answer entirely. The rationale — OSD reads the page
  image, so it stays right where a corrupt text layer misleads langdetect
  — holds, but OSD is a guess from a single sampled page and is wrong
  often enough to matter. On a 1,787-document siphonophore build from
  2026-04-18 (a pre-v1.0 tree, and the only corpus at hand),
  **188 papers had a p>=0.99 Latin-script detection overruled by a
  non-Latin OSD verdict** — Bigelow 1914 read as Cyrillic, Alvarino 1976b
  as Greek, Broch 1928 as Japanese.

  How many of those were actually damaged depends on which packs the host
  had. Resolution returns `[]` when the OSD-named pack isn't installed,
  and the fallback union rescues the paper. So the **129 Cyrillic verdicts
  are the certainly-harmed set** — `rus` is in the stock
  `ocr_languages_default`, so it is always there to win — while the 59
  exotic verdicts (jpn 24, ara 13, tha 6, han 5, ell 2, ben 2, hin 2) only
  bite where that pack happens to be installed. Of the 129, the English
  ones are rescued by the `eng` suffix anyway, leaving **40 papers whose
  own language pack was displaced** (fra 19, spa 14, dan 2, deu 2, hrv 1,
  ita 1, nld 1). None tripped `gibberish_after_ocr`.

  Re-OCRing 10 affected papers both ways (3 body pages each, same source
  PDF, only `-l` differing, on a host with all 126 packs installed — the
  worst case) measures what a wrong pack costs when it does win:
  **wrong-script glyphs -68%** (752 -> 242) and **correct diacritics
  +227%** (256 -> 836), with word count flat at +1%. Mean gibberish score
  moved -2.6%, which is why no gate fired. The failure is not "text with
  the accents stripped": a pack in the wrong script rewrites whole words,
  and on four of the ten papers *not one* of the language's diacritics
  survived. Chun 1888 lost a clause to hallucinated Thai numerals —
  `Sie besassen eine völlig runde Schwimmglocke mit relativ sehr kleinem`
  came out as `เ 1111 mit relativ sehr kleinem`. Car & Hadži 1914
  (Croatian, OSD said Cyrillic) rendered `dalje redovna opažanja` as
  `dalje гейоупа opazanja`.

  What this does *not* establish is the state of any shipped corpuscle.
  The tree measured here predates v1.0 by four months, its file-type mix
  is nothing like the current one (1,379 born-digital / 254
  broken-text-layer / 154 scanned, against 734 / 14 / 1,021 reported for
  the v1.0-era build in #185), and classification has changed underneath
  it: Chun 1888 is `broken_text_layer` with an OSD verdict of Thai in that
  tree and `scanned` with no OSD verdict at all today. Its April text
  carries 1,176 umlauts, so `tha` was not installed and the fallback
  rescued it. Whether papers in the v1.0 corpuscle are affected is
  unmeasured — that corpus was not available here.

  The two signals are now unioned, OSD first. Tesseract takes a multi-pack
  `-l` without complaint, so a wrong OSD verdict costs one surplus pack
  instead of the right one. An OSD verdict of `Latin` is unchanged: it has
  no entry in the script→pack map, so it never formed a disagreement in
  the first place.

- **`tesseract_packs` in `scan_detection.json` now records what ocrmypdf
  is actually given.** It held targeted resolution alone while the caller
  appended `eng` on top, so a paper OCR'd with `rus+eng` was filed as
  `['rus']` — despite a comment claiming the field mirrored the real
  invocation. Operators grep this field to find bad OCR, and the
  half-truth is what made the OSD bug above look worse than it was: the
  original diagnosis blamed a missing `eng` fallback that had been there
  all along. The value now comes from the same function `prepare_pdf`
  calls, so the two cannot drift again. `tesseract_pack_available` still
  reports *targeted* resolution, so it keeps its meaning now that the
  pack list is never empty.

- **`TestCitationGraph::test_references_match_corpus_papers` no longer
  fails on a small corpus.** It reads "zero in-corpus citations" as
  evidence of broken reference parsing, which requires that a match was
  likely to begin with. On the 4-paper demo it isn't — no demo paper
  cites another (Pugh 2001 and Dunn 2005 both cite *Pugh 1975/1989*
  against a corpus holding *Pugh 2001*) — so two papers failed on every
  clean local build, and an operator validating a fresh install got a red
  suite with no indication it was expected. The check now stands down
  below 25 documents and is unchanged above it.

- **`resume_scenario` is now actually deselected from a bare `pytest`.**
  Both `pytest.ini` and `tests/test_resume_scenario.py` documented it as
  deselected by default, but nothing implemented that — no `addopts`, no
  `collection_modifyitems` hook — so a plain `pytest` built the demo
  corpuscle twice mid-run and took ~8 minutes instead of ~2. Now an
  `addopts` line does what the marker description always claimed. `pytest
  -m resume_scenario` still opts in, since `-m` is last-wins.

- **The README install block now installs `pngquant` on macOS arm64.**
  It listed `bash tools/install_ocr_extras.sh` as the step that provides
  `pngquant`, but conda-forge has no osx-arm64 build, so on Apple Silicon
  the script prints a Homebrew command and exits without installing
  anything. An operator following the README verbatim finished the
  install believing they had it, and every scanned paper thereafter came
  out at `--optimize 1` — 90 MB instead of 35 MB on a 45-page Russian
  scan. INSTALL.md had this right all along; the README's one-command
  block did not.

## [1.0.0] - 2026-08-05

Through most of this cycle `dev` carried a plain `0.6.0` — the version
actually released on 2026-06-04 — because the post-release step that
reintroduces a pre-release suffix had been skipped. Everything built from
`dev` between that release and this one therefore stamped itself v0.6.0,
whatever it actually contained. A `1.0.0rc1` suffix restored the
distinction late in the cycle, and this release drops it; see
[CONTRIBUTING.md](CONTRIBUTING.md) for the convention.

### Theme — v1.0 installability

1.0 is the version strangers install, so the cycle's organizing
principle is that **a green CI badge must mean a fresh install works**.
The v0.6 MCP surface freeze holds — no new tools. See
[dev_docs/PLAN.md](dev_docs/PLAN.md) for the wave plan.

Two things ran the cycle's fixes past its tests. A pre-release UX pass
worked through the README as a new user would, against a smoke corpus of
35 papers spanning 13 languages, 4 Fraktur, 5 Russian and 3 image-only —
which is where most of the OCR and first-run work below came from. Then a
full production run on Bouchet: **1,769 of 1,769 documents, 0 stage
failures, 261,093 chunks embedded, 934 of 934 eligible figures through the
vision pass, and no job in the chain ending TIMEOUT.** The first attempt
at that run is what surfaced the silent SLURM-chain failure fixed below.

### Changed (breaking)

- **`mcp` pinned to `2.0.0`, `mcpsrv/` migrated to the 2.x API**
  ([#156](https://github.com/caseywdunn/corpus/issues/156)). `mcp` was
  unpinned and PyPI's latest became `2.0.0`, which removed
  `mcp.server.fastmcp` — so a fresh `conda env create` produced 18
  test-collection errors and a `corpus serve` that could not start,
  while CI stayed green on a cached pre-2.0 env. `FastMCP` is now
  `MCPServer` (`mcp.server`), `Image` moves to `mcp.server.mcpserver`,
  and the private `_mcp_server` backing attribute is `_lowlevel_server`.
  The 38-tool surface, the `@mcp.tool()` schema generation, the
  `ToolManager.call_tool` instrumentation monkeypatch, and the SSE
  transport are all unchanged — verified end-to-end over a real
  `--transport sse` server (bearer auth, initialize handshake, 38 tools,
  every smoke layer). **Operators pinning `mcp` themselves must move to
  `2.0.0`;** 2.0 additionally pulls `httpx2`, `mcp-types`,
  `opentelemetry-api`, `pyjwt`, and `python-multipart`.
- **Every pip dependency is now pinned exactly**
  ([#158](https://github.com/caseywdunn/corpus/issues/158)): `ocrmypdf`,
  `lancedb`, `langdetect`, `uvicorn`, `anthropic`, `pytesseract`,
  `qwen-vl-utils`, and `accelerate` join the already-pinned ML stack
  (#98) and `mcp` (#156). CI re-resolves `environment.yaml` on every run,
  so an unpinned dep made CI itself nondeterministic — any push could
  break for reasons unrelated to the push. `lancedb` (#71) and `mcp`
  (#156) had each already done it. Bumping one is now a deliberate,
  reviewable commit across `environment.yaml`, `requirements.txt`, and
  `pyproject.toml`.
- **OCR now runs on everything that is not digitally native, not just on
  what reads badly.** `detect_scan_type` only ever inspected the
  *content* of a PDF's text layer, so `born_digital` meant "reads
  plausibly", not "was produced digitally" — a scan carrying any
  third-party OCR layer was trusted and never re-OCR'd. On the 32-paper
  smoke corpus that was 29 of 32 papers, including all four Fraktur ones,
  whose layers are visibly corrupt but score *below* the gibberish
  threshold. The entire multilingual OCR apparatus, `deu_latf` included,
  ran on 4 documents out of 32.

  Detection now asks whether pages *are* scans, independently of the text
  layer: `_scanned_page_fraction` measures how many sampled pages carry a
  single image covering ≥50% of the page. The populations separate
  perfectly on the reference corpus — every scan 0.50–1.00, every
  born-digital paper exactly 0.00 — and pages are sampled across the
  whole document rather than the first N, because Kawamura 1911a is a
  born-digital English typescript for pages 0–7 and a Japanese scan from
  page 12, and front-only sampling reported 0% raster and skipped OCR on
  13 pages of Japanese. Mixed volumes get `--redo-ocr`, which replaces
  OCR text while leaving genuine digital text alone, so a bound-in
  typescript isn't rasterized to fix a scanned half; uniformly-scanned
  documents get `--force-ocr`. Never `--skip-text`, which would preserve
  the layer we just rejected.

  Language selection had the same circularity — it was read off the layer
  we distrust, which sent Olfers 1824 (German Fraktur) to the Catalan
  pack. It now comes from OCRing a 5-page sample, per page and unioned,
  because bilingual originals-plus-translations are routine in this
  literature (verified on `jpn+eng`, `eng+fra`, `eng+rus`, `eng+spa`
  pairs). Pages are sampled from the middle 15–85% — the front is covers
  and plates, the tail is references — and Tesseract OSD picks each
  page's script so it is OCR'd with only that script's packs. Probe pages
  whose own OCR output is gibberish are discarded: OSD is not infallible,
  and on Bernstein 1934 it read a Cyrillic table as Latin, producing text
  langdetect called Catalan at p=0.86. Confidence cannot catch that;
  gibberish score can (0.61, against 0.12–0.39 for that document's
  genuine pages). Two results are recorded in the config so they are not
  re-litigated: probe DPI must stay 300 (at 200, Kawamura lost its
  Japanese and Linnaeus 1735 went from correctly finding nothing to a
  confident wrong `ca`), and langdetect ships no Latin profile, so Latin
  can only ever be mis-identified — a confidence floor routes those to a
  union containing `lat`.

  **This re-runs extraction on every existing corpuscle and changes
  stored text.** On the smoke corpus: 35 papers / 3,323 chunks / 420
  figures, 31 documents OCR'd where 4 were before, Olfers 1824 going from
  `!Oir llt/ne €lHbfl1rr.` to readable German and Eschscholtz 1825 from 3
  extracted taxa to 6. `zero_references_unexpected` rose 8 → 12, and that
  is *correct*: Grobid had been fabricating references out of mojibake (De
  Haan 1827's 23 "references" had surnames `Den Nam Wui I .`, `Jdt`,
  `Ta`, `On`), verified by running Grobid against the original and the
  re-OCR'd PDF side by side.

### Added

- **T3 — scheduled clean-room CI lane**
  ([`.github/workflows/clean-room.yml`](.github/workflows/clean-room.yml),
  [#158](https://github.com/caseywdunn/corpus/issues/158)). Weekly plus
  `workflow_dispatch`. Runs the documented install path end to end —
  cold conda solve, the real `docker-compose.yml`, demo `corpus run`,
  `corpus_required`, SSE round-trip — with the HuggingFace model cache
  **disabled**, so a genuine first-run model download is exercised the
  way a new user experiences it (the path that 429'd in #140).
  `dev_docs/ec2_smoke.sh` is relabelled **T3-bare** and remains the
  manual pre-release check for the bare-host apt/miniforge bootstrap.
- **`corpus prefetch`** ([#159](https://github.com/caseywdunn/corpus/issues/159)).
  Downloads the three models the pipeline otherwise fetches on first use —
  docling's page-layout model, docling's TableFormer, and the ~4.3 GB
  BGE-M3 embedding model — with retry and backoff, because HuggingFace
  429s anonymous traffic and a shared institutional NAT looks like abuse
  from the other side. `--include-vision` adds the ~16 GB local VLM.
  Prints the cache directory and warns when `HF_HOME`/`HF_HUB_CACHE` are
  unset. INSTALL.md documents the offline pattern: prefetch where there is
  internet, then run with `HF_HUB_OFFLINE=1` where there isn't.
- **`corpus check` gained three pre-flight probes** (#159,
  [#160](https://github.com/caseywdunn/corpus/issues/160)) — the OCR
  toolchain (`tesseract` / `ocrmypdf` / `ghostscript` on PATH; a failure,
  exit 3), the Tesseract language packs against
  `ocr.ocr_languages_default` (a warning naming the missing codes and the
  `tools/install_tessdata.sh` invocation that fixes them), and the model
  cache (a warning, never a network call — safe on an air-gapped node).
  Previously the only `shutil.which` calls lived in `pipeline/scan.py`,
  i.e. checked mid-run; and skipping `install_tessdata.sh` — a *required*
  post-install step — silently OCR'd Cyrillic and Fraktur against the
  English pack.
- **Install documentation for clusters and offline hosts**
  ([#153](https://github.com/caseywdunn/corpus/issues/153), from a user's
  own install report). INSTALL.md now covers redirecting conda/pip/HF
  caches out of a small home directory, the `HOME` override some sites
  need, running Grobid under Singularity in a batch job, and the fact that
  `grobid.url` must name the allocated compute node rather than
  `localhost`.
- **T1-compose — `docker-compose.yml` is exercised on every push**
  ([#161](https://github.com/caseywdunn/corpus/issues/161)). T1 starts
  Grobid as a GHA `services:` container, so the compose file every
  non-HPC user actually runs was covered by no test — which is how both
  #146 (wrong default image, crash-looping on Apple Silicon) and #157
  (broken healthcheck) shipped. The new job boots the real file, waits
  for `/api/isalive`, and requires `docker inspect` to report `healthy`,
  which is the standing regression guard for #157.

  It earned its keep immediately: on its first run it caught that
  `docker-compose.yml` was missing `-XX:-UseContainerSupport`, without
  which this image's JVM aborts at startup on some cgroup-v2 hosts
  (upstream JDK bug, #72 — a `CgroupV2Subsystem.getInstance` NPE) and
  Grobid never binds a port. A local cgroup-v2 host boots fine without it
  while GHA's cgroup-v2 runners fail every time, so the flag now ships on
  by default; `-Xmx4g` pins the heap regardless, so it costs nothing.
  Users on an affected host previously saw only "Grobid never comes up".
- **`tools/install_ocr_extras.sh`** — installs `pngquant`, mirroring
  `tools/install_tessdata.sh`'s contract (idempotent, exits 0 either
  way). Without `pngquant`, `pipeline/scan.py` silently degrades
  `ocrmypdf` from `--optimize 2` to `--optimize 1`, which was documented
  as "OCR output will be larger, not wrong" — a bad undersell on scanned
  material. Measured on Beklemishev 1969, a 45-page Russian scan: source
  PDF 1.1 MB, `--optimize 1` **90 MB**, `--optimize 2` **35 MB**, 61%
  smaller. Since v1.0 re-OCRs every scan rather than trusting its text
  layer, most of a historical corpus now takes that path — 31 of 35
  papers in the smoke corpus, where 4 did before — so on a full library
  this is the difference between fitting a disk quota and not.

  `pngquant` cannot simply go in `environment.yaml`: conda-forge has
  linux-64 and osx-64 builds but no osx-arm64, and conda env files have no
  platform conditionals, so listing it would break `conda env create` on
  Apple Silicon, a supported target. `jbig2enc` is on conda-forge for no
  platform at all. Both verified with `conda search --platform` rather
  than by trusting the existing comment. The script installs `pngquant`
  from conda-forge where a build exists — **no root required, so it works
  on an HPC account** — prints the brew/apt line where it doesn't, and for
  `jbig2enc` detects an HPC-ish environment and offers `module avail` plus
  a source build instead of a useless `sudo apt-get`. It also notes that
  `conda install` wants a login node, since compute nodes often have no
  outbound network. `corpus check` promotes a missing `pngquant` from info
  to **warn**, quantifies the cost, and names the script; `jbig2enc` stays
  informational, because it compresses bitonal images only while
  `pngquant` covers the colour and greyscale scans that dominate this
  material.
- **Language is visible on the served surface.** Nothing in the tool
  surface exposed language, on a corpus whose entire premise is
  multilingual historical literature — a client asked "what languages are
  in this corpus?" had to infer it from titles. `list_papers` now returns
  `language` / `languages` and accepts a `language=` filter (a bilingual
  volume answers to both its codes), and `corpus_summary` gains
  `by_language` and `n_papers_ocred`. Additive only: the 38-tool freeze
  holds, and no existing field changes shape.
- **An API-stability policy**
  ([dev_docs/API_STABILITY.md](dev_docs/API_STABILITY.md)), deferred out of
  v0.6 because a policy written over a moving surface documents intentions
  rather than commitments. It says what the freeze covers (the 38 tool
  names, parameter names and meanings, response field names, the
  pagination convention, the error vocabulary, profile semantics), what it
  does not, which changes are additive versus breaking, and the
  deprecation path — one full minor plus 90 days, announced in the
  docstring, because that is the tool description MCP clients actually
  read and MCP has no out-of-band warning channel.

  Two exclusions are worth reading before depending on this: on-disk
  artifacts carry their own `ARTIFACT_SCHEMA_VERSION` and are not a public
  API, and **extracted content is explicitly not covered** — a stable API
  does not promise stable data, as this release's re-OCR change
  demonstrates. Pin a built corpuscle, not a version, when you need
  reproducibility.
- **The distribution channel is stated rather than implied**
  (INSTALL.md §"How corpus is distributed"). git clone + conda is
  canonical and the only supported channel; pin a release with
  `git clone --branch v1.0.0`. **There is no PyPI package by decision, not
  by omission:** pip cannot install tesseract, ghostscript, pngquant,
  pandoc or Grobid, so `pip install corpus` would import cleanly and then
  fail on the first scanned PDF — the exact failure class this cycle
  exists to remove. Zenodo carries the citable DOI, which is what a PyPI
  presence would have been standing in for.
- **The OCR budget scales with the document.** `stage_timeouts.ocr` is now
  a *floor* rather than the whole budget, and two keys join it:
  `stage_timeouts.ocr_per_page` (default 30 s) sets the effective
  per-document cap to `max(ocr, ocr_per_page × pages)`, and
  `ocr.tesseract_page_timeout` (default 900 s) caps a single page. See
  **Fixed** below for what each of them was losing.

### Changed (breaking, MCP surface)

- **Figure licensing: the gate decides, not the client**
  ([#154](https://github.com/caseywdunn/corpus/issues/154)). Three defects
  fixed together, because two of them pulled in opposite directions.

  **Advisory metadata no longer leaks into permissive use.** `get_figure`
  and `get_figure_url` injected `publishable` / `license_source`
  regardless of profile, so a model *just authorized* to display a figure
  read `"publishable": false` beside it and withheld it. The default
  profile is `report`, the server refuses nothing under it, and figures
  were being withheld anyway. `license` / `license_url` / `attribution`
  are still present in every profile — captions need them — but the
  clearance *determination* now appears only under a strict profile
  (`manuscript` / `presentation`) or on an explicit
  `get_figure(..., include_licensing=True)`. **`publishable` is gone from
  every response shape**, replaced where it appears by
  `publication_clearance`.

  **"Unknown" is no longer conflated with "restricted."**
  `publishable=0` meant *both* "the rightsholder forbade this" and "we
  could not establish public domain". In the served corpus that collapse
  was total: 55,177 works (86%) were `publishable=0,
  license_source=unknown` and **not one** was asserted
  `all-rights-reserved`. A new `publication_clearance` reports five
  states — `public_domain`, `licensed_open`, `restricted`,
  `undetermined`, `no_record` — and refusal messages name the state that
  caused them, spelling out that `no_record` is an absence of evidence
  rather than a prohibition. Derived from existing columns, so **no
  authority-DB rebuild is required**. Unrecognized license strings now
  warn at build time instead of being silently NULLed, which is what made
  a typo'd `license = {CC-BY 4.0}` (space, not hyphen) block figures as
  firmly as an explicit refusal.

  **`get_figure_roi_image` no longer bypasses the gate.** It took no
  `profile` and never consulted the licensing check, so a client refused
  by `get_figure_image` under `manuscript` could obtain the same pixels
  through it — including the whole uncropped figure, via the no-pixel-ROI
  fallback. It now accepts `profile` and refuses before touching disk.
  `tests/test_freeze_contract.py` gained a check that *every*
  pixel-returning tool accepts `profile`, since the v0.6 plan had warned
  about exactly this class and the warning was honored for
  `figure_http.py` while this tool was missed.

  Also on the HTTP route: `?profile=` is now validated (an unknown value
  400s instead of silently falling through to the server default, which
  `resolve_profile`'s own contract said callers must prevent), and the
  query string is parsed with `parse_qs` so a percent-encoded `label`
  resolves instead of missing its crop. `mcpsrv/default_instructions.md`
  gained a licensing section telling the served model to display what the
  server returns and to pass `profile="manuscript"` when it is
  publication-bound.

- **Lexicon figure retrieval expands synonyms**
  ([#143](https://github.com/caseywdunn/corpus/issues/143)).
  `get_figures_for_lexicon_term` and `get_figure_dossier_for_term` did a
  case-insensitive substring count of the *single string the caller
  passed*, so a query for a lexicon term silently missed every figure
  whose caption used a different surface form of the same concept —
  `wing` didn't find `ala`, `ala` didn't find `wing`, and neither found
  `forewing`. Ingestion had always been synonym-aware (each mention
  records the `matched_text` found and the `canonical` it resolves to);
  retrieval just ignored that layer.

  A query is now resolved to its canonical term and matched against every
  surface form observed in the corpus. Rows report `matched_surfaces` and
  the resolved `canonical`; an unrecognized term degrades to a literal
  search and says `resolved: false` rather than silently looking
  synonym-aware. Measured on the demo corpus — figures returned, before →
  after: `feeding polyps` 0 → 13, `swimming bells` 0 → 16, `siphons`
  0 → 13, `gastrozooid` 11 → 13.

  Because the surface forms come from what extraction actually matched
  rather than from the lexicon YAML (which a distilled bundle doesn't
  ship), this also works **across languages**: the demo's Russian paper
  contributes `нектофор` → `nectophore` and `гастрозоид` →
  `gastrozooid`, so an English query now reaches Cyrillic captions
  (`нектофор` 2 → 16). The corresponding limitation is that a declared
  synonym appearing in no paper's text is not expanded.

### Fixed

- **OCR no longer loses whole pages silently.** OCR was dropping entire
  pages, nondeterministically, and reporting success. Linnaeus 1735, same
  command and same input across two runs — per-page character counts:

      run A  [43, 414, 5027, 6420, 10575, 3796, 6012, 0, 9758, 8462, 6124, 0, 0]
      run B  [43, 414, 5027,    0,     0, 3796, 6012, 0,    0,    0,    0, 0, 0]

  56,631 characters became 15,292. Both runs logged "OCR completed
  successfully", both exited 0, and `corpus status` showed 100% on every
  stage. Reproducing with stderr captured named the mechanism:
  `[tesseract] took too long to OCR - skipping`. That is ocrmypdf's
  `--tesseract-timeout`, whose documented behaviour is to give up on a
  page but *copy the preprocessed page into the final output* — a blank
  page and a clean exit. The pipeline never set it, so it took ocrmypdf's
  default, far too tight for a dense 300-dpi historical scan with seven
  language packs loaded, and load-dependent enough not to reproduce
  reliably. Three fixes: `--tesseract-timeout` is now set explicitly from
  `ocr.tesseract_page_timeout` (default 900 s/page), which takes Linnaeus
  to 92,616 characters across all 13 pages — the *better* of the two runs
  above was still missing 39% of the document, and pages 12–13 were blank
  in both; ocrmypdf's stderr is logged **on success**, not only on
  failure, which is the single line that made this invisible, since every
  message that matters arrives on a clean exit; and a post-OCR check
  counts and names pages that came out with no text, because blank pages
  and plates are legitimate but a *run* of them is not. Exposure scales
  with the re-OCR change above — 31 documents now take this path where 4
  did — and on a cluster with array jobs competing for cores it would have
  been worse and equally silent.
- **A monograph no longer fails OCR for being long.**
  `stage_timeouts.ocr` was a flat 1800 s cap on one document's whole
  `ocrmypdf` call, and on timeout the stage fails outright rather than
  degrading. Measured throughput varies ~15× with how many Tesseract
  language packs are loaded, not just with page count — Totton 1965a, 314
  pages under `eng` alone, ran 1.3 s/page; Linnaeus 1735, 13 pages under a
  7-pack union, ran 20.0 s/page. Extrapolated to the largest document in
  the full siphonophore library (delle Chiaje 1830–31, 1,549 pages) that
  is 34 minutes at best and ~8.6 hours at worst, so the flat 30-minute cap
  failed even the best case and a full-library run would have lost its
  biggest monographs. Raising the flat number is the wrong fix: a budget
  big enough for 1,549 pages lets a genuinely hung 3-page paper burn hours
  before failing, which is the thing the timeout exists to prevent. So
  `ocr` became a floor and the effective cap is
  `max(ocr, ocr_per_page × pages)` — 30 min for a 2-page paper, 157 min
  for Totton, 12.9 h for delle Chiaje. The chosen budget is logged on
  every OCR call so it can be seen rather than inferred. Note that the
  per-page rate is already the worst (7-pack) observed case, so it is
  deliberately *not* scaled again by pack count; an earlier version did,
  and double-counted its way to a 51-hour budget.
- **CJK text is no longer scored for gibberish, and an untrusted language
  is no longer served.** Two defects, both from Latin-prose assumptions.
  `_gibberish_score` excluded CJK *tokens*, which was not enough:
  Yamamori 2014 still flagged `gibberish_after_ocr` despite OCRing
  correctly under `jpn`, because the quality gate scores `text.json`, and
  what remains after excluding CJK is `['##', 'Li', '\_', "『'", '=']` —
  page numbers, figure labels and markup fragments, which score as
  garbage however good the OCR was. On a document that is 55–63% CJK the
  measure has no validity in either direction, so it now returns `0.0`
  above a 0.30 CJK share rather than producing a misleading number. That
  is safe only because `_scanned_page_fraction` now answers "is this a
  scan?" independently, so the heuristic is no longer the only thing
  between a corrupt text layer and the corpus. (This was the fifth defect
  traced to a Latin-prose assumption in that one function.) Separately,
  `_scan_facts` fell back to a paper's `detected_language` even when
  `language_trusted` was `False` — laundering a value the pipeline had
  explicitly rejected into the served API. Linnaeus 1735 is Latin, its
  corrupt text layer reads as Catalan, and the new `language` field
  reported Catalan; it now reports no language, which is the honest
  answer. Verified on a full 35-paper rebuild: corpus text 1,943,374 →
  2,206,893 characters (+13.6%), Linnaeus alone 38,349 → 301,833; bogus
  author-initial panels 3 → 0; `section_class` coverage 11.1% → 16.3%; one
  remaining gibberish flag, on a genuinely garbled plate-only volume.
- **`tools/install_tessdata.sh` installed almost nothing on a fresh
  env.** conda-forge's tesseract 5.5.2 now bundles 158 language packs, and
  they are `tessdata_fast` builds — verified byte-exact against upstream,
  `rus` at 3.9 MB against 15.3 MB for `best`. The script's "file exists →
  skip" check therefore fetched only `deu_latf` and left the low-accuracy
  models in place, under a banner reading `Source: tessdata_best`. It now
  tracks what it fetched in a marker file and replaces anything else, and
  gained `eng`, `--force` and `--help`. The claim that tesseract "ships
  only English LSTM data" was false and is corrected in README.md,
  INSTALL.md and `environment.yaml`.
- **`corpus prefetch` missed docling's HybridChunker tokenizer** while
  asserting "every required model is already cached". Under the documented
  `HF_HUB_OFFLINE=1` recipe, chunking then fell back to the naive
  character chunker — 16 chunks became 1 for the same paper, exit 0 either
  way. `prefetch_docling` now exercises the chunker and counts the
  tokenizer in `DOCLING_REPOS`, so prefetch stops short-circuiting, and
  the fallback logs at `ERROR` with the remedy. `corpus prefetch` manages
  **four** models, not three; README.md, INSTALL.md and the module
  docstring all said three.
- **`get_chunks_for_topic` described `score` as "cosine similarity" while
  returning `_distance`.** That docstring is the tool description every
  MCP client reads, so thresholding or sorting on it inverted the ranking.
- **`resolve_reference` failed on any multi-author query.** Everything
  before the year was captured as one surname, so `Totton Bargmann 1965`
  searched for an author literally named "Totton Bargmann" and matched
  nothing — and typing both surnames is the natural way to look up a
  two-author work, so it failed on first use. It now tries the whole blob
  first (multi-word surnames like `van Soest`, `De Haan`, `Lo Bianco` are
  real in this literature), then each surname in turn; `not_found` reports
  `authors_tried`.
- **`get_citation_graph` truncated the commonest question it is asked.**
  It applied its 50-edge per-node cap even at `depth=1`, where the runaway
  [#87](https://github.com/caseywdunn/corpus/issues/87) guards against
  cannot occur — one node is expanded and `max_total_edges` already bounds
  the walk. So "show me this paper's bibliography" silently returned a
  fraction of a hub work's references. Now unbounded at depth 1, 50
  beyond. The `work_id` tiebreak also degenerated to alphabetical whenever
  citation counts tie, which in a small corpus is nearly always (one list
  came back cut off at "Jacobs 1937"); ties now keep document order. The
  docstring's motivating example was also wrong — it cited "Totton &
  Bargmann 1965's 155 references" from an MCP client's report, repeated
  unverified; the real figure is 213 in `references.json` and 210 citation
  rows. It now also warns that a hub work's depth-1 graph runs ~55 kB,
  past what some MCP clients pass through in one tool result — a
  *transport* limit, not corpus-side truncation, so `truncated: false` can
  be accurate while the client still fails to deliver the payload. Pass an
  explicit `max_edges_per_node` when a bounded response is wanted.
- **BibTeX escapes no longer reach served metadata**
  ([#177](https://github.com/caseywdunn/corpus/issues/177)). A `.bib` field
  is LaTeX source, not text, and corpus was ingesting it as text — so five
  titles in a 702-paper build were served reading
  `… (R.Br.) A.Braun \& Vatke`, and any citation a client composed carried
  the backslash. The `.bib` was correct; `\&` is the required escaping for
  a literal ampersand.

  `bib/parser.py` now decodes a value at the single boundary every field
  already passed through: escaped specials (`\&`, `\%`, `\_`, `\#`, `\$`),
  group-protection braces (`{DNA}`, `{V}iburnum` — these were already
  dropped), and accent commands, which matter here because this corpus's
  authors are largely European — `M{\"u}ller` → `Müller`, `Ma\'nko` →
  `Mańko`, `Fran\c{c}ois` → `François`, `\O rsted` → `Ørsted` (a space
  *terminates* a LaTeX letter command rather than being content). Output is
  NFC-normalized, so a served name compares equal to one typed directly
  instead of differing by normalization form. Applied to `author` and
  `journal` as well as `title`, since `format_citations` composes from all
  three.

  `corpus bib export` gained the inverse escaping, which is required rather
  than symmetric decoration: emitting a bare `&` produces a `.bib` that
  breaks the moment it reaches LaTeX, and a bare `%` is worse — it comments
  out the rest of the line, so a file corpus exported would silently lose
  fields when any real BibTeX tool re-read it. Export → import now
  round-trips a title unchanged.

  One deliberate limitation, documented in the tests: a *literal* brace
  does not survive, because brace removal is total and
  [#141](https://github.com/caseywdunn/corpus/issues/141) depends on an
  escaped brace leaving nothing behind. The two behaviors conflict, and
  OCR-emitted stray braces are far commoner in this material than titles
  that genuinely contain one.
- **`scan_detection.json` was never in the served bundle whitelist**, so
  the `scan_file_type` field the index has always exposed was silently
  `null` whenever anyone served a distilled bundle rather than the build
  tree. Found while surfacing language on the served surface.
- **Caption parsing invented figure panels out of people's initials.**
  `(A. Agassiz)` is a species authority and `Photo credit to C. Munro` is
  a credit line; both became panels — 17 bogus records across 32 papers,
  including `{'label': 'C', 'description': 'Munro'}`. An `X.` match is now
  read as an initial when a capitalised word follows *and* a credit phrase
  precedes, or an opening paren immediately precedes, or a closing paren or
  year follows the surname. A caption that genuinely opens a panel with a
  capitalised taxon name still parses.
- **`section_class`: fixed the genuine misses.** The 89% null rate is
  mostly *correct* rather than the English-only vocabulary first assumed —
  German `ZUSAMMENFASSUNG` and French `INDEX BIBLIOGRAPHIQUE` classify
  fine, and the nulls are taxonomic names, paper-specific headings,
  running heads and OCR noise, because descriptive taxonomic literature
  does not use IMRaD sections. The real misses: Russian
  `МАТЕРИАЛ И МЕТОДИКА` (singular, missed by the stricter
  `материалы и методы`), morphology/anatomy → `description` in three
  languages, and `CONCLUSIONS ON THE METHODS OF DEVELOPMENT` being
  labelled `methods` because first-match-wins ordered `methods` first.
- **A fresh corpuscle's `--dry-run` no longer prints six ERRORs.**
  `corpus run --dry-run` reported failure seconds after `corpus check`
  reported the host ready, because five guards treated "not built yet" as
  a failure; each is now dry-run-aware and reports the real plan. The
  [#139](https://github.com/caseywdunn/corpus/issues/139) taxonomy
  precondition is unchanged for real runs — a dry-run writes nothing, so
  it cannot produce the empty `taxa.json` that guard exists for. Also in
  the same pass: `--dry-run` created `output/`, `documents/` and
  `vector_db/` while printing "No files written", which additionally made
  the second dry-run behave differently from the first; `~` in config
  paths was not expanded, yielding `<corpuscle>/~/data/…` in the error
  message; the config template told you to uncomment *every* line of a
  taxonomy block that then yields an invalid `worms`+`path` mix; and
  `corpus status` printed a `--filter-gate` hint that does nothing without
  `--list-hashes`.
- **One agreed citation form** across `CITATION.cff`, its packaged copy
  and README.md (Church, Mańko, Zapata, Dunn 2026), so `--cite`,
  `--cite=bibtex` and `bundle_manifest.json` all match.
- **QC visualizations no longer die on large scans.** Pillow's
  decompression-bomb ceiling was killing visualization rendering on 4
  papers once re-OCR started rendering scans. It is raised around that
  render and restored in a `finally`, since it is process-wide state — the
  guard exists to stop a hostile upload exhausting memory, and here the
  image is one we just rendered from the operator's own PDF.
- **Every invocation no longer opens with two lines of torch internals.**
  torch's CUDA driver/runtime `UserWarning` printed ahead of any corpus
  output on every invocation *including `--help`*, while `corpus check`
  already reports GPU status properly; it is now silenced at the two sites
  that probe for an accelerator. Likewise transformers' "Token indices
  sequence length is longer than…" during chunking, which `HybridChunker`
  provokes deliberately while measuring where to split, and which reads
  like a crash.
- **The SLURM pipeline chain no longer fails silently.** The 2026-08-02
  siphonophore build stalled with stage 1 at 97.9% and stages 2–4 never
  run, from two independent faults. Stage 1 had no walltime margin: 8
  tasks × 256 PDFs against a 24 h wall meant tasks 2 and 3 drew the
  OCR-heavy scans and hit the wall at 252/256 and 225/256, while the tasks
  that finished cleared it by as little as 55 min. And `afterok` on a job
  array is all-or-nothing, so two TIMEOUT tasks had SLURM cancel Pass 3b,
  Embed and Finalize with no log, no mail and no queue entry — the build
  looked finished and nothing said otherwise. Fixes: `BATCH_SIZE` default
  256 → 64 (~3× margin at the worst observed rate, and it spreads the
  OCR-heavy tail instead of concentrating it), plus
  `--mail-type=FAIL,TIME_LIMIT` and `--open-mode=append`; `NUM_BATCHES` is
  derived from the corpuscle's own `input_pdfs` rather than defaulting to
  1, so the array can no longer be too small for the library; a
  chain-watchdog job on `afterany` always reports the array outcome and
  names the cancelled downstream jobs, and the dependency state is printed
  at submit time. Two ordering bugs went with them: `finalize` now depends
  on `afterok:EMBED:PASS3B` rather than embed alone — they are siblings, so
  `bundle` could previously start while Pass 3b was still rewriting
  `figures.json` and Pass 3c still renaming split-panel PNGs, both of which
  `mcpsrv/bundle.py` copies into `_serve/` — and Grobid moved to
  `week`/48 h, because it is submitted *before* stage 1, so an equal 24 h
  wall on `day` guaranteed Grobid died first and papers processed after
  that got placeholder metadata that implicit resume will not retry. It
  costs nothing: the `afterany` cancel job still tears Grobid down when
  stage 1 ends. `NUM_PASS3B_BATCHES` also no longer defaults to
  `NUM_BATCHES` — Pass 3b sends only multi-panel figures to the VLM, 934 of
  21,789 records on this corpus, so one GPU task suffices and matching
  stage 1's array would just queue against the 16-GPU per-user cap.
- **The `corpus_required` metadata checks compared typography, not
  meaning** ([#167](https://github.com/caseywdunn/corpus/issues/167)). All
  five failures came from two causes, neither a defect in the pipeline or
  the data. Whitespace: the checks compared BibTeX strings to docling's
  extracted text literally, so `claudanielis, a` vs `claudanielis , a` and
  `Einige histologische` vs `Einige  histologische` failed on titles that
  are plainly present. A `_norm()` now collapses whitespace and drops
  space-before-punctuation. Schneider 1891 only started failing when v1.0
  began re-OCRing scans that carry a text layer — and that re-OCR is an
  improvement (`Prof. J. Victor (Jarus` → `Prof. J. Victor Carus`), so the
  test was punishing a fix, and would equally have failed on any docling
  upgrade that shifted tokenization. Cross-script metadata: Stepanjants
  1970 records an English *translation* as its title while the body reads
  `СИФОНОФОРЫ РАЙОНА…`, is "Stepanjants" in the bib and "Степаньянц" on
  the page. Latin metadata for foreign-language papers is normal practice,
  no extraction work will make one string contain the other, and checking
  it tests transliteration, which neither these tests nor the pipeline
  do — those three checks now skip with the reason stated. The trigger is
  a *share*, not a majority: Stepanjants 1970 is 21,084 Latin characters
  against 15,075 Cyrillic, because Russian taxonomic papers are dense with
  Latin binomials and journal names, so a majority test would call it a
  Latin document and never fire. ≥20% non-Latin fires, and skip messages
  quote the measured figure (`body is 42% cyrillic`) so the reasoning is
  checkable rather than asserted.

  **The three `--deselect` flags are gone from CI.** T1, T2 and the T3
  clean-room lane each skipped these functions by name, so the repaired
  tests ran nowhere — and skipping whole functions also hid every paper
  the checks *do* work on, meaning a genuine extraction regression would
  not have been caught. The exclusions now live inside the tests, where
  they are per-paper and state their reason.
- **A job still loading Grobid's models no longer looks like a failed
  one.** SLURM reports `RUNNING` when the container starts, but Grobid
  binds `:8070` another ~10–60 s later, after loading its Wapiti models.
  Every documented readiness probe used `curl -fsS`, whose `-S`
  un-silences errors, so that normal window printed `curl: (7) …
  Connection refused` — and the runbook's `until` loop printed it once
  per retry. One Grobid job was cancelled 108 s in, mid-model-load, on
  the strength of that message. `dev_docs/BOUCHET.md` §6 now uses the
  bounded, quiet poll `slurm/batch_pipeline.sh` has always run, so the
  manual path and the orchestrator behave alike, and says which log to
  read if it genuinely times out. The same wait reached `README.md`
  (which probed a bare `curl` immediately after `docker compose up -d`,
  with no wait at all), `dev_docs/PLATFORM_SMOKE.md`, and the usage
  headers of both `slurm/` scripts, which had stopped at `export
  GROBID_URL=…`. The two `PLATFORM_SMOKE` Singularity legs also gained
  the tmp and logs binds they were missing — without them Grobid answers
  HTTP 500 to every request, or crashes on startup.
- **Grobid's per-job scratch directories no longer leak, or decide the
  job's exit status.** `slurm/batch_grobid.sh` cleaned up with `trap 'rm
  -rf …' EXIT`. Grobid still holds `grobid-service.log` open when
  SIGTERM lands, so NFS silly-renames it to `.nfsXXXX` and the `rmdir`
  fails with `Directory not empty`; as the last command of an `EXIT`
  trap under `set -euo pipefail` that also returned 1 as the job's exit
  status. Seven stale directories had accumulated under `cache/` since
  April. Retrying does not help — the file takes longer to clear than
  `KillWait` allows — so cleanup now removes the contents, always
  succeeds, and each submit sweeps directories left by earlier jobs.
  Note that `-empty` alone is not a safe liveness signal: Grobid writes
  `grobid_tmp` per request and leaves it empty when idle, so that rule
  is scoped to `grobid_logs` (where a live instance always holds
  `grobid-service.log`), with an age-based backstop for both trees.
  Leftover empty `cache/grobid_logs/<jobid>/` directories are now
  expected between submits rather than a symptom.
- **`corpus check` and `corpus run` now agree about which Grobid they
  mean.** `run` has honored `$GROBID_URL` over the config's `grobid.url`
  since [#138](https://github.com/caseywdunn/corpus/issues/138) — on HPC
  Grobid lands on a compute node whose hostname isn't known until submit
  time, so a static `config.yaml` cannot name it. `check` did not, so it
  probed `localhost:8070` and hard-failed on exactly the setup it exists
  to approve, then printed `docker compose up -d grobid` as the remedy —
  unactionable on a cluster with no Docker. `check` now applies the same
  override, and all three copies of that remediation string (both CLI
  gates and the orchestrator's warn-only path) name the cluster route as
  well: `sbatch slurm/batch_grobid.sh` then `export
  GROBID_URL=http://<node>:8070`.
- **`corpus check` validates `bib` and `lexicon` paths.** Neither was
  stat'd, so a typo'd `bib:` passed pre-flight and then exited the run
  partway through `pipeline.main`, while a typo'd `lexicon:` silently
  produced no category artifacts at all. Both are resolved exactly as
  `corpus run` resolves them: unset is informational (a supported
  configuration), a set-but-missing path is a precondition failure.
- **The bib parser no longer discards the rest of the file on one bad
  entry** ([#141](https://github.com/caseywdunn/corpus/issues/141)).
  Importing a 19,834-entry export parsed only 2,258 of them, with a single
  WARNING and a summary that looked entirely plausible. Three fixes:
  (a) **the root cause** — neither brace scanner honored backslash
  escapes, so Grobid OCR output like `author = {Des, Ej\{aims}` counted
  the escaped brace as an opening one and the entry never closed. The
  quoted-string scanner had always handled backslashes; only the brace
  path was missing it. (b) A genuinely malformed entry now costs *one
  entry*: the scan recovers at the next top-level `@` instead of stopping.
  (c) The summary reconciles parsed entries against the number of `@`
  markers in the file and warns on a shortfall, which is what makes a
  truncation visible at all. An unbalanced brace inside a *value* — which
  used to swallow the entry's remaining fields with no message
  whatsoever — now logs, naming the entry and the field.
- **Curated `.bib` entries stay out of the warning tier**
  ([#142](https://github.com/caseywdunn/corpus/issues/142),
  [#152](https://github.com/caseywdunn/corpus/issues/152), completing
  [#100](https://github.com/caseywdunn/corpus/issues/100)).
  `import_bibtex` early-`continue`d on its no-field-diff branch before
  ever calling `apply_entry`, making #100's "stamp `bib_imported_at` even
  with no diff" logic unreachable dead code — the branch worked and was
  even unit-tested, but nothing reached it. A full round-trip re-import
  therefore stamped only entries that happened to have edits (20 of
  19,834 in the reported corpus), so `format_citations` kept emitting
  "generated via reconciliation, check if correct" on hand-curated works;
  and since an authority-DB rebuild clears the stamp and forces a
  re-import, curation was defeated on every rebuild. Also collapsed a
  vestigial `return 0 if no_match == 0 else 0` whose branches were
  identical.
- **A configured-but-missing taxonomy fails instead of silently skipping**
  ([#139](https://github.com/caseywdunn/corpus/issues/139)). The first
  full production run on Bouchet put 1763 papers through with empty
  `taxa.json` and an empty `taxon_mentions.sqlite`, because
  `taxonomy.source: worms` needs outbound internet that compute nodes
  don't have and the per-paper stage merely warned. `corpus run` now
  passes `--require-taxonomy` to the extract stage whenever the corpuscle
  configures `taxonomy.source`, and that stage refuses to run rather than
  produce empty annotations, with a message naming the login-node
  pre-build and the `dwca` export. A corpuscle with no `taxonomy:` block,
  or an explicit `--no-taxa`, is unaffected — that remains a supported
  configuration. README §Taxonomy documents the internet requirement and
  the export→dwca workflow.
- **`corpus prefetch` was broken on its main path.** It called
  `emb.encode(...)`, but the embedding backend's method is `embed` — so a
  genuine cold prefetch downloaded all 4.8 GB, *then* failed with an
  `AttributeError`, then retried it six times with backoff before
  reporting "failed after 6 attempts" as if the network were at fault. Two
  fixes: the correct method, and `_with_retry` no longer retries
  programming errors (`AttributeError`/`TypeError`/`NameError`/
  `ImportError`), which it cannot fix and which retrying only disguises.
  A test now binds the call to the real `EmbeddingBackend` API, since
  every existing test mocked the backend and the clean-room lane exercises
  the cold download through `corpus run` rather than `corpus prefetch` —
  so nothing caught it.
- **A correct first run no longer ends in warnings.** Building a corpuscle
  with no `taxonomy:` block — a documented, supported setup — emitted two
  `WARNING`s about `taxonomy.sqlite` and `instructions.md` being absent
  from the bundle, and `corpus status` marked the taxonomy `✗` alongside
  real artifacts. Both are optional by design (and since #139 a
  *configured* but missing taxonomy fails the run outright, so absence
  here means "not wanted"). They are now an `INFO` line and a `–  (optional
  — no taxonomy: block configured)` marker. The `✗` is reserved for gaps
  that matter.
- **Run logs no longer name scripts that don't exist.** Ten loggers were
  named after root-level scripts removed in v0.3 (`package_for_serve`,
  `build_taxon_mentions`, `ingest_taxonomy`, …), so every line of a run log
  pointed at a filename a user could not find. Renamed to `corpus.*`.
- **Documentation swept against the code.** A review after the cycle
  landed found docs describing removed or renamed things:
  `dev_docs/LICENSING.md` still taught the `publishable` wire field that
  #154 removed; `dev_docs/MCP_TOOLS.md` described four figure tools as
  they were before #154/#143; and `dev_docs/OVERVIEW.md` handed newcomers
  a file-inventory table listing eleven root-level scripts that were
  folded into packages back in v0.3 (#60), two of them as broken links.
  All corrected, and every relative link in the docs now resolves.
- **Documentation swept again after the OCR and prefetch changes**, since
  four of the five defects that pass found were created by those changes.
  `install_ocr_extras.sh` had reached the README install block only, so
  every other install recipe still omitted it — INSTALL.md's cluster
  block, BOUCHET.md §2, PLATFORM_SMOKE.md, and the clean-install
  walkthrough — and the cluster ones matter most, because that is where
  the disk difference lands; BOUCHET.md was also missing `pip install -e .`
  entirely, and still carried the false "ships only English LSTM data"
  tesseract claim corrected everywhere else, in the one document someone
  follows while standing on the cluster. README said "born-digital PDFs
  with a clean text layer skip OCR entirely", now actively misleading: the
  classification is geometric, so a scan arriving with third-party OCR
  baked in is still a scan and gets re-OCR'd. OVERVIEW.md described scan
  detection without its primary signal, and had none of `--redo-ocr`, the
  OCR language probe, or the scaling timeouts. Checked mechanically and
  found clean: relative links (0 broken), documented config keys (all 7
  valid against the pydantic schema), documented CLI verbs against the
  real subcommand set, and the demo paper counts — including the README's
  new claim that the demo OCRs three of its four papers, verified against
  `scan_detection.json` (Pugh, Stepanjants and Schneider are
  `raster_page_images`; Marrus alone is `born_digital`).
- **`docker-compose.yml`'s Singularity hint was missing `--pwd
  /opt/grobid`.** Docker honours the image `WORKDIR`; Singularity starts
  in the host cwd instead, so the line as written dies before the service
  starts — `[FATAL tini] exec ./grobid-service/bin/grobid-service failed:
  No such file or directory`. Hit for real on Bouchet.
  `slurm/batch_grobid.sh` and PLATFORM_SMOKE.md both passed `--pwd`
  already; only this comment was wrong, and it is the one a reader is most
  likely to copy. BOUCHET.md had described the hand-rolled failure as a
  service that "starts and then fails every request" — that is the symptom
  *after* `--pwd` is added, while the first failure is that it never
  starts. Both are now spelled out with the error text.
- **Stopped hard-coding the siphonophore library's size.** The library
  grows as papers are added, so every exact count baked into the docs is
  wrong by the next `git lfs pull`; BOUCHET.md alone asserted three
  different figures. Replaced with phrasing that scales, and the rule is
  recorded in AGENTS.md. Two corrections worth stating rather than
  burying: the header now says how to get the real count, because
  `ls … | wc -l` returns 17 — the library is nested in letter
  subdirectories — so it takes `find … -iname '*.pdf' | wc -l`, or just
  `corpus check`, which reports it. And `NUM_BATCHES` is derived from
  corpus size rather than being a fixed tuning constant, so it cannot be
  left alone: `pipeline/main.py` slices `all_hashes[start:end]`, so a
  surplus array task gets an empty slice and exits immediately while a
  shortfall leaves papers unprocessed **with no error** — asymmetric, so
  the docs say to round up and explain why.
- **BOUCHET.md no longer calls the sample build a "dry run."** It used the
  phrase for two opposite things within one screen: `corpus run --dry-run`
  prints the plan without writing artifacts, while eight lines later "Dry
  run (20–50 papers) before the full corpus" headed a section that creates
  a corpuscle, submits real CPU and GPU jobs, and writes the artifacts you
  then inspect. A reader coming off step 7 could reasonably take the
  section below it to write nothing too, when it is in fact the last
  checkpoint before committing the full library. Renamed to "smoke test",
  which its own first sentence already called it.
- **`dev_docs/UX_REVIEW_v1.0.md` removed.** The pre-release review is
  retained outside the repo; its actionable findings are filed as
  [#164](https://github.com/caseywdunn/corpus/issues/164)–[#173](https://github.com/caseywdunn/corpus/issues/173),
  each carrying its own reproducing evidence, and the fixes it drove are
  in this section. Three of its findings were **withdrawn** after
  measuring against the artifacts rather than trusting an MCP client's
  inference, and three of those four would have been actively harmful to
  "fix": reconciliation is conservative by design and its near-misses are
  correct rejections of same-author/same-year different works; figure
  licensing derives public-domain status from year, so `license: null` on a
  pre-1931 work is fully cleared rather than unknown; and `section_class`'s
  high null rate is mostly correct, because descriptive taxonomic
  literature does not use IMRaD headings.
- **Dead script names removed from user-visible messages.** A served MCP
  error payload told users to run `python embed_chunks.py <output_dir>`
  (removed in v0.3) — it now says `corpus run --only embed`. Likewise the
  serve-time taxon-mention warning, twelve `corpus run` help strings, and
  the orchestrator's logger name, which labelled every run log line
  `update_corpus` after a script of that name.
- **`corpus prefetch` no longer warns about `HF_HOME` when it has nothing
  to download.** The nag now fires only when models will actually be
  written.
- **README's disk estimate includes the models.** It said "several times
  the size of the original PDFs", which omits the ~5 GB of models that
  dominate for a small corpus.
- **`corpus check` no longer greener than `corpus run`.** It reported
  "ready" for a zero-PDF `input_pdfs`, which `corpus run` refuses outright —
  and run's refusal points the user at "`corpus check` for the full
  pre-flight surface", which was then less complete than the thing it
  deferred to. Now a precondition failure with the same wording and the same
  `--skip-checks` escape hatch. Found by re-running the T4 operator
  walkthrough after this cycle's CLI changes.
- **Grobid's compose healthcheck actually works now**
  ([#157](https://github.com/caseywdunn/corpus/issues/157)). It shelled
  out to `curl`, which the `lfoppiano/grobid:0.8.1` image does not ship
  (nor `wget`, `nc`, or `python3` — verified by `docker exec`), so
  `corpus-grobid` reported `(unhealthy)` forever while serving perfectly.
  It now does a real HTTP GET of `/api/isalive` through bash's
  `/dev/tcp`, greps the response for `true`, and names `bash` explicitly
  because `CMD-SHELL` runs `/bin/sh`, which in this image is dash and has
  no `/dev/tcp`. `docker inspect` reports `healthy`.
- **docling's picture classifier is no longer downloaded or run**
  ([#140](https://github.com/caseywdunn/corpus/issues/140)).
  `do_picture_classification` was `True`, which fetched
  `DocumentFigureClassifier-v2.5` from HuggingFace and ran it on every
  figure — but nothing read the result: `figure_type` comes from our own
  `classify_figure()` heuristic and the served vocabulary is corpus's own
  (`figure` / `plate` / `subpanel`). One fewer model download on a new
  user's first run, and one fewer HuggingFace 429 surface (this model was
  the first fetch to fail in the #140 CI outage). The layout model and
  TableFormer are still downloaded. Existing corpuscles re-run the
  figure/extract stage on the next `corpus run`.
- **Test dependencies are no longer runtime dependencies**
  ([#162](https://github.com/caseywdunn/corpus/issues/162)). `pytest`,
  `pyflakes`, and `ipykernel` moved from `[project].dependencies` to a
  `[project.optional-dependencies].dev` extra, so installing corpus no
  longer drags in a test runner and a Jupyter kernel. Development clones
  want `pip install -e ".[dev]"`; the conda env supplies them
  unconditionally as before.
- **`requires-python` is bounded to `>=3.12,<3.13`**
  ([#163](https://github.com/caseywdunn/corpus/issues/163)). It was
  unbounded while `environment.yaml` pins 3.12, CI tests only 3.12, and
  the #98 known-good ML set was verified against nothing else — so pip
  was free to attempt an untested 3.13+ install via the documented
  pip-only path.
- **Shell completions offer every verb.** The three hand-maintained
  completion snippets had drifted: `taxonomy` was missing from all of
  them. Adding `prefetch` surfaced it, and a new test now compares each
  snippet against the live argparse subcommand list so it can't drift
  again.
- **Corrected a wrong claim about CI caching** that had propagated into
  #156, #158, and `CONTRIBUTING.md`. `setup-miniconda` does not cache the
  solved conda environment and this repo adds no `actions/cache` for it,
  so T0/T1/T2 have always re-resolved `environment.yaml` from scratch.
  The reason `mcp` 2.0.0 went unnoticed for 24 days was not a stale
  cache: it was that CI only runs on push and nobody pushed to `dev`
  between 2026-07-06 and 2026-07-30. The first push after the release
  failed immediately.

## [0.6.0] - 2026-06-04

### Theme — v0.6 road-to-1.0

v0.6 is the **API-freeze** cycle: no new feature tools, instead a
one-time pass to finalize the public MCP tool surface, fix known
correctness bugs, and harden ops — because after 1.0 a change to any
tool default, signature, or response shape is a breaking change. The
served surface is now **38 MCP tools** with a uniform error payload and
a consistent pagination convention; this is the surface 1.0 freezes.

### Changed (breaking)

- **MCP surface frozen at 38 tools.** The redundant singular tools were
  removed in favor of their batched plurals
  ([#88](https://github.com/caseywdunn/corpus/issues/88) §2.3):
  `format_citation` → `format_citations`, `get_paper` → `get_papers`,
  `get_chunk` → `get_chunks` (pass a single-element list for the
  one-item case). `mcpsrv/default_instructions.md` routes citations
  through `format_citations`.
- **Removed `translate_chunk`** ([#124](https://github.com/caseywdunn/corpus/issues/124)).
  It was the only server-side LLM-call tool; the MCP server is now a
  purely deterministic retrieval layer. Translation is an analysis step
  an MCP client (itself an LLM) does on chunk text it already retrieved
  via `get_chunks` — it doesn't need a server-side Claude call, key, or
  cache. Removed pre-1.0 so it isn't frozen into the stable surface.
- **Breaking response-shape defaults**
  ([#88](https://github.com/caseywdunn/corpus/issues/88) Part 1):
  `lexicon_matrix` returns per-term **summaries by default**
  (`detail=True` for the full grid — fixes a 71 MB → 684 KB runaway);
  the older figure tools (`get_figures_for_taxon`,
  `get_figures_for_lexicon_term`) return a **caption preview** by
  default (`full_caption=True` for the verbatim caption,
  [#85](https://github.com/caseywdunn/corpus/issues/85));
  `get_chunks_for_topic` gains `with_cites=True` in-text cite refs.
- **Uniform tool error payload** (freeze gate,
  [#88](https://github.com/caseywdunn/corpus/issues/88) §3). Every tool
  error return now carries a human `error` message **plus a machine
  `code`** (`not_found` / `ambiguous` / `invalid_argument` /
  `not_configured` / `no_results` / `unavailable` / `empty_item` /
  `forbidden`). Clients branching on the citation `error` token
  (`"not_found"` / `"ambiguous"`) must branch on `code` instead.
- **Pagination naming** reconciled to `limit` everywhere, except the
  multi-section dossier/graph tools, which keep descriptive `max_*`
  caps (one call, several independently-bounded sections — the
  documented exception). A `tests/test_freeze_contract.py` meta-test
  enforces both the 38-tool set and the naming rule.
- **Output-type profiles replace `--allow-unpublishable`**
  ([#101](https://github.com/caseywdunn/corpus/issues/101)). Figure +
  citation gating is now driven by a per-call `profile=` arg
  (`report` / `manuscript` / `presentation`); the server
  `--default-profile` sets the fallback for calls that omit it (default
  `report`, permissive — startup warns). New `list_output_profiles` /
  `get_active_profile` discovery tools. Figure responses gain
  `license` / `attribution` fields.
- **Figure panel detection: `--figure-panels` + OCR floor default-on**
  ([#102](https://github.com/caseywdunn/corpus/issues/102)). The
  `vision:` config block became `figures:`; `--vision-backend` +
  `--content-aware-figures` collapsed into one
  `--figure-panels {ocr,vision-local,vision-claude,off}` (default
  `ocr`, a CPU-only floor). A legacy `vision:` block now fails config
  validation with the migration mapping. Existing corpuscles re-run the
  figure stage once on the next `corpus run`. `corpus run` accepts
  `--figure-panels` as a per-run override of `figures.panel_detection`;
  `--no-vision` is now a deprecated alias for `--figure-panels ocr`.
- **Figure resolution: native per-figure DPI by default**
  ([#121](https://github.com/caseywdunn/corpus/issues/121)). Figures
  were rasterized at docling's 72 dpi default and looked grainy in
  print. Extraction now defaults to `figures.resolution_mode: native`:
  each figure is rendered at its **source's native pixel density** (a
  600-dpi scan figure stays 600 dpi — resolution varies per figure),
  with `figures.vector_dpi` (default 300) for resolution-less vector
  figures and an optional `figures.max_dpi` cap. `resolution_mode:
  fixed` keeps a uniform `figures.images_scale` (default 2.0 → 144 dpi)
  instead. Affects future ingests only; lift an existing bundle with
  `backfill_figure_dpi.py --native`.

### Added

- **`format_citations`** batched citation formatter and **`search_taxon`
  `parent_chain`** ancestry walk
  ([#88](https://github.com/caseywdunn/corpus/issues/88)).
- **Breadth + edge caps on `get_citation_graph`**
  ([#87](https://github.com/caseywdunn/corpus/issues/87)):
  `max_edges_per_node` / `max_total_edges` + a `truncated` flag.
- **`/healthz` capability report + refuse-to-serve on degraded
  capability** ([#91](https://github.com/caseywdunn/corpus/issues/91)).
  `/healthz` returns JSON capability flags and **503** when a backing
  index is degraded; `get_chunks_for_topic` raises a hard error (not
  empty rows) on a degraded semantic index.
- **Central per-invocation run log**
  ([#90](https://github.com/caseywdunn/corpus/issues/90)):
  `<output_dir>/runs/<timestamp>/run.log` archives argv, resolved
  config, dependency-stack versions, and stage success/failure counts
  (top-level `run.log` kept as "latest").
- **`corpus debug-pdf`** single-PDF debug runner with per-stage tracing
  ([#92](https://github.com/caseywdunn/corpus/issues/92)).
- **Per-tool instrumentation shim** in `mcpsrv/app.py` (call/error/
  latency counters) feeding `/healthz` and the run log.
- **`backfill_figure_dpi.py`** re-renders an existing bundle's figures
  in place from their stored bbox + source PDF, at a fixed `--scale` or
  `--native` (per-figure source DPI), without re-running docling
  ([#121](https://github.com/caseywdunn/corpus/issues/121)).

### Fixed

- **Bib-provenance preserved through import + reconciliation**
  ([#100](https://github.com/caseywdunn/corpus/issues/100)): a
  reference in the user-edited `.bib` stays authoritative `bib`
  provenance after an unchanged re-import and after reconciliation, so
  curated references no longer warn spuriously.
- **`build_taxon_mentions` freshness is fingerprint-aware**
  ([#95](https://github.com/caseywdunn/corpus/issues/95)): re-ingest is
  gated on the taxonomy sha recorded in each paper's taxa-stage
  fingerprint, not just `taxa.json` mtime (unreliable across HPC nodes).
- **HuggingFace implicit-token warning** silenced via a shared
  `HF_HUB_DISABLE_IMPLICIT_TOKEN` setter
  ([#97](https://github.com/caseywdunn/corpus/issues/97)).
- **WoRMS `isMarine=0` gap documented**
  ([#96](https://github.com/caseywdunn/corpus/issues/96)).
- **Diacritic author lookup**
  ([#122](https://github.com/caseywdunn/corpus/issues/122)):
  `get_works_by_author("Müller")` returned 0 while `"Muller"` returned
  the real set — the `surname_normalized` column is diacritic-stripped
  at build but was queried with a plain `.strip().lower()`. Every author
  site now queries with the same `normalize_for_key` normalizer (`Müller
  ≡ Muller`; `Mueller` stays distinct). Serve-side fix — no rebuild.

### Migration (v0.5 → v0.6)

| Old | New |
| --- | --- |
| `format_citation(...)` | `format_citations(queries=[...] / work_ids=[...] / paper_hashes=[...])` |
| `get_paper(hash)` | `get_papers(hashes=[hash])` |
| `get_chunk(hash, chunk_id)` | `get_chunks(hash, chunk_ids=[chunk_id])` |
| `translate_chunk(hash, chunk_id, target_language)` | removed — translate the text from `get_chunks` client-side |
| citation error token `{"error": "not_found"}` | `{"error": <message>, "code": "not_found"}` (branch on `code`) |
| `lexicon_matrix()` full grid | summary by default; `lexicon_matrix(detail=True)` for the grid |
| `get_figures_for_taxon` full `caption_text` | caption preview; `full_caption=True` for verbatim |
| server flag `--allow-unpublishable` | server `--default-profile` + per-call `profile=` |
| config `vision.backend: {none,local,claude}` | config `figures.panel_detection: {off,ocr,vision-local,vision-claude}` (default `ocr`) |
| CLI `--vision-backend` / `--content-aware-figures` | CLI `--figure-panels {ocr,vision-local,vision-claude,off}` (on both `corpus run` and `pipeline.main`); `--no-vision` → alias for `--figure-panels ocr` |
| figures rendered at a fixed 72 dpi | `figures.resolution_mode: native` (default) — per-figure source DPI; `fixed` + `figures.images_scale` for uniform DPI; `backfill_figure_dpi.py --native` lifts existing bundles |

## [0.5.0] - 2026-05-29

### Theme

v0.5 is the **served-bundle quality** cycle. v0.4 made the pipeline that
*produces* a corpus reliable; v0.5 turns to the surface a downstream LLM
client *consumes* — the MCP tool surface — and closes the two
trust-breaking gaps external evaluators surfaced there:

1. **Citation grounding** ([#79](https://github.com/caseywdunn/corpus/issues/79)) —
   amalgamated citations the model assembled in its own context window
   (right title, wrong author/journal). The fix moves citation
   formatting server-side so the model pastes a single authoritative
   string instead of recombining structured fields.
2. **Cache-friendly dossier tools** ([#76](https://github.com/caseywdunn/corpus/issues/76)) —
   the `for each X, call get_Y` access pattern burned ~1.95 M
   cache_read tokens on enumerative prompts. Pre-aggregated dossier
   tools return structured indices + IDs + headings, with full text
   pulled selectively (~85 % token / ~92 % cache_read reduction on the
   worked p02 example).

Non-breaking on the pipeline side — v0.4.x corpuscles build as-is. The
MCP surface gains tools (`format_citation`, six dossier tools) and
additive token-saving flags but removes or renames nothing; clients
that ignore the new surface keep working.

### Added

- **`format_citation` MCP tool + provenance cascade**
  ([#79](https://github.com/caseywdunn/corpus/issues/79)). Formats a
  citation server-side from the authoritative `biblio_authority.sqlite`
  record and returns `formatted` / `inline` strings the client pastes
  verbatim, plus a provenance-tier `warning` footnote. Stops the
  client-side recombination where amalgamation happened. Shipped across
  schema + importer stamp (`219a4a2`), `bib/format.py` author-year
  formatter (`4d16064`), `BiblioAuthority.provenance()` tier classifier
  (`fc3a7e9`), and the tool wiring (`e8f075a`). `vancouver` / `bibtex`
  styles + a `bibliography.citation_style` config field are deferred to
  follow-ups until a second style exists.
- **Citation-grounding instructions rewrite + quality harness**
  (part of [#79](https://github.com/caseywdunn/corpus/issues/79)).
  `mcpsrv/default_instructions.md` rewrote the bibliography section as
  an explicit routing rule (call `format_citation`, paste `formatted` /
  `inline` and preserve `warning` verbatim, handle `not_found` /
  `ambiguous` without fabricating) — `1d8b485`. A real Claude-API
  round-trip harness (`tests/test_prompt_quality.py`) scores expected
  citation emission + parenthetical trace-back + forbidden-hallucination
  windows, gated behind `RUN_PROMPT_QUALITY=1` + `ANTHROPIC_API_KEY`;
  scoring functions have always-run in-process unit tests. The new
  `prompts:` ground-truth YAML block is documented in
  `dev_docs/TESTING.md`.
- **Dossier tool suite — six cache-friendly aggregate tools**
  ([#76](https://github.com/caseywdunn/corpus/issues/76)).
  `corpus_summary` (`ed6a229`); `get_taxon_dossier` + the `get_chunks`
  workhorse pair (`a007da7`); `get_taxon_lexicon_slice` (`db767e7`);
  `lexicon_matrix` + `get_lexicon_term_dossier` (`b7195d7`); the
  figure-dossier pair (`05e775c`); `get_papers` +
  `get_taxon_subtree_dossier` (`1d5e6ed`). MCP tool count 27 → 38;
  100 unit tests across the six groups; documented in
  `dev_docs/MCP_TOOLS.md`.
- **`get_chunks_for_topic` `with_text=False` mode**
  ([#82](https://github.com/caseywdunn/corpus/issues/82)). Metadata-only
  semantic-search response (~80 chars/row vs ~600) — scan IDs + score,
  then drill in via `get_chunks(paper_hash, chunk_ids=[...])`. Projected
  ~45 % cut on typical k=10 scan-then-focus workflows.
  Backwards-compatible (default stays `True`). `40d56d1`.
- **`with_text=False` on `get_chunks_by_section` + `get_chunks_for_taxon`**
  ([#84](https://github.com/caseywdunn/corpus/issues/84)). Mirrors the
  #82 pattern on the two remaining chunk surfaces. The new flag drops
  the chunk `text` field and emits `len_chars` instead, so a caller
  can scan IDs + section + headings cheaply and drill in with
  `get_chunks(paper_hash, chunk_ids=[...])`. Default behaviour
  unchanged — additive flag, no contract break.
- **Success summary on a clean `corpus run`**
  ([#57](https://github.com/caseywdunn/corpus/issues/57)). New
  `_write_run_log` helper writes the `corpus status --report` content to
  `<output_dir>/run.log` at the tail of `corpus run`; the trailing
  success message points at the file. Overwrites the prior log; dry-run
  short-circuits before writing. `d365792`.

### Changed

- **Long `@mcp.tool()` docstrings trimmed + prompt-cache docs**
  ([#81](https://github.com/caseywdunn/corpus/issues/81)). Trimmed the
  eight longest tool docstrings (~25 kc → ~19.6 kc, ~1.4 k tokens saved
  per uncached session, ~20 % catalog reduction) and added a
  `dev_docs/MCP_TOOLS.md` section showing Anthropic API clients where to
  place `cache_control` breakpoints for the ~5 k-token static prefix.
  `0dbc0ca`.
- **MCP tools now reject `limit < 1` instead of treating `limit=0` as
  unlimited** ([#86](https://github.com/caseywdunn/corpus/issues/86)).
  `get_chunks_by_section`, `get_chunks_for_taxon`, `get_taxon_mentions`,
  `get_figures_for_taxon`, and `get_figures_for_lexicon_term` used to
  return the entire candidate set when called with `limit=0` (the
  `rows[:limit] if limit else rows` pattern). A client passing
  `limit=0` reasonably expects "zero rows" or an error — getting
  back a 2000-row response was a footgun. The tools now return
  `[{"error": "limit must be >= 1 …"}]` for `limit < 1`, and silently
  clamp absurdly large limits to a per-tool ceiling
  (`mcpsrv.app.MAX_LIMIT = 500`). Existing callers passing positive
  limits keep working unchanged.
- **`taxonomy.root_id` accepts the DwC-A LSID string format**
  ([#78](https://github.com/caseywdunn/corpus/issues/78)). Schema
  relaxed to `Optional[Union[int, str]]`: `worms` keeps the integer
  AphiaID path; `dwc` / `dwca` now accept either a bare integer taxonID
  or an LSID string, matching what the CLI already accepted. Downstream
  code already wrapped with `str(...)`, so no callers changed.
  Regression tests pin both shapes. `7be13f8`.
- **Install docs reframed around arm64-native conda**
  ([#77](https://github.com/caseywdunn/corpus/issues/77)). README +
  INSTALL.md: the requirement is an arm64-native conda (Anaconda,
  Miniconda, or Miniforge all qualify), not Miniforge specifically. New
  `conda info | grep platform` pre-gate; `corpus check` now hard-fails
  on Rosetta'd Python on macOS. `CONDA_SUBDIR=osx-arm64` documented as
  the alternative. Installation section restructured (supported
  platforms first, macOS items consolidated, Grobid startup moved to its
  own section as a run-time concern). `b92e355` + follow-ups.

### Fixed

- **Lexicon YAML loader silently degraded on malformed input.**
  `load_lexicon` now rejects list-shaped term entries with a
  `ValueError` naming the category, term, offending type, and corrected
  shape; `main.py` narrowed its exception swallow to file / YAML errors
  so schema violations propagate instead of degrading to no-op
  annotation. `f4b36b0`.
- **Code-review batch — release-blocker + 8 follow-ups** (`09758a0`).
  Addressed the findings from the pre-release Codex review pass across
  the new MCP surface.
- **Pinned the docling/torch/transformers/sentence-transformers stack**
  ([#98](https://github.com/caseywdunn/corpus/issues/98)). The
  extraction/ML stack was unpinned in `environment.yaml`,
  `requirements.txt`, and `pyproject.toml`, so every build resolved the
  latest PyPI wheels. `docling` 2.95/2.96 silently break extraction on
  macOS arm64 (MPS) — empty output, no crash — which surfaced as a
  reproducible T2 failure with no code change. Pinned to the last-green
  set (`docling==2.94.0`, `transformers==5.8.1`, `torch==2.12.0`,
  `sentence-transformers==5.5.0`). A deliberate forward-move is tracked
  in the issue.
- **`corpus run` now fails on a zero-yield extraction**
  ([#99](https://github.com/caseywdunn/corpus/issues/99)). When every
  processed document produced zero chunks (e.g. a silently-failing
  extraction backend), the run previously logged "all steps succeeded"
  and packaged an empty served bundle. The extract stage now exits
  non-zero on a corpus-wide zero-chunk result and warns on individual
  empty documents, so the failure is loud instead of shipped.

## [0.4.0] - 2026-05-13

### Theme

v0.4 is the **operational hardening** cycle. v0.3 collapsed the CLI
surface and the per-corpuscle config; v0.4 takes the resulting machine
through CI it didn't have, a platform-portability pass that fixes the
silent failure modes external HPC users hit, and the install-onboarding
papercuts that first-time operators raised.

Three through-lines:

1. **CI tiers** ([#75](https://github.com/caseywdunn/corpus/issues/75)) —
   GitHub Actions on every push + every PR. T0 (lint + unit) takes ~3 min;
   T1 (Linux + Grobid) and T2 (macOS arm64) each build the 4 + 1 demo
   and exercise implicit resume inline (~13 / 9 min). The classes of
   regression that previously surfaced weeks later in PLATFORM_SMOKE.md
   now block PRs within ~15 min.
2. **Silent-failure cleanup on HPC + Apple Silicon** —
   [#70](https://github.com/caseywdunn/corpus/issues/70),
   [#71](https://github.com/caseywdunn/corpus/issues/71),
   [#72](https://github.com/caseywdunn/corpus/issues/72), plus a
   macOS-specific KMP / Rosetta pass. All four had the same shape:
   the pipeline kept running but quietly lost output (papers, the
   served bundle, or the entire resume mechanism).
3. **Install-onboarding fixes from external testers** —
   [#73](https://github.com/caseywdunn/corpus/issues/73) (install
   ordering + config template indentation).

No breaking changes for end users; v0.3.x corpuscles work as-is.

### Added

- **Tiered CI on GitHub Actions**
  ([#75](https://github.com/caseywdunn/corpus/issues/75)).
  Two workflows; four functional tiers; both fire on every push and
  every PR (no branch filter — topic branches get the same signal
  the daily-cadence `dev` line does).
  - **T0 — lint + unit** (`.github/workflows/lint.yml`). pyflakes
    gate via `tests/test_no_undefined_names.py` + the ~314 unit
    tests that don't need a built corpus. ~3 min wall time.
  - **T1 — Linux integration** (`.github/workflows/integration.yml`,
    `ubuntu-24.04`). Grobid as a service container, `corpus run`
    on the 4-paper demo, bundle-manifest checks, audit-clean check,
    MCP / SSE round-trip via `tools/smoke_test_sse.py`, then the
    4 + 1 implicit-resume scenario inline (copy
    `tests/fixtures/round2_paper/Siebert_etal2011.pdf` into demo/,
    re-run, assert `skipped == 4 ∧ embedded ≥ 1 ∧ failed == 0` —
    canonical regression check for #71). ~13 min wall.
  - **T2 — macOS arm64** (`macos-15`). Same shape as T1 with
    `grobid.disable: true` (Docker Desktop isn't on GHA macOS
    runners); reference-extraction tests skipped (Linux T1 covers
    the Grobid-XML path). Catches macOS-arm64-specific regressions
    including platform-specific resume behavior. ~9 min wall.
  - Two manual tiers in dev_docs: T3 clean-room EC2
    (`dev_docs/ec2_smoke.sh`, ~25 min, pre-release ritual) and T4
    operator walkthrough (`dev_docs/clean_install_walkthrough.sh`,
    full CLI surface coverage).
  - **Pytest markers**: `corpus_required` (~73 ground-truth tests
    that need a built demo; T0 deselects, T1/T2 select after the
    build) and `resume_scenario` (the standalone pytest port of
    the 4 + 1 scenario, kept for local-dev iteration; CI does the
    equivalent via shell assertions inline rather than running the
    pytest twice).
  - **Status badges**: README scoped to `main` (release-state
    signal users care about); CONTRIBUTING scoped to `dev` (live
    development-line signal contributors want).
  - **Release ritual gated on green CI**
    ([CONTRIBUTING.md](CONTRIBUTING.md)). Two new checkpoints: dev
    CI green before bumping the version, main CI green before
    tagging. The dev → main merge is the last point where a regression
    can be caught before it gets stamped into `bundle_manifest.json`.

- **`corpus taxonomy export | ingest` — DwC-A round-trip**
  (commit `23448f0`). The reverse of the existing ingest is now a
  first-class verb: dump a built `taxonomy.sqlite` back out as a
  Darwin Core Archive `.zip`. Use cases: share a snapshot without
  forcing the recipient to walk WoRMS again; commit a fixture into
  a downstream repo so CI exercises the `dwca` ingest path without
  network calls. Round-trip property:
  `corpus taxonomy ingest --source dwca --input <export.zip>`
  recovers the same `taxa` row set as the source SQLite.

- **Demo ships a pre-built Siphonophorae DwC-A** (commit `447437f`).
  `demo/taxonomy.zip` is the order-level WoRMS snapshot
  (~1,200 taxa, ~30 KB) baked via the new `corpus taxonomy export`.
  The demo's `config.yaml` switches to
  `source: dwca, path: ./taxonomy.zip`. First `corpus run` on the
  demo ingests the taxonomy from a local file in seconds — no
  rate-limited WoRMS REST walk. Re-export with the documented
  `pipeline.taxonomy_ingest` + `corpus taxonomy export` two-step
  in the demo's `config.yaml` comment block.

- **pyflakes static-analysis gate**
  (`tests/test_no_undefined_names.py`, commit `721ed5d`). Catches
  undefined-name (guaranteed runtime `NameError`) findings on the
  source tree, runs in seconds, gates T0. Surfaced and fixed nine
  pre-existing missing-imports + bare-typing bugs as part of the
  initial adoption.

- **EC2 clean-room smoke script** (`dev_docs/ec2_smoke.sh`, commit
  `4c86104`). One-shot platform-portability check on a fresh Ubuntu
  EC2 host: apt deps + miniforge + conda env + `pip install -e .`
  + Grobid via Docker + demo `corpus run` + bundle distillation +
  SSE round-trip. Exits 0 iff every success criterion in
  `PLATFORM_SMOKE.md` passes. ~25 min wall, ~$2–3 EC2 cost. The
  release ritual now runs this as the clean-cache counterpart to
  GHA's warm-cache T1.

- **Platform-portability smoke runbook**
  (`dev_docs/PLATFORM_SMOKE.md`, commit `a245d0e`, narrowed in
  `5da9ff4`). Pre-release manual smoke against macOS arm64 +
  linux-x86_64 Bouchet with explicit success criteria, mirrored
  programmatically by `ec2_smoke.sh`. Retitled in v0.4 to "manual
  fallback / release-time verification" with a banner pointing at
  the CI tiers as the authoritative coverage.

- **README "Using the MCP server" section.** Covers the
  text-only-chat vs. report-generation split: chat for exploration,
  report path (`pandoc` / `pdflatex` for PDFs with figures) for
  morphological-diversity summaries, character matrices, and CSV /
  TSV / JSONL / BibTeX exports. New example query: *Write a PDF
  report with LaTeX showing all nectophore images for* Nanomia.

- **World Flora Online** added to the DwC-A source table in the
  README, with the Zenodo DOI for the recent snapshot — fills the
  plant-taxonomy gap left by WoRMS / Catalogue of Life pointers.

### Changed

- **Demo slimmed from 11 → 4 + 1 papers** (commit `d597f1f`).
  Four papers in `demo/` (born-digital English Dunn-etal-2005,
  born-digital English Pugh-2001, scanned German Fraktur
  Schneider-1891, scanned Russian Stepanjants-1970) for the
  standard build path; one paper (Siebert-etal-2011) held back in
  `tests/fixtures/round2_paper/` for the 4 + 1 implicit-resume
  scenario. Round-1 wall time dropped from ~25 min to ~2 min on
  macOS arm64 (and ~5 min in CI), which is what made the demo
  viable as a CI fixture in the first place. Still exercises every
  OCR + extraction code path the 11-paper version did (Fraktur,
  Cyrillic, born-digital).

- **README Installation section reordered**
  ([#73](https://github.com/caseywdunn/corpus/issues/73), commit
  `4e473f0`). New "Prerequisites" subsection up front lists Docker
  and miniforge before any conda command. New "Clone and install"
  block leads with `git clone` (previously omitted — operators
  had to infer they needed to clone the repo before running
  `conda env create`) and folds in `tools/install_tessdata.sh` as
  the fourth step. Drops the now-redundant "Docker is a
  prerequisite" paragraph that previously appeared after the
  install commands.

- **Config template — explicit indentation note above the
  `taxonomy:` block** (`pipeline/config.template.yaml`,
  [#73](https://github.com/caseywdunn/corpus/issues/73)). The
  2-space nesting tripped beroe's smoke install: stripping the
  leading `# ` from `#   source: worms` without preserving the
  2-space indent put `source:` at column 0 and the YAML loader
  raised "block mapping" at parse time. A new paragraph above the
  block spells out the indentation contract.

- **CI tier structure documented in CONTRIBUTING, not README**
  (commits `1cfed8b`, `b67f3c4`). User-facing README keeps the
  status badges (scoped to `main`); CONTRIBUTING gets the full
  tier table + local equivalents (right next to "Tests" and "What
  to run before opening a PR"), plus its own badges scoped to
  `dev` for the live signal contributors want.

### Fixed

- **`find_all_pdfs` re-ingesting `processed.pdf` from prior runs**
  ([#75](https://github.com/caseywdunn/corpus/issues/75), commit
  `5373623`). `pipeline/io.py:find_all_pdfs` did a bare
  `input_dir.rglob("*.pdf")` with no output-directory exclusion.
  When `input_pdfs` is an ancestor of `output_dir` (the demo's
  canonical `input_pdfs: .` with `output_dir: ./output`), every
  subsequent `corpus run` walked the per-paper
  `documents/<HASH>/processed.pdf` artifacts. Those have different
  SHA-256s than the originals because OCR overlays a text layer,
  so the pipeline counted them as new documents and the corpus
  doubled on every re-run. New `exclude_under` keyword param
  plumbed through from `pipeline.main`; `corpus check` already
  used the equivalent helper. Unit test in
  `tests/test_find_all_pdfs.py`. T1/T2 CI now exercises round-2
  against the demo's literal config, so a regression here blocks
  PRs.

- **LanceDB `list_tables()` return-shape break**
  ([#71](https://github.com/caseywdunn/corpus/issues/71), commit
  `f81a4b8`). LanceDB 0.30.x changed `list_tables()` from a list to
  a generator-like view. The old `args.table_name in table_names`
  membership check consumed the generator on first read and
  returned False on the second, so `pipeline.embed` then tried to
  re-`create_table` the existing table and crashed mid-build.
  Wrapped in a `lancedb_table_names()` helper that materializes to
  a list once per stage invocation. `corpus run` against any
  pre-existing build no longer aborts at embed; implicit resume
  works again across the LanceDB version bump. Regression check
  inline in T1/T2.

- **Bundle distillation absolute-path audit on per-category
  lexicon JSONs** ([#70](https://github.com/caseywdunn/corpus/issues/70),
  commit `d925cbe`). The §10 audit (no absolute paths in any
  served artifact) previously only scrubbed `summary.json` /
  `figures.json` / `taxa.json`. Multi-category lexicons emit one
  `<hash>/<category>.json` per top-level lexicon category
  (anatomy, biogeography, methods, …), each carrying the absolute
  path of the source `lexicon.yaml` in its `input_fingerprint`.
  Renamed `_scrub_taxa` → `_scrub_input_fingerprint_path` and
  applied to `taxa.json` plus every per-category JSON via the
  existing `_per_paper_lexicon_outputs(hash_dir)` helper.
  Regression test in `test_package_for_serve.py`; CI greps the
  `Path scrub: rewrote N files; audit clean.` log line.

- **docling crashes on HPC nodes without `g++` in PATH**
  ([#72](https://github.com/caseywdunn/corpus/issues/72), commit
  `d925cbe`). `torch._inductor` JIT-compiles certain ops via the
  system `g++` at runtime when triggered by docling's layout /
  table code paths. HPC compute nodes typically don't have `g++`
  on the default PATH unless a GCC module is loaded, so inductor
  bubbles `OSError: [Errno 14] Bad address: 'g++'` up as a
  docling crash and the affected paper is silently lost from the
  build (1 / 11 in the demo on Bouchet; extrapolates badly to
  production-scale corpora). Verified upstream that `module load
  GCC` alone doesn't help — Bouchet's GCC exposes `g++` but not
  `cc1plus`, so inductor still fails. Fixed by setting
  `TORCH_COMPILE_DISABLE=1` in `pipeline/__init__.py` before any
  submodule loads torch / docling. Docling doesn't depend on
  inductor for correctness; the (small) compile-time perf gain
  only matters for hot tensor ops we don't hit in the pipeline.
  Operators who want inductor JIT can override with
  `TORCH_COMPILE_DISABLE=0` in their env.

- **macOS Apple Silicon — duplicate libomp abort at first
  `import torch`** (commit `d925cbe`). pip-installed torch ships
  its own `libomp.dylib`, scikit-learn ships another, conda-forge's
  `llvm-openmp` ships a third. Without an override, the first
  `import torch` after numpy aborts with the OpenMP duplicate-
  library message and the whole CLI dies before any user code
  runs. `pipeline/__init__.py` sets `KMP_DUPLICATE_LIB_OK=TRUE`
  on darwin (matching what `tools/run_mcp_server.sh` already did
  for the serve path since v0.2).

- **`instructions.md` shepherded into `output_dir` + served
  bundle** (commit `a04aff6`). The README and bundled template
  document `instructions.md` as a corpuscle-root file (next to
  `config.yaml`), but every downstream reader
  (`mcpsrv.main`'s `InitializeResult.instructions` default,
  `mcpsrv.bundle`'s served-bundle whitelist) looks under
  `<output_dir>/`. Without forwarding, the user's corpus-specific
  nudges never reached the served bundle. `corpus run` now copies
  the corpuscle-root `instructions.md` into `output_dir` at the
  end of the pipeline, alongside the bundle build. Same commit
  also fixes a separate "stage count off-by-one" in
  `corpus check`'s ok-line count.

- **`corpus check` silently dropped its own status lines + bib
  staged-overrides NameError** (commit `3393a48`). The validation
  pass produced no output until the final summary, so operators
  couldn't see which checks were passing vs. being silently
  skipped — and on a non-zero exit, the failure message was the
  only signal they got. Status lines now stream as each check
  runs. Same commit: fixed a `NameError` on `db_path` in
  `bib.authority`'s staged-bib path that surfaced as `Could not
  apply staged bib overrides: name "db_path" is not defined` on
  every fresh build (the staged-bib code path had never been
  exercised — caught by the new pyflakes gate).

- **`mcpsrv/tools/chunks.py` missing imports** (commit `dabde66`).
  `translate_chunk` referenced `json` (used to format a tool
  response) and `EmbeddingError` (caught in the embedding-error
  handler of `get_chunks_for_topic`) without importing either.
  The test suite never exercised the affected branches; pyflakes
  caught both during gate adoption.

- **Grobid JDK cgroup-v2 NPE on modern Ubuntu hosts** (commit
  `40ae330`). `lfoppiano/grobid:0.8.1`'s bundled JVM hits a known
  NullPointerException in `CgroupV2Subsystem.getInstance` on
  cgroup-v2 hosts when container-aware sizing is on. Ubuntu 24.04
  GHA runners default to cgroup v2. The EC2 smoke script and
  `integration.yml`'s Grobid service set
  `JAVA_OPTS='-XX:-UseContainerSupport -Xmx4g -Xms1g'` to skip the
  broken codepath; `-Xmx` / `-Xms` are set explicitly so we don't
  need the auto-detection that triggers the NPE.

- **`tools/smoke_test_sse.py` repaired for v0.3 CLI** (commit
  `aea1456`). The SSE round-trip tool referenced the long-removed
  `mcp_server.py` v0.2 entry point. T1 in the new CI matrix
  exercises this every push.

- **EC2 smoke: `pngquant` + `jbig2enc` best-effort** (commit
  `b29aae7`). `ocrmypdf`'s optional native helpers aren't on every
  Ubuntu AMI's default `universe` channel, but the runtime
  degrades gracefully without them (`pipeline.scan` drops
  `--optimize 2 → 1` when `pngquant` isn't on PATH). The smoke
  script now warns and proceeds rather than aborting the whole
  run if `apt install` fails for the helpers.

### Carried from 0.3.1

- See the `## [0.3.1] - 2026-05-11` section below for the
  `get_figure_url` fix (cherry-picked into dev so the v0.4 cycle
  inherits it).

## [0.3.1] - 2026-05-11

### Fixed

- **MCP figure-byte regression from v0.1 — figures now embeddable
  into operator-generated PDFs again**
  ([#69](https://github.com/caseywdunn/corpus/issues/69)).
  `get_figure_image` (introduced in v0.2.0) returns bytes through
  the MCP SDK's `Image()` content channel, which clients render
  inline for the human reader but do *not* expose to the model as
  raw or base64 bytes the model can re-emit through `Write`/`Bash`.
  v0.2 also scrubbed `get_figure`'s `image_path` from absolute →
  relative-to-`documents/`, removing the v0.1 fallback (local stdio
  clients used to `Read` the absolute path directly). Net effect:
  no operator-side path to materialize a figure PNG on disk for
  pandoc / LaTeX / PDF assembly.

  Adds **`get_figure_url`**, a sibling tool that returns an HTTP URL
  the caller fetches via `Bash: curl -H "$auth_header" -o
  /tmp/fig.png "$url"`. Bytes flow over HTTP outside the MCP
  JSON-RPC channel, so they don't burn context tokens regardless of
  figure size (~50 tokens per response vs. ~700K for a 2 MB
  base64-encoded scan). Same shape on local stdio and SSE/AWS
  deploys.

  - **SSE mode**: route mounts alongside `/sse` on the same
    uvicorn endpoint, gated by the existing bearer-token middleware.
  - **stdio mode**: a daemon-thread uvicorn side-car binds an
    ephemeral 127.0.0.1 port at startup and mints a one-shot
    bearer token. The `get_figure_url` tool returns
    `{url, auth_header, mime_type, publishable, license,
    license_source, fetch_hint}`.
  - **#51 publishable gate** is honored identically to
    `get_figure_image` — refuses URLs for unpublishable figures
    unless the server is started with `--allow-unpublishable`.
  - `get_figure` and `get_figure_image` are unchanged; the new tool
    is purely additive.

## [0.3.0] - 2026-05-11

### Breaking changes

- **The 13 root-level Python scripts are gone.** `process_corpus.py`,
  `update_corpus.py`, `embed_chunks.py`, `mcp_server.py`,
  `corpus_status.py`, `build_biblio_authority.py`,
  `build_taxon_mentions.py`, `backfill_intext_citations.py`,
  `reconcile_corpus_to_biblio.py`, `ingest_taxonomy.py`,
  `package_for_serve.py`, `bib_export.py`, and `bib_import.py` are
  deleted in v0.3 ([#60](https://github.com/caseywdunn/corpus/issues/60)).
  Their logic now lives as private modules under `pipeline/`, `bib/`,
  or `mcpsrv/` (e.g. `pipeline.embed`, `pipeline.status`,
  `pipeline.taxonomy_ingest`, `bib.authority`, `bib.reconcile`,
  `mcpsrv.bundle`). Operator-facing entry point is the unified
  `corpus` binary (`corpus run`, `corpus status`, `corpus serve`,
  `corpus bib export|import`, `corpus init`).
- **The `--resume` flag is gone.** `corpus run` is always idempotent;
  re-runs only do work whose inputs have changed (per-stage state +
  fingerprints from v0.2). Use `--force-rebuild` for the rare
  clean-rebuild case
  ([#60](https://github.com/caseywdunn/corpus/issues/60)).
- **The repo-root `config.yaml` is gone.** Per-corpuscle
  `config.yaml` (scaffolded by `corpus init`) is now the single
  source of truth for *this* corpuscle's inputs + system tuning;
  built-in defaults backstop missing keys
  ([#59](https://github.com/caseywdunn/corpus/issues/59)).
- **Operators upgrading from v0.2 must rebuild their corpuscles
  from scratch.** No migration tool, no thin shims. Recipe:
  `pip install -e .` (or `pip install git+...@v0.3.0`),
  `cd <corpuscle>`, `corpus init && $EDITOR config.yaml`,
  `corpus run`.

### Added

- **Unified `corpus` CLI**
  ([#60](https://github.com/caseywdunn/corpus/issues/60)) — one
  binary, cargo/git/gh/kubectl-style verbs (`run`, `status`,
  `serve`, `init`, `bib export|import`, plus stubs for `check` (#62)
  and `completion` (#61)). Global `--config` / `-c PATH` (env:
  `CORPUS_CONFIG`) is git-style pre-verb. The new `corpus run`
  reads the per-corpuscle config, validates against the
  pydantic schema (#59), resolves relative paths against the
  config file's parent (#61), and dispatches the existing
  orchestrator with implicit resume.
- **Per-corpuscle pydantic config schema + bundled template + `corpus init`**
  ([#59](https://github.com/caseywdunn/corpus/issues/59)) — see
  earlier entry; full `config.yaml` surface (input_pdfs,
  output_dir, bib, lexicon, taxonomy, vision, grobid, bibliography,
  licensing) plus carried-over OCR / chunking / quality-gate blocks.
  Field-level validation errors point at the exact key + value.
- **Shared rich console layer**
  ([#63](https://github.com/caseywdunn/corpus/issues/63)) — single
  `pipeline/console.py` Console; emoji status symbols + braille
  spinners + bars on TTY, clean ASCII fallback on SLURM `.out` /
  CI logs / journalctl. Diagnostic logging stays plain text.
- **`pyproject.toml` + `corpus` console_scripts entry point**
  ([#58](https://github.com/caseywdunn/corpus/issues/58)) — project
  becomes pip-installable; `corpus` lands on PATH via
  `[project.scripts]` after `pip install -e .` (development) or
  `pip install git+https://github.com/caseywdunn/corpus.git@vX.Y.Z`
  (deploys). Package version resolves from
  `pipeline.version.__version__` via
  `[tool.setuptools.dynamic]`, so `pip show corpus`,
  `corpus --version`, and the bundle manifest never drift. The
  unified subcommand surface (run, check, status, serve, init, bib,
  completion) lands across the rest of v0.3
  ([#60](https://github.com/caseywdunn/corpus/issues/60) and
  siblings); this entry covers packaging only.

### Changed

- [INSTALL.md](INSTALL.md) and [README.md](README.md) install
  snippets now use `pip install -e .` after `conda env create`.
  [DEPLOY.md](DEPLOY.md) §"On-host setup" likewise
  swaps `pip install -r requirements.txt` for `pip install -e <repo>`.
  `requirements.txt` is retained for the AWS deploy parity per
  [CONTRIBUTING.md](CONTRIBUTING.md) §"Dependencies — two files, on
  purpose"; both manifests must stay in sync.

#### UX polish from the clean-install walkthrough

A second pass over operator-facing output, surfaced by running
[dev_docs/clean_install_walkthrough.sh](dev_docs/clean_install_walkthrough.sh)
on a fresh EC2 box from zero PDFs to a serving MCP endpoint. Each
change below traces back to a confusing moment in that exercise.

- **Per-paper roadsigns in `corpus run` output.** Banner block at
  the top of each paper (`[N/M] Filename.pdf (shorthash)`); every
  per-stage line inside `pipeline.runner` is prefixed with the same
  `[<stem>]` tag via a `LoggerAdapter`. Docling's per-document
  banners (`Processing document processed.pdf` /
  `Finished converting document processed.pdf in N sec`) get the
  same `[Filename.pdf (shorthash)]` annotation via a scoped
  `logging.Filter`, since every PDF is staged under the literal
  `processed.pdf` name before docling sees it.
- **Unified `Filename.pdf (shorthash)` convention.** Every operator-
  facing site that pairs a paper's filename with its hash now uses
  the same single-space `(...)` form: `corpus run` banners,
  `corpus status --sort-by`, `corpus status --propose-skips`, the
  dry-run bucket lists. The `hash X` / `[for X, hash Y]` variants
  are gone.
- **`corpus run --dry-run` names which papers fall in each bucket.**
  Three buckets — `would full-process`, `would partial-process`,
  `would skip` — each lists its members (capped at 20 per bucket
  with `... and N more`). Operators on a 200-paper refresh can
  now confirm the delta before committing CPU.
- **`corpus run --dry-run` ends with a dry-run-aware success line.**
  No more "run complete. Try `corpus serve` next." after a plan-
  only invocation that wrote nothing.
- **`corpus status --report` legibility.** Stage-completion bars now
  carry a `0% / 50% / 100%` scale axis above them, so a wall of
  100% bars actually reads as 100%. Quality flags get a per-gate
  explainer (one paragraph of operator-facing prose pulled from a
  module-level `_GATE_INFO` dict) plus severity and a follow-up
  command (`corpus status --filter-gate <name>`); the old output
  was just `<gate_id>  <count>`. Affected papers in `--sort-by` /
  `--propose-skips` outputs are surfaced by filename, not bare hash.
- **`corpus serve --check` validates both bundles.** Refactored
  into `_check_one_bundle(label, path, ...)` and dispatched twice
  when both `<output_dir>/` (build) and `<output_dir>/_serve/`
  (distilled) exist. Each line is prefixed `[build bundle]` /
  `[distilled bundle]` so the operator can attribute a failure
  to the right tree. Missing `_serve/` is now an explicit warn line
  instead of a silent gap. The `bundle_manifest.json` check now
  varies by bundle type (required-for-distilled, expected-absent-
  for-build) instead of the previous one-size-fits-all warning that
  fired on every healthy local build.
- **`corpus --help` and per-verb help.** Each subcommand parser now
  has an operator-facing `description=` paragraph; the passthrough
  verbs (`status`, `serve`) also have `epilog=` blocks that list the
  forwarded flags (`--transport`, `--port`, `--auth-token-file`,
  `--report`, `--json`, `--sort-by`, …) — argparse can't introspect
  the downstream module, so they used to be invisible from `--help`.
- **`corpus -q` quiets the whole subprocess tree.** Previously,
  `-q` set the parent CLI's log level to WARNING but the actual
  chatter came from `pipeline.orchestrator` and its sub-
  subprocesses, each of which configured its own logger at INFO.
  Verbosity now propagates via a `CORPUS_LOG_LEVEL` env var that
  every package `__init__.py` honors at import time
  (`logging.basicConfig` is a no-op after the first call, so
  configuring at package-import-time wins before any module's own
  setup runs — no per-module refactor needed). `print_status` also
  gates `ok`/`info` lines on the root logger level, so quiet mode
  silences the CLI's own progress sigils too. Result on a 12-paper
  exercise: `corpus -q run --dry-run` now produces 1 log line
  instead of ~30.
- **`taxonomy_ingest` no longer looks hung.** The WoRMS walker
  now logs every 25 records with a running rate (`worms: 250
  records walked (2.3/s)`), and the opener flags the rate-limit
  + expected duration (`large subtrees can take 10+ minutes`)
  instead of emitting a single "Walking WoRMS from AphiaID X ..."
  line and blocking for minutes with no further output.
- **`corpus status --report` quality-flag follow-up.** Affected
  papers can be listed with the printed-inline
  `corpus status --filter-gate <name>` command.
- **Log format unified.** Three modules
  (`pipeline/log.py::setup_root_logging`, `pipeline/embed.py`,
  `pipeline/status.py`) had drifted from
  `%(asctime)s %(levelname)s %(name)s: %(message)s` to a
  timestampless variant, so subprocess output lost its timestamps
  partway through a single `corpus run`. All sites now consistent.
- **`print_status` escapes rich markup in caller messages.**
  Labels like `[build bundle]` were silently dropped on a TTY
  because `rich` interpreted them as style tags; now escaped via
  `rich.markup.escape()`. The off-TTY (plain `print()`) branch
  was already safe.
- **`config.template.yaml` taxonomy block ships commented out.**
  A new corpuscle from `corpus init` no longer hard-codes the
  Siphonophorae example as the default `taxonomy.source` /
  `root_id`; the operator uncomments and picks. The README's
  taxonomy section was expanded with a per-source comparison
  table (`worms` / `dwca` / `dwc`) and the worked example was
  demoted to a code-comment example.
- **Docker called out as a prerequisite in the README.**
  Linked the official install docs + the `apt install docker.io`
  one-liner; the prior copy went straight into
  `docker compose up -d grobid` with no setup advice.
- **`corpus bib export/import` and `corpus serve` no longer print
  a `runpy` warning.** `bib/__init__.py` and `mcpsrv/__init__.py`
  used to eagerly import their `__main__`-runnable submodules,
  which trips a runpy warning when those modules are then invoked
  as `python -m bib.export` / `python -m mcpsrv.main`. Switched
  the at-risk imports to lazy module-level `__getattr__`.
- **Filenames alongside hashes throughout.** A new
  `filename_by_hash` map in the `corpus status` rollup, plus
  `_label_for_paper()` helper, surfaces the original filename
  next to the short hash in every human-readable rendering.
  `--list-hashes` stays bare (it's the xargs-friendly programmatic
  interface).

### Fixed

- **WoRMS AphiaID 1267 was Cnidaria, not Siphonophorae.** The
  README, demo config, and config template all asserted that
  `root_id: 1267` was the order Siphonophorae. It is the phylum
  Cnidaria. Anyone running the demo or copying the example was
  kicking off a BFS walk of all of Cnidaria (~25k+ taxa, ~8–10
  hours at the 0.3s WoRMS rate-limit) instead of Siphonophorae
  (~700–1000 taxa, ~10–15 min). Real AphiaID is `1371`; verified
  against the WoRMS REST API.
- **lancedb `Connection.table_names()` deprecation.** Swept six
  call sites (`pipeline/embed.py`, `pipeline/io.py`,
  `mcpsrv/bundle.py`, `mcpsrv/indexes.py`) to `list_tables()`.
  Verified with `-W error::DeprecationWarning` on the embed /
  bundle / indexes / io test suites.

### Added

- **[dev_docs/clean_install_walkthrough.sh](dev_docs/clean_install_walkthrough.sh)** —
  copy-paste UX walkthrough: fresh conda env → editable install →
  tessdata + grobid → `corpus init` → run → status → SSE serve
  smoke-test with bearer auth → add more PDFs → re-run idempotently
  → tour `--help` / `--cite` / `--version` / every `--dry-run`
  variant / `bib export-import` / `completion` / the underlying
  `python -m` entries. Reference for re-running after CLI-
  affecting changes.

## [0.2.0] - 2026-05-08

A hardening + iteration release. v0.2 closes out v0.1's deferred items
(vision pass at corpus scale, expanded lexicon + taxonomy support,
Bouchet batch-script fixes), adds a real update lifecycle (per-stage
resume, input fingerprints, `update_corpus.py` orchestrator), and lays
in a robustness + observability layer (structured `stage_failures[]`,
silent-failure quality gates, `corpus_status` rollup, network circuit
breakers). Three substantive bugs surfaced during release validation
and ship in this version: a SLURM-array slicing race that corrupted
per-doc state files
([#55](https://github.com/caseywdunn/corpus/issues/55)), an
`input_fingerprint` gap that silently skipped lexicon refreshes on
`--resume` ([#56](https://github.com/caseywdunn/corpus/issues/56)),
and a broken `tesseract-data-*` block in `environment.yaml` that
broke fresh installs
([#52](https://github.com/caseywdunn/corpus/issues/52)). The deploy
stack also moved from single-EC2 + Let's Encrypt to a shared ALB +
ACM pattern that supports multiple organism corpuscles on one
hostname-routed load balancer. Geographic extraction
([#13](https://github.com/caseywdunn/corpus/issues/13)) and trait
extraction ([#14](https://github.com/caseywdunn/corpus/issues/14))
are pushed to later releases.

### Added

#### Update lifecycle

- **`update_corpus.py` orchestrator**
  ([#32](https://github.com/caseywdunn/corpus/issues/32)) — one command
  runs the pipeline + every post-pipeline script in dependency order
  with `--resume`. Forwards the full Stage 1 flag surface so callers
  don't lose pipeline knobs. Makes "add papers and update everything"
  a one-liner.
- **Granular per-stage resume**
  ([#28](https://github.com/caseywdunn/corpus/issues/28)) — replaces
  the all-or-nothing `summary.json` completion marker with per-stage
  status tracked in `pipeline_state.json`. Stages whose artifact is
  already present and whose inputs haven't drifted are skipped on
  `--resume`.
- **Input fingerprints on annotation artifacts**
  ([#29](https://github.com/caseywdunn/corpus/issues/29)) —
  `lexicon_sha256` and `taxonomy_snapshot_date` are stamped into
  `chunks.json` / `taxa.json` / `anatomy.json` so staleness is
  detectable without re-reading the source files.
- **`--re-annotate-stale` flag**
  ([#33](https://github.com/caseywdunn/corpus/issues/33)) — re-runs
  only the lexicon categories whose fingerprint changed. Subsumed
  by per-stage resume + per-category fingerprints.
- **Idempotency contracts on post-pipeline scripts**
  ([#30](https://github.com/caseywdunn/corpus/issues/30)) —
  `build_biblio_authority`, `build_taxon_mentions`,
  `backfill_intext_citations`, and `reconcile_corpus_to_biblio`
  audited + tested for re-run safety.
- **`--audit-orphans`**
  ([#31](https://github.com/caseywdunn/corpus/issues/31)) —
  read-only listing of `documents/<HASH>/` directories (and LanceDB
  rows) whose source PDF is no longer in the input set. Deletion
  stays manual.

#### Robustness + observability

- **Structured `stage_failures[]` + per-stage timing**
  ([#34](https://github.com/caseywdunn/corpus/issues/34)) — reason
  codes (`timeout`, `crash`, `external_unavailable`,
  `unsupported_format`, `corrupted`, `quality_gate`) replace v0.1's
  free-text `errors[]`. Every downstream tool reads from the new
  schema.
- **Silent-failure quality gates**
  ([#36](https://github.com/caseywdunn/corpus/issues/36)) — flag
  empty text layers, gibberish OCR output, all-black figures,
  zero-reference papers, and collapsed-extraction chunks. Surfaces
  what v0.1 was happily writing without complaint.
- **`corpus_status`**
  ([#40](https://github.com/caseywdunn/corpus/issues/40)) — single
  command rolls up stage completion, failure breakdown by reason,
  quality flags, and staleness across the whole build. Dashboard
  for everything else in this section and the update lifecycle.
- **Huge-document hard cap**
  ([#35](https://github.com/caseywdunn/corpus/issues/35)) —
  page-count gate with a structured `too_large` flag for monographs
  that would blow past reasonable wallclock. Haeckel 1888
  *Challenger Siphonophorae* is the canary.
- **External-service flakiness layer**
  ([#37](https://github.com/caseywdunn/corpus/issues/37)) — shared
  retry + backoff + circuit breaker + cache for Grobid, BHL,
  CrossRef, OpenAlex via `pipeline/external.py`. `STRICT_NETWORK`
  env var fails-fast for release builds.
- **Standardized `--dry-run`**
  ([#41](https://github.com/caseywdunn/corpus/issues/41)) — across
  all pipeline + post-pipeline CLIs, replacing the prior 4-of-10
  inconsistent state.

#### Vision + figures

- **SLURM array support for the vision pass**
  ([#27](https://github.com/caseywdunn/corpus/issues/27)) —
  `batch_pass3b.sh` now parallelizes the same way
  `batch_pipeline.sh` does. Prerequisite for running Pass 3b
  (Qwen2.5-VL-7B) at corpus scale
  ([#11](https://github.com/caseywdunn/corpus/issues/11)).
- **Archaic plate prefixes + Roman-numeral support** in figure
  extraction ([#16](https://github.com/caseywdunn/corpus/issues/16))
  — recovers figure numbers from 19th-c. and early-20th-c. caption
  conventions (`Tab. III.`, `Pl. XII.`) the v0.1 heuristics didn't
  match.

#### Indices + features

- **BibTeX round-trip curation**
  ([#26](https://github.com/caseywdunn/corpus/issues/26)) —
  `bib_export` serializes `biblio_authority.sqlite` to BibTeX;
  `bib_import` brings hand edits back into the authority database.
  Closes the loop on user curation of corpus bibliography.
- **Multi-category lexicons**
  ([#24](https://github.com/caseywdunn/corpus/issues/24)) —
  `--lexicon CATEGORY:PATH` accepts any number of YAML files
  (anatomy, biogeography, …); top-level keys define categories.
  Annotations carry a `category` field so per-category resume works.
- **Plant-source taxonomy support**
  ([#23](https://github.com/caseywdunn/corpus/issues/23)) —
  `ingest_taxonomy.py` accepts non-WoRMS Darwin Core archives,
  including `taxa.tsv` Taxon-core layouts
  ([#22](https://github.com/caseywdunn/corpus/issues/22)).
- **Expanded default Tesseract language pack coverage**
  ([#46](https://github.com/caseywdunn/corpus/issues/46)) — the OCR
  fallback union covers more 19th-c. European scripts.

#### MCP server

- **Corpuscle-aware server name**
  ([#17](https://github.com/caseywdunn/corpus/issues/17)) — the
  deployed MCP server identifies itself as `corpus:<corpuscle>`
  rather than the bare `corpus`, with `__version__` from
  `pipeline/version.py` surfaced via `bundle_info`.
- **Inline figure bytes from `get_figure_image`** — returns PNG
  content directly so MCP clients don't need filesystem access to
  the bundle.
- **Two-layer client-side instructions.** The
  `InitializeResult.instructions` blob now joins a packaged default
  ([`mcpsrv/default_instructions.md`](mcpsrv/default_instructions.md))
  with an optional per-corpuscle override prepended on top. Default
  guidance (defer to corpus taxonomy + bibliography, preserve
  historical terminology) ships with the server and reaches every
  client with no operator action; corpus-specific nudges land first
  via `<corpuscle>/instructions.md` when present.
  [`templates/instructions.md`](templates/instructions.md) is the
  starting scaffold for the corpuscle layer.
- **`/healthz` liveness probe** on the SSE transport. Returns
  `200 ok` without requiring the bearer token — mounted ahead of
  the auth middleware in `mcpsrv.transport._HealthzASGI` so uptime
  monitors and reverse-proxy readiness checks don't need the shared
  secret. Exposed through `deploy/nginx.conf` as `location = /healthz`.

#### Internal

- **`pipeline/` package**
  ([#42](https://github.com/caseywdunn/corpus/issues/42),
  [#44](https://github.com/caseywdunn/corpus/issues/44),
  [#45](https://github.com/caseywdunn/corpus/issues/45)) — the
  ~1700-line `process_corpus.py` is split into focused submodules
  (`pipeline.annotate`, `pipeline.chunking`, `pipeline.metadata`,
  `pipeline.figure_passes`, `pipeline.scan`, `pipeline.extract`,
  `pipeline.runner`, `pipeline.main`). `process_corpus.py` is now a
  thin CLI shim. No behavior change.
- **`mcpsrv/` package**
  ([#15](https://github.com/caseywdunn/corpus/issues/15)) — the
  ~2050-line `mcp_server.py` is split per concern (papers /
  taxonomy / bibliography / figures / transports). No behavior
  change.
- **`bib/` package**
  ([#43](https://github.com/caseywdunn/corpus/issues/43)) — bundles
  `bib_metadata`, `bib_export`, and `bib_import` into a single
  namespace.
- **Single-source `__version__`** in `pipeline/version.py`, stamped
  into bundle manifests + MCP `bundle_info`. `main` / `dev` / `vN`
  branching model adopted; `dev` carries a PEP 440 pre-release
  suffix (`0.2.0.dev0` → `0.2.0`).
- **`schema_version` field** on persistent artifacts so future
  schema changes can be detected without parsing.
- **`pngquant` added** to the build environment
  ([#25](https://github.com/caseywdunn/corpus/issues/25)) so
  bundled PNGs ship pre-compressed.

### Changed

- **Repo layout cleanup.** Library modules `figures.py`, `taxa.py`,
  `grobid_client.py`, `embeddings.py`, `vision.py`, `external.py`,
  and `version.py` moved from the repo root into the `pipeline/`
  package — they were never run directly, only imported. Import
  sites: `from figures import …` → `from pipeline.figures import …`
  (and likewise for the other six). One-off utilities
  (`dedup_ghost_works.py`, `unify_doi_corpus_key.py`) moved to
  `tools/`. `PLAN.md` moved to `dev_docs/`. The repo root now holds
  only user-facing CLI entry points.
- **Anatomy-only naming scrubbed throughout.** The lexicon system
  is no longer hard-coded to anatomy; field and file names use
  `lexicon` / `category` instead. Documentation refreshed in lock
  step.
- **`CLAUDE.md` → `AGENTS.md`** — generic agent-guidance filename
  for compatibility with non-Claude assistants.
- **MCP hot paths trimmed for token cost** — `list_papers` row
  shape reduced; previously unbounded MCP list returns are now
  bounded.
- **Per-PDF subprocess gated by platform**
  ([#18](https://github.com/caseywdunn/corpus/issues/18)) — the
  docling subprocess wrapper, originally added to contain a Linux
  segfault, no longer runs on macOS where it was pure overhead.
- **Deploy architecture: ALB + ACM, multi-organism support.** Was
  single-EC2 + nginx + Let's Encrypt; now each EC2 sits behind a
  shared Application Load Balancer (`siphonophores-mcp-alb`) that
  terminates TLS with an ACM wildcard cert. nginx on the instance
  drops to plain `:80` (no certbot, no TLS). Adds a default
  `/health` endpoint nginx serves directly for ALB target-group
  health checks (separate from the `/healthz` MCP probe under
  bearer auth). One listener rule per organism (host-header match
  → target group), so adding `cnidaria.siphonophores.org` etc. is
  a target-group + listener-rule + DNS CNAME, no new TLS infra.
  CloudFormation `BucketName` parameter now optional via a new
  `CreateBucket` flag so re-deploys can attach to an existing
  bundle bucket. Full runbook rewritten: [DEPLOY.md](DEPLOY.md).
- **`INSTALL.md` promoted to repo root.** Moved from `dev_docs/`
  since the content (jbig2enc, OCR language packs, Grobid, pip
  fallback) is user-facing and belongs alongside README,
  CONTRIBUTING, and CHANGELOG. `dev_docs/` keeps its
  maintainer-only docs (PLAN, DEPLOY, BOUCHET, MCP_TOOLS, TESTING,
  OVERVIEW). All cross-references updated.

### Fixed

- **Bouchet batch scripts**
  ([#20](https://github.com/caseywdunn/corpus/issues/20),
  [#21](https://github.com/caseywdunn/corpus/issues/21)) —
  native-library load problems, CUDA / torch wheel selection, and a
  fail-loud preflight so configuration drift surfaces before the
  array launches instead of after each task fails individually.
- **`deploy/stack.yaml`**
  ([#6](https://github.com/caseywdunn/corpus/issues/6)) — removed
  stale default-VPC assumptions that broke deploys on accounts
  without a default VPC.
- **Docling extraction fails loudly** instead of writing placeholder
  text when the parser produces nothing usable.
- **Missing PyMuPDF surfaces as a stage failure** instead of a
  silent skip.
- **`--skip-pipeline` rejected with `--re-annotate-stale`** —
  previously these silently combined into a no-op.
- **`rapidfuzz` required at import** in `unify_doi_corpus_key` and
  `dedup_ghost_works` — failure surfaces at startup, not partway
  through a run.
- **`ingest_to_vector_db` import** wired into `pipeline.main` — was
  silently dropped during the package extraction.
- **`package_for_serve` discovers lexicon outputs dynamically** —
  no more hardcoded list that drifted from the multi-category
  lexicon work; also bundles `instructions.md` and warns on missing
  top-level files.
- **Three issues filed during release validation**
  ([#47](https://github.com/caseywdunn/corpus/issues/47),
  [#48](https://github.com/caseywdunn/corpus/issues/48),
  [#49](https://github.com/caseywdunn/corpus/issues/49)) — path
  handling, cross-paper leakage in a query path, and an incorrect
  rank field. Closed before release.
- **Stage 1 SLURM array tasks could process the same PDF
  concurrently and corrupt per-doc state files**
  ([#55](https://github.com/caseywdunn/corpus/issues/55)). The
  pre-resume filter at the top of `pipeline.main.main()` ran
  *before* batch slicing, so each array task's slice depended on
  disk state at task-start time. Tasks starting later saw fewer
  remaining hashes; their slice indices then landed on different
  members of the list, producing overlapping batches. The two
  writers raced on a shared `pipeline_state.json.tmp` filename,
  leaving either an interleaved-bytes payload (31 corrupted
  summaries in the production run that surfaced this) or a
  `FileNotFoundError` on the rename. Fixed by slicing on the
  unfiltered hash list (`_slice_hashes_for_batch` in
  `pipeline.main`); resume skipping is now done per-doc inside the
  loop. Belt-and-braces: `_save_pipeline_state` and
  `create_summary_json` now use per-writer tmp filenames
  (`.tmp.<pid>.<ns>`) so any future regression that re-introduces
  concurrent writers on the same hash gets last-write-wins
  atomicity instead of corruption. Regression tests in
  `tests/test_batch_slicing.py`.
- **`--resume` outer fast-path ignored `input_fingerprint`, so
  editing a lexicon or taxonomy never invalidated already-recorded
  stages** ([#56](https://github.com/caseywdunn/corpus/issues/56)).
  The per-doc fast-path skip in `pipeline.main` called
  `_all_stage_artifacts_complete` without forwarding the live
  taxonomy/lexicon fingerprints, so a doc whose
  `taxa_and_lexicon_extraction` stage was recorded under (e.g.) a
  taxonomy-only fingerprint was silently skipped on a re-run that
  loaded a real lexicon — even though the inner per-stage gate
  would have correctly flagged it stale, that gate was unreachable.
  Effective result: the only way to force a lexicon refresh was
  bumping `PIPELINE_VERSION`. Production occurrence: 12,669 docs
  in job `11108546` carried a fingerprint without `lexicons`
  because the source YAML had loaded as `{}`; the resume run
  silently passed over all of them. Fixed by adding
  `expected_fingerprints` to `_all_stage_artifacts_complete` and a
  shared `_expected_fingerprints_for_run` helper that both
  `main.py`'s outer gate and `runner.py`'s inner gate now use,
  keeping the two staleness notions in lockstep. Regression tests
  in `test_per_stage_resume.py`.
- **`environment.yaml` references to `tesseract-data-<code>`
  packages** ([#52](https://github.com/caseywdunn/corpus/issues/52),
  closes [#9](https://github.com/caseywdunn/corpus/issues/9)) —
  the 12 per-language entries added in #46 listed package names
  that don't exist on conda-forge, so a fresh `conda env create`
  failed with `PackagesNotFoundError`. Existing dev envs survived
  because incremental updates silently skipped the missing
  packages. Replaced with a new
  [`tools/install_tessdata.sh`](tools/install_tessdata.sh) helper
  that downloads the default fallback set (deu, fra, rus, lat,
  spa, por, chi_sim, chi_tra, jpn, ell, kor, grc, deu_latf)
  directly from `tessdata_best`. Idempotent; takes a custom
  language list as positional args; honors `TESSDATA_DIR` for
  non-conda installs. Subsumes the v0.1 deu_latf-on-Bouchet manual
  install (#9) — `deu_latf` is now part of the default download.

### Removed

- **`bib_metadata.py` shim** — the deprecated re-export shim is
  removed. Migrate any leftover `from bib_metadata import …`
  imports to `from bib import …` (or `from bib.parser import …`
  for private helpers).

### Deferred / known limitations

- **Streamable HTTP transport with OAuth**
  ([#5](https://github.com/caseywdunn/corpus/issues/5)) — deferred
  indefinitely. SSE + bearer-token works for the ~20-collaborator
  deploy target; revisit only if wider distribution materializes
  or MCP clients drop SSE support.
- **Geographic mention layer**
  ([#13](https://github.com/caseywdunn/corpus/issues/13)) — pushed
  to v2.0+; the mention-layer surface is likely to be reworked at
  the major-version boundary.
- **Trait extraction + identification keys**
  ([#14](https://github.com/caseywdunn/corpus/issues/14)) — v0.3
  candidate. Substantial enough to warrant its own plan section
  when picked up.
- **Embedding model migration path**
  ([#38](https://github.com/caseywdunn/corpus/issues/38)) —
  design-only for v0.2; implement when a model swap is actually
  needed (BGE-M3 v2 or a different family).
- **MCP server lazy index loading**
  ([#39](https://github.com/caseywdunn/corpus/issues/39)) —
  premature at 1.8K papers; documented as a known scaling cliff
  for 10K+ corpuses.

## [0.1.0] - 2026-05-01

First tagged release. The pipeline ingests a corpus of scientific-literature
PDFs, extracts text/figures/references, builds a set of structured indices
(taxonomy, anatomy, bibliography, citations, figures, embeddings), and serves
them to MCP clients. End-to-end on ~1,800 siphonophore papers spanning
late-18th-century printed monographs through born-digital 2025 articles.

### Added

#### Pipeline

- **Two-stage pipeline** with hash-addressed per-paper artifacts. Stage 1
  (`process_corpus.py`, CPU-only) does scan classification, OCR, docling
  extraction, Grobid metadata, chunking, and figure detection. Stage 2
  (`embed_chunks.py`, GPU) does BGE-M3 embeddings into LanceDB. Each unique
  PDF is identified by the first 12 hex chars of its SHA-256; all artifacts
  live under `<output_dir>/documents/<HASH>/`. Presence of `summary.json` is
  the completion marker that `--resume` checks.
- **Three-class scan detection** (`born_digital` / `scanned` /
  `broken_text_layer`) using a langdetect + gibberish-score + Tesseract-OSD
  visual-script cross-check. Routes to ocrmypdf with the appropriate mode
  (`--skip-text` vs `--force-ocr`) and language packs.
- **Language-aware OCR** via ocrmypdf with per-document Tesseract pack
  selection. Supports English, German (incl. `deu_latf` Fraktur for 19th-c.
  literature), French, Russian, Latin. Falls back to a default union when no
  language is detectable.
- **Three-pass figure pipeline.**
  - **Pass 1+2** — docling Picture extraction with classification, dedupe,
    semantic filenames (`fig_3_a.png`), and PyMuPDF fallback when docling
    finds nothing.
  - **Pass 2.5** — caption-panel parser (regex over `A.` / `(A)` / `A–C.`
    conventions, multilingual prefixes) + missing-figures scan that catches
    `Fig. N` mentions in body text whose figure docling didn't extract.
  - **Pass 3a** — Tesseract ROI pass with `image_to_data` + caption
    cross-check, opt-in via `--content-aware-figures`.
  - **Pass 3b** — vision-LLM ROI pass via Claude or local Qwen2.5-VL-7B
    (`vision.LocalVLMBackend`). Same result schema as 3a. Selected via
    `--vision-backend {claude,local}`.
  - **Pass 3c** — compound-figure rename. When a single docling extraction
    contains multiple `Fig. N` labels, partition ROIs by parent figure,
    rename to range notation (`fig_3-4.png`), record `previous_filenames[]`,
    and emit a standalone figure record per recovered sub-figure.
- **Figure+caption as a first-class joint object** in `figures.json`:
  caption text + bbox, panel descriptions, ROIs, cross-references back to
  chunks. Per-paper `figures_report.html` provides a thumbnail + caption
  contact sheet for human QC.
- **Per-page QC visualizations** (`visualizations/page_N_visualization.png`)
  overlay text-cell boxes (red) and figure bboxes (yellow) on each page.
- **Per-paper logs** at `documents/<HASH>/pipeline.log` via
  `per_pdf_file_log` context manager. All `logging` calls inside the
  per-paper block are captured to that file in addition to the root handler.

#### Indices

- **WoRMS taxonomic backbone** ingested via `ingest_taxonomy.py` from any
  Darwin Core taxonomy snapshot. Stored as `taxonomy.sqlite` alongside
  per-paper artifacts; consumed by every taxon-keyed query.
- **Bibliographic authority database** (`build_biblio_authority.py` →
  `biblio_authority.sqlite`). Unifies corpus papers, cited references, and
  taxonomic authority strings into a single graph. Three-phase build:
  seed-from-corpus → ingest-cited-references-with-cascading-match → link
  taxonomic authorities. GUID priority: DOI → BHL Part/Item → normalized
  citation key (`corpus:haeckel|1888|report on the siphonophorae col`).
  ~65,000 references across the corpus deduplicated into ~15–25K unique
  works.
- **BHL enrichment** for historical references — fuzzy author+year+title
  match against the Biodiversity Heritage Library API, gated behind
  `--enrich-bhl`. Two-stage cascade with query cache and resume support.
- **Taxon mention database** (`build_taxon_mentions.py` →
  `taxon_mentions.sqlite`). gnfinder-driven detection over chunks, resolved
  against the configured taxonomy snapshot. Supports abbreviated-form
  expansion (`A. elegans` → `Agalma elegans` based on most recent genus
  context). Backs all taxon-keyed MCP tools.
- **In-text citation graph** (`backfill_intext_citations.py` →
  `intext_citations.json` per paper). Parses Grobid TEI body
  `<ref type="bibr">` elements, joins each to a chunk offset, and resolves
  the target attribute to a `work_id` in the authority database.
- **Anatomy-term index** — exact-match + stemming pass over a curated
  YAML lexicon (nectophore, bract, palpon, gonophore, pneumatophore, …).
  Per-paper offsets in `anatomy.json`.
- **Vector index** — BGE-M3 dense embeddings (1024-dim, 8K context,
  multilingual) stored in LanceDB. Replaces the prior OpenAI-backed
  embedding path. Auto device-detect (cuda → mps → cpu) via
  `embeddings.detect_device()`.

#### MCP server

- **26-tool MCP server** (`mcp_server.py`) over the indices, built on
  FastMCP with both stdio and SSE-over-HTTP transports. Eager in-memory
  index at startup. Tool surface covers papers, taxonomy, anatomy,
  bibliography, citations, and figures. Full list in
  [`dev_docs/MCP_TOOLS.md`](dev_docs/MCP_TOOLS.md).
- **SSE transport with bearer-token auth** for remote-client connections
  (Claude Desktop Custom Connectors, Claude Code, claude.ai). Single shared
  token in SSM Parameter Store / `/etc/corpus/token`; per-user tokens
  deferred until revocation becomes a real need.
- **`bundle_info` tool** reports manifest fields (`bundle_version`,
  `pipeline_git_sha`, `embedding_model`, `taxonomy_snapshot_date`, paper /
  figure / chunk counts). Lets clients detect stale endpoints.
- **`<corpuscle>/instructions.md` served as `InitializeResult.instructions`** — a
  per-corpus prompt/orientation file that ships to MCP clients on connect.
- **Corpuscle pattern** — per-instance state (lexicon, taxonomy snapshot,
  authority overrides) lives in one user-controlled directory; the pipeline
  is multi-corpus by configuration.

#### Operational

- **SLURM job-array parallelization** for Stage 1 via
  `slurm/batch_pipeline.sh`. `--batch-index` / `--batch-size` deterministic
  slicing of the sorted hash list. Companion `batch_process_corpus.sh` (day
  partition), `batch_pass3b.sh` (gpu_h200, Qwen2.5-VL), `batch_embed.sh`
  (gpu, BGE-M3), `batch_grobid.sh`, `batch_biblio.sh`. All resumable; can be
  chained via `--dependency=afterok`.
- **Subprocess isolation for segfaults** — Stage 1 wraps the docling parse
  in a subprocess so a segfault in one PDF doesn't kill the batch.
- **Test suite** (`tests/`) — ground-truth tests for 3 curated papers and
  corpus-wide structural/consistency checks across all 1,787 papers. See
  [`dev_docs/TESTING.md`](dev_docs/TESTING.md).
- **AWS deployment runbook** — `deploy/stack.yaml` (CloudFormation),
  `deploy/nginx.conf`, systemd unit, `update.sh`, and a CLI-only runbook in
  [`DEPLOY.md`](DEPLOY.md). Two-bundle model: build bundle
  (Bouchet, ~10 GB) vs. served bundle (S3, ~3 GB). `package_for_serve.py`
  walks `documents/<HASH>/`, copies whitelisted files only, and writes a
  versioned `bundle_manifest.json`.

### Changed

- **Switched embeddings from OpenAI → local BGE-M3.** Eliminates network /
  rate-limit / silent zero-vector failure modes. Vector dim 1536 → 1024;
  any pre-existing LanceDB index must be rebuilt. OpenAI backend removed
  entirely from `embed_chunks.py` and `requirements.txt`.
- **Replaced naive 1000-char chunker with docling's `HybridChunker`.**
  Tokenizer-aware, respects section/heading structure.
- **Replaced `extract_metadata` stub with Grobid client.** Grobid runs as a
  Docker / Singularity service; `process_corpus.py` calls
  `/api/processFulltextDocument` and persists the full TEI-XML alongside
  `metadata.json` for re-parsing without reinvoking Grobid.
- **8-char hash prefix → 12-char.** Prevents collisions across thousands of
  PDFs, verified at directory creation.
- **`config.yaml` is now actually loaded.** Previously a stub; `load_config`
  in `process_corpus.py` parses it, dead config sections were pruned.
- **Embedding failures now raise**, not silently insert zero vectors. Same
  contract on the local backend as the prior OpenAI path was supposed to
  provide.
- **Repo layout reorganized** — top-level CLIs at root, batch scripts in
  `slurm/`, helper utilities in `tools/`. Snakemake legacy pipeline deleted
  (process_corpus.py is the only entry point).
- **Demo corpus walkthrough** — `README.md` rewritten for biological users;
  `demo/` ships a 5-paper example corpus + `instructions.md` noting that
  Velella and Porpita are not siphonophores.

### Fixed

- **Visual-script-mismatch false positives** in scan detection. Tightened
  the gibberish-score threshold from 0.25 → 0.40, recovering ~23 papers
  that were being incorrectly forced through OCR despite having clean Latin
  text layers (e.g., `Rossietal2008.pdf`).
- **Filename-year fallback** when Grobid emits an empty `<date/>`. Falls
  back to a 4-digit year extracted from the filename rather than dropping
  the year entirely.
- **Figure dedup** — exclude furniture from groups; order panels by
  position so coequal-panel figures get stable A/B/C labels.
- **Subprocess crash handling** — process group isolation + memory bumps
  so a docling segfault on one PDF doesn't take out the rest of the batch.
- **`sync_to_s3` empty-array expansion fix** for macOS bash 3.2.
- **Ubuntu 24.04 user-data fix** — no awscli apt package, skip the upgrade
  step that fails on first boot.
- **EBS resize wait loop** — modern AWS reports `optimizing` not
  `optimized`; the runbook polled the wrong state.

### Removed

- **Snakefile + Snakemake-only scripts.** `process_corpus.py` is the sole
  pipeline entry point.
- **OpenAI embedding backend.** Replaced by local BGE-M3 (see Changed).
- **5 orphan scripts** with zero external references (commit `80a1a8d`).
- **Hard-coded taxonomic-group specificity** in filenames and config —
  the pipeline is now corpus-agnostic; siphonophore-specific data lives in
  the corpuscle.

### Deferred / known limitations

- **Figure-number extraction on historical scans** — ~538 of 1,787 papers
  (~30%) have `Fig. N` references in body text but no extracted figure
  numbers. Mostly 19th-c. and early-20th-c. scans whose caption formatting
  doesn't match the heuristics that drive `parse_figure_number`. Figures
  are still extracted; only the numbering is missing. Tracked in
  [#16](https://github.com/caseywdunn/corpus/issues/16) for v0.1.x.
- **MCP server identity is not corpuscle-aware** — the deployed server
  identifies itself as `corpus` rather than after the corpuscle it serves,
  and is not version-stamped at the server level (the bundle manifest is
  versioned, but the server name isn't). Affects multi-corpus deployments.
  Tracked in [#17](https://github.com/caseywdunn/corpus/issues/17) for v0.1.x.
- **Fraktur Tesseract pack on Bouchet** — 19th-c. German scans (Goldfuss
  1820, Pagenstecher 1869, Brandt 1837, Donitz 1871, etc.) currently OCR
  to whitespace because the `deu_latf` pack isn't installed in the build
  environment. Tracked in [#9](https://github.com/caseywdunn/corpus/issues/9)
  for v0.1.x; install + reprocess the affected papers to recover them.
- **Vision-pass coverage at corpus scale** — `batch_pass3b.sh` (Qwen2.5-VL)
  is wired and tested but did not run as part of the v0.1 rebuild. Pass 3c
  (compound-figure split) and the bulk of the 6,841 pre-existing
  `missing_figures[]` records are unresolved. Tracked in
  [#11](https://github.com/caseywdunn/corpus/issues/11) for v0.1.x.
- **Geographic extraction** (§12 Layer 3 in dev_docs/PLAN.md) — not yet implemented.
  Tracked in [#13](https://github.com/caseywdunn/corpus/issues/13).
