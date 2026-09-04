# PLAN.md — Corpus pipeline (v1.3)

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
against something other than itself.

**v1.2 (2026-08-29) is the one v1.3 grows directly out of.** Its theme
paragraph is in the CHANGELOG; what matters here is the inheritance. The
gold set is now a standing instrument rather than a one-cycle project — any
extraction change can be scored against it before and after, which is the
condition v1.3's work inherits and should keep using. Its unfinished
fingerprint-based regression reference (#187) is carried below. The
library-inspection work (#217) moves with the skills-and-usage cycle to v1.4.

v1.2 also left a standing lesson worth stating once, because it recurred
three times in that cycle and once more at its own release: **a measurement
is not a result until you have looked at what it is measuring.** The
plate-legend fix shipped on a recall gain and quietly cost precision;
#254's OCR page-blanking was invisible to every quality gate the pipeline
had, because each one measured the document rather than the machine that
built it; and a stray file sat in the repo root through four green CI runs
and two release PRs because nothing looks there.

**v1.3 is the evidence-integrity and auditability cycle.** A preliminary
acceptance pass over the newest reference corpuscle made wrong figure-caption
associations the most conspicuous errors in actual use. The same build exposed
whole-document OCR loss, append-only embeddings and order-dependent reference
reconciliation. Those are one product problem: corpus can return plausible
evidence without being able to show that it is the evidence printed on the
page, or that a re-run would return the same current corpus.

The v1.3 outcome is therefore: **corpus either preserves and links source
evidence correctly, or exposes its uncertainty; a re-run cannot leave stale
evidence; and an operator can inspect the result against the original page.**
The skills-and-usage work previously assigned to v1.3 moves to v1.4 with a
narrower scope. It should consume trustworthy captions, references and update
semantics rather than design around their current failures.

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
targets the **gold corpuscle** (built in v1.2) — big enough for a rate to mean
something, small enough to rebuild per release. That comparison was run
by hand against a viburnum rebuild during v1.1 and is what proved the
3.4 GB → 2.3 GB drop was #184's re-encoding rather than lost content.

The `--deselect` question that paragraph raised is answered: **none**.
[#167](https://github.com/caseywdunn/corpus/issues/167) removed all three
flags, T1 now runs bare `-m corpus_required`, and T2 additionally ignores
`test_reference_extraction.py` because Grobid is disabled there. Both
workflow files state the reasoning inline.

## v1.3 — evidence integrity and auditability

This cycle has one admission rule: work belongs here only when it changes
whether evidence is preserved, associated with the right source, reproducible
after an update, or inspectable when it is uncertain. A new workflow, export
surface or convenience feature does not qualify. The four workstreams below
are ordered: the page-level view is an instrument for the caption work; the
caption and reference decisions are build outputs; update correctness makes
those outputs reproducible; the execution-plane contract keeps them out of
the running server.

### 1. Source-page and figure evidence — the primary fidelity gate

The newest reference-corpuscle acceptance pass makes this the first-order
user problem. A missing caption is visible; a confidently returned caption
from the wrong figure is plausible misinformation. v1.2 measured caption
binding at 0.574 recall / 0.878 precision, with the weakest layouts in the
oldest material; [FIGURE_PARSING.md](FIGURE_PARSING.md) is the stable account
of what those numbers mean.

- [ ] **Repair and re-measure caption association**
  ([#195](https://github.com/caseywdunn/corpus/issues/195)). Bind on validated
  figure-number evidence where available, distinguish a prose caption from a
  bare label, and record the chosen candidate, rejected candidates, source,
  page distance and confidence. A weak association must be exposed as weak in
  the stored artifact and MCP result, not presented as an ordinary caption.
  Full regression replay now moves the default served surface from 0.551 /
  0.890 to 0.580 / 0.910 over the same persisted Docling inputs; keep this open
  until the remaining panel-split acceptance check and clean rebuild are done.
- [x] **Resolve the grouped-plate failure represented by `Vanhoeffen1906`**
  ([#203](https://github.com/caseywdunn/corpus/issues/203)). Bind enumerating
  caption blocks, parse lists of figure numbers separately from lettered
  panels, admit those plates to ROI splitting, and cross-check the proposed
  split against `missing_figures`. Scope is the measured layout family, not a
  general segmentation rewrite.
- [x] **Make vision output completeness a hard fact**
  ([#269](https://github.com/caseywdunn/corpus/issues/269)). A response clipped
  by the token budget is incomplete, never a successful panel result; size the
  budget from emitted structure or fail/retry explicitly.
- [x] **Build the optional, re-runnable per-document page report**
  ([#274](https://github.com/caseywdunn/corpus/issues/274)): original processed
  page, cell/figure/ROI overlay, selectable parsed text, and page-level
  statistics/provenance in one HTML document. It is build/operator-side QC,
  reads `docling_doc.json`, and never ships in the served bundle. Decide in
  the same change whether the existing unconditional visualization PNGs become
  opt-in or on-demand; do not replace them with another corpus-wide raster set.
- [x] **Close the completed figure-measurement issue**
  ([#194](https://github.com/caseywdunn/corpus/issues/194)) after confirming
  the scorer and recorded result still match the current gold set. Remaining
  caption behavior belongs to #195/#203 rather than another measurement task.
- [x] **Close the whole-document silent-loss paths**:
  `no_text_layer` / `vendor_boilerplate_only` preserving the wrong layer
  ([#264](https://github.com/caseywdunn/corpus/issues/264)); visual-script
  comparison hidden behind the gibberish threshold
  ([#266](https://github.com/caseywdunn/corpus/issues/266), superseding
  [#172](https://github.com/caseywdunn/corpus/issues/172)'s broader symptom);
  image placeholders
  satisfying `empty_text`
  ([#267](https://github.com/caseywdunn/corpus/issues/267)); and OCRmyPDF
  exiting zero with no recovered text
  ([#268](https://github.com/caseywdunn/corpus/issues/268)). Provide
  [#186](https://github.com/caseywdunn/corpus/issues/186)'s per-document OCR
  override as the explicit escape hatch for a classifier that still gets a
  document wrong.

**Acceptance:** each conspicuous failure from the newest reference corpuscle
has a small regression fixture; caption recall and precision are reported
before and after by era and layout; precision does not trade away silently for
recall; and every document that recovers no meaningful text fails an error-level
quality gate regardless of process exit code or markup placeholders.

### 2. Reference evidence — deterministic and non-destructive

The full proposal in #240 mixes a necessary data-model repair with speculative
model-assisted blocking. v1.3 takes the deterministic foundation and leaves
the latter out until the new observation set can measure whether it is needed.

- [x] **Repair the inputs before clustering:** stop recording journal names
  as work titles ([#226](https://github.com/caseywdunn/corpus/issues/226)) and
  stop a malformed DOI from short-circuiting stronger title/author evidence
  ([#239](https://github.com/caseywdunn/corpus/issues/239)).
- [x] **Separate reference observations from canonical works**
  ([#240](https://github.com/caseywdunn/corpus/issues/240), deterministic
  core only). Keep one immutable/re-derivable observation per citing-paper
  reference; map it to a canonical work in a separate relation carrying the
  method, score and producer version. Reconciliation becomes a pure,
  re-runnable decision over observations rather than a destructive merge.
- [x] **Centralize the works/citations schema used by tests**
  ([#237](https://github.com/caseywdunn/corpus/issues/237)) before migrating
  it, so five hand-written test schemas cannot validate five different models.
- [x] **Preserve the frozen MCP wire surface.** Existing reference tools read
  a compatibility view/materialization with their parameter names, defaults
  and response fields unchanged. Better data is expected; an accidental API
  migration inside a data-model change is not.
- [ ] **Re-measure `get_missing_references`**
  ([#155](https://github.com/caseywdunn/corpus/issues/155)) after the migration,
  using the raw observation set to explain every canonicalization decision.

**Deferred from #240:** embedding-based candidate blocking and local-model
adjudication. Neither enters v1.3 without a measurement from the deterministic
system showing which unresolved population it would fix, a versioned model
input, and an auditable persisted verdict.

**Acceptance:** reconciliation from a clean build and reconciliation after an
incremental paper addition produce the same observation-to-work mapping;
running it twice is a no-op; no raw observation is destroyed; and every mapping
can answer which evidence and rule produced it.

### 3. Truthful updates and release evidence

- [ ] **Define and implement the corpuscle update contract**
  ([#265](https://github.com/caseywdunn/corpus/issues/265)): additions,
  removals, changed PDF bytes, same-hash derived-content changes, bib edits,
  config edits and version upgrades each state what invalidates, what is
  replaced and what is pruned.
- [ ] **Extend input fingerprints**
  ([#174](https://github.com/caseywdunn/corpus/issues/174)) to every input that
  can change an artifact, including the relevant bib fields, filename and
  resolved configuration. A key belongs in a stage fingerprint only when that
  stage actually consumes it.
- [ ] **Make embedding replacement atomic per document and delete the Stage 1
  fake completion marker**
  ([#271](https://github.com/caseywdunn/corpus/issues/271)). Re-embedding one
  hash replaces its rows; it never skips stale text or appends a second
  generation. A marker is written only after the table commit and carries the
  chunk fingerprint, model, dimension and committed row count. Bundle metadata
  validates all markers rather than sampling an arbitrary one.
- [ ] **Land the fingerprint-based release reference**
  ([#187](https://github.com/caseywdunn/corpus/issues/187)). Diff pipeline
  output, quality flags and manifest facts for the fixed gold corpuscle; test
  counts remain a CI activity signal, not a data regression reference.

**Acceptance:** for every supported change class, an incremental run and a
clean rebuild have the same current document set, artifact fingerprints,
cross-paper mappings and vector rows. `corpus status` and the bundle manifest
never infer completion from a placeholder or from a sampled marker.

### 4. Execution planes and a thin server

Write the normative model in [OVERVIEW.md](OVERVIEW.md), not in this ephemeral
plan. Use **execution planes** because pipeline "stage" already has a different
meaning:

| Plane | Owns | May compute | Must not do |
| --- | --- | --- | --- |
| Library curation | PDFs, BibTeX, curator directives | inspect, validate, propose/review source edits | write corpus artifacts or duplicate generic PDF/OCR knowledge |
| Build/materialization | OCR, extracted text/figures, associations, databases, embeddings, served bundle | expensive batch/GPU/API work; deterministic resumable transforms | depend on `tools/` or `skills/`; leave ambiguous partial state |
| Serve/query | an immutable bundle plus disposable caches | bounded lookup, filtering, authorization, formatting and compatible query embedding | OCR, reconciliation, corpus-wide mutation, external enrichment or general LLM calls |
| Client/agent | user intent and deliverables | synthesis, translation, workflow orchestration and presentation | serve as the enforcement point for licensing, provenance or access control |

- [x] **Document the ownership and data flow** in OVERVIEW; summarize it in
  README and the contributor invariants in AGENTS.md. Data flows library →
  build → immutable bundle → bounded server response → client output. Client
  feedback becomes an explicit reviewed library edit, not server mutation.
- [ ] **Audit the served path against that contract.** In particular: move
  on-demand figure crops out of the bundle or materialize them at build time;
  make figure-download URLs work behind the reference reverse proxy without
  returning the shared MCP bearer token in a model-visible response; and apply
  one bounded-input policy to every MCP collection/list parameter. File the
  implementation issues before changing the frozen surface.
- [ ] **Strengthen the freeze gate.** Snapshot tool signatures/defaults and
  representative response schemas in addition to the existing tool-name,
  error-shape and licensing checks.

The query embedder is the deliberate exception to "no models in the server":
the query vector must be compatible with the stored index. It stays lazy,
version-checked and bounded to the active request. Thin means logically
read-only and operationally bounded, not computation-free.

### v1.3 release gate

- [ ] Every user-visible caption failure selected from the newest reference
  corpuscle is fixed or returned with explicit uncertainty, and the gold
  caption-binding report has been inspected rather than merely generated.
- [ ] The page report makes the original page, parsed text, chosen caption and
  competing evidence reviewable without manually joining four artifacts.
- [ ] Clean and incremental builds agree for vector rows and reference
  mappings on the fixed regression corpuscle.
- [ ] The served bundle can be mounted read-only; all MCP calls still work,
  including remote whole-figure and panel download through the deployed proxy.
- [ ] T0, T1/T2, T3, T3-bare where platform behavior changed, and the relevant
  T5 fidelity scorers pass under CONTRIBUTING.md's release ritual.

## v1.4 — skills and usage

The v1.3 work makes this cycle possible; it is not part of it. v1.4 is the
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
not part of v1.4. If the clade-monograph acceptance run exposes another
evidence-integrity defect, fix it as a defect; do not expand the skills cycle
into another pipeline redesign.

## Unscheduled

Not claimed by v1.3 or v1.4. This is a selected orientation list, not a second
issue tracker; GitHub issues are authoritative. Dependencies that matter are
stated inline. Split because the two halves get picked up for different
reasons: a known defect is picked up when it bites someone, an unbuilt feature
when something makes it worth building.

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

- [ ] **`pipeline/CITATION.cff` is a hand-maintained copy**
  ([#173](https://github.com/caseywdunn/corpus/issues/173)) with nothing
  enforcing sync — currently byte-identical, so this is a guard rather
  than a fix. `tests/test_citation_cff.py` (v1.1.1) now validates both
  files independently, which is adjacent but not the same check.
- [ ] **README is ambiguous about where `instructions.md` lives**
  ([#171](https://github.com/caseywdunn/corpus/issues/171)).
- [ ] **The local VLM loads in float32 on MPS**
  ([#258](https://github.com/caseywdunn/corpus/issues/258)). Validate the dtype
  fix on Apple Silicon together with any move beyond the current docling pin.
- [ ] **Migrate the lint gate to Ruff**
  ([#259](https://github.com/caseywdunn/corpus/issues/259)). Preserve F821 as
  its own hard gate while making the existing `# noqa` suppressions real, then
  resolve the remaining findings rather than suppressing them wholesale.
- [ ] **Report the BHL enrichment hit rate**
  ([#260](https://github.com/caseywdunn/corpus/issues/260)); the existing
  counters are initialized but never updated.
- [ ] **Close the repaired taxonomy-ingest issue** after verifying the merged
  idempotence and in-place repair tests
  ([#262](https://github.com/caseywdunn/corpus/issues/262)).
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

## Reference: target queries

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`.

| # | Pattern | Status entering v1.3 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer ([#13](https://github.com/caseywdunn/corpus/issues/13), deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; v1.3 repairs caption/reference evidence and makes uncertainty inspectable before v1.4 adds the clade-monograph skill ([#179](https://github.com/caseywdunn/corpus/issues/179)) |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred ([#14](https://github.com/caseywdunn/corpus/issues/14)) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Indices in place; a corpus-scale vision run landed in v1.0 and figure detection became measurable against truth in v1.2. On the current gold corpuscle the raw artifact is 0.936 recall / 0.876 precision, while the stricter default MCP type filter is 0.875 / 0.985. Caption binding remains the conspicuous evidence error and is v1.3's primary fidelity gate |
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
