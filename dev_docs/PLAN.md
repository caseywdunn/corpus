# PLAN.md — Corpus pipeline (v1.5)

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
against something other than itself; **v1.3** made evidence auditable —
captions carry provenance, references are separated from canonical works, a
re-run cannot leave stale evidence, and an operator can inspect a parse
against the original page ([v1.3.0](https://github.com/caseywdunn/corpus/blob/v1.3.0/dev_docs/PLAN.md),
Zenodo [10.5281/zenodo.22647843](https://doi.org/10.5281/zenodo.22647843)); and
**v1.4** went after the results that were wrong without saying so and the build
hazards behind them, taking the tracker from 34 open issues to 13
([v1.4.0](https://github.com/caseywdunn/corpus/blob/v1.4.0/dev_docs/PLAN.md)).

**v1.4 (2026-09-09) is the one v1.5 grows directly out of.** Its theme
paragraph is in the CHANGELOG; what matters here is the inheritance. The
defects that were wrong without saying so are closed: a text layer of
unmappable glyph indices no longer counts as clean, OCR packs are chosen from
the page images and corroborated before they override, abbreviated binomials
resolve, and the MCP tools report what they returned rather than only whether
they truncated. Alongside them the build hazards — a fixed Grobid port, a GPU
allocation quietly becoming a CPU one, an unbounded build, silent long stages.

Three instruments carry forward and should keep being used: the gold set from
v1.2, the build-comparison harness (`tools/qc/build_reference.py`), and — new
in v1.4 — the **local measurement substrate**, which sizes a change against all
1,775 reference documents in minutes on a laptop. Nine of v1.4's twenty-two
items had their shape changed by measurement before they were written, six
against the issue's own proposal, and six needed deletion or nothing rather
than new code. **Measure before writing** is now the cheapest step, not the
expensive one.

Two things v1.4 established that v1.5 should not re-derive:

- **A comparison matrix with no self-comparison in it measures nothing.** The
  first local-VLM report compared every dtype against every other and against
  production's ROIs, and everything disagreed with everything — uninterpretable
  until two runs of one dtype established what the model's disagreement with
  *itself* looks like. That control is now the first thing to run, not the last.
- **Pass 3b geometry is not reproducible**, run to run on one accelerator or
  across accelerators. Measured, recorded in the standing gates, and the reason
  small ROI deltas between builds are not a defect signal on their own.

What v1.4 did *not* do is build anything a user talks to. That is v1.5.

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

## v1.4 — shipped 2026-09-09

Closed and preserved in [the tag's copy of this
file](https://github.com/caseywdunn/corpus/blob/v1.4.0/dev_docs/PLAN.md); the
theme paragraph is in [CHANGELOG.md](../CHANGELOG.md). Twenty-two items, the
tracker from 34 open issues to 13.

Two carried forward rather than closed, both surfaced by the release rebuild
and both real:

- **Stage 1's 256 GB is sized for the wrong thing.** One document OOMed at
  268 GB inside a full 8-worker task and peaked at 38 GB when re-run alone, so
  the bound that matters is worker concurrency, not per-document memory.
  #182's `docling` knobs exist and are unset; nothing yet says what to set them
  to. Raising `--mem` treats the symptom.
- **Discovery-mode figure counts are unexplained.** One plate went 56 -> 35
  `discovery_materialized` between the v1.3.0 and v1.4.0 builds on
  byte-identical input images. ROI *coordinates* are known-nondeterministic and
  excluded from the update contract; open-ended discovery *counts* are a
  different path and were never measured. Settle it by running vision twice on
  that document before treating it as noise.

Also worth an hour: CI lanes install via `apt-get update` on GitHub's runner
image, which ships a Google Chrome repo corpus never uses. A hash-sum mismatch
there took out T3 and T1/T2 on the release commit three times on 2026-09-09.
Dropping `/etc/apt/sources.list.d/google-chrome.list` before the update removes
the whole class.

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
  ([#283](https://github.com/caseywdunn/corpus/issues/283), superseding the
  [#98](https://github.com/caseywdunn/corpus/issues/98) follow-up). Still
  `docling==2.94.0`; current is **2.126.0**, so the pin is 32 minor versions
  behind. #98 asked whether 2.95 or 2.96 broke MPS extraction — **that box is
  dropped, not carried.** Thirty versions on it is archaeology: whatever broke
  has almost certainly been rewritten, and the answer would not say whether
  2.126 works. The question with a consumer is "does current docling extract
  correctly on arm64?", which is one test rather than a bisect. Needs
  Apple-Silicon hardware; budget for API churn rather than a version bump,
  since nobody has read docling's changelog across that range.

  Two things make this a smaller bet than when #98 was written. #99's guard
  now treats a corpus-wide zero-chunk result as a hard error, so the silent
  empty-bundle failure that made 2.96 dangerous fails loudly. And **v1.2's
  fidelity harness ([#193](https://github.com/caseywdunn/corpus/issues/193))
  gives this a criterion it never had** — "better or worse" against the gold
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

- [x] **The local VLM loads in float32 on MPS**
  ([#258](https://github.com/caseywdunn/corpus/issues/258)). Closed in v1.4:
  measured on an M2 Max and the default stays float32. Half precision keeps
  every panel but moves the boxes (mean IoU 0.65–0.75 against float32, worst
  panel 0.0), while two float32 runs agree at 1.0 — greedy decoding over digit
  tokens, so a flipped digit moves a coordinate hundreds of pixels. The memory
  premise stands unfixed: float32 measured 37.69 GB on Metal, so a 32 GB Mac
  still cannot run the default, and `figures.vision_dtype` is the fallback
  there with documented worse geometry.
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
