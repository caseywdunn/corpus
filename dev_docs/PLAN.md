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
mean a fresh install works, validated by 1,769 of 1,769 documents on
Bouchet with 0 stage failures; **v1.1.1** made `CITATION.cff` schema-valid,
which is what had been silently blocking Zenodo archival since v0.3.0; and
**v1.2** built the gold transcription set and scored the extractor against
it — the first signal this project has had that measures the pipeline
against something other than itself.

**v1.2 (2026-08-29) is the one v1.3 grows directly out of.** Its theme
paragraph is in the CHANGELOG; what matters here is the inheritance. The
gold set is now a standing instrument rather than a one-cycle project — any
extraction change can be scored against it before and after, which is the
condition v1.3's work inherits and should keep using. Two of v1.2's items
did not finish and are carried below with their reasons: #187's drift half
and #217.

v1.2 also left a standing lesson worth stating once, because it recurred
three times in that cycle and once more at its own release: **a measurement
is not a result until you have looked at what it is measuring.** The
plate-legend fix shipped on a recall gain and quietly cost precision;
#254's OCR page-blanking was invisible to every quality gate the pipeline
had, because each one measured the document rather than the machine that
built it; and a stray file sat in the repo root through four green CI runs
and two release PRs because nothing looks there.

**v1.3 is the skills-and-usage cycle.** Where v1.2 asked whether the corpus
is *right*, v1.3 asks whether anyone can *use* it. Scope below, and as in
every cycle since v1.1 it is a candidate list rather than a commitment —
the material is expected to redirect it.

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

## v1.3 — skills and usage

The packaging and recipes that turn 38 frozen tools into an answer to a
question someone actually has. Additive by construction: skills live beside
the server, not inside the API contract, so nothing here reopens the surface
v0.6 froze.

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
  v1.2's gold-corpuscle work, where the accuracy half landed and the drift
  half did not.
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
- [ ] **Extend per-stage fingerprints to bib entries, filenames and
  config keys** ([#174](https://github.com/caseywdunn/corpus/issues/174)),
  so editing the input bib stops being a silent no-op. The per-document bib
  directives v1.2 shipped rest on it — #188 depends on this or on a
  documented caveat.

**Housekeeping**

- [ ] **`pipeline/CITATION.cff` is a hand-maintained copy**
  ([#173](https://github.com/caseywdunn/corpus/issues/173)) with nothing
  enforcing sync — currently byte-identical, so this is a guard rather
  than a fix. `tests/test_citation_cff.py` (v1.1.1) now validates both
  files independently, which is adjacent but not the same check.
- [ ] **README is ambiguous about where `instructions.md` lives**
  ([#171](https://github.com/caseywdunn/corpus/issues/171)).
- [ ] **Nothing in CI looks at the repo root.** A stray 9-byte file named
  `6` — a `%PDF-1.4` fragment from a shell redirect typo in `ec63c92` —
  sat next to `README.md` through four green CI runs and a full release
  PR, and was caught by eye at the v1.2.0 tag boundary, one merge from
  being in a citable Zenodo archive permanently. T0 lints `pipeline/`,
  `mcpsrv/`, `bib/` and `tools/` for undefined names; no check has an
  opinion about the top level. An allowlist test over `git ls-files |
  grep -v /` is cheap and would have failed on the offending commit.
  Same shape as the `tools/` pyflakes gap below: the check that would
  have caught it did not exist because nobody had been bitten yet.
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
  **Untracked** — file an issue if picked up. v1.2's figure-fidelity
  scoring ([#194](https://github.com/caseywdunn/corpus/issues/194)) is what
  would size it.
- **Vision pass corpus-scale validation.**
  [#11](https://github.com/caseywdunn/corpus/issues/11) is **closed** as
  carried-out-in-code: `corpus run` invokes the vision pass whenever
  `figures.panel_detection` selects a vision backend and the host
  capability check passes. The corpus-scale run happened in v1.0 — 934 of
  934 eligible figures through the vision pass on Bouchet. What remains
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

| # | Pattern | Status entering v1.3 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer ([#13](https://github.com/caseywdunn/corpus/issues/13), deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; citation-trust gap closed by [#79](https://github.com/caseywdunn/corpus/issues/79) in v0.5, provenance preservation by [#142](https://github.com/caseywdunn/corpus/issues/142) / [#152](https://github.com/caseywdunn/corpus/issues/152) in v1.0; synthesis recipe scoped for v1.3 as the clade-monograph skill ([#179](https://github.com/caseywdunn/corpus/issues/179)) |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred ([#14](https://github.com/caseywdunn/corpus/issues/14)) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Indices in place; the corpus-scale vision run landed with v1.0's Bouchet production run (934/934 eligible figures); figure coverage became measurable against truth in v1.2 (0.923 recall at 0.967 precision on the served surface) |
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
