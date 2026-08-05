# PLAN.md — Corpus pipeline (v1.1)

v0.1 (2026-05-01) shipped the full extraction → annotation → indexing
→ MCP-serving stack. v0.2 hardened internals (per-stage resume, input
fingerprints, structured failure schema, quality gates). v0.3 collapsed
the user surface into one CLI plus a per-corpuscle `config.yaml`. v0.4
closed the silent-failure modes external testers and clean-room HPC
hosts surfaced, and gated the next cycle on tiered CI. v0.5 was the
served-bundle quality cycle — LLM-side citation amalgamation
([#79](https://github.com/caseywdunn/corpus/issues/79)), the dossier
suite ([#76](https://github.com/caseywdunn/corpus/issues/76)), and a
pinned ML stack after an unpinned docling release silently broke macOS
arm64 extraction ([#98](https://github.com/caseywdunn/corpus/issues/98),
[#99](https://github.com/caseywdunn/corpus/issues/99)).

**v0.6 (2026-06-04) was the API-freeze cycle and it landed.** The MCP
surface is frozen at **38 tools** with a uniform error payload and a
consistent pagination convention: the redundant singular tools were
removed in favor of batched plurals, `lexicon_matrix` went
summary-by-default, `get_citation_graph` gained breadth caps, the
output-type profile ontology
([#101](https://github.com/caseywdunn/corpus/issues/101)) shipped, and
`tests/test_freeze_contract.py` now guards the surface mechanically.
Ops work landed alongside: `/healthz` capability reporting
([#91](https://github.com/caseywdunn/corpus/issues/91)), per-invocation
run logs ([#90](https://github.com/caseywdunn/corpus/issues/90)),
`corpus debug-pdf` ([#92](https://github.com/caseywdunn/corpus/issues/92)),
and the `--figure-panels` selector
([#102](https://github.com/caseywdunn/corpus/issues/102)). The v0.6
punch list is preserved in the
[v0.6.0 tag's history](https://github.com/caseywdunn/corpus/blob/v0.6.0/dev_docs/PLAN.md).

**v1.0 (2026-08-05) was the installability cycle and it landed.** The
organizing principle was that **a green CI badge must mean a fresh
install works** — and that a *silent* badge means nothing at all. Both
halves shipped: **T3**, a clean-room lane that installs from
`environment.yaml` on a weekly clock with the HuggingFace cache disabled
so a genuine first-run model download happens, and **T1-compose**, which
boots the real `docker-compose.yml` on every push. Nothing had exercised
the compose file before, which is how a default image that crash-looped
on Apple Silicon ([#146](https://github.com/caseywdunn/corpus/issues/146))
and a healthcheck probing for an absent `curl`
([#157](https://github.com/caseywdunn/corpus/issues/157)) both reached
users. `mcp` had gone unpinned into 2.0.0 and broken every clean install
for 24 days ([#156](https://github.com/caseywdunn/corpus/issues/156));
the remaining nine pip deps are now pinned exactly
([#158](https://github.com/caseywdunn/corpus/issues/158)).

The secondary through-line was **no silent wrongness** — code that
produced plausible-looking bad output rather than errors. A bib parser
that discarded all but 2,258 of 19,834 entries on one escaped brace
([#141](https://github.com/caseywdunn/corpus/issues/141)), curated
references stuck in the warning tier
([#142](https://github.com/caseywdunn/corpus/issues/142),
[#152](https://github.com/caseywdunn/corpus/issues/152)), a taxonomy
stage that skipped into an empty `taxa.json` on compute nodes with no
internet ([#139](https://github.com/caseywdunn/corpus/issues/139)), and
OCR falling back to English on Cyrillic pages
([#160](https://github.com/caseywdunn/corpus/issues/160)).
`corpus prefetch` ([#159](https://github.com/caseywdunn/corpus/issues/159))
closed the first-run model download for offline and rate-limited hosts.
On the served surface, figure licensing stopped letting the client decide
([#154](https://github.com/caseywdunn/corpus/issues/154)) and lexicon
figure retrieval learned its own synonyms
([#143](https://github.com/caseywdunn/corpus/issues/143)).

Two external install reports drove much of it —
[#145](https://github.com/caseywdunn/corpus/issues/145) (macOS Apple
Silicon) and [#153](https://github.com/caseywdunn/corpus/issues/153)
(HPC) — and the v0.6 38-tool freeze held throughout: no tool was added,
and where 1.0 changed served behavior it was to make an existing tool
honest. The formal API-stability policy deferred out of v0.6 shipped as
[API_STABILITY.md](API_STABILITY.md), and distribution was settled
explicitly rather than by omission: git + conda stays canonical, with
Zenodo for citability. Validation was a full production run on Bouchet —
**1,769 of 1,769 documents, 0 stage failures, 261,093 chunks embedded,
934 of 934 eligible figures through the vision pass, and no job in the
chain ending TIMEOUT.** The v1.0 punch list is preserved in the
[v1.0.0 tag's history](https://github.com/caseywdunn/corpus/blob/v1.0.0/dev_docs/PLAN.md).

**v1.1 is the post-1.0 correctness cycle.** Its contents were triaged at
the 1.0 release rather than absorbed into it: the pre-release UX review
(#164–#173), the viburnum production build (#174–#177, of which only #177
was pulled into 1.0), and two items carried across several cycles
([#98](https://github.com/caseywdunn/corpus/issues/98),
[#155](https://github.com/caseywdunn/corpus/issues/155)). Nothing in it
blocks an install, which is the bar 1.0 set.

The freeze permits all of it: every item is a bug fix or an additive
change under [API_STABILITY.md](API_STABILITY.md), so none waits for a
major. **Scope is not yet fixed** — §1 is a candidate list in
dependency-free groups, not a wave plan. The cycle's shape should come
from whichever group turns out to matter to real corpora; the viburnum
items (#175, #165, #176) are the ones with a build behind them.

Doc map unchanged: architectural background in
[OVERVIEW.md](OVERVIEW.md); per-feature history in
[CHANGELOG.md](../CHANGELOG.md); the API contract in
[API_STABILITY.md](API_STABILITY.md); HPC operations in
[BOUCHET.md](BOUCHET.md); deployment in [DEPLOY.md](../DEPLOY.md);
platform-portability criteria in
[PLATFORM_SMOKE.md](PLATFORM_SMOKE.md). Open work is tracked in
[GitHub issues](https://github.com/caseywdunn/corpus/issues).

## 0. Standing gates

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

**Baselines need re-measuring.** The last recorded ones were taken at
`1dbb69a` — T0 at 623 passed / 2 skipped, `corpus_required` at 163
passed / 14 skipped / 4 xfailed, SSE smoke at 3 layers / 38 tools — but
that predates both the release commit and
[#167](https://github.com/caseywdunn/corpus/issues/167), which changed
what the corpus-wide comparisons assert and therefore which `--deselect`
flags T1 needs. Re-measure on `v1.0.0` before using any of it as a
regression reference (§1, Housekeeping).

## 1. v1.1 candidate list

No wave ordering — none of these blocks another, and none is a
precondition for verifying the rest, which is what made v1.0's waves
necessary.

### Served-surface correctness

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
- [ ] **Reference reconciliation**
  ([#155](https://github.com/caseywdunn/corpus/issues/155)).
  `get_missing_references` is dominated by unreconciled citation-string
  variants of works already in the corpus —
  `resolve_reference "Bigelow 1911 The Siphonophorae"` returns 53
  matches, one canonical and ~50 variants. Doing it properly means DOI
  normalization, block-and-cluster canonicalization, alternate-DOI
  aliasing, junk filtering, and probably an auditable LLM adjudication
  pass: a cycle of its own. **The cheap slice** is the tool-level
  mitigation — collapse candidates by normalized (author, year, title)
  before ranking, so the same work stops appearing as N rows. That makes
  the tool honest without fixing reconciliation, and it is additive.

### Extraction quality

- [ ] **Tesseract OSD overrides a p=1.0 language detection**
  ([#176](https://github.com/caseywdunn/corpus/issues/176)) with no union
  and no `eng` fallback, and there is no per-document override — 6 of 123
  OCR'd papers in the viburnum build.
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
  then advance deliberately. Needs Apple-Silicon hardware.

### Operator surface

- [ ] **Progress heartbeat during long per-document stages**
  ([#170](https://github.com/caseywdunn/corpus/issues/170)).
- [ ] **`--filter-gate` is silently ignored without `--list-hashes`**
  ([#169](https://github.com/caseywdunn/corpus/issues/169)) — the hint
  was fixed in 1.0, the underlying flag was not.
- [ ] **Surface the naive-chunker fallback in `corpus status`**
  ([#168](https://github.com/caseywdunn/corpus/issues/168)).
- [ ] **Extend per-stage fingerprints to bib entries, filenames and
  config keys** ([#174](https://github.com/caseywdunn/corpus/issues/174)),
  so editing the input bib stops being a silent no-op.

### Housekeeping

- [ ] **`pipeline/CITATION.cff` is a hand-maintained copy**
  ([#173](https://github.com/caseywdunn/corpus/issues/173)) with nothing
  enforcing sync — currently byte-identical, so this is a guard rather
  than a fix.
- [ ] **README is ambiguous about where `instructions.md` lives**
  ([#171](https://github.com/caseywdunn/corpus/issues/171)).
- [ ] **Add `tools/` to the pyflakes gate** in
  `tests/test_no_undefined_names.py`. It lints `pipeline/`, `mcpsrv/` and
  `bib/` only, so the Python scripts under `tools/` never get the
  NameError check [#75](https://github.com/caseywdunn/corpus/issues/75)
  built it for.
- [ ] **Re-measure the test baselines** on `v1.0.0` and record which
  `--deselect` flags T1 still needs after #167 (see §0).
- [ ] **16 source comments cite PLAN.md sections that no longer exist** —
  `§10` (×8), `§9` (×4), `§7` (×1), and `§3` (×3, where the current §3 is
  a different topic), across `mcpsrv/`, `pipeline/`, `tests/` and
  BOUCHET.md. They point at a v0.1-era numbering this file has been
  pruned past several times. Either cite the stable doc (OVERVIEW.md,
  DEPLOY.md) or drop the pointer; a roadmap that gets pruned every cycle
  is the wrong anchor for a code comment.

### Features, additive by construction

- [ ] **A `skills/` plugin directory**
  ([#178](https://github.com/caseywdunn/corpus/issues/178)), **a
  clade-monograph skill**
  ([#179](https://github.com/caseywdunn/corpus/issues/179)), and **a
  README quick start**
  ([#180](https://github.com/caseywdunn/corpus/issues/180)).

### External contribution

- [ ] **[PR #144](https://github.com/caseywdunn/corpus/pull/144)** from
  @beroe — original-description linking against in-corpus works
  (423 → 505 of 598 ctenophore taxa). Mergeable and substantive, but it
  is a 4-commit external change that has never run CI, and it touches
  `bib/authority.py`, which #154 rewrote during the 1.0 cycle. It wants
  a review and a CI run, not a fast merge. Note it overlaps #175: both
  are `parse_authority`, from opposite ends.

## 2. Target queries (evergreen reference)

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`.

| # | Pattern | Status entering v1.1 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer ([#13](https://github.com/caseywdunn/corpus/issues/13), deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; citation-trust gap closed by [#79](https://github.com/caseywdunn/corpus/issues/79) in v0.5, provenance preservation by [#142](https://github.com/caseywdunn/corpus/issues/142) / [#152](https://github.com/caseywdunn/corpus/issues/152) in v1.0; synthesis recipe not yet scoped |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred ([#14](https://github.com/caseywdunn/corpus/issues/14)) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Indices in place; the corpus-scale vision run landed with v1.0's Bouchet production run (934/934 eligible figures), the figure-coverage audit is still undone (untracked — see §4) |
| Q5 | "Summarize `<author X>`'s comments about `<author Y>`." | Indices in place |
| Q6 | "Summarize `<topic>` across the corpus." | Indices in place; cache cost addressed by dossier tools [#76](https://github.com/caseywdunn/corpus/issues/76) in v0.5 |
| Q7 | "Plot species described per decade." | Indices in place |
| Q8 | "Summarize what is known about `<anatomy>`." | Indices in place; figure-retrieval synonym blindness fixed by [#143](https://github.com/caseywdunn/corpus/issues/143) in v1.0, though translations still miss inflected forms ([#165](https://github.com/caseywdunn/corpus/issues/165), §1); figure-coverage audit still undone (untracked — see §4) |

## 3. Versioning + release ritual

`__version__` in [pipeline/version.py](../pipeline/version.py) is the
single source of truth and is stamped into every persistent artifact
(bundle manifest, MCP `bundle_info`).
[CONTRIBUTING.md](../CONTRIBUTING.md) covers the branching model and
release ritual. Note **step 8** — pruning this document at release time.
It was skipped at v0.6.0, which is why this file described a finished
cycle in the present tense for two months; it was carried out at v1.0.0.

## 4. Carryover to v1.2+ (not scheduled for v1.1)

Longer-horizon than §1 — either net-new features (safe to add after 1.0
without breaking the frozen surface) or work awaiting motivation. Carried
so they stay visible.

- **Container distribution image.** The only channel that can ship the
  full native toolchain in one artifact: Docker is already a
  prerequisite for Grobid, `docker-compose.yml` could bring up grobid +
  corpus together, and #153's HPC user already runs Apptainer, which
  pulls straight from a Docker registry. Costs: bind-mounting the PDF
  directory, GPU passthrough for the local VLM, and a large image with
  torch in it. 1.0's install path is verified continuously now, which was
  the precondition — worth its own issue.
- **`export_to_disk` + `suggest_command`**
  ([#88](https://github.com/caseywdunn/corpus/issues/88) Part 2) and the
  **`corpus export` CLI** ([#93](https://github.com/caseywdunn/corpus/issues/93))
  — bulk-export surface for "the LLM isn't the consumer" workflows.
  Additive, so they don't constrain the freeze.
- **`verify_claim`** ([#123](https://github.com/caseywdunn/corpus/issues/123)).
  Per-claim ledger anchoring as a thin similarity-only wrapper over
  `get_chunks_for_topic`. New tool — post-freeze by construction.
- **Drift detection** ([#80](https://github.com/caseywdunn/corpus/issues/80)).
  Pre-run diff of resolved config + per-input content SHAs.
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
  **Untracked** — file an issue if picked up.
- **Vision pass corpus-scale validation.**
  [#11](https://github.com/caseywdunn/corpus/issues/11) is **closed** as
  carried-out-in-code: `corpus run` invokes the vision pass whenever
  `figures.panel_detection` selects a vision backend and the host
  capability check passes. The corpus-scale run happened in v1.0 — 934 of
  934 eligible figures through the vision pass on Bouchet. What remains is
  the **figure-coverage audit**: count figures with `pass3c_status` set,
  sum `missing_figures[]` lengths, and establish what "eligible" excluded.
  **Untracked** — release-validation work, not coding work.
- **Evaluate Cloud Run vs the EC2+ALB stack**
  ([#89](https://github.com/caseywdunn/corpus/issues/89)). Deployment
  decision, not part of the MCP API contract.

## 5. Out of scope (longer horizon)

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
