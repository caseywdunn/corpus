# PLAN.md — Corpus pipeline (v0.6)

v0.1 (2026-05-01) shipped the full extraction → annotation → indexing
→ MCP-serving stack. v0.2 (2026-05-08) hardened internals: per-stage
resume, input fingerprints, structured failure schema, quality gates,
and an ALB-fronted multi-organism deploy. v0.3 (2026-05-11) collapsed
the user surface: one CLI installed once, a per-corpuscle
`config.yaml`, automatic database management, capability-aware
defaults. v0.4 (2026-05-17) closed the silent-failure modes external
testers and clean-room HPC hosts surfaced, and gated the next cycle
on tiered CI. v0.5 (2026-05-29) was the served-bundle quality cycle:
it closed the two trust-breaking gaps external evaluators found in the
MCP surface — LLM-side citation amalgamation (`format_citation` +
provenance cascade, [#79](https://github.com/caseywdunn/corpus/issues/79))
and the ~2 M cache_read enumerative access pattern (the six-tool
dossier suite, [#76](https://github.com/caseywdunn/corpus/issues/76)) —
plus token-efficiency follow-ups
([#81](https://github.com/caseywdunn/corpus/issues/81),
[#82](https://github.com/caseywdunn/corpus/issues/82),
[#84](https://github.com/caseywdunn/corpus/issues/84),
[#86](https://github.com/caseywdunn/corpus/issues/86)). The release
itself surfaced one more silent-failure mode and its cause: an
unpinned docling/torch stack let a same-day docling release silently
break macOS arm64 extraction, which shipped an empty bundle reporting
success. v0.5 closed both — pinning the ML stack to a known-good set
([#98](https://github.com/caseywdunn/corpus/issues/98)) and making
`corpus run` fail loudly on a zero-yield extraction
([#99](https://github.com/caseywdunn/corpus/issues/99)).

**v0.6 theme is not yet set.** The cycle opens with a candidate pool
(below) drawn from v0.5 carryover and the issues filed during the v0.5
cycle, but the through-line and punch list still need to be chosen.
Treat §1 as a backlog to prioritize, not a committed scope.

Doc map unchanged: architectural background in
[OVERVIEW.md](OVERVIEW.md); per-feature history in
[CHANGELOG.md](../CHANGELOG.md); HPC operations in
[BOUCHET.md](BOUCHET.md); deployment in [DEPLOY.md](DEPLOY.md); the
full v0.5 plan archived in the
[v0.5.0 tag's history](https://github.com/caseywdunn/corpus/blob/v0.5.0/dev_docs/PLAN.md).
Open work is tracked in
[GitHub issues](https://github.com/caseywdunn/corpus/issues).

## 1. v0.6 candidate pool (unprioritized)

Filed and visible, not yet scoped into a punch list. Grouped by theme
so a v0.6 through-line is easier to pick.

### Carryover from v0.5

- **Drift detection**
  ([#80](https://github.com/caseywdunn/corpus/issues/80)).
  Diff the resolved config + per-input content SHAs across runs and
  report the diff before the orchestrator starts, so an operator can
  see *why* a re-run is doing more work than expected.
- **Column-store shape for `lexicon_matrix`**
  ([#83](https://github.com/caseywdunn/corpus/issues/83)).
  Optional parallel-arrays + `row_schema` row shape instead of
  repeated JSON keys per row (~20–30% token saving at 100-row scale).
  Held pending a prompt-suite analysis showing large-matrix usage
  matters.
- **Figure-number extraction on old/scanned papers**
  ([#16](https://github.com/caseywdunn/corpus/issues/16)).
  ~538 of 1,787 papers have unparsed figure numbers; modest heuristic
  work.
- **Vision pass corpus-scale validation**
  ([#11](https://github.com/caseywdunn/corpus/issues/11)).
  Operational run + figure-coverage audit on Bouchet against the v0.5
  build. Release-validation, not coding.
- **Move the docling pin forward**
  ([#98](https://github.com/caseywdunn/corpus/issues/98) follow-up).
  v0.5 pinned `docling==2.94.0` (+ torch/transformers/sentence-
  transformers) because docling 2.95/2.96 silently break macOS arm64
  (MPS) extraction. Reproduce on an arm64 Mac, determine whether it's
  a docling API change to adapt to or an upstream bug to report, then
  advance the pin deliberately.

### MCP surface refinements (filed during v0.5)

- **Caption-preview default on older figure tools**
  ([#85](https://github.com/caseywdunn/corpus/issues/85)).
- **Breadth + edge caps on `get_citation_graph` traversal**
  ([#87](https://github.com/caseywdunn/corpus/issues/87)).
- **MCP tooling evaluation + consideration**
  ([#88](https://github.com/caseywdunn/corpus/issues/88)).

### Operations + observability (filed during v0.5)

- **MCP server health checks + structured logging**
  ([#91](https://github.com/caseywdunn/corpus/issues/91)).
- **Central per-invocation log for non-dry-run runs**
  ([#90](https://github.com/caseywdunn/corpus/issues/90)).
  Pairs with the v0.5 `run.log` work ([#57]).
- **Evaluate Cloud Run vs the EC2+ALB stack**
  ([#89](https://github.com/caseywdunn/corpus/issues/89)).

### Pipeline + CLI (filed during v0.5)

- **`corpus export` CLI — streaming bulk exports outside the MCP
  server** ([#93](https://github.com/caseywdunn/corpus/issues/93)).
- **Standalone single-PDF debug runner**
  ([#92](https://github.com/caseywdunn/corpus/issues/92)).
- **Flip the figure-publishability gate default** to permissive,
  opt-in strict ([#94](https://github.com/caseywdunn/corpus/issues/94)).
- **`build_taxon_mentions` skips papers regardless of `taxa.json`
  freshness** ([#95](https://github.com/caseywdunn/corpus/issues/95)).
- **WoRMS DwC-A excludes `isMarine=0` records** — doc gap that bites
  freshwater/terrestrial users
  ([#96](https://github.com/caseywdunn/corpus/issues/96)).
- **HuggingFace token warning**
  ([#97](https://github.com/caseywdunn/corpus/issues/97)).

## 2. Target queries (evergreen reference)

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`.

| # | Pattern | Status entering v0.6 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer ([#13](https://github.com/caseywdunn/corpus/issues/13), deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; citation-trust gap closed by [#79](https://github.com/caseywdunn/corpus/issues/79) in v0.5; synthesis recipe not yet scoped |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred ([#14](https://github.com/caseywdunn/corpus/issues/14), v0.6+ candidate) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Vision pass on by default since v0.3; full-corpus run pending [#11](https://github.com/caseywdunn/corpus/issues/11) |
| Q5 | "Summarize `<author X>`'s comments about `<author Y>`." | Indices in place |
| Q6 | "Summarize `<topic>` across the corpus." | Indices in place; cache cost addressed by dossier tools [#76](https://github.com/caseywdunn/corpus/issues/76) in v0.5 |
| Q7 | "Plot species described per decade." | Indices in place |
| Q8 | "Summarize what is known about `<anatomy>`." | Vision pass on by default since v0.3; full-corpus run pending [#11](https://github.com/caseywdunn/corpus/issues/11) |

## 3. Versioning + release ritual

`__version__` in [pipeline/version.py](../pipeline/version.py) is the
single source of truth and is stamped into every persistent artifact
(bundle manifest, MCP `bundle_info`).
[CONTRIBUTING.md](../CONTRIBUTING.md) covers the branching model and
release ritual.

## 4. Out of scope (carried from v0.5)

- [#14](https://github.com/caseywdunn/corpus/issues/14) — Trait
  extraction + identification keys (Q3). Substantial enough to
  warrant its own plan section when picked up; v0.6+ candidate.
- [#13](https://github.com/caseywdunn/corpus/issues/13) — Geographic
  mention layer. Deferred to v2.0+; the mention-layer surface is
  likely to be reworked at the major-version boundary.
- [#38](https://github.com/caseywdunn/corpus/issues/38) — Embedding
  model migration path. Design-only; implement when a model swap is
  actually needed (BGE-M3 v2 or a switch to a different family).
- [#39](https://github.com/caseywdunn/corpus/issues/39) — MCP server
  lazy index loading. Premature at 1.8K papers; documented as a
  known scaling cliff; revisit when a corpuscle pushes 10K+.
- [#5](https://github.com/caseywdunn/corpus/issues/5) — Streamable
  HTTP transport with OAuth. Deferred indefinitely. SSE +
  bearer-token works for the ~20-collaborator deploy target.
- Multi-region failover, autoscaling, or blue/green deploys for the
  AWS served bundle. Single-instance per corpuscle is fine until
  it isn't.
- Authentication beyond bearer tokens / OAuth (Cognito,
  institutional SSO).
- Mirror to Cloudflare R2 / Backblaze B2 for cost. Defer until S3
  egress shows up on a bill.
- A thin HTML/web UI on top of the MCP server. Out of scope until
  the MCP-only experience has actual non-Claude-Desktop users.
