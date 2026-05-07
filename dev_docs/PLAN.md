# PLAN.md — Corpus pipeline (v0.2)

v0.1 (2026-05-01) shipped the full extraction → annotation → indexing →
MCP-serving stack on a ~1,800-paper siphonophore corpus. v0.2 hardens
v0.1's rough edges and fills in the deferred items: vision pass at
corpus scale, expanded taxonomy + lexicon coverage, a Streamable-HTTP
transport with OAuth, batch-script cleanup, granular per-stage resume,
robustness + quality gates against silent-failure modes, and a
`corpus_status` view of the build. Geographic extraction is pushed to
v2.0+ (see §4).

This document is scoped to v0.2 work. Architectural background and
pipeline internals live in [OVERVIEW.md](OVERVIEW.md);
per-feature history in [CHANGELOG.md](../CHANGELOG.md); HPC operations in
[BOUCHET.md](BOUCHET.md); deployment in
[DEPLOY.md](DEPLOY.md). Open work is tracked in
[GitHub issues](https://github.com/caseywdunn/corpus/issues); the
v0.2-scoped subset is listed below.

## 1. v0.2 punch list

### Implementation order

Items below are grouped thematically, but implemented in this order so
data shapes land before consumers and small wins land before
architectural lifts:

1. **Quick wins** — #6, #22, #25, #31, #41. Small fixes, a read-only
   audit, and `--dry-run` flags. Low risk; lands first to build
   momentum.
2. **Robustness layer** — #34 → #36 → #40, then #35, #37. Structured
   `stage_failures[]`, quality gates, `corpus_status.py`. Establishes
   the failure / status data model before anything restructures the
   pipeline.
3. **Update lifecycle** — #28 → #29, then #30, #32, #33. Granular
   per-stage resume + input fingerprints. #34's structured-failure
   schema lands first so #28 inherits it instead of inventing a
   parallel one.
4. **Slot in anywhere** — #11 + #27 (vision at corpus scale), #15
   (mcp_server refactor), batch-script cleanup (#20, #21), feature
   work (#23, #24, #26), #12, #16, #17. Independent enough to
   land in any session.

The thematic groupings below remain useful reference; read this list
when planning a session.

### Hardening — close out v0.1 rough edges

- [#9](https://github.com/caseywdunn/corpus/issues/9) — Install
  `deu_latf` Fraktur pack on Bouchet. ~6–8 19th-c. German scans
  (Goldfuss 1820, Pagenstecher 1869, Brandt 1837, Dönitz 1871, …)
  currently OCR to whitespace; reprocess them after the pack lands.
- [#16](https://github.com/caseywdunn/corpus/issues/16) — Figure-number
  extraction on old/scanned papers. ~538/1787 papers have `Fig. N`
  mentions in body text but no extracted figure number. Improve
  `parse_figure_number` heuristics or fall back to running-text
  mention positions when caption parsing fails.
- [#17](https://github.com/caseywdunn/corpus/issues/17) — MCP server
  identity. Partially addressed in v0.2-dev: `__version__` from
  [pipeline/version.py](../pipeline/version.py) now surfaces via the `bundle_info` tool.
  Remaining: server name should be corpuscle-aware (e.g.,
  `corpus:siphonophores` rather than `corpus`).
- [#6](https://github.com/caseywdunn/corpus/issues/6) — `deploy/stack.yaml`
  default-VPC assumption. Fails on accounts without a default VPC;
  parameterize or pre-flight.
- [#22](https://github.com/caseywdunn/corpus/issues/22) —
  `ingest_taxonomy.py` DwC field-name fix (`taxon` vs. `taxa`).
- [#25](https://github.com/caseywdunn/corpus/issues/25) — `pngquant`
  missing from `environment.yaml`.
- Bouchet batch-script cleanup
  ([#20](https://github.com/caseywdunn/corpus/issues/20),
  [#21](https://github.com/caseywdunn/corpus/issues/21)) —
  library-load problems across batch scripts, CUDA / torch fixes.

### Vision + figures

- [#11](https://github.com/caseywdunn/corpus/issues/11) — **Vision pass
  at corpus scale.** Pipeline is wired (`batch_pass3b.sh`,
  Qwen2.5-VL-7B local backend, Claude remote backend) and tested on
  the demo set; v0.1 didn't run it on the full corpus. Resolves the
  bulk of the 6,841 pre-existing `missing_figures[]` records and the
  compound-figure splits (Pass 3c).
- [#27](https://github.com/caseywdunn/corpus/issues/27) — Add SLURM
  array processing to `batch_pass3b.sh` so it parallelizes the same
  way `batch_pipeline.sh` does. Prerequisite for #11 at any reasonable
  wallclock.

### Indices + features

- [#26](https://github.com/caseywdunn/corpus/issues/26) — Round-trip
  bibliography: MCP export + shell import.
- [#23](https://github.com/caseywdunn/corpus/issues/23) — Expand
  taxonomy ingest beyond marine taxa (currently WoRMS-shaped
  assumptions); plant DwC archives as the test case.
- [#24](https://github.com/caseywdunn/corpus/issues/24) — Expand
  anatomy lexicon.

### Update lifecycle — make re-runs cheap and correct

The user-curatable inputs (anatomy lexicon, taxonomy snapshot, bib
file, input PDF set) all change between runs. v0.1 has no incremental
path for any of them: the only safe response to a lexicon edit or a
taxonomy bump is "rebuild every paper from Docling onward." This
section makes update flows first-class so iteration on inputs is
fast and obviously-correct.

- [#28](https://github.com/caseywdunn/corpus/issues/28) — **Granular
  per-stage resume.** Linchpin of the rest of this section. Replaces
  the all-or-nothing `summary.json` completion marker with per-stage
  status and renames the numbered passes to describe-what-they-do
  (`extract.py`, `parse_metadata.py`, …). Touches every batch script
  and the resume logic.
- [#29](https://github.com/caseywdunn/corpus/issues/29) — **Input
  fingerprints on annotation artifacts.** Stamp `lexicon_sha256` /
  `taxonomy_snapshot_date` into `chunks.json` / `taxa.json` /
  `anatomy.json`. Without this there is no signal for which papers'
  annotations are stale.
- [#30](https://github.com/caseywdunn/corpus/issues/30) — **Idempotency
  audit + tests** for the four post-pipeline scripts
  (`build_biblio_authority.py`, `build_taxon_mentions.py`,
  `backfill_intext_citations.py`, `reconcile_corpus_to_biblio.py`).
- [#31](https://github.com/caseywdunn/corpus/issues/31) — **Orphan
  detection.** `--audit-orphans` lists `documents/<HASH>/` (and
  LanceDB rows) whose source PDF is no longer in the input set.
  Read-only; deletion stays manual.
- [#32](https://github.com/caseywdunn/corpus/issues/32) —
  **`update_corpus.py` wrapper.** One command that runs the pipeline
  and post-pipeline scripts in dependency order with `--resume`.
  Makes "add papers and update everything" a one-liner.
- [#33](https://github.com/caseywdunn/corpus/issues/33) — **Lexicon
  round-trip.** Edit `lexicon.yaml` (multi-category, top-level keys
  `anatomy:` / `biogeography:` / …) and re-annotate only the affected
  categories on the next `--resume`. Subsumed by per-stage resume +
  per-category fingerprints in `pipeline_state.json` (#28 / #29).

### Robustness + observability

The pipeline is resilient (it doesn't crash on bad inputs) but it
isn't *diagnosable* — failures are free-text strings in `errors[]`,
silent-failure modes (gibberish OCR, empty extractions) ship without
a flag, and there's no single command to ask "what's the state of
this build." This section closes those gaps.

- [#34](https://github.com/caseywdunn/corpus/issues/34) — **Per-stage
  wallclock caps + structured `stage_failures[]`** with reason codes
  (`timeout`, `crash`, `external_unavailable`, `unsupported_format`,
  `corrupted`, `quality_gate`). Replaces free-text `errors[]`.
- [#35](https://github.com/caseywdunn/corpus/issues/35) — **Huge-document
  path.** Page-count gate + chunked OCR for 500+ page monographs
  (Haeckel 1888 *Challenger Siphonophorae* is the canary).
- [#36](https://github.com/caseywdunn/corpus/issues/36) — **Quality
  gates** that catch silent failures: empty text, gibberish OCR
  output, all-black figures, suspicious zero-reference papers,
  collapsed-extraction chunks.
- [#37](https://github.com/caseywdunn/corpus/issues/37) — **External-service
  flakiness.** Standardized retry + backoff + caching for Grobid,
  BHL, CrossRef, OpenAlex; circuit breakers; `--strict-network` for
  release builds.
- [#40](https://github.com/caseywdunn/corpus/issues/40) —
  **`corpus_status.py`.** Single command that reports stage
  completion, failure breakdown by reason, quality flags, and
  staleness. The dashboard for everything else in this section and
  in §1 Update lifecycle.
- [#41](https://github.com/caseywdunn/corpus/issues/41) —
  **Standardize `--dry-run`** across all pipeline + post-pipeline
  CLIs. Currently inconsistent (4 of 10 scripts have it).

### Internal

- [#15](https://github.com/caseywdunn/corpus/issues/15) — Refactor
  `mcp_server.py`. The 2050-line single file has outgrown a single
  module; split per-concern (papers / taxonomy / bibliography /
  figures / transports). No behavior change.
- [#12](https://github.com/caseywdunn/corpus/issues/12) — Stage 2
  parallelism. `batch_embed.sh` already parallelizes; assess whether
  further work is needed (per-document chunking concurrency, BGE-M3
  batch sizing).

## 2. Target queries (evergreen reference)

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`. v0.1 made Q2, Q5–Q7 fully answerable; Q1, Q4, Q8
remain partial pending v0.2 work; Q3 is gated on trait extraction
deferred post-v0.2 ([#14](https://github.com/caseywdunn/corpus/issues/14)).

| # | Pattern | Status after v0.1 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer (#13, deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; synthesis recipe not yet scoped |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred (#14, post-v0.2) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Needs vision pass at corpus scale (#11) |
| Q5 | "Summarize `<author X>`'s comments about `<author Y>`." | Indices in place |
| Q6 | "Summarize `<topic>` across the corpus." | Indices in place |
| Q7 | "Plot species described per decade." | Indices in place |
| Q8 | "Summarize what is known about `<anatomy>`." | Needs vision-pass figure tagging at corpus scale (#11) |

## 3. Versioning + release ritual

`__version__` in [pipeline/version.py](../pipeline/version.py) is the single source of
truth and is stamped into every persistent artifact (bundle manifest,
MCP `bundle_info`). [CONTRIBUTING.md](../CONTRIBUTING.md) covers the
branching model and release ritual; the short version: `dev` carries
a PEP 440 pre-release suffix (`0.2.0.dev0` → `0.2.0a1` → `0.2.0`),
the release commit drops the suffix, the next commit on `dev`
reintroduces one for the next target.

## 4. Out of scope for v0.2

- [#13](https://github.com/caseywdunn/corpus/issues/13) — Geographic
  mention layer (NER + GeoNames + locality table; the third layer in
  the `*_mentions` family alongside biblio + taxon). Deferred to
  v2.0+; the mention-layer surface is likely to be reworked at the
  major-version boundary.
- [#14](https://github.com/caseywdunn/corpus/issues/14) — Trait
  extraction + identification keys (Q3). Substantial enough to
  warrant its own plan section when picked up; v0.3 candidate.
- [#38](https://github.com/caseywdunn/corpus/issues/38) — Embedding
  model migration path (versioned LanceDB tables + dual-read).
  Design-only for v0.2; implement when a model swap is actually
  needed (BGE-M3 v2 or a switch to a different family).
- [#39](https://github.com/caseywdunn/corpus/issues/39) — MCP server
  lazy index loading. Premature at 1.8K papers; document as a known
  scaling cliff and revisit when a corpuscle pushes 10K+.
- [#5](https://github.com/caseywdunn/corpus/issues/5) — Streamable
  HTTP transport with OAuth. Deferred indefinitely. SSE +
  bearer-token works for the ~20-collaborator deploy target; revisit
  only if wider distribution actually materializes or if MCP clients
  drop SSE support.
- Multi-region failover, autoscaling, or blue/green deploys for the
  AWS served bundle. Single-instance is fine until it isn't.
- Authentication beyond bearer tokens / OAuth (Cognito, institutional
  SSO).
- Mirror to Cloudflare R2 / Backblaze B2 for cost. Defer until S3
  egress shows up on a bill.
- A thin HTML/web UI on top of the MCP server. Out of scope until the
  MCP-only experience has actual non-Claude-Desktop users.
