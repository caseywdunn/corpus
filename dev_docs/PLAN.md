# PLAN.md — Corpus pipeline (v0.4)

v0.1 (2026-05-01) shipped the full extraction → annotation → indexing
→ MCP-serving stack. v0.2 (2026-05-08) hardened internals: per-stage
resume, input fingerprints, structured failure schema, quality gates,
and an ALB-fronted multi-organism deploy. v0.3 (2026-05-11) collapsed
the user surface: one CLI installed once, a per-corpuscle
`config.yaml`, automatic database management, capability-aware
defaults.

**v0.4 is the operational hardening cycle.** v0.3 left the machine in
a usable shape, but the first real exposure to external testers and a
clean-room HPC environment surfaced several silent-failure modes
(papers lost mid-build, doubled corpora on re-run, the served bundle
refusing to render). v0.4 closes those gaps and gates the next cycle
on CI that didn't exist before.

The framing is intentionally narrower than v0.3's: no breaking
changes, no new surface area, just the corner cases that turn
"usable" into "doesn't quietly fail on the second user." v0.3.x
corpuscles work as-is.

Doc map unchanged: architectural background in
[OVERVIEW.md](OVERVIEW.md); per-feature history in
[CHANGELOG.md](../CHANGELOG.md); HPC operations in
[BOUCHET.md](BOUCHET.md); deployment in [DEPLOY.md](DEPLOY.md); the
full v0.3 plan archived in the
[v0.3.1 tag's history](https://github.com/caseywdunn/corpus/blob/v0.3.1/dev_docs/PLAN.md).
Open work is tracked in
[GitHub issues](https://github.com/caseywdunn/corpus/issues).

## 1. v0.4 punch list

### Tiered CI on GitHub Actions

The class of regression that took weeks to surface in v0.3 (bundle
audit, LanceDB shape break, `instructions.md` shepherding) had no
automated detection path. v0.4 adds GHA tiers that run on every
push + every PR so future regressions block PRs instead of riding
into the next pre-release smoke. See
[CONTRIBUTING.md](../CONTRIBUTING.md) §"Continuous integration" for
the full tier table.

- [x] **T0 — lint + unit, every push, every branch**
  ([#75](https://github.com/caseywdunn/corpus/issues/75)).
  pyflakes gate + the ~314 unit tests with no corpus dependency,
  ~3 min wall.
- [x] **T1 + T2 — demo build + serve + resume scenario, on every
  push + every PR**
  ([#75](https://github.com/caseywdunn/corpus/issues/75)).
  T1 (Linux + Grobid) and T2 (macOS arm64, `grobid.disable: true`)
  each build the 4-paper demo, run `corpus_required` tests, then
  the 4 + 1 implicit-resume scenario inline. ~13 + 9 min wall;
  matrix bounded by T1. Two markers added (`corpus_required` and
  `resume_scenario`) so local pytest can carve out the same
  subset.
- [x] **Status badges + release ritual gated on CI green**
  ([#75](https://github.com/caseywdunn/corpus/issues/75)).
  README badge → `main` (release state), CONTRIBUTING badge →
  `dev` (active development line). Release ritual now has two
  explicit checkpoints: dev CI green before bumping the version,
  main CI green before tagging.

### Silent-failure cleanup on HPC + Apple Silicon

Four issues with the same shape: pipeline kept running but quietly
lost output. All have regression checks now in T1/T2.

- [x] **`find_all_pdfs` re-ingesting `processed.pdf` from prior runs**
  ([#75](https://github.com/caseywdunn/corpus/issues/75)).
  `input_pdfs: .` + `output_dir: ./output` (the demo's layout)
  doubled the corpus on every re-run. New `exclude_under` keyword
  on `find_all_pdfs`; unit test in `tests/test_find_all_pdfs.py`.
- [x] **LanceDB `list_tables()` return-shape break**
  ([#71](https://github.com/caseywdunn/corpus/issues/71)).
  0.30.x changed list → generator-view, broke implicit resume
  silently. Wrapped in a one-shot list materializer; regression
  inline in T1/T2.
- [x] **Bundle absolute-path audit on per-category lexicon JSONs**
  ([#70](https://github.com/caseywdunn/corpus/issues/70)).
  `_scrub_taxa` only knew `taxa.json`; multi-category lexicons
  emit one `<hash>/<category>.json` apiece, all leaking
  `lexicon.yaml`'s absolute path. Generalized to
  `_scrub_input_fingerprint_path` over every per-category artifact;
  CI greps the `audit clean.` log line.
- [x] **`TORCH_COMPILE_DISABLE=1` for docling/inductor g++ crash**
  ([#72](https://github.com/caseywdunn/corpus/issues/72)).
  `torch._inductor` JIT-shells `g++` at runtime; HPC nodes
  without GCC on PATH crash silently per affected paper.
  External tester confirmed `module load GCC` alone doesn't fix
  it (Bouchet exposes `g++` but not `cc1plus`). Set in
  `pipeline/__init__.py` before any torch import.
- [x] **macOS KMP duplicate-libomp abort**.
  Three OpenMP runtimes (torch, scikit-learn, conda-forge) collide
  at `import torch` on Apple Silicon; first call dies before any
  user code. `KMP_DUPLICATE_LIB_OK=TRUE` set on darwin in
  `pipeline/__init__.py` (matches what `tools/run_mcp_server.sh`
  already did for the serve path since v0.2).
- [x] **Grobid JDK cgroup-v2 NPE on modern Ubuntu hosts**.
  `lfoppiano/grobid:0.8.1`'s JVM hits a known NPE in cgroup-v2
  container-aware sizing on Ubuntu 24.04 (the GHA runner default).
  `JAVA_OPTS=-XX:-UseContainerSupport` + explicit `-Xmx` /
  `-Xms` in CI's Grobid service and in the EC2 smoke script.
- [x] **`instructions.md` shepherded into `output_dir` + served
  bundle**.
  Documented as a corpuscle-root file but every downstream reader
  looked under `<output_dir>/`. `corpus run` now copies it at the
  end of the pipeline.

### Demo as CI fixture

The 11-paper demo took ~25 min on a warm macOS arm64 laptop —
unviable as a per-PR CI fixture. Slimming + a local taxonomy
snapshot drop it to ~2 min on the same hardware.

- [x] **Demo slimmed from 11 → 4 + 1 papers**.
  Four PDFs in `demo/` (born-digital × 2, Fraktur, Cyrillic)
  cover every OCR + extraction path the 11-paper version did.
  One PDF held back in `tests/fixtures/round2_paper/` for the
  4 + 1 implicit-resume scenario.
- [x] **Pre-built Siphonophorae DwC-A in `demo/`**.
  `demo/taxonomy.zip` is the order-level WoRMS snapshot baked
  via the new export verb; demo's `config.yaml` switches to
  `source: dwca`. First `corpus run` ingests the taxonomy from a
  local file in seconds — no rate-limited WoRMS REST walk.
- [x] **`corpus taxonomy export | ingest` — DwC-A round-trip**.
  First-class verb to dump a built `taxonomy.sqlite` back out as
  a Darwin Core Archive `.zip`. Enables the demo-as-fixture move
  + lets operators share snapshots without forcing recipients to
  walk WoRMS again.

### Install onboarding

External-tester feedback ([#73](https://github.com/caseywdunn/corpus/issues/73))
caught two friction points the maintainer-eye missed.

- [x] **README Installation section reordered**.
  New Prerequisites subsection (Docker, miniforge) ahead of any
  conda command; new Clone-and-install block starts with
  `git clone` and includes `tools/install_tessdata.sh` as the
  fourth canonical step.
- [x] **Config template — indentation note above the `taxonomy:`
  block**.
  beroe stripped the leading `# ` from `#   source: worms` and
  hit a YAML "block mapping" parse error. Template now spells
  out the 2-space contract.

### Quick fixes

- [x] **pyflakes static-analysis gate**
  (`tests/test_no_undefined_names.py`).
  Catches undefined-name (guaranteed runtime `NameError`)
  findings; surfaced and fixed nine pre-existing imports during
  adoption.
- [x] **`corpus check` silent status lines + bib NameError**.
  Validation lines now stream as each check runs; the
  `bib.authority` staged-bib `db_path` NameError that the
  pyflakes gate caught is fixed.
- [x] **`mcpsrv/tools/chunks.py` missing `json` + `EmbeddingError`
  imports**.
- [x] **`tools/smoke_test_sse.py` repaired for v0.3 CLI**.
- [x] **EC2 smoke `pngquant` + `jbig2enc` best-effort**.
  AMIs without `universe` channel no longer abort the smoke;
  `pipeline.scan` already degrades gracefully when these helpers
  aren't present.
- [x] **EC2 clean-room smoke script + programmatic gate**
  (`dev_docs/ec2_smoke.sh`).
  One-shot platform-portability check from a bare Ubuntu EC2
  instance; exits 0 iff every PLATFORM_SMOKE.md criterion passes.

### Docs

- [x] **README — "Using the MCP server" section** between
  remote-deploy and the additional-docs index. Covers the
  text-only-chat vs. report-generation split, with character
  matrices + CSV / TSV / JSONL / BibTeX exports as the report-
  path examples.
- [x] **README — linked first mentions** of Yale HPC (Bouchet
  docs), AWS, Darwin Core, World Flora Online.
- [x] **CI tier table in CONTRIBUTING**, not README — the table
  is contributor-facing, README keeps the badges only.
- [x] **PLATFORM_SMOKE retitled** to "manual fallback /
  release-time verification" with a banner pointing at the CI
  tiers as the authoritative coverage. Same content but the
  framing is now "what CI doesn't catch" rather than "the
  primary check."

## 2. Carryover to v0.5

Originally scoped to v0.3, didn't fit in v0.4 either, no fresh
issue yet for any of them. Carried over so they stay visible.

- **Drift detection.** Hash the resolved config (input paths +
  per-input content SHA) into a corpuscle-side state file so a
  re-run that's doing more than expected can show *why*.
- **Success summary on clean `corpus run`** (half-shipped in v0.3 —
  the next-step pointer landed, the
  [#57](https://github.com/caseywdunn/corpus/issues/57)
  `--report` block + `<output_dir>/run.log` line did not).
- **Lexicon YAML loader bug.** A non-mapping-of-mappings lexicon
  silently degrades to no-op annotation with
  `WARNING ... 'list' object has no attribute 'get'`. Either
  tighten the schema check + emit an actionable error, or accept
  both shapes.
- **Figure-number extraction on old/scanned papers**
  ([#16](https://github.com/caseywdunn/corpus/issues/16)).
  ~538 of 1,787 papers have unparsed figure numbers; modest
  heuristic work.
- **Vision pass corpus-scale validation**
  ([#11](https://github.com/caseywdunn/corpus/issues/11)).
  Operational run + figure-coverage audit on Bouchet against the
  v0.4 build. Release-validation, not coding.

## 3. Target queries (evergreen reference)

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`.

| # | Pattern | Status after v0.4 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer ([#13](https://github.com/caseywdunn/corpus/issues/13), deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; synthesis recipe not yet scoped |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred ([#14](https://github.com/caseywdunn/corpus/issues/14), v0.5+ candidate) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Vision pass on by default since v0.3; full-corpus run pending [#11](https://github.com/caseywdunn/corpus/issues/11) |
| Q5 | "Summarize `<author X>`'s comments about `<author Y>`." | Indices in place |
| Q6 | "Summarize `<topic>` across the corpus." | Indices in place |
| Q7 | "Plot species described per decade." | Indices in place |
| Q8 | "Summarize what is known about `<anatomy>`." | Vision pass on by default since v0.3; full-corpus run pending [#11](https://github.com/caseywdunn/corpus/issues/11) |

## 4. Versioning + release ritual

`__version__` in [pipeline/version.py](../pipeline/version.py) is the
single source of truth and is stamped into every persistent artifact
(bundle manifest, MCP `bundle_info`).
[CONTRIBUTING.md](../CONTRIBUTING.md) covers the branching model and
release ritual.

**v0.4 release note:** non-breaking release; v0.3.x corpuscles work
as-is. The CHANGELOG `[Unreleased]` → `[0.4.0]` entry is in place;
the release ritual's two new CI gates (dev green before version
bump, main green before tag) are exercised first time on this
release.

## 5. Out of scope for v0.4

- [#14](https://github.com/caseywdunn/corpus/issues/14) — Trait
  extraction + identification keys (Q3). Substantial enough to
  warrant its own plan section when picked up; v0.5+ candidate.
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
