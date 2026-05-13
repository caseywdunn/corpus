# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
  