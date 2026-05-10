# PLAN.md — Corpus pipeline (v0.3)

v0.1 (2026-05-01) shipped the full extraction → annotation → indexing
→ MCP-serving stack. v0.2 (2026-05-08) hardened internals: per-stage
resume, input fingerprints, structured failure schema, quality gates,
and an ALB-fronted multi-organism deploy.

**v0.3 is about the user surface.** The pipeline today is internally
correct but operationally wide: thirteen CLI entry points at the repo
root, an opt-in `--resume` flag everywhere, manual orchestration of
cross-paper databases, taxonomy ingest separate from the main
pipeline, and per-corpuscle inputs (PDF dir, bib path, lexicon path,
taxonomy source) carried as long command lines. v0.3 collapses that
surface to **one CLI installed once, a `config.yaml` per corpuscle
directory, automatic database management, and capability-aware
defaults**. The framing of v0.2 was "we built the right machine";
v0.3 is "we make the machine usable by someone who isn't us."

This document is scoped to v0.3 work. Architectural background and
pipeline internals live in [OVERVIEW.md](OVERVIEW.md); per-feature
history in [CHANGELOG.md](../CHANGELOG.md); HPC operations in
[BOUCHET.md](BOUCHET.md); deployment in [DEPLOY.md](DEPLOY.md). Open
work is tracked in
[GitHub issues](https://github.com/caseywdunn/corpus/issues); the
v0.3-scoped subset is listed below.

## 1. v0.3 punch list

Tracked under the
[v0.3.0 milestone](https://github.com/caseywdunn/corpus/milestone/1).
Checkbox status below mirrors per-issue acceptance criteria; tick a
box here when its issue closes so the scannable progress lives in
one file.

### One CLI, code installed once

The unifying simplification: install `corpus` once
(`pip install -e .` from the source clone for development, or
`pip install git+https://github.com/caseywdunn/corpus.git@vX.Y.Z`
for deploys), then operate on any number of corpuscle directories.
A corpuscle is just a directory containing `config.yaml`; `corpus`
reads it from cwd by default and from `--config <path>` otherwise.
Multiple corpora = multiple corpuscle directories, all served by
the same installed binary — no duplicated source trees drifting
out of sync. Stage 1, Stage 2, post-pipeline cross-paper
databases, and bundle distillation all run from one entry point.

- [ ] **Single user-facing CLI: `corpus`.**
  ([#60](https://github.com/caseywdunn/corpus/issues/60))
  All operator interactions
  go through one binary, exposed as a setuptools entry point so
  the invocation is just `corpus <verb>`, not `python corpus.py`.
  Subcommand structure follows the cargo / git / gh / kubectl
  pattern — one tool, several verbs:
  - `corpus run` — the pipeline. Subsumes today's
    `update_corpus.py`, `process_corpus.py`, `embed_chunks.py`,
    `build_biblio_authority.py`, `build_taxon_mentions.py`,
    `backfill_intext_citations.py`,
    `reconcile_corpus_to_biblio.py`, `ingest_taxonomy.py`, and
    `package_for_serve.py`. Carries `--dry-run`,
    `--force-rebuild` / `--force-rebuild-<stage>`, `--no-vision`,
    `--no-bundle`, `--no-prune` / `--force-prune`, `--enrich-bhl`.
  - `corpus check` — environment + config pre-flight (detail
    below). Distinct from `corpus run --dry-run` (which prints the
    *plan*); answers "can the next run actually succeed on this
    host?" without consulting the output tree.
  - `corpus status` — coverage + error report against the built
    corpuscle (detail below). Carries the
    [#54](https://github.com/caseywdunn/corpus/issues/54)
    triage flags (`--sort-by`, `--tail`, `--propose-skips`,
    `--skipped`, `--re-process-flagged`), the
    [#57](https://github.com/caseywdunn/corpus/issues/57)
    `--report` view, and `--json` for scriptable consumption
    (CI / dashboards want structured; rich text is for humans).
  - `corpus serve` — replaces today's `mcp_server.py`. Reads the
    same `config.yaml` for bundle path, instructions, and
    bearer-token file.
  - `corpus bib export` / `corpus bib import` — replace today's
    `bib_export.py` / `bib_import.py`. Independently invocable;
    operators still run them between or before `corpus run`
    (including before the first run, see below).
  - `corpus init` — scaffolds `config.yaml` in cwd from the
    bundled template (which ships inside the installed `corpus`
    package, not in the corpuscle directory). The discoverable
    on-ramp for a new corpuscle: `mkdir my-corpus && cd my-corpus
    && corpus init && $EDITOR config.yaml`.

  Every current root-level Python script is either renamed,
  absorbed as a subcommand, or demoted to an internal module
  (kept importable for tests and ad-hoc debugging, dropped as a
  user-facing CLI). The repo root ends up with one operator
  binary and zero ambiguity about which script does what.
- [ ] **v0.2 → v0.3 is a clean break.**
  ([#60](https://github.com/caseywdunn/corpus/issues/60))
  No migration tool, no thin
  shims printing "use `corpus run` instead" — the old
  root-level CLIs (`process_corpus.py`, `update_corpus.py`,
  `embed_chunks.py`, `mcp_server.py`, etc.) are simply gone.
  Operators upgrading from v0.2 rebuild their corpuscles from
  scratch under the new `config.yaml`-driven flow:
  `corpus init && $EDITOR config.yaml && corpus run`. The v0.3
  CHANGELOG entry needs an explicit `### Breaking changes`
  subsection enumerating the dropped/renamed CLIs and the
  dropped `--resume` flag so operators see the change before
  they `git pull`.
- [x] **Packaging: `pyproject.toml` + `corpus` console_scripts entry
  point.**
  ([#58](https://github.com/caseywdunn/corpus/issues/58))
  A minimal `pyproject.toml` at the repo root
  (setuptools backend — boring, well-understood, no lock-file
  churn) makes the project pip-installable; `corpus` lands on
  PATH via `[project.scripts]` after `pip install -e .`
  (development) or
  `pip install git+https://github.com/caseywdunn/corpus.git@v0.3.0`
  (AWS deploy, replacing today's `git clone &&
  pip install -r requirements.txt`). Existing `pipeline/`,
  `mcpsrv/`, and `bib/` packages stay where they are. Out of
  scope for v0.3: PyPI publishing (the name `corpus` is taken;
  would need a different distribution name and a real
  release-artifact policy), conda-forge, and a `src/corpus/`
  namespace restructure (nice-to-have but not blocking — the
  global `pipeline` / `mcpsrv` / `bib` namespace pollution is
  annoying but contained). Defer all three until there's an
  audience that isn't cloning.
- [x] **Version handling: single source in `pipeline/version.py`,
  `pyproject.toml` reads it dynamically.**
  ([#58](https://github.com/caseywdunn/corpus/issues/58))
  Once `pyproject.toml`
  exists there are two version surfaces (Python imports +
  package metadata) that must not drift. Pattern: keep
  `pipeline/version.py` as the source of truth (existing
  imports + bundle-manifest stamping unchanged), and have
  `pyproject.toml` resolve the package version from it via
  `[tool.setuptools.dynamic] version = {attr =
  "pipeline.version.__version__"}`. Release ritual in
  [CONTRIBUTING.md](../CONTRIBUTING.md) §Releasing stays exactly
  as-is — bump the `.py` file, commit, tag, push — and `pip
  show corpus`, `corpus --version`, and the bundle manifest all
  agree without a second bump step. Rejected alternatives:
  reading the package version via `importlib.metadata` at
  runtime (breaks on uninstalled source trees, brittle for
  `python -m pytest` without `pip install -e .` first) and
  `setuptools-scm` (overkill for a one-developer-plus-
  collaborators project; produces gnarly `0.3.0.dev12+gabcd123`
  strings on untagged commits).
- [ ] **Global `--config` option, pre-verb (git-style).**
  ([#60](https://github.com/caseywdunn/corpus/issues/60); short
  alias `-c` and `CORPUS_CONFIG` env var on
  [#61](https://github.com/caseywdunn/corpus/issues/61))
  All
  subcommands resolve their config the same way: `corpus
  --config <path> <verb> [args]`. Default is `./config.yaml`
  relative to cwd, which is what the corpuscle-as-directory
  model wants ninety-five percent of the time; `--config` is the
  documented escape hatch for scripts that don't `cd` first,
  multi-config experiments, and the demo (which becomes plain
  `corpus --config demo/config.yaml run` rather than a special
  `--demo` flag). Discoverable via `corpus --help` rather than
  scattered across each subcommand's flag set.
- [ ] **Operator UX: shared `rich` console layer.**
  ([#63](https://github.com/caseywdunn/corpus/issues/63))
  A single
  `pipeline/console.py` Console instance backs every `corpus`
  subcommand. Emoji status symbols (✓/✗/⚠ with ASCII fallback),
  braille spinners on stages without a known total, and progress
  bars on stages where it is (per-paper Stage 1, embed batches,
  the finalize tail). `console.is_terminal` gates everything: a
  TTY operator gets the rich rendering, while SLURM `.out`
  capture, CI logs, and `journalctl` readers see clean ASCII —
  uniform call sites, no per-site `if tty:` branches. Pattern
  matches sharkmer's `indicatif` usage: one `show_progress` bool,
  hidden Progress when not a TTY. Critically, the diagnostic
  `logging` stream stays plain text — rich is the operator-facing
  layer (interactive run, success summary), not the structured
  event log that lands in journals and grep pipelines. The
  [#57](https://github.com/caseywdunn/corpus/issues/57) report is
  the most visible payoff: rich's `Table` + colored ✓/✗ + bars
  render the per-stage / lexicon-coverage / artifact-presence
  sections natively, and degrade to plain text in SLURM `.out`
  with no extra code. ~1.5 MB install dep; small relative to the
  torch + transformers footprint, and `rich` is already in the
  transitive dep tree via the `mcp` and `anthropic` SDKs.
- [ ] **Implicit resume on `corpus run`.**
  ([#60](https://github.com/caseywdunn/corpus/issues/60))
  Drop the `--resume` flag.
  The pipeline is always idempotent: re-runs do only the work
  whose inputs have changed (per-stage state, content
  fingerprints, PDF set diff — all infrastructure already in
  place from v0.2). Explicit `--force-rebuild` (and per-stage
  `--force-rebuild-<stage>`) covers the rare clean-rebuild case.
  Ctrl-C mid-run is part of this contract: per-stage state is
  flushed atomically as each stage completes (per #55), so the
  next `corpus run` picks up exactly where the previous left off.
- [ ] **Strict `corpus run --dry-run`.**
  ([#60](https://github.com/caseywdunn/corpus/issues/60))
  Prints the plan — what would
  run, what would skip, what would be deleted — without touching
  disk. Today's per-stage `--dry-run` (#41) honors this stage by
  stage; under `corpus run`, audit every write path so the
  dry-run guarantee holds across all stages, post-pipeline
  scripts, and the bundle distillation step.
- [ ] **`corpus check` detail.**
  ([#62](https://github.com/caseywdunn/corpus/issues/62))
  Validates `config.yaml` against a
  pydantic schema (field-level errors point at the exact key and
  value — `vision.backend must be one of {local, claude, none},
  got 'claud'`), probes GPU availability, pings the configured
  Grobid endpoint, checks `ANTHROPIC_API_KEY` presence when the
  vision backend is `claude`, sanity-checks output disk space,
  and prints a green/yellow/red report. The same pydantic schema
  is loaded by every other subcommand at startup, so config
  errors surface uniformly, not as a `KeyError` halfway through a
  Stage 1 run. Cheap to write once the config-driven CLI exists,
  and saves long debug sessions on Bouchet when something is
  misconfigured before SLURM kicks off. Doesn't read the output
  tree — distinct from `corpus status` (postmortem against a
  built corpuscle).
- [ ] **Success summary + next-step pointer.**
  ([#60](https://github.com/caseywdunn/corpus/issues/60); closes the
  report half of [#57](https://github.com/caseywdunn/corpus/issues/57))
  On clean completion,
  `corpus run` prints the same scannable report `corpus status
  --report` produces, plus the `corpus serve` command to launch a
  local MCP server against the new corpuscle. The report shape is
  the one specified in
  [#57](https://github.com/caseywdunn/corpus/issues/57): per-stage
  pass rate (with terminal bars), lexicon coverage per category
  with fingerprint match (denominator = all papers, breakdown
  into `ok` / `stale` / `missing`), cross-paper artifact presence
  (`taxonomy.sqlite`, `biblio_authority.sqlite`,
  `taxon_mentions.sqlite`, `vector_db/lancedb` ✓/✗), and a flat
  tally of `stage_failures.reason_code × stage` and
  `quality_flags.gate`. The reference script in #57 is the
  starting point. The same payload is also written to
  `<output_dir>/run.log` (one line of JSON per finished run plus
  a human-readable rendering) so operators reconstructing a run
  three weeks later don't have to scroll back through SLURM
  `.out` files or `journalctl`. Complements the existing
  per-paper `documents/<HASH>/pipeline.log` (detailed per-stage
  trace, overwritten on each re-run): `run.log` is the summary
  *across* runs, `pipeline.log` is the detail *within* one run.
  Closes [#57](https://github.com/caseywdunn/corpus/issues/57).
- [ ] **Bundle distillation in line.**
  ([#60](https://github.com/caseywdunn/corpus/issues/60))
  `corpus run` walks
  `package_for_serve.py` so a successful run produces a
  ready-to-ship served bundle alongside the build bundle.
  `--no-bundle` for users who only want the build artifacts.
- [ ] **`corpus bib export` / `corpus bib import` work before first
  `corpus run`.**
  ([#67](https://github.com/caseywdunn/corpus/issues/67))
  Today `bib_import` writes into
  `biblio_authority.sqlite`, which only exists after the pipeline
  has run. Under v0.3, `corpus bib import` either creates an
  empty authority DB and stages the edits, or persists the edits
  to a sidecar file the next `corpus run` picks up — the design
  call is part of the issue. Either way, the operator can curate
  the bibliography ahead of the first run, and the curation lands
  on the resulting authority DB without a second
  `corpus bib import` step.

### CLI conventions

The universal-CLI affordances every operator looks for first.
None of them are deep architectural decisions; their absence
just makes the tool feel rough. All tracked under
[#61](https://github.com/caseywdunn/corpus/issues/61).

- [ ] **`corpus --version` / `-V`** prints the package version (from
  `pipeline.version.__version__`), and the git SHA when the
  current install was sourced from a development clone. Matches
  every other CLI on the system.
- [ ] **`corpus --cite`** prints the citation. Default is a one-line
  plain-text rendering for the README / email use case;
  `corpus --cite=bibtex` emits a paste-ready BibTeX entry (so
  `corpus --cite=bibtex >> refs.bib` just works). Source of
  truth is a new `CITATION.cff` at the repo root — GitHub's
  standard citation-file format, which also gives the GitHub
  page a "Cite this repository" button automatically. The
  resolved citation is stamped into `bundle_manifest.json`
  alongside `bundle_version` / `pipeline_git_sha` and surfaced
  through the MCP `bundle_info` tool, so downstream LLM clients
  can cite the tool that gave them the literature without the
  user having to dig for it. Single source, three surfaces
  (CLI, GitHub, MCP), zero drift.
- [ ] **`corpus completion bash|zsh|fish`** generates a shell
  completion script the operator sources from their dotfiles.
  Matches `gh`, `kubectl`, `cargo`, `aws`, `docker`. Operators
  tab-complete subcommand names + `--config` paths constantly;
  absence is a sharp edge.
- [ ] **`CORPUS_CONFIG` env var** as an alternative to
  `--config PATH`. Matches `KUBECONFIG`, `AWS_PROFILE`,
  `RUSTUP_HOME`. Useful in CI / SLURM scripts where exporting
  once at the top is cleaner than threading the flag through
  every invocation. Precedence: `--config` flag > env var >
  `./config.yaml` in cwd.
- [ ] **`-c` short alias for `--config`.** Small, but a real
  ergonomic ask once an operator types it ten times a day.
- [ ] **Defined exit codes.** `0` success, `1` generic failure, `2`
  config error (missing/invalid `config.yaml`, schema mismatch,
  unresolvable input path), `3` precondition not met (`corpus
  check` fail, missing bundle, port in use). Lets CI and SLURM
  `--dependency=afterok:` chains branch on the *kind* of failure
  rather than just success/failure.
- [ ] **Verbosity levels.** `-v` / `-vv` raise the log threshold
  (INFO → DEBUG → all child loggers chatty); `-q` / `--quiet`
  drops it to WARNING. Replaces today's boolean `--verbose`.
  SLURM `.out` capture often wants `-q`; debugging operators
  want `-vv`.
- [ ] **Path resolution.** Relative paths in `config.yaml`
  (`input_pdfs:`, `bib:`, `lexicon:`, `output_dir:`,
  `taxonomy.path:`) resolve against the directory containing the
  config file, not cwd. Means `cd demo && corpus run` and
  `corpus --config demo/config.yaml run` from anywhere give
  identical results — the demo's `input_pdfs: ./*.pdf` always
  finds `demo/<files>.pdf`. cwd matters only for "where does
  `corpus` find `config.yaml` when `--config` and
  `CORPUS_CONFIG` are both absent."

### Per-corpuscle `config.yaml`

The current `config.yaml` mixes system-wide tuning (OCR language
packs, quality-gate thresholds) with no corpuscle-specific surface;
all per-corpuscle inputs are CLI flags. v0.3 turns `config.yaml`
into the single source of truth for *this* corpuscle's inputs and
settings — one file per corpuscle directory, separate from the
installed code. All tracked under
[#59](https://github.com/caseywdunn/corpus/issues/59).

- [ ] **`config.template.yaml` ships inside the installed `corpus`
  package; `config.yaml` is per-corpuscle.** First-run check: if
  `corpus run` finds no `config.yaml` in cwd (and `--config`
  isn't set), it errors with `no config.yaml in cwd; run "corpus
  init" to scaffold one, or pass --config PATH`. `corpus init`
  copies the bundled template into cwd, after which the operator
  edits the input paths. No silent auto-create — operators should
  make a deliberate choice about where their PDFs / bib / lexicon
  live. (`config.yaml` may or may not be tracked in git: the
  corpuscle directory might be its own repo, in which case
  committing `config.yaml` is the operator's choice, not the
  tool's.)
- [ ] **Move per-corpuscle inputs into config.** Today's CLI flags
  fold in:
  - `input_pdfs:` (replaces the positional input arg)
  - `output_dir:` (defaults to `./output`)
  - `bib:` (replaces `--bib`)
  - `lexicon:` (replaces `--lexicon`)
  - `taxonomy:` block — `source: worms | dwc | dwca`,
    `root_id:` (for WoRMS), `path:` (for local DwC archives)
  - `vision:` block — `backend: local | claude | none`,
    `model:` overrides
  - `grobid:` block — `url:`, optional `disable: true`
  - `bibliography:` block — `enrich_bhl: false` (default; opt-in
    to BHL enrichment per the auto-build bullet below)
  - `licensing:` block — `pd_cutoff_years: 95` (default;
    US-specific public-domain cutoff per
    [#51](https://github.com/caseywdunn/corpus/issues/51));
    other policy fields land as the #51 design lands
- [ ] **Drift detection.** Hash the resolved config (input paths +
  per-input content SHA) into a corpuscle-side state file. A
  mismatch on the next run logs which keys drifted and which stages
  it invalidates, so an operator can see *why* a re-run is doing
  more than they expected.
- [ ] **Demo is a regular corpuscle.** `demo/` is just a corpuscle
  directory that happens to be tracked in the corpus source repo:
  it has its own `config.yaml`, `lexicon.yaml`,
  `siphonophores.bib`, `instructions.md`, and the 11 PDFs.
  Running the demo is identical to running any other corpuscle:
  `cd demo && corpus run`. (Or `corpus --config demo/config.yaml
  run` from anywhere, but `cd` is the documented path so the
  demo experience matches the real-corpus experience.)
  `demo/output/` lands in `.gitignore` so the source repo stays
  clean. Doubles as the smoke-test fixture (see Tests bullet
  below).
- [ ] **Tests pin the demo corpuscle as fixture.** `CORPUS_OUTPUT_DIR`
  env var still works as the explicit override; default fixture
  resolves to `demo/output/` after a clean `corpus run` against
  `demo/`. The fixture-fallback logic in `tests/conftest.py` is
  already close; tighten the precedence so a missing env var
  resolves deterministically to the demo output rather than the
  legacy `output/` ambiguity.

### Auto-build, auto-detect, auto-clean

Today's pipeline fails open: configure a taxonomy source but skip
`ingest_taxonomy.py` and the resulting MCP server returns "DB not
loaded" for every taxonomy query, with no easy way to find out
why. v0.3 closes that gap by making the unified entry point
inspect the config and act.

- [ ] **Auto-build cross-paper databases.**
  ([#64](https://github.com/caseywdunn/corpus/issues/64))
  If `taxonomy:` is configured
  in `config.yaml` but `taxonomy.sqlite` doesn't exist, build it.
  Same for `biblio_authority.sqlite` and `taxon_mentions.sqlite`.
  Hash the input source for each (bib SHA-256, taxonomy snapshot
  SHA-256, lexicon SHA-256) into the corpuscle-side state; rebuild
  the database when the input hash changes. `--force-rebuild` and
  per-database `--force-rebuild-<db>` flags cover the ops escape.
  `bibliography:` block in `config.yaml` carries `enrich_bhl: false`
  by default — slow, rate-limited, network-dependent, opt-in via
  config or `--enrich-bhl` CLI override.
- [ ] **Vision pass: opt-out, capability-aware.**
  ([#65](https://github.com/caseywdunn/corpus/issues/65))
  Run the vision pass
  by default. If the configured backend isn't usable on this host
  (no GPU for `local`, no `ANTHROPIC_API_KEY` for `claude`), skip
  with a clear log message and a one-line nudge for what to do
  about it. `--no-vision` for explicit skip. Closes the long tail
  of [#11](https://github.com/caseywdunn/corpus/issues/11) for the
  default-config path.
- [ ] **Orphan cleanup is the default.**
  ([#66](https://github.com/caseywdunn/corpus/issues/66); extends
  [#31](https://github.com/caseywdunn/corpus/issues/31))
  If a PDF disappears from the
  configured input dir, its `documents/<HASH>/` directory and
  matching LanceDB rows are deleted on the next run. Today's
  `--audit-orphans` (#31) is read-only; v0.3 makes deletion the
  default, with `--no-prune` for the read-only mode and a safety
  rail (refuse to prune if more than N% of documents would be
  affected — likely a config mistake or unmounted volume) gated
  behind `--force-prune`.

### Long-run resilience

A real corpus runs for hours. v0.2's `update_corpus.py` chains the
post-pipeline cross-paper builds as a foreground subprocess loop —
fine for the demo, fragile when a laptop closes or an SSH session
drops mid-run. v0.3's unified CLI keeps the same shape; the
fragility is the same. The fix is operational, not structural: a
SLURM analogue that chains the cross-paper tail, and a documented
detached-run path for laptop runs.

- [ ] **`slurm/batch_finalize.sh` for the cross-paper tail.**
  ([#57](https://github.com/caseywdunn/corpus/issues/57))
  Today's
  `slurm/batch_pipeline.sh` chains Grobid → Stage 1 → Pass 3b +
  Embed and stops; the four cross-paper builds run on the login
  node manually after the array completes. Under v0.3 the
  finalize phase of `corpus run` lands in `batch_finalize.sh`,
  chained `--dependency=afterok:` on the embed job. Single CPU
  job, no GPU. Mirrors the existing `slurm/batch_pipeline.sh`
  ergonomics and accepts the same opt-in env vars
  (`ENRICH_BHL=1`, `SERVE_BUNDLE_DIR=...`).
- [ ] **Documented detached-run recipe for laptops.**
  ([#57](https://github.com/caseywdunn/corpus/issues/57))
  A `nohup` / `tmux` recipe in the README under "First time run"
  so detaching is the documented path, not a workaround.

### Curation surface — what ships, what's reusable

`package_for_serve.py` whitelist-copies every paper with a complete
`summary.json` into the served bundle, with no operator control over
exclusions and no per-figure rights metadata for downstream reuse.
v0.3 closes both gaps so operators can shape what's exposed and
clients (especially manuscript-authoring agents) can know what's
safe to embed.

- [ ] [#54](https://github.com/caseywdunn/corpus/issues/54) — **PDF QC
  workflow: skip flag, content-vs-technical gates, worst-first
  triage.** A `serve = false` BibTeX field (curated via the #26
  round-trip) propagates to `works.serve` in
  `biblio_authority.sqlite`; `package_for_serve.py` honors it as
  the only deploy-time exclusion gate. The pipeline runs
  unchanged so excluded papers stay in `corpus status` rollups
  for comparison against new candidates. New `corpus status`
  flags — `--sort-by <metric> --tail N`, `--propose-skips`,
  `--skipped`, `--re-process-flagged <gate>` — give operators a
  worst-first triage list. Adds a small set of QC metrics
  (language coverage, citation resolution rate, page-text
  histogram, dictionary hit rate, figure-to-page ratio,
  head/foot pollution score, Grobid header confidence) with
  matching `quality_flags` gates. Implementation order in the
  issue.
- [ ] [#51](https://github.com/caseywdunn/corpus/issues/51) — **Figure
  licensing + publishable gate.** `license` / `licenseurl` BibTeX
  fields (SPDX short identifiers + a small custom vocabulary:
  `public-domain`, `all-rights-reserved`,
  `publisher-permission`, `unknown`) propagate through the
  authority DB into per-figure metadata. A `publishable` boolean
  is derived at build time from license + age cutoff (default
  `pd_cutoff_years=95`, US-specific, configurable per corpuscle).
  `get_figure_image` refuses to return image bytes when
  `publishable=false`, with the reason in an error field;
  `--allow-unpublishable` MCP server flag covers local-only
  rights-holder cases. Both raw fields and the derived flag are
  exposed so clients with non-default jurisdictions / use cases
  can re-derive. Adds a `dev_docs/LICENSING.md` documenting the
  policy and known limitations (reprint chains, museum-photo
  copyright, etc.). Implementation order in the issue;
  orthogonal to #54 (skip gate is "appears in MCP at all";
  publishable gate is "image is reusable in derived
  publication").

### Documentation surface

The README drifts long because operator material (remote deploy,
on-host setup) lives next to user-onboarding material (install,
demo walkthrough). v0.3 splits them, and rolls the surface
collapse from §1 ("One CLI, code installed once") through into
the docs. All tracked under
[#68](https://github.com/caseywdunn/corpus/issues/68).

- [ ] **DEPLOY.md → repo root.** Sits alongside INSTALL.md (already
  at root post-v0.2) and CONTRIBUTING.md as the operator-facing
  top-level docs. `dev_docs/` keeps maintainer-only docs (PLAN,
  BOUCHET, MCP_TOOLS, TESTING, OVERVIEW).
- [ ] **Migrate "Deploying MCP server remotely" out of README.md.**
  The ~60-line block at README §"Deploying MCP server remotely"
  (bearer-token generation, SSE startup, smoke test, three-way
  client-config matrix) belongs in `DEPLOY.md` next to the AWS
  runbook. README keeps a short pointer + "Deploying MCP server
  locally" (which is the path most users actually take).
- [ ] **README walkthrough slim-down** post-unified-CLI. The "First
  time run" / "Adding and updating documents" / "Vision pass"
  sections each collapse to 1-3 lines once `corpus run` is the
  only entry point and resume is implicit. Rewrite them once the
  CLI work lands so the README stops describing scripts that no
  longer exist.

### Operator clarity

- [ ] [#47](https://github.com/caseywdunn/corpus/issues/47) — **MCP
  server failures should be self-diagnosing.** Today's startup
  errors send operators to the systemd journal + nginx logs.
  Surface the most common failure modes (missing bundle, stale
  token, port already in use, missing cross-paper DB, wrong
  bundle layout) with actionable messages at startup, and add
  `corpus serve --check` for a pre-flight that exits with a
  diagnosis without binding the port. (Distinct from the
  pipeline-side `corpus check`: this one validates the
  serve-time surface — bundle present, token readable, port
  free, cross-paper DBs loadable.)

### Quick wins

- [ ] [#50](https://github.com/caseywdunn/corpus/issues/50) — `pandoc`
  in `environment.yaml`.

### Carryover

- [ ] [#11](https://github.com/caseywdunn/corpus/issues/11) — Vision
  pass at corpus scale. Mostly addressed by the opt-out + capability
  detection above ([#65](https://github.com/caseywdunn/corpus/issues/65));
  remaining work is verifying a full-corpus run
  end-to-end on Bouchet after the unified CLI lands.
- [ ] [#16](https://github.com/caseywdunn/corpus/issues/16) — Figure-
  number extraction on old/scanned papers. Carried from v0.2;
  ~538 of 1,787 papers have unparsed figure numbers. Modest
  heuristic work, fits the v0.3 cycle if there's room.

## 2. Target queries (evergreen reference)

The eight target query patterns the corpus is designed to serve.
Generic shapes; concrete instantiations live in the corpuscle's
`instructions.md`. v0.2 made Q4 + Q8 fully answerable once the
vision-tagged figures land at corpus scale (still pending — see
[#11](https://github.com/caseywdunn/corpus/issues/11) under
Carryover); Q3 remains gated on trait extraction
([#14](https://github.com/caseywdunn/corpus/issues/14), v0.4
candidate).

| # | Pattern | Status after v0.2 |
| --- | --- | --- |
| Q1 | "List all collection locations of `<species>`." | Partial — needs geographic mention layer (#13, deferred to v2.0+) |
| Q2 | "Compose a monographic review of `<genus>`." | Indices in place; synthesis recipe not yet scoped |
| Q3 | "Make a key to identify species in `<genus>`." | Trait extraction deferred (#14, v0.4 candidate) |
| Q4 | "List all valid species + one-paragraph summary + diagnostic figures." | Vision pass in default config in v0.3; full-corpus run pending #11 |
| Q5 | "Summarize `<author X>`'s comments about `<author Y>`." | Indices in place |
| Q6 | "Summarize `<topic>` across the corpus." | Indices in place |
| Q7 | "Plot species described per decade." | Indices in place |
| Q8 | "Summarize what is known about `<anatomy>`." | Vision pass in default config in v0.3; full-corpus run pending #11 |

## 3. Versioning + release ritual

`__version__` in [pipeline/version.py](../pipeline/version.py) is the
single source of truth and is stamped into every persistent artifact
(bundle manifest, MCP `bundle_info`).
[CONTRIBUTING.md](../CONTRIBUTING.md) covers the branching model and
release ritual; the short version: `dev` carries a PEP 440 pre-release
suffix (`0.3.0.dev0` → `0.3.0a1` → `0.3.0`), the release commit drops
the suffix, the next commit on `dev` reintroduces one for the next
target.

**v0.3 release note:** this cycle is a breaking-change release
(see §1 "v0.2 → v0.3 is a clean break" bullet). The CHANGELOG
entry must lead with an explicit `### Breaking changes`
subsection enumerating the dropped/renamed CLIs and the dropped
`--resume` flag; operators upgrading from v0.2 must rebuild
their corpuscles from scratch.

## 4. Out of scope for v0.3

- [#14](https://github.com/caseywdunn/corpus/issues/14) — Trait
  extraction + identification keys (Q3). Substantial enough to
  warrant its own plan section when picked up; v0.4 candidate.
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
