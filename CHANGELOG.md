# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`/healthz` liveness probe** on the SSE transport. Returns `200 ok`
  without requiring the bearer token — mounted ahead of the auth
  middleware in `mcpsrv.transport._HealthzASGI` so uptime monitors and
  reverse-proxy readiness checks don't need the shared secret. Exposed
  through `deploy/nginx.conf` as `location = /healthz`.

### Fixed

- **Stage 1 SLURM array tasks could process the same PDF concurrently
  and corrupt per-doc state files**
  ([#55](https://github.com/caseywdunn/corpus/issues/55)). The
  pre-resume filter at the top of `pipeline.main.main()` ran *before*
  batch slicing, so each array task's slice depended on disk state at
  task-start time. Tasks starting later saw fewer remaining hashes;
  their slice indices then landed on different members of the list,
  producing overlapping batches. The two writers raced on a shared
  `pipeline_state.json.tmp` filename, leaving either an interleaved-
  bytes payload (31 corrupted summaries in the production run that
  surfaced this) or a `FileNotFoundError` on the rename. Fixed by
  slicing on the unfiltered hash list (`_slice_hashes_for_batch` in
  `pipeline.main`); resume skipping is now done per-doc inside the
  loop. Belt-and-braces: `_save_pipeline_state` and
  `create_summary_json` now use per-writer tmp filenames
  (`.tmp.<pid>.<ns>`), so any future regression that re-introduces
  concurrent writers on the same hash gets last-write-wins atomicity
  instead of corruption. Regression test in
  `tests/test_batch_slicing.py`.
- **README §"Deploying MCP server remotely"** referenced an EC2 + ALB +
  CloudFront pattern that doesn't match the actual `deploy/stack.yaml`
  (single EC2 + nginx + Let's Encrypt). Corrected.
- **`environment.yaml` references to `tesseract-data-<code>` packages**
  ([#52](https://github.com/caseywdunn/corpus/issues/52)) — the 12
  per-language entries added in v0.2.0 (#46) listed package names that
  don't exist on conda-forge, so a fresh `conda env create` failed with
  `PackagesNotFoundError`. Existing dev envs survived because
  incremental updates silently skipped the missing packages. Replaced
  with a new `tools/install_tessdata.sh` helper that downloads the
  default fallback set (deu, fra, rus, lat, spa, por, chi_sim, chi_tra,
  jpn, ell, kor, grc, deu_latf) directly from `tessdata_best`.
  Idempotent; takes a custom language list as positional args; honors
  `TESSDATA_DIR` for non-conda installs.

## [0.2.0] - 2026-05-07

A hardening + iteration release. v0.2 closes out v0.1's deferred items
(vision pass at corpus scale, expanded lexicon + taxonomy support,
Bouchet batch-script fixes), adds a real update lifecycle (per-stage
resume, input fingerprints, `update_corpus.py` orchestrator), and lays
in a robustness + observability layer (structured `stage_failures[]`,
silent-failure quality gates, `corpus_status` rollup, network circuit
breakers). Geographic extraction
([#13](https://github.com/caseywdunn/corpus/issues/13)) and trait
extraction ([#14](https://github.com/caseywdunn/corpus/issues/14)) are
pushed to later releases.

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
- **`deu_latf` Fraktur Tesseract pack**
  ([#9](https://github.com/caseywdunn/corpus/issues/9)) — must be
  installed manually (see [`dev_docs/INSTALL.md`](dev_docs/INSTALL.md))
  before reprocessing 19th-c. German scans (Goldfuss 1820,
  Pagenstecher 1869, Brandt 1837, Dönitz 1871). Carried over from
  v0.1.
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
  [`dev_docs/DEPLOY.md`](dev_docs/DEPLOY.md). Two-bundle model: build bundle
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
  