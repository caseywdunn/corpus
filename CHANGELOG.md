# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-04-XX

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

- **OCR-step failure on ~33 docs** — ocrmypdf produces 0-text outputs on a
  subset of `force_ocr` and `skip_text` invocations. Diagnostic ticket with
  flag matrix in [#10]. Affects three real Cyrillic-mojibake papers
  (Stepanjants 1966 et al.) plus various very-old scans. Workaround
  pending.
- **Vendor-watermark misclassification** — three papers (Karplus2014,
  Browne1905, Fleming1828) classified `born_digital` because the only
  text layer is a per-page ProQuest/BHL banner. Tracked in [#8].
- **Fraktur pack not yet installed on Bouchet** — 19th-c. German scans
  (Goldfuss 1820, Pagenstecher 1869, Brandt 1837, Donitz 1871, …) currently
  OCR to whitespace. Tracked in [#9].
- **Geographic extraction (Q1 in PLAN.md §8)** — not yet implemented.
  