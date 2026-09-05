# Architecture Overview

The corpus pipeline processes a collection of scientific literature PDFs into a searchable, queryable knowledge base. It handles born-digital papers, scanned historical plates, multilingual text, and 19th-century Fraktur typefaces.

## Repo layout

The top-level entry-point scripts are thin shims; the implementation is grouped into packages:

| Package | Role |
|---|---|
| `pipeline/` | Stage 1 + Pass 3b/3c orchestrator and shared library modules. Split into `scan.py` (OCR), `extract.py` (docling), `metadata.py` (Grobid + bib), `chunking.py`, `annotate.py` (taxa + lexicons), `figure_passes.py`, `runner.py` (per-paper orchestrator), `main.py` (CLI), and supporting `config.py` / `io.py` / `log.py` / `stages.py`. Shared library modules: `figures.py`, `taxa.py`, `grobid_client.py`, `embeddings.py`, `vision.py`, `external.py`, `version.py`. |
| `mcpsrv/` | MCP server. `app.py` defines the MCPServer instance; `tools/{papers,taxonomy,bibliography,figures,chunks,lexicon,profiles}.py` register the `@mcp.tool()`-decorated functions (see [MCP_TOOLS.md](MCP_TOOLS.md) for the catalog); `transport.py` handles stdio + SSE; `indexes.py` is the eager in-memory index. |
| `bib/` | BibTeX import / export round-trip plus shared metadata helpers (`parser.py`, `importer.py`, `export.py`). |
| `slurm/` | SLURM batch scripts (Bouchet). |
| `deploy/` | CloudFormation, nginx config, systemd unit, sync + update shell scripts. |
| `tests/` | Ground-truth + corpus-wide consistency tests. |

There are no top-level entry-point scripts: `corpus` (from `[project.scripts]`) is the single CLI, and it dispatches into the packages above. The root-level scripts this doc used to list — `process_corpus.py`, `mcp_server.py`, and friends — were folded into packages in v0.3 (#60).

## Execution planes and data ownership

An **execution plane** says where a kind of decision belongs. It is broader
than a pipeline stage: stages are resumable steps inside one build, while the
planes span work before the build, the build itself, the running service and
the client consuming that service.

| Plane | Owns | May compute | Must not do |
|---|---|---|---|
| **Library curation** | Source PDFs, BibTeX and curator directives such as `ocrlang` / `keeppages` | Inspect sources, validate metadata and propose or review source edits | Write derived corpus artifacts; carry private copies of generic PDF/OCR rules |
| **Build/materialization** | OCR, extracted text and figures, associations, databases, embeddings and the served bundle | Expensive batch/GPU/API work and deterministic, resumable transforms | Depend on `tools/` or `skills/`; publish ambiguous partial state as complete |
| **Serve/query** | An immutable served bundle plus disposable caches | Bounded lookup, filtering, authorization, formatting and compatible query embedding | OCR, reconciliation, corpus-wide mutation, external enrichment or general LLM calls |
| **Client/agent** | User intent, workflow state and deliverables | Synthesis, translation, orchestration and presentation | Become the enforcement point for licensing, provenance or access control |

The normal data flow is one way:

```text
library -> build/materialization -> immutable bundle -> bounded response -> client output
```

Feedback travels back as an explicit, reviewed edit to the library or build
configuration; a query does not mutate the corpus. Generic knowledge about
PDFs, OCR or artifact formats belongs in `pipeline/`, even when first needed
by a library-curation workflow. Read-only inspection commands may live in
`tools/`, and collection-specific judgments belong in the library itself or
in a skill. `pipeline/` therefore never imports from `tools/` or `skills/`.

Thin does not mean computation-free. Semantic search must embed the active
query with the same versioned model as the stored vectors. That operation is
lazy and bounded to one request; it does not alter the corpus. Authorization
and provenance checks also remain at the server boundary because a client
nudge is not enforcement.

Panel crops now use a process-private disposable cache outside the bundle
(#275), capped at 128 entries and 128 MiB. Keys include source-image bytes and
ROI coordinates; replacement pixels cannot reuse an old crop. Cache eviction
is safe because image responses hold their own bytes, and HTTP panel requests
can regenerate crops directly. The temporary cache follows the host's `TMPDIR`
and is removed at normal process exit. No figure request writes to the bundle.
Remote figure URLs use five-minute, process-local capabilities scoped to one
figure, optional panel and output profile (#276). They never expose the MCP
bearer credential or authorize MCP routes. The download boundary rechecks
licensing. Reverse proxies use an explicitly configured public base, never
untrusted forwarded headers; deployment details are in [DEPLOY.md](../DEPLOY.md).

Every explicit MCP list selector has the same budget (#277): at most 500 items,
4,096 characters per item and 65,536 characters per list. Tools reject excess
with `invalid_argument` in their existing error shape **before** consulting the
corpus. They never truncate selectors silently. Omitted/empty selectors retain
their documented semantics; this input budget is not a new output-pagination
contract. Existing implicit whole-paper/all-corpus queries and aggregate
response-size limits need separate review, not an unnoticed wire change.
The freeze tests pin all tool signatures/defaults, SDK input schemas and
representative response shapes alongside the licensing and error invariants.

## Content-addressed storage

Every unique PDF is identified by the first 12 hex characters of its SHA-256 hash. All artifacts for that PDF live under `<output_dir>/documents/<HASH>/`. Identical PDFs found at multiple input paths are processed once; every discovered path is recorded in `summary.json`. Implicit resume checks per-stage records in `pipeline_state.json`, including producer version and input fingerprints; the presence of `summary.json` alone does not prove completion. Embedding has its own committed-row evidence, described below.

```
output/
  documents/
    a1b2c3d4e5f6/
      processed.pdf            # OCR'd copy (if scanned)
      docling_doc.json         # full docling parse
      scan_detection.json      # born_digital | scanned | broken_text_layer
      text.json                # extracted text (markdown)
      figures.json             # figure metadata, captions, classifications
      figures/                 # extracted figure PNGs
      page_report.html         # optional on-demand page/caption audit; build-only
      metadata.json            # Grobid bibliographic metadata + references
      chunks.json              # text chunks for embedding
      taxa.json                # taxon mentions (incl. input fingerprint)
      <category>.json          # one per category in --lexicon (anatomy.json,
                               # biogeography.json, …)
      intext_citations.json    # in-text <ref type="bibr"> resolved to work_ids
      pipeline.log             # per-document processing log
      summary.json             # provenance + per-stage status + quality flags
    ...
  vector_db/
    lancedb/                   # LanceDB vector index
    <HASH>_embedded.done       # per-hash embedding completion markers
  taxonomy.sqlite              # Darwin Core snapshot
  biblio_authority.sqlite      # deduplicated works + citation graph
  taxon_mentions.sqlite        # cross-paper taxon index
```

## Pipeline stages

The pipeline is split into independent stages so that CPU-heavy work, GPU work, and API-cost work can run on different hardware and be resumed independently. Three terms recur in this doc and in the code — they are *not* synonyms:

- **Stage** — a top-level, independently-resumable phase the orchestrator and SLURM scripts schedule, often on distinct hardware: **Stage 1** (CPU per-paper processing), the optional **vision pass** (GPU/API), and **Stage 2** (embedding, GPU).
- **Step** — one operation in the ordered sequence *inside* a stage. Stage 1 has the six steps tabled below; the figure lifecycle has eleven (see [Figure pipeline](#figure-pipeline)).
- **Pass** — a *figure-pipeline-internal* numbering (Pass 2.5, 3a, 3b, 3c) that predates the stage scheme and survives in code field names (`pass3_status`, `pass3_backend`) and the `slurm/batch_pass3b.sh` job. Passes 2.5 / 3a / 3c are CPU work invoked from the Stage 1 runner; only **Pass 3b** (vision) needs a GPU/API, which is why it is the one pass that appears as its own subsection below and can be scheduled separately. All four passes are detailed in the Figure pipeline section.

### Stage 1: CPU processing (`pipeline/`)

Orchestrated by `pipeline.runner.run_pdf_processing_pipeline`, six steps per PDF:

| Step | What happens | Output |
|---|---|---|
| **Scan detection** (`pipeline/scan.py`) | **Page geometry first** (v1.0): the fraction of sampled pages carrying a single full-page image decides whether a document is a scan, independent of what its text layer claims — a scan bundled with third-party OCR is still a scan. Then the text-layer heuristics for everything else: volume gate, langdetect, gibberish-score, Tesseract OSD visual-script cross-check. A non-Latin `ocrlang` pin against a Latin-only layer triggers that cross-check regardless of gibberish score. Classifies as `born_digital`, `scanned`, or `broken_text_layer`; per-document `ocrmode` can force and select the OCR route while preserving the detected verdict for audit. | `scan_detection.json` |
| **PDF preparation** (`pipeline/scan.py`) | Copies born-digital PDFs as-is unless `ocrmode` explicitly overrides detection. Everything else is OCR'd — `--force-ocr` for a uniform scan or a rejected/no-content layer, `--redo-ocr` for a mixed volume so a bound-in digital typescript isn't rasterized. Language comes from OCRing a 5-page sample rather than the (distrusted) text layer, per page and unioned so bilingual originals-plus-translations get both packs; failing that, the fallback union (eng/deu/fra/rus/lat/spa/por/chi_sim/chi_tra/jpn/ell/kor + deu_latf Fraktur, via `ocr.ocr_languages_default`). Per-page and per-document timeouts both scale — see `ocr.tesseract_page_timeout` and `stage_timeouts.ocr_per_page`. A zero exit with every output page textless is persisted as `ocr_no_text_recovered`, not accepted as successful transcription. | `processed.pdf` |
| **Text + figure extraction** (`pipeline/extract.py`) | Docling parses the PDF into structured text and figure regions. Figures go through a classification/caption pipeline (see [Figure pipeline](#figure-pipeline) below). Falls back to raw PyMuPDF image extraction when docling finds nothing. | `text.json`, `figures.json`, `figures/*.png` |
| **Metadata extraction** (`pipeline/metadata.py`, `bib/`) | Grobid extracts title, authors, year, DOI, abstract, section structure, and parsed references. `--bib` overrides the header from a curated BibTeX. Falls back to placeholder when Grobid is unavailable. | `metadata.json` |
| **Chunking** (`pipeline/chunking.py`) | Splits extracted text via docling's `HybridChunker` (tokenizer-aware, respects section/heading structure), with section-class labels. | `chunks.json` |
| **Annotation** (`pipeline/annotate.py`) | Per-chunk taxon mentions (against the DwC taxonomy snapshot) and lexicon matches (one pass per category in `--lexicon`). Each output file stamps a per-category `input_fingerprint`; the stage-completion record in `pipeline_state.json` mirrors it so a per-category resume detects exactly which categories changed. | `taxa.json`, `<category>.json` |

Stage 1 supports SLURM job-array parallelization via `--batch-index` / `--batch-size`. Each array task deterministically processes a slice of the sorted hash list. See [BOUCHET.md](BOUCHET.md) for operational details.

### Vision pass (optional, GPU/API): figure panels + compound figures

This is **Pass 3b** of the figure pipeline — the only figure pass that needs a GPU or a paid API, which is why it is called out at stage level here. It runs a vision-language model (Qwen2.5-VL-7B-Instruct locally, or Claude via API; selected by `--figure-panels {vision-local|vision-claude}`) over extracted figures to detect multi-panel structure, compound plates spanning several figure numbers, and per-panel ROI boxes, writing the results back into `figures.json`. **Pass 3c** (compound-figure resolution + range-notation renaming, e.g. `fig_3-4.png`) is a CPU follow-on that runs automatically in the same invocation whenever Pass 3b — or the OCR-based Pass 3a — flagged a compound.

Both vision backends size their output-token budget from the caption-declared
panel count, treating the configured token value as a floor. If a response
reaches that cap, the backend retries that figure once at twice the budget.
Completion is a data-integrity boundary: a second capped response, one that
cannot be parsed, or one that omits either required output list is recorded as
`vision_backend_failed` with `pass3_error`; it is never treated as an empty
detection. Only an explicit complete response with empty `panels` and
`embedded_figures` lists becomes `no_labels_found` (#269).

By default both run inline inside the Stage 1 per-paper runner; on HPC the GPU-bound Pass 3b is commonly split into its own scheduled job (`slurm/batch_pass3b.sh`) after the CPU stage completes. See [Figure pipeline](#figure-pipeline) for the full pass-by-pass behavior and what is lost when the vision pass is skipped.

### Stage 2: Embedding (`pipeline/embed.py`)

Reads `chunks.json` per hash and produces vector embeddings stored in LanceDB.

- **Backend**: local sentence-transformers with [BGE-M3](https://huggingface.co/BAAI/bge-m3) (1024-dim, multilingual), runs on CUDA/MPS/CPU via `embeddings.detect_device()`. Embedding failures raise rather than silently inserting zero vectors.
- Each document is replaced in **one LanceDB transaction**. A fresh internal
  `embedding_generation` labels the incoming rows; a merge inserts that
  generation and deletes every previous row for the same PDF hash. Shorter
  chunk lists and legacy duplicates are handled without touching other papers.
  An explicitly empty chunk list removes old rows. Missing or malformed chunk
  artifacts fail instead of being mistaken for an empty document.
- Only Stage 2 writes `<HASH>_embedded.done`, after the transaction succeeds,
  using an atomic file replacement. The marker records producer/schema
  versions, input fingerprint, model, dimension, table, row count, generation
  and embedding-space producer identity (repository revision or local bytes).
  The fingerprint covers the ordered chunk text and every metadata field stored
  in the index: paths, bibliographic fields, page count, chunk IDs, section and
  headings. Run timestamps and unrelated annotation fields do not invalidate it.
- `--resume` checks the input fingerprint **and** a projected census of committed
  row counts/generations. Matching model/dimension alone is insufficient. A
  crash after the table commit but before the marker replacement leaves evidence
  that fails this check, even if the old and new chunk counts are identical.
  Retry safely replaces the document again; it never appends a duplicate copy.
- `--dry-run --resume` uses the same evidence without loading the model or
  writing rows. `--rebuild` drops and re-embeds the **whole** table, including
  when two different models have the same dimension; combining it with
  `--pdf-hash` is rejected. Legacy markers lack proof and incur a one-time
  re-embedding on the next resume.
- Bundling validates every document's marker against current inputs and rows,
  rejects orphan rows and mixed models, and derives model/dimension from the
  validated set before touching the served destination. `corpus status --report`
  uses the same check instead of treating an index directory as completion.
  A build with neither an index nor current-document markers can still be served
  without semantic search; an incomplete index is not that mode.

The portable `embedding_producer.json` sidecar in a new served bundle carries
the identity verified across all document receipts. Query embedding stays lazy
but loads the build's model and pins its recorded repository revision. A
same-dimension model override is not compatibility: mismatched IDs or loaded
producers disable semantic search with an explicit degraded status. Custom
local models require `--embedding-model` on the serving host and matching
content hashes; private build paths are not shipped. Legacy bundles still use
their manifest model with a weaker-proof warning when revision evidence is
absent. Build-directory serving checks every document's available identity
receipt at startup; it does not sample one marker on each query.

Changing the embedding producer requires a whole-table `--rebuild`, even at
the same dimension. Legacy producer receipts migrate only when all current
documents are included, never through a partial update that could mix spaces.
The identity records the normalized SentenceTransformer recipe; device/package
floating-point reproducibility is measured by exact vector comparisons, not
promised by a model-name match. Custom embedding backends should expose their
own `producer`; otherwise the receipt explicitly has declared-model-ID proof.

### Corpuscle update contract

The build plane owns update decisions (#265). A PDF hash identifies source
bytes, **not** the validity of all derived data. The target invariant is that
an incremental run and a clean rebuild yield the same current document set,
logical vector rows and reference mappings. Internal transaction IDs and run
timestamps are provenance, not semantic differences.

| Change | Required invalidation / replacement | Current coverage |
|---|---|---|
| Add a PDF | Build the new hash; refresh cross-paper materializations and bundle | Existing implicit-add path; Stage 2 preserves other hashes |
| Remove a PDF | Retire its derived directory and vector rows; remove current cross-paper and served references | Strict source inventory, recoverable retirement, current bibliography/taxon-index pruning and fresh bundle replacement; real-vector clean/incremental deletion regression tested |
| Change PDF bytes | Treat as an addition plus removal of the old hash if no other source path retains it | Same rules as above |
| Same hash, changed chunks or embedded metadata | Replace that document's vector rows; publish a new receipt after commit | Implemented and tested in Stage 2 (#271) |
| Rename or change BibTeX | Refresh consumed paths/metadata and their descendants; leave unrelated artifacts alone | Resolved per-paper entries (including absence) and basenames fingerprint metadata; source-path inventory refreshes even when extraction skips; Stage 2 notices the changed metadata/paths (#174) |
| Change configuration or curator directives | Fingerprint resolved values in their consumers and descendants; replace subordinate caches too | Stage 1 settings below, OCR directives and annotation inputs are tracked; broader provenance/whole-build acceptance remains open (#174) |
| Upgrade producer or embedding model | Invalidate incompatible receipts; rebuild the whole vector table when switching models | Stage/embedding producer checks and model guard exist; whole-corpus release comparison remains pending (#187) |

These are **single-writer build directories**. Parallel extraction workers may
own disjoint hashes, but do not overlap embedding/bundling with edits to their
inputs, another embedding writer, or an update of the producer checkout. Use a
separate build and served directory for release candidates. Per-document vector
transactions do not make the complete corpus build atomic. Bundling assembles
and audits a fresh sibling directory before replacing the destination; a
failed copy leaves the previous destination untouched. Replacement is an
offline operation, not a live-server hot-swap protocol: stop serving during
publication and restart afterward. Rebuilding a bundle copies its complete
served content; it does not reuse files based on timestamps.

Source retirement requires a complete, readable inventory and excludes the
build subtree even when the source root contains it. `--force-prune` bypasses
the percentage safety rail, **not** unreadable/missing-source errors or failed
vector cleanup. Retired document artifacts and their embedding receipts remain
under `<output_dir>/.retired/documents-*`; historical reference observations
remain in the bibliography DB but no longer form current edges. Taxon-index
receipts compare content digests and retire missing document members, never
infer unchanged evidence from timestamps. A malformed replacement fails the
post phase without deleting the previously indexed evidence.

Previous bundles are retained beside the destination as
`.<bundle-name>-previous-*/bundle`, and failed staging trees as
`.<bundle-name>-staging-*`. Budget space for a complete additional bundle plus
retained generations; inspect and explicitly remove old generations when no
longer needed. These archives are never part of the active served bundle.
The remaining rows above are release gates, not claims that all update classes
already converge. Metadata fingerprints use canonical JSON of the resolved
BibTeX entry, never a digest of the whole bibliography. Filename is a separate
consumed input because fallback titles/years and provenance use it even without
a bib match. Migrating a metadata receipt alone refreshes metadata, not OCR,
Docling or chunking. Separately, **legacy builds without configuration receipts
need a one-time Stage 1 refresh**: their original settings cannot be proven.
Use `corpus run --dry-run` to inspect the work before launching that migration.
`corpus status` with a config also inventories current PDF bytes and compares
resolved per-paper BibTeX, OCR/page directives, source paths/copies, lexicon
categories and the built taxonomy snapshot with receipts. Its JSON separates
`configuration_drift` from `source_input_drift`, including added/removed hashes.
Configured unreadable sources fail the audit; missing `input_pdfs` explicitly
means the source inventory was not checked. This is read-only build/operator
work, not an MCP query. It does not contact services or prove custom model
weights unchanged or infer previous CLI overrides. It checks configured
taxonomy-source receipts without refreshing the source. These limits are
included in the report's scope.
Metadata/annotation-only reruns also preserve the existing figure/vision layer;
the inexpensive figure report is refreshed because it consumes the header.

Annotation has an explicit output-set receipt (`annotation_outputs.json`),
even when neither taxonomy nor lexicon is enabled. Both resume gates verify
the receipt and its file digests. Removing a category, emptying its term list,
or deliberately disabling annotation moves superseded files into per-document
`annotation_history/`; a verified absent taxonomy output also retires its
taxon-index rows. A missing configured lexicon/taxonomy source is an error,
not a request to delete annotations. Legacy annotation stages migrate their
receipts once without redoing OCR or figure extraction. Neither receipts nor
annotation history ship in the served bundle.

Taxonomy ingestion produces a complete snapshot, not an accumulating union.
Its `meta.input_fingerprint` identifies source kind, root selection, parser
receipt version and consumed DwC bytes (the selected Taxon core for DwC-A).
Unchanged inputs reuse the snapshot without writes. Changes—including removing
a root restriction—build a fresh SQLite before replacing the old one; failed
ingestion leaves the old database untouched, and superseded snapshots remain
under `.retired/taxonomy-*.sqlite`. Legacy snapshots without proof rebuild once.
Full `corpus run` notices source drift; phase-split extraction requires a matching
pre-built snapshot. WoRMS is deliberately pinned between explicit
`corpus taxonomy ingest --rebuild` refreshes; no status request polls the API.
Taxonomy-to-work `authority_match` links are likewise derived over the current
snapshot and bibliography on each post pass. Replaced authorship strings and
newly acquired works replace old links; an unchanged result writes nothing.
Unreferenced taxonomy-only stubs are retired, while cited works and links with
other provenance remain intact. This does not broaden the existing zoological
authority-parsing policy (#175).

#### Stage 1 configuration and cache ownership

`pipeline/build_inputs.py` records readable setting values under each stage's
`input_fingerprint.config`; both resume gates use the same map. Removing a
setting compares against its built-in default. Explicit CLI panel/model
settings take precedence over the config file. Scheduler placement and output
paths are not extraction inputs.

| Setting family | First consumer; downstream invalidation |
|---|---|
| OCR detection defaults, thresholds and language-probe controls | Scan detection; prepared PDF and everything derived from it |
| OCR optimization, workers and timeouts | PDF preparation; extraction, metadata, chunks, figures and annotations |
| Figure raster settings and extraction accelerator | Docling extraction; chunks, figures and annotations, not Grobid metadata |
| `chunking.max_tokens` | Fallback chunker; annotation and chunk/figure links. Docling's HybridChunker still uses its own default tokenizer/limit; this setting does not configure it |
| Grobid disablement, header/citation consolidation, request timeout and producer identity | Metadata and its quality checks; not OCR or figure detection |
| Panel mode, resolved vision model and producer identity | Figure materialization and links; not OCR, stored text/chunks or metadata |
| Huge-document and quality-gate thresholds | Their respective checks |

Producer and dependent completion receipts are cleared **before** producer
writes, including forced reruns. Failed or interrupted work cannot retain an
old success receipt. Logs and `summary.json`'s `processing_summary.rerun_reasons`
name changed inputs.

Grobid TEI reuse requires a matching prepared-PDF digest, consolidation settings,
raw-citation request policy,
pipeline version, reported service version, optional `grobid.producer_id` and
TEI payload digest. Unverified/stale TEI and its receipt
move into the build-only `metadata_cache_history/` directory; they cannot be
resurrected as active citations if Grobid is unavailable. BibTeX-only edits
reuse a verified TEI. No history is added to the served bundle.

Fulltext requests explicitly set `includeRawCitations=1`. The reference `raw`
field and observation's raw citation retain Grobid's extracted source string,
independently of its best-effort title/author parsing. This is OCR/extraction
evidence, not a claim that the string is a perfect transcription. A request
policy receipt invalidates older Stage 1 completions and TEI caches once;
future no-op runs reuse the verified result. Empty parsed bibliographies can
still coexist with reference text in the source and must not be described as
proof that a paper cites nothing.

Panel-mode/model changes reconstruct fresh base figures from `processed.pdf`
before applying caption, ROI and compound passes, including when switching
panels off. This pays one fresh figure extraction, not another OCR run, and
does not keep a permanent duplicate raster cache. Full Docling reruns also
replace the old figure tree and discard an obsolete Docling sidecar if none is
produced. Staging preserves the old generation on ordinary publication errors;
receipt invalidation makes interrupted publication retryable. This is not a
whole-build transaction. Chunk/figure links are replaced rather than appended.

The independently scheduled vision overlay has its own receipt. Success
preserves the CPU-floor receipt so later CPU resume does not overwrite vision
results; failure invalidates the shared figure layer for recovery. It refreshes
links, source-page mappings, the figure report and quality flags too.

`corpus status` compares Stage 1 config settings when it resolves the output
from config (or when an explicit `--config` is supplied). Text output is bounded;
`--json` includes all configuration differences. CLI/CPU-floor overrides may
legitimately differ. The separate source-input audit also checks PDF inventory,
renames, current BibTeX, annotation inputs and taxonomy source receipts. Neither
audit contacts remote services or instantiates a model. Whole-build
clean/incremental parity remains a separate acceptance check; a zero-difference
status report does not close it. All of this belongs to the build plane, not
the MCP request path.

#### Vision producer evidence

Vision fingerprints resolve the actual default model instead of recording
`null`. They include implementation bytes (prompts and generation handling),
relevant installed package versions and generation settings. For ordinary
HuggingFace models, the identity is the repository revision, with the loaded
model's reported revision taking precedence over the offline cache preview.
For a custom local model directory, all non-hidden regular files are
content-hashed, including weights; changing bytes with the same size/mtime is
detected. This can make status/dry-run expensive for large local models, but
does not load them into an inference engine. Missing cached revisions are
explicitly unverified, not assumed current. Results retain `vision_producer`
in `figures.json` as well as their stage receipts.

A repository revision is not a fresh checksum of every cached weight file:
do not edit HuggingFace's immutable cache in place; use a custom local model
directory for modified weights. Remote vision services expose a model ID,
not attestable weights. Pin dated IDs and set/change `figures.producer_id` for
custom deployments or aliases whose identity changes out of band. Grobid's
analogous operator assertion is `grobid.producer_id` below. These proof levels
are explicit limits, not a promise to detect undisclosed server-side changes.
Model directories/services must remain stable during a run. Runtime hardware
and numerical reproducibility still require the independent build comparison.

#### Grobid capability and recovery

A successful local metadata write is not proof of successful Grobid extraction.
`metadata.json.grobid` and the metadata stage receipt distinguish intentional
disablement, unavailable capability and complete extraction. A curated BibTeX
header still wins, but does not hide a missing reference list. Status counts
recorded `disabled`, `unavailable`, `request_failed`, `parse_failed`, `extracted`
and `cached` outcomes independently of stage timing success. Legacy evidence is
`unknown`. These are persisted outcomes, not live health checks.

| Transition | Next build behavior |
|---|---|
| Disabled → enabled and reachable | Refresh metadata/references only, preserving unrelated extraction and figures |
| Unavailable at startup → reachable | Retry fallback papers through both resume gates |
| Individual request/parse failure → reachable | Retry the affected paper; successful papers retain their receipts |
| Successful extraction → temporary outage | Preserve complete receipts and verified TEI; BibTeX edits can still reuse that TEI |
| Enabled → deliberately disabled | Archive active TEI; replace Grobid-derived header/references/in-text data with BibTeX/filename fallback |
| Service version or declared producer identity changes | Refresh metadata and reject the old TEI cache |

Repeated offline builds do not churn fallback artifacts. Malformed/error-page
responses never become successful TEI cache entries; failed parses retire the
active cache so retries cannot loop over the same bad payload. Archived TEI is
retained under `metadata_cache_history/`, not automatically restored or bundled.
Default operation permits fallback; `--strict-network` returns failure for
startup unavailability or transient request failures instead of accepting them
as completed extraction.

The build checks liveness and the reported version once at startup. Grobid's
[`/api/version`](https://grobid.readthedocs.io/en/latest/Grobid-service/) reports
service version/revision, not a complete custom-model/configuration identity.
For those deployments, set `grobid.producer_id` to an image/model/configuration
digest or other immutable identifier and change it with the producer. This is
an **operator declaration**, not remotely verified evidence. Unreported versions
are explicitly unknown; same-version custom model changes cannot be detected
without updating the declared ID. Keep the service deployment fixed throughout
a run, just as the pipeline checkout must stay fixed.

Service URLs are placement, not model identity: a moved service with the same
version/declared ID does not force another extraction. URL resolution is explicit
CLI flag → `GROBID_URL` → config → default. Disablement wins over a URL. The
configured Grobid request timeout is applied to the client and fingerprinted.
Dry runs and status do not probe the service; a dry-run recovery estimate cannot
know whether an unavailable service has since returned.

Legacy metadata receipts without this capability evidence need one metadata
refresh, not an OCR/Docling rerun. TEI lacking the new provenance is archived;
have Grobid available during migration to regenerate its reference evidence.

## Figure pipeline

A corpuscle can span centuries of literature — 19th-century engraved plates with facing-page captions, mid-century half-tone figures, and born-digital vector panels — so figure extraction is the part of the pipeline most exposed to layout variation. Quality, segmentation, classification, and caption association are an explicit long-term optimization target; new document layouts will keep surfacing new failure modes. This section documents the full lifecycle, the code that implements each step, where resolution is set (and where it is *not* lost), and what degrades when an optional pass is skipped. The two leaf modules are `pipeline/figures.py` (pure functions: classify, parse, dedupe, caption, ROIs, linking) and `pipeline/figure_passes.py` (the Pass 2.5 / 3a / 3b orchestration); both are driven per-paper from `pipeline/runner.py`, with extraction itself in `pipeline/extract.py`.

### Lifecycle (per PDF)

The eleven steps below are all driven per-PDF from `pipeline/runner.py`. The "Stage · when" column ties each step back to the stage vocabulary above: steps 1–7 and 11 run during Stage 1 (CPU) — every paper, except the conditional PyMuPDF fallback (step 5). Steps 8–10 are the figure **passes** selected by `figures.panel_detection` (Pass 3a `ocr` is the default floor; Pass 3b is opt-in via the `vision-*` modes; Pass 3c auto-follows a compound); of those, only Pass 3b leaves the CPU stage (GPU/API).

| # | Step (figure pass) | Code | Stage · when | Writes |
|---|---|---|---|---|
| 1 | **Docling extraction** — iterate `document.pictures`, render each to PNG, capture page + bbox | `pipeline/extract.py:165-257` (render/save at `:184-189`, `:232`) | Stage 1 · always | `figures/*.png`, base fields |
| 2 | **Caption association** — structural evidence plus ranked geometric candidates | `pipeline/figures.py:extract_caption_info` | Stage 1 · always | selected caption fields plus status, confidence and bounded candidate evidence |
| 3 | **Figure-number parsing** — multilingual, incl. Roman numerals and bounded OCR repair | `pipeline/figures.py:parse_figure_number`; structural-link-only recovery in `parse_structural_caption_number` | Stage 1 · always | `figure_number`, `figure_number_source` |
| 4 | **Classification** — figure / plate / subpanel / graphical_element / unclassified | `pipeline/figures.py:classify_figure` `:425-471` | Stage 1 · always | `figure_type` |
| 5 | **PyMuPDF fallback** — raw embedded-image extraction when docling finds nothing | `pipeline/extract.py:273-354` (gated `:277`, `fitz.Pixmap` `:310`) | Stage 1 · only if docling yields 0 figures and the PDF is not a scan | `width`, `height`, `extraction_method: pymupdf` |
| 6 | **Deduplication + panel grouping** — bbox-overlap merge, whole-figure-vs-subpanel split, panel-letter assignment in reading order | `pipeline/figures.py:dedupe_figures` `:494-599` (overlap `:394-413`, reading order `:474-491`) | Stage 1 · always | `panel_letter`, refined `figure_type` |
| 7 | **Pass 2.5** — caption panels, grouped-plate figure targets + missing-figure detection | `pipeline/figure_passes.py`; `detect_missing_figures` in `pipeline/figures.py` | Stage 1 · always | `panels_from_caption`, `plate_figures_from_caption`, top-level `missing_figures[]` |
| 8 | **Pass 3a** — OCR panel / known grouped-plate ROIs | `pipeline/figure_passes.py`; `detect_figure_rois` in `pipeline/figures.py` | Stage 1 · default floor (`--figure-panels ocr`) | `rois[]` (`source: ocr:tesseract`), `pass3_status`, `pass3_target_kind` |
| 9 | **Pass 3b** — vision-LLM panel, grouped-plate, compound + bare-plate number detection (replaces 3a when selected) | `pipeline/figure_passes.py`; backends in `pipeline/vision.py` | Vision pass · `--figure-panels vision-local\|vision-claude`; GPU/API, separable on HPC | `rois[]` (`source: vision:…`), `parent_figure_index`, `pass3_backend`, audited plate-number candidates + logical records |
| 10 | **Pass 3c** — compound-figure resolution | `pipeline/figures.py:1298-1502` (trigger `runner.py:282-293`) | Stage 1 · auto when a 3a/3b status ends in `_compound` | renames PNG to `fig_3-4.png`, new `image_shared_with` sub-figure records |
| 11 | **Chunk-figure linking** | `pipeline/figures.py:link_chunks_to_figures` `:1510-1581` | Stage 1 · always | `figure_refs` (on chunks), `referenced_in_chunks` (on figures) |

`figures.json` intentionally contains both **physical detections** and
**logical evidence records**. A grouped historical plate is one physical image
but can support many independently retrievable figure identities. Older
caption-derived children name the host with `shares_image_with`; Pass 3c and
vision-discovered children use `image_shared_with`. Both mean "reuse this
host image", not "another image was detected". Physical-detection QC excludes
those children and collapses typed panel siblings; caption-binding QC retains
them because each child is a distinct source identity. Consumers must not use
raw entry counts as physical-figure counts.

Figure records are written to `<HASH>/figures.json`. Generate the optional
self-contained page audit when a decision needs visual review:

```bash
python -m pipeline.page_report <output>/documents/<HASH> --pages 12-13,17
```

It places the processed page, toggled PDF-word / figure / ROI / caption-
candidate overlays, selectable Docling text, caption decisions and page-level
statistics in `<HASH>/page_report.html`. It can be rerun for any page subset,
is not copied into the served bundle, and replaces the old unconditional
`visualizations/*.png` raster set. Reports over more than 200 selected pages
require an explicit `--max-pages 0` override. Bbox coordinate systems differ
by path and are tagged per record (`bbox_coord_system`:
`pdf_pts_bottom_left` for Docling, `pdf_pts_top_left` for PyMuPDF).

### Figure resolution — where it is set, and what is (not) downscaled

The saved figure's resolution is fixed **entirely at extraction time**; nothing downstream — dedup, the vision passes, or the MCP server — ever resizes or re-encodes the stored PNG.

- **Docling path (the common case), `resolution_mode: native` (default, #121).** Docling detects + classifies figures (rendering its own crops at `figures.images_scale`), then a PyMuPDF pass (`pipeline.figures.render_figures`, called from `extract.py`) **re-renders each figure's bbox at its source's native pixel density** and overwrites the saved PNG. So a figure backed by a 600-dpi scan stays 600 dpi and a 150-dpi one stays 150 — resolution tracks the source and **varies per figure**, rather than a single fixed DPI. A **vector** figure has no native resolution, so it renders at `figures.vector_dpi` (default **300**); `figures.max_dpi` optionally caps dense full-page scans. The bbox region (not the raw embedded xref) is rendered, so vector annotations + composite raster/vector + multi-panel plates survive at native fidelity.
- **Docling path, `resolution_mode: fixed`.** Skips the native pass; the saved PNG is docling's render at **`72 × images_scale` dpi** (default `2.0` → 144 dpi; `1.0` → 72 dpi was the grainy pre-v0.6 behavior). A single uniform DPI for every figure — predictable size, but down-renders high-res scans. Both modes only affect future ingests; lift an existing bundle in place with `tools/backfill_figure_dpi.py [--native | --scale S]` (re-renders from stored bbox + `processed.pdf`, no docling re-run).
- **PyMuPDF fallback.** `fitz.Pixmap(doc, xref)` (`pipeline/extract.py:310`) pulls the embedded image at its **native stored resolution** — no render, no scaling. So, paradoxically, the fallback path can yield *higher*-resolution figures than the primary docling path; but it only fires when docling extracts zero figures and the PDF is not a scan (`:277`). The native pass skips these (already native).
- **No serve-time downscaling.** Whole-figure requests return the original PNG
  bytes. Panel requests crop the original resolution through the shared
  disposable cache, without thumbnailing. Crop inputs above 128 MiB compressed
  or 64 million pixels are rejected explicitly to bound query-time resource
  use, never silently reduced. `get_figure_roi_image.image_path` identifies a
  disposable crop logically; it is not a new file in the bundle. Retrieve its
  bytes through `get_figure_image` or `get_figure_url`.
- **The 2000px clamp is model-input only.** Pass 3b downsamples the image handed to the vision model to ≤2000px on the long side (`pipeline/vision.py:226-243`) purely to bound API cost; it does **not** touch the saved figure.
- **Metadata gap to be aware of.** Only the PyMuPDF path records pixel `width`/`height`; the docling path records no dimensions and **no dpi** in `figures.json`, so stored resolution can't currently be audited from metadata alone — worth closing when figure-quality work starts.

### Caption association

Captions are the highest-value annotation per figure and the hardest in
historical layouts. `extract_caption_info` records the result as evidence,
not just as a string:

1. **Structural candidate.** Resolve Docling's picture-to-caption link. Text
   is read by provenance span rather than only from the whole `TextItem`,
   because Docling can merge running prose from one page with a figure label
   on the next. Tightly adjacent caption openers and bodies are joined even
   when Docling labels a panel continuation as ordinary text or a list item;
   each added fragment must extend a valid letter-panel declaration.
2. **Geometric candidates.** Inspect provenance-level figure labels on the
   picture page and immediately following page. Same-page candidates outrank
   facing-page candidates; structural evidence outranks a proximity tie. A
   same-page candidate naming a different figure can override a structural
   link only when it is within the caption-component gap and at least 24 PDF
   points closer; the displaced link remains in the rejected evidence.
3. **Ownership check.** A next-page heuristic is rejected when a substantial
   picture on that page is a better local owner for the candidate. Historical
   facing-page captions remain possible, but are marked low confidence.
4. **Grouped and facing-page legends.** A caption block naming several figure
   numbers is split into per-number entries; numeric lists/ranges are separate
   from lettered panels. Exact-count duplicate assignments are reconciled into
   a one-to-one mapping. An entry matching a bare-label figure on the preceding
   page enriches that record instead of being attached to the current page's
   plate. The exact number match is `bound` / `medium`, even though its page
   distance is one. When several logical records share one plate image, Pass
   2.5 stores their numbers as `plate_figures_from_caption`, records which were
   independently cited as missing before expansion, and admits the host image
   to one numeric ROI pass. The resulting regions are distributed back to the
   individual records; they never enter `panels_from_caption` as fake letters.
5. **Number provenance.** Ordinary labels remain start-anchored. A Docling
   structural link may additionally recover a measured damaged label after a
   short species heading (for example `Physalia physalis Fic. 6`); the same
   embedded text is never admitted by proximity search. The raw caption is
   preserved, and `figure_number_source` distinguishes `caption_start`,
   `docling_caption_link_embedded_ocr_label`, and the plate-legend
   reconciliation variants. OCR digit repairs are intentionally narrow
   (`ro`/`ror` in label context); arbitrary letter runs are not numbers.

Pass 2.5 recognizes period, parenthesized, comma and range panel styles through
letter L. Small noncontiguous sets are rejected as probable initials; a strong
multi-label set can retain printed gaps such as A–H, K, L. Period-style marker
scanning stops at an abbreviation glossary, preventing OCR-split keys such as
`C. rad.lat` from becoming invented panels.

The selected text stays complete. `caption_candidates` retains at most five
chosen/rejected records with candidate text capped at 600 characters, bbox,
source, page distance, geometric distance, confidence and rejection reason.
The summary fields are:

- `caption_status`: `bound`, `uncertain` or `unbound`;
- `caption_confidence`: `high`, `medium`, `low` or null;
- `caption_kind`: `prose_caption`, `bare_label`, `unlabelled_caption` or null;
- `caption_page_distance`: caption page minus picture page;
- `caption_source`: the producing rule, including `docling_caption_link`,
  `heuristic_proximity`, `plate_legend`, the plate/facing-page reconciliation
  variants (including `preceding_page_plate_legend`), the compound-split
  variants, or null.

`figure_number_source` is separate from `caption_source`: the first records
how the identifier was read or reconciled, while the second records why the
caption belongs to this picture. Keeping both prevents an OCR repair from
looking like source text and prevents a strong caption link from obscuring a
weaker number inference.

An unopposed cross-page heuristic is `uncertain` / `low`, never an ordinary
caption. No candidate yields `unbound`, with empty caption text. The figure MCP
responses expose these summary fields; for an older bundle the server
derives only the compatible summary from its already-stored caption, source
and pages (with `caption_kind: unknown`)—it does not rerun association at
query time.

### Selecting a pass — and what is lost when one is skipped

Panel detection is selected by **`figures.panel_detection`** in the corpuscle's `config.yaml` (`ocr` / `vision-local` / `vision-claude` / `off`; #102), overridable per run with **`corpus run --figure-panels <mode>`**. The deprecated `corpus run --no-vision` is an alias for `--figure-panels ocr`. Under the hood `corpus run` translates this into the `--figure-panels` flag on `pipeline.main` / the orchestrator (which is where it appears in `slurm/` jobs).

Since #102 the default is **`ocr`** — Pass 3a is the CPU floor and runs on a plain `corpus run`. So the default output now includes OCR `rois`; only `panel_detection: off` reproduces the pre-#102 "no panel geometry" behavior.

- **With `off` (no Pass 3a or 3b → no `rois`):** there is no geometric segmentation of multi-panel figures or known grouped plates. Region retrieval can only fall back to the whole image plus the caption-derived `panels_from_caption` / `plate_figures_from_caption` descriptions.
- **Pass 3a (OCR, the `ocr` default) vs Pass 3b (vision):** 3a runs Tesseract on a 3×-upscaled image. It searches an ordinary figure only for its caption-declared A/B/C labels; on a known grouped plate it instead searches for the exact caption-declared number allow-list. This avoids treating arbitrary chart numbers as figures. It has low recall on line art and engravings. 3b (Qwen2.5-VL locally or Claude via API, the `vision-local` / `vision-claude` modes) handles both target kinds more reliably and is the **only** pass that discovers an *unexpected compound* — several logical figures in an image that caption analysis had not already expanded. It also admits a bare, confidently bound `plate` record to unconditioned number discovery. That path requires at least two distinct Arabic number+bbox candidates at confidence ≥0.80; retains accepted and rejected candidates on the host; and emits one logical, image-sharing figure record per accepted region with `caption_status: unbound`. Self-reported confidence is not sufficient when the output contradicts itself: a fabricated alphabetic panel grid whose boxes are copied onto numeric candidates is retained as rejected evidence and materializes nothing. It never turns a model description into source caption text. Ordinary figures, prose-captioned plates, uncertain plate links and the OCR pass remain allow-list-only. The two passes are mutually exclusive.
- **Without Pass 3b → no unexpected-compound recovery in Pass 3c:** a previously unknown `fig_3-4.png`-style image stays a single record with an ambiguous `figure_number`, and `missing_figures[]` entries that 3c would have matched remain unresolved. A caption-enumerated historical plate is different: Pass 2.5 has already created its logical records, and the default OCR pass may locate their regions without invoking 3c.

Historical legends can sit on the same page as a plate, on the following page,
or on the preceding facing leaf. The preceding-leaf rule is deliberately not
mere adjacency: the legend page must carry an explicit `PLATE N` heading and
the following page must contain exactly one matching `PLATE N` picture.
Measured OCR-damaged `Fic.`/`Fics.` openers are admitted only inside that
plate-number context. Plate identity and child-figure identity are separate
namespaces, so Plate X and Figure 10 can share an image without deduplication.

The separately schedulable `corpus run --only vision` path has artifact parity
with inline vision: after Pass 3b it runs compound materialization (Pass 3c)
and rebuilds bidirectional chunk/figure links. The running MCP server performs
none of this reconstruction; it reads the materialized records and crops their
persisted regions.

Because these passes are GPU/API-cost-bearing, the corpus-scale validation of figure coverage is tracked separately ([#11](https://github.com/caseywdunn/corpus/issues/11)), and figure-number recovery on old/scanned papers is an open gap ([#16](https://github.com/caseywdunn/corpus/issues/16)). When a new layout exposes a new failure mode, capture it as a ground-truth fixture under `tests/` so the fix is guarded against regression.

## Taxonomic annotation

When a Darwin Core taxonomy SQLite snapshot is available (`<corpuscle>/taxonomy.sqlite`, built by `pipeline/taxonomy_ingest.py`), the pipeline annotates chunks with recognized taxon names. The snapshot can be built from any DwC source — the lab default for siphonophores is WoRMS pruned to Siphonophorae (`--source worms --root-id 1371`), but `ingest_taxonomy.py <corpuscle> --source dwc --input <Taxon.tsv>` ingests any downloaded DwC export, optionally pruned to a subgraph via `--root-id <taxonID>`. Schema follows the Darwin Core Taxon class (`taxonID`, `scientificName`, `parentNameUsageID`, `acceptedNameUsageID`, …).

> **WoRMS is marine-only (`isMarine=0` records are excluded) (#96).** The WoRMS DwC-A backbone — and the `--source worms` REST walk — only carry taxa flagged marine. A corpuscle whose literature spans freshwater or terrestrial groups (or marine taxa with non-marine relatives) will hit *silent* resolution failures: those names simply don't resolve to a `taxonID`, so they're dropped from `taxa.json` with no error, and `search_taxon` returns `not_found`. For a non-marine or mixed clade, build the snapshot from a GBIF or WFO DwC export (`--source dwc --input Taxon.tsv`) instead of WoRMS. This is a data-source limitation, not a pipeline bug — there is no flag that makes WoRMS emit non-marine records.

## Lexicon-driven annotation

A user-supplied lexicon tags chunks with category-specific terms. The YAML is two-level — top-level keys are categories, each value is a flat term map:

```yaml
anatomy:
  nectophore:
    synonyms: [nectophores, swimming bell]
biogeography:
  pelagic:
    synonyms: [open water]
```

Pass it with `--lexicon path/to/lexicon.yaml`; each category emits its own `<hash>/<category>.json` (so `anatomy.json`, `biogeography.json`, …). See [demo/lexicon.yaml](../demo/lexicon.yaml) for a worked siphonophore example.

Lexicons are inputs you maintain alongside the literature, not part of the tool. Each category's content is fingerprinted independently (SHA-256 over the canonical JSON of just that section) and recorded both inside the `<category>.json` artifact and in the per-paper `pipeline_state.json` completion record. On `--resume`, editing one section re-runs `taxa_and_lexicon_extraction` against the new fingerprint; sections whose hash didn't change stay cached.

## Cross-paper databases

Built from per-paper artifacts after Stage 1 finishes. All four are independently rebuildable without touching OCR / extraction.

| Database | Builder | What it stores |
|---|---|---|
| `taxonomy.sqlite` | `pipeline/taxonomy_ingest.py` | Darwin Core taxon backbone + synonymy. |
| `biblio_authority.sqlite` | `bib/authority.py` (+ `bib/reconcile.py`) | Append-only reference observations, auditable observation-to-work decisions, canonical works, and the materialized citation graph; work identities resolve to DOI / BHL Part / normalized citation key. |
| `taxon_mentions.sqlite` | `pipeline/taxon_mentions.py` | Cross-paper taxon-name index from gnfinder + abbreviated-form expansion (`A. elegans` → `Agalma elegans`). |
| `intext_citations.json` (per-paper) | `pipeline/intext_citations.py` | TEI body `<ref type="bibr">` elements joined to chunk offsets and resolved to `work_id` in the bibliographic authority. |

`corpus run` chains the pipeline, embeddings, and these four post-pipeline steps in dependency order (see `pipeline/orchestrator.py`); `corpus status` reports what completed, what is stale, and what carries quality flags.

### Reference evidence and canonical works

`work_documents` is the current PDF-to-work membership relation. Multiple
documents may share a canonical work (for example, separately scanned volumes
with one DOI). `works.corpus_hash` is only a deterministic representative for
legacy/frozen responses, preferring a served member; it is not an inventory.
Paper-hash lookups use indexed membership queries, and every member can supply
its own reference observations and citation edges.

License, serving exclusions, OCR/page directives and curation notes remain
document-local. BibTeX export emits a separate, uniquely keyed entry for each
member, carrying `corpus_hash`; imports use that key to avoid applying one
scan's policy to all scans of the work. A work-only import with document-local
fields is rejected when multiple documents would be affected. Reconciliation
and maintenance merges move memberships before removing a canonical row.

Authority seeding compares the content of consumed metadata fields, not just
mtime. A DOI edit moves the affected membership without orphaning the old
graph; reference materialization then rebuilds against the current membership
and canonical metadata fingerprint. Removing a document artifact directory
deactivates its membership and reference set but retains historical reference
observations. This does not by itself remove an artifact directory whose PDF
was deleted upstream; source-inventory retirement belongs to the update
contract, not authority inference.

`references.json` is evidence, not a canonical-work table. The authority
builder content-addresses every bibliography occurrence into
`reference_observations`; its raw citation, parsed fields and source-set
membership are append-only. `reference_current_sets` is the replaceable pointer
that says which observed source set is current for each paper. Reprocessing or
removing a paper changes that pointer, not its historical evidence.

`observation_work` is the independently re-derivable verdict. Every current
observation maps to one canonical `works` row and carries the match method,
score and rule-producer version. If any current reference set or the producer
version changes, the builder discards derived verdicts and uncurated ghosts and
recomputes the complete current mapping in a stable content-address order. An
unchanged run does not update timestamps or rows. This makes a clean build and
an incremental addition converge on the same current map while retaining the
conservative rule: weak evidence creates a separate work rather than silently
misrouting a citation.

The deterministic resolver normalizes DOI wrappers and percent encoding, and
tests narrowly shaped DOI corruption only when title evidence independently
agrees. To escape a damaged or reordered first-author block, it may match an
in-corpus work across all blocks only when the year and complete normalized
author-surname set agree, both titles are substantive, exactly one candidate
passes, and the titles satisfy the established exact or dual fuzzy thresholds.
A conflicting DOI is preserved in the observation and named in the mapping
method; it is never silently rewritten. Title/year agreement without author-set
agreement is review evidence, not permission to merge (#155, #225).

The existing `citations` table is a compatibility materialization of that map,
not a second source of truth. The frozen MCP bibliography tools continue to
read it unchanged; clustering, reconciliation and evidence mutation therefore
remain build-plane work rather than leaking into the running server.
`bib.reconcile` may still choose the surviving canonical row for a corpus PDF,
but it records the scored decision in `work_reconciliation_decisions` and never
deletes the underlying reference observation. Maintenance merges redirect
`observation_work` before removing a superseded canonical row.

## MCP server (`mcpsrv/`)

Exposes the processed corpus as an MCP (Model Context Protocol) server that LLM clients can query. The server is a read-only view over per-paper artifacts; it does not store data of its own. The server entry point `corpus serve` is a thin shim — the implementation lives in `mcpsrv/`, with the `@mcp.tool()`-decorated functions split across `mcpsrv/tools/{papers,taxonomy,bibliography,figures,chunks,lexicon,profiles}.py`. See [MCP_TOOLS.md](MCP_TOOLS.md) for the full tool surface and count.

## Steering the client session

The server often needs to shape *how* the client uses the corpus — report vs. manuscript structure, cross-validation of synthesized claims, citation discipline — not just answer queries. MCP is **client-driven**: a server can only add to the model's context as a *response* to a client-initiated request (`initialize`, `tools/call`, `prompts/get`, `resources/read`), and can **never unilaterally push a turn** mid-session. Elicitation, sampling, and notifications do *not* inject free-form prompting into the model's context. So all steering rides on three response channels, in increasing specificity:

1. **`instructions`** — `InitializeResult.instructions`, sourced from `<corpuscle>/instructions.md`, injected once at session start. The always-on baseline: tool routing, citation rules, "respect the active output profile." Established.
2. **Tool-result guidance** — short guidance text appended to a tool's response payload, which lands in context the moment that tool is called. This is the idiomatic *just-in-time* nudge channel; it beats front-loading everything into `instructions` (which the model reads once and drifts from). `format_citations`'s provenance warning is the established example. The intended extension is to make this **output-profile-aware** (tracked in [#101](https://github.com/caseywdunn/corpus/issues/101)): e.g. under a `manuscript` profile, figure tools remind the model to verify publishability and emit the attribution string, and synthesis retrievals nudge cross-validation against multiple sources.
3. **MCP Prompts** (`prompts/list` / `prompts/get`, surfaced in Claude Code as `/mcp__corpus__<name>`) — user-invoked structural scaffolds (a manuscript skeleton, a monographic-review recipe, a cross-validation checklist). Zero context cost until invoked. Forward-looking, tracked in #101.

Two rules govern the idiom:

- **Soft steering** (document structure, cross-validation, house style) goes in the nudge channels above, keyed to the client-selected output profile. It is advisory — the model may ignore it.
- **Hard requirements** (figure publishability, citation provenance) are enforced **server-side at the tool boundary**, never delivered as a nudge the model can ignore. This is the #79/#100 lesson: a trust-critical gate must be code, not a prompt.

Tool-result nudges cost tokens on every call and pull against the served-bundle payload-trimming work (#76, #81–86), so keep them terse and conditional rather than emitting on every call.

## Key files

| File / package | Role |
|---|---|
| `pipeline/` | Stage 1 + Pass 3b/3c orchestrator (split per-stage; CLI in `pipeline/main.py`) |
| `pipeline/main.py` | Stage 1 CLI: per-paper extraction, OCR, metadata, chunking, annotation |
| `pipeline/orchestrator.py` | Runs pipeline + post-pipeline steps in dependency order (backs `corpus run`) |
| `pipeline/status.py` | Rollup of stage completion, quality flags, staleness (backs `corpus status`) |
| `pipeline/embed.py` | Stage 2 embedding (BGE-M3 → LanceDB) |
| `pipeline/figures.py` | Figure extraction, classification, caption parsing, chunk-figure linking |
| `pipeline/vision.py` | Vision-LLM backends (local Qwen2.5-VL, Claude API) |
| `pipeline/taxa.py` | Taxonomy DB access, taxon-mention extraction, lexicon loaders |
| `pipeline/grobid_client.py` | Grobid TEI parsing for headers, references, in-text citations |
| `pipeline/embeddings.py` | BGE-M3 embedding backends (local, HF, etc.) |
| `pipeline/external.py` | Shared retry + circuit breaker + `--strict-network` mode |
| `pipeline/version.py` | Single-source `__version__` stamped into every artifact |
| `mcpsrv/` | MCP server implementation (MCPServer, eager index, stdio + SSE transports) |
| `pipeline/prefetch.py` | Model prefetch + offline cache inspection (backs `corpus prefetch`) |
| `bib/` | BibTeX parser, importer, exporter (round-trip curation) |
| `config.yaml` | Pipeline configuration (loaded by `pipeline.config.load_config`) |
| `slurm/batch_pipeline.sh` | SLURM orchestrator: chains Grobid, Stage 1 array, cleanup, Pass 3b, Embed |
| `slurm/batch_process_corpus.sh` | SLURM Stage 1 batch script (supports job arrays) |
| `slurm/batch_pass3b.sh` | SLURM Pass 3b batch script (GPU) |
| `slurm/batch_embed.sh` | SLURM embedding batch script (GPU) |
| `slurm/bouchet_paths.sh` | Shared path definitions for all batch scripts |
| `pipeline/taxonomy_ingest.py` | Build Darwin Core taxonomy SQLite from a DwC file, archive, or the WoRMS API |
| `bib/authority.py` | Bibliographic authority DB (deduplicated works + citation graph) |
| `pipeline/taxon_mentions.py` | Cross-paper taxon mentions SQLite |
| `pipeline/intext_citations.py` | TEI body → `intext_citations.json` per paper |
| `bib/reconcile.py` | Merge ghost cited-references onto corpus papers |
| `mcpsrv/bundle.py` | Whitelist + manifest the served bundle for S3 / EC2 deploy |
| `demo/lexicon.yaml` | Example multi-category lexicon (siphonophore anatomy under `anatomy:`) |
