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

Two current serve-path details still need to be brought into this contract:
panel crops are cached under the mounted bundle instead of a separate
disposable cache, and remote figure URLs expose a bearer credential in a
model-visible response. They are known architecture gaps, not precedents for
adding more serve-time materialization.

## Content-addressed storage

Every unique PDF is identified by the first 12 hex characters of its SHA-256 hash. All artifacts for that PDF live under `<output_dir>/documents/<HASH>/`. Identical PDFs found at multiple input paths are processed once; every discovered path is recorded in `summary.json`. The presence of `summary.json` is the completion marker that `--resume` checks; per-stage resume (#28) additionally tracks each stage's artifact independently so a lexicon edit only re-runs the annotation pass.

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
- Embeddings are batched. Per-hash completion is marked by `<HASH>_embedded.done`.
- `--resume` skips already-embedded hashes. `--rebuild` drops the table and re-embeds (use when switching models).

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
- **No serve-time downscaling.** `get_figure_image` returns the PNG byte-for-byte (`mcpsrv/tools/figures.py:707-708`); the HTTP route does `target.read_bytes()` with no processing (`mcpsrv/figure_http.py:158`). Panel crops are cut from the full-resolution PNG on demand and cached (`tools/figures.py:710-732`), inheriting source resolution. There is no max-dimension cap, thumbnailing, or re-encode anywhere on the serve path.
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
