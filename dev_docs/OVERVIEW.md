# Architecture Overview

The corpus pipeline processes a collection of scientific literature PDFs into a searchable, queryable knowledge base. It handles born-digital papers, scanned historical plates, multilingual text, and 19th-century Fraktur typefaces.

## Repo layout

The top-level entry-point scripts are thin shims; the implementation is grouped into packages:

| Package | Role |
|---|---|
| `pipeline/` | Stage 1 + Pass 3b/3c orchestrator and shared library modules. Split into `scan.py` (OCR), `extract.py` (docling), `metadata.py` (Grobid + bib), `chunking.py`, `annotate.py` (taxa + lexicons), `figure_passes.py`, `runner.py` (per-paper orchestrator), `main.py` (CLI), and supporting `config.py` / `io.py` / `log.py` / `stages.py`. Shared library modules: `figures.py`, `taxa.py`, `grobid_client.py`, `embeddings.py`, `vision.py`, `external.py`, `version.py`. |
| `mcpsrv/` | MCP server. `app.py` defines the FastMCP instance; `tools/{papers,taxonomy,bibliography,figures,chunks}.py` register the 27 `@mcp.tool()` decorated functions; `transport.py` handles stdio + SSE; `indexes.py` is the eager in-memory index. |
| `bib/` | BibTeX import / export round-trip plus shared metadata helpers (`parser.py`, `importer.py`, `export.py`). |
| `slurm/` | SLURM batch scripts (Bouchet). |
| `deploy/` | CloudFormation, nginx config, systemd unit, sync + update shell scripts. |
| `tests/` | Ground-truth + corpus-wide consistency tests. |

The top-level [`process_corpus.py`](../process_corpus.py) and [`mcp_server.py`](../mcp_server.py) are CLI shims that re-export the packages above so existing invocations keep working.

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
      visualizations/          # per-page QC overlays (red=words, yellow=figures)
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

The pipeline is split into independent stages so that CPU-heavy work, GPU work, and API-cost work can run on different hardware and be resumed independently.

### Stage 1: CPU processing (`pipeline/`)

Orchestrated by `pipeline.runner.run_pdf_processing_pipeline`, six steps per PDF:

| Step | What happens | Output |
|---|---|---|
| **Scan detection** (`pipeline/scan.py`) | Multi-stage heuristic — text-volume gate, langdetect on the existing text layer, gibberish-score, and (for borderline cases) a Tesseract OSD visual-script cross-check — classifies as `born_digital`, `scanned`, or `broken_text_layer` | `scan_detection.json` |
| **PDF preparation** (`pipeline/scan.py`) | Copies born-digital PDFs as-is; runs `ocrmypdf` on scanned/broken PDFs with language-appropriate Tesseract packs (the default fallback union covers eng/deu/fra/rus/lat/spa/por/chi_sim/chi_tra/jpn/ell/kor + deu_latf Fraktur — configurable via `ocr.ocr_languages_default`) | `processed.pdf` |
| **Text + figure extraction** (`pipeline/extract.py`) | Docling parses the PDF into structured text and figure regions. Figures go through a classification/caption pipeline (see [Figure pipeline](#figure-pipeline) below). Falls back to raw PyMuPDF image extraction when docling finds nothing. | `text.json`, `figures.json`, `figures/*.png`, `visualizations/*.png` |
| **Metadata extraction** (`pipeline/metadata.py`, `bib/`) | Grobid extracts title, authors, year, DOI, abstract, section structure, and parsed references. `--bib` overrides the header from a curated BibTeX. Falls back to placeholder when Grobid is unavailable. | `metadata.json` |
| **Chunking** (`pipeline/chunking.py`) | Splits extracted text via docling's `HybridChunker` (tokenizer-aware, respects section/heading structure), with section-class labels. | `chunks.json` |
| **Annotation** (`pipeline/annotate.py`) | Per-chunk taxon mentions (against the DwC taxonomy snapshot) and lexicon matches (one pass per category in `--lexicon`). Each output file stamps a per-category `input_fingerprint`; the stage-completion record in `pipeline_state.json` mirrors it so a per-category resume detects exactly which categories changed. | `taxa.json`, `<category>.json` |

Stage 1 supports SLURM job-array parallelization via `--batch-index` / `--batch-size`. Each array task deterministically processes a slice of the sorted hash list. See [BOUCHET.md](BOUCHET.md) for operational details.

### Pass 3b: Vision analysis (`process_corpus.py --vision-backend`)

A GPU pass that runs a vision-language model (Qwen2.5-VL-7B-Instruct locally, or Claude API) over extracted figures to detect:
- **Panel structure**: whether a figure is a multi-panel composite
- **Compound figures**: plate-style composites spanning multiple figure numbers
- ROI bounding boxes for individual panels within composites

Results are written back into `figures.json` per hash.

### Pass 3c: Compound figure renaming

Automatically renames figure files for detected compound figures using range notation (e.g., `fig_3-4.png` for a plate containing figures 3 and 4). Runs as part of the same invocation as Pass 3b.

### Stage 2: Embedding (`embed_chunks.py`)

Reads `chunks.json` per hash and produces vector embeddings stored in LanceDB.

- **Backend**: local sentence-transformers with [BGE-M3](https://huggingface.co/BAAI/bge-m3) (1024-dim, multilingual), runs on CUDA/MPS/CPU via `embeddings.detect_device()`. Embedding failures raise rather than silently inserting zero vectors.
- Embeddings are batched. Per-hash completion is marked by `<HASH>_embedded.done`.
- `--resume` skips already-embedded hashes. `--rebuild` drops the table and re-embeds (use when switching models).

## Figure pipeline

A corpuscle can span centuries of literature — 19th-century engraved plates with facing-page captions, mid-century half-tone figures, and born-digital vector panels — so figure extraction is the part of the pipeline most exposed to layout variation. Quality, segmentation, classification, and caption association are an explicit long-term optimization target; new document layouts will keep surfacing new failure modes. This section documents the full lifecycle, the code that implements each step, where resolution is set (and where it is *not* lost), and what degrades when an optional pass is skipped. The two leaf modules are `pipeline/figures.py` (pure functions: classify, parse, dedupe, caption, ROIs, linking) and `pipeline/figure_passes.py` (the Pass 2.5 / 3a / 3b orchestration); both are driven per-paper from `pipeline/runner.py`, with extraction itself in `pipeline/extract.py`.

### Lifecycle (per PDF)

| # | Step | Code | Always run? | Writes |
|---|---|---|---|---|
| 1 | **Docling extraction** — iterate `document.pictures`, render each to PNG, capture page + bbox | `pipeline/extract.py:165-257` (render/save at `:184-189`, `:232`) | yes | `figures/*.png`, base fields |
| 2 | **Caption association** — structural link then proximity heuristic | `pipeline/figures.py:extract_caption_info` `:655-752` | yes | `caption_text`, `caption_source`, `caption_page`, `caption_bbox` |
| 3 | **Figure-number parsing** — multilingual, incl. Roman numerals | `pipeline/figures.py:parse_figure_number` `:334-359` (`_FIGURE_PREFIX` `:45-50`, roman `:84-100`) | yes | `figure_number` |
| 4 | **Classification** — figure / plate / subpanel / graphical_element / unclassified | `pipeline/figures.py:classify_figure` `:425-471` | yes | `figure_type` |
| 5 | **PyMuPDF fallback** — raw embedded-image extraction when docling finds nothing | `pipeline/extract.py:273-354` (gated `:277`, `fitz.Pixmap` `:310`) | only if docling yields 0 figures and the PDF is not a scan | `width`, `height`, `extraction_method: pymupdf` |
| 6 | **Deduplication + panel grouping** — bbox-overlap merge, whole-figure-vs-subpanel split, panel-letter assignment in reading order | `pipeline/figures.py:dedupe_figures` `:494-599` (overlap `:394-413`, reading order `:474-491`) | yes | `panel_letter`, refined `figure_type` |
| 7 | **Pass 2.5 — caption panels + missing-figure detection** | `pipeline/figure_passes.py:36-95`; `detect_missing_figures` `pipeline/figures.py:269-331` | yes | `panels_from_caption`, top-level `missing_figures[]` |
| 8 | **Pass 3a — OCR panel ROIs** *(optional)* | `pipeline/figure_passes.py:98-151`; `detect_figure_rois` `:1083-1138` | `--content-aware-figures` | `rois[]` (`source: ocr:tesseract`), `pass3_status` |
| 9 | **Pass 3b — vision-LLM panel + compound detection** *(optional, supersedes 3a)* | `pipeline/figure_passes.py:154-213`; backends in `pipeline/vision.py` | `--vision-backend {claude\|local}` | `rois[]` (`source: vision:…`), `parent_figure_index`, `pass3_backend` |
| 10 | **Pass 3c — compound-figure resolution** *(auto when 3a/3b ran)* | `pipeline/figures.py:1298-1502` (trigger `runner.py:282-293`) | when a 3a/3b status ends in `_compound` | renames PNG to `fig_3-4.png`, new `image_shared_with` sub-figure records |
| 11 | **Chunk-figure linking** | `pipeline/figures.py:link_chunks_to_figures` `:1510-1581` | yes | `figure_refs` (on chunks), `referenced_in_chunks` (on figures) |

Figure records are written to `<HASH>/figures.json`; the per-page QC overlays in `<HASH>/visualizations/` (yellow figure bboxes, red word boxes) are the primary *visual* regression check, since most figure-quality properties resist unit assertions. Bbox coordinate systems differ by path and are tagged per record (`bbox_coord_system`: `pdf_pts_bottom_left` for docling, `pdf_pts_top_left` for PyMuPDF).

### Figure resolution — where it is set, and what is (not) downscaled

The saved figure's resolution is fixed **entirely at extraction time**; nothing downstream — dedup, the vision passes, or the MCP server — ever resizes or re-encodes the stored PNG.

- **Docling path (the common case).** The pipeline constructs `PdfPipelineOptions` with `generate_picture_images=True` but **leaves `images_scale` unset** (`pipeline/extract.py:63-70`), so docling renders each picture crop from a page raster at its *default* scale (docling's base render is 72 dpi at scale 1.0). Critically, this is a **rasterization at a fixed scale, not preservation of the source's native embedded resolution** — a born-digital paper with a high-resolution embedded figure is still captured at the default render scale. **This is the single biggest figure-quality lever and the most likely first optimization**: raising `images_scale` (e.g. to 2–4×) trades disk/CPU for sharper crops. Validate the exact dpi against the pinned docling (`docling==2.94.0`) before tuning, since the constant lives in docling, not this repo.
- **PyMuPDF fallback.** `fitz.Pixmap(doc, xref)` (`pipeline/extract.py:310`) pulls the embedded image at its **native stored resolution** — no render, no scaling. So, paradoxically, the fallback path can yield *higher*-resolution figures than the primary docling path; but it only fires when docling extracts zero figures and the PDF is not a scan (`:277`).
- **No serve-time downscaling.** `get_figure_image` returns the PNG byte-for-byte (`mcpsrv/tools/figures.py:707-708`); the HTTP route does `target.read_bytes()` with no processing (`mcpsrv/figure_http.py:158`). Panel crops are cut from the full-resolution PNG on demand and cached (`tools/figures.py:710-732`), inheriting source resolution. There is no max-dimension cap, thumbnailing, or re-encode anywhere on the serve path.
- **The 2000px clamp is model-input only.** Pass 3b downsamples the image handed to the vision model to ≤2000px on the long side (`pipeline/vision.py:226-243`) purely to bound API cost; it does **not** touch the saved figure.
- **Metadata gap to be aware of.** Only the PyMuPDF path records pixel `width`/`height`; the docling path records no dimensions and **no dpi** in `figures.json`, so stored resolution can't currently be audited from metadata alone — worth closing when figure-quality work starts.

### Caption association

Captions are the highest-value annotation per figure and the hardest in historical layouts. `extract_caption_info` (`pipeline/figures.py:655-752`) tries two paths in order:

1. **Docling structural link** (`:686-699`): `picture.captions[0].resolve(document)` — when docling's layout model already bound a caption to the picture, its text + provenance are taken directly and tagged `caption_source: docling_caption_link`.
2. **Proximity heuristic** (`:704-752`): scans every `document.texts` item whose text *begins* with a figure label (`_FIGURE_NUMBER_IN_CAPTION_RE`, `:67-70`), restricted to the figure's own page or the **immediately following page** — facing-page captions are routine in plate monographs where the caption leaf faces the plate. Candidates are scored by vertical gap between figure and caption (in PDF points), with a large `_CROSS_PAGE_PENALTY` (1000, `:712`) biasing toward same-page matches; the smallest-gap candidate wins and is tagged `caption_source: heuristic_proximity`. There is no previous-page fallback and no absolute distance threshold (closest wins).

When neither path finds text, the figure keeps `caption_source: null` — a useful signal for figure-quality triage.

### Optional passes — what is lost when they are skipped

All of steps 8–10 are **off by default**; a plain `corpus run` produces figures, captions, numbers, classification, dedup/panel-letters, `missing_figures`, and chunk links, but no pixel-level panel geometry and no compound splitting.

- **Without Pass 3a or 3b (no `rois`):** there is no geometric segmentation of multi-panel figures. `get_figure_image(label="B")` can only fall back to the whole figure; panel-level retrieval degrades to the caption-derived `panels_from_caption` descriptions (text, not crops).
- **Pass 3a (OCR) vs Pass 3b (vision):** 3a (Tesseract on a 3×-upscaled crop) detects only printed letter labels and has low recall on line-art/engraved plates; 3b (Qwen2.5-VL locally or Claude Haiku via API, selected by `--vision-backend`) is substantially more reliable and is the **only** pass that detects *compound* figures — a single extracted image that actually contains several numbered figures (common in plate-heavy monographs). 3b supersedes 3a when both are requested.
- **Without Pass 3b → no Pass 3c:** compound plates are never split. A `fig_3-4.png`-style image stays a single record with an ambiguous `figure_number`, its embedded sub-figures remain invisible to `figure_number`-keyed queries, and the `missing_figures[]` entries that 3c would have matched to recovered sub-figures (`image_shared_with` links, real numbers + captions) are left unresolved. The `missing_figures` list itself is still produced (Pass 2.5) as a coverage signal.

Because these passes are GPU/API-cost-bearing, the corpus-scale validation of figure coverage is tracked separately ([#11](https://github.com/caseywdunn/corpus/issues/11)), and figure-number recovery on old/scanned papers is an open gap ([#16](https://github.com/caseywdunn/corpus/issues/16)). When a new layout exposes a new failure mode, capture it as a ground-truth fixture under `tests/` so the fix is guarded against regression.

## Taxonomic annotation

When a Darwin Core taxonomy SQLite snapshot is available (`<corpuscle>/taxonomy.sqlite`, built by `ingest_taxonomy.py`), the pipeline annotates chunks with recognized taxon names. The snapshot can be built from any DwC source — the lab default for siphonophores is WoRMS pruned to Siphonophorae (`--source worms --root-id 1371`), but `ingest_taxonomy.py <corpuscle> --source dwc --input <Taxon.tsv>` ingests any downloaded DwC export, optionally pruned to a subgraph via `--root-id <taxonID>`. Schema follows the Darwin Core Taxon class (`taxonID`, `scientificName`, `parentNameUsageID`, `acceptedNameUsageID`, …).

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
| `taxonomy.sqlite` | `ingest_taxonomy.py` | Darwin Core taxon backbone + synonymy. |
| `biblio_authority.sqlite` | `build_biblio_authority.py` (+ `reconcile_corpus_to_biblio.py`) | Deduplicated works graph: corpus papers, cited references, taxonomic-authority strings; resolved to DOI / BHL Part / normalized citation key. |
| `taxon_mentions.sqlite` | `build_taxon_mentions.py` | Cross-paper taxon-name index from gnfinder + abbreviated-form expansion (`A. elegans` → `Agalma elegans`). |
| `intext_citations.json` (per-paper) | `backfill_intext_citations.py` | TEI body `<ref type="bibr">` elements joined to chunk offsets and resolved to `work_id` in the bibliographic authority. |

[`update_corpus.py`](../update_corpus.py) chains the pipeline, embeddings, and these four post-pipeline scripts in dependency order; [`corpus_status.py`](../corpus_status.py) reports stage completion, quality flags, and stale annotations.

## MCP server (`mcpsrv/`)

Exposes the processed corpus as an MCP (Model Context Protocol) server that LLM clients can query. The server is a read-only view over per-paper artifacts; it does not store data of its own. The server entry point [`mcp_server.py`](../mcp_server.py) is a thin shim — the implementation lives in `mcpsrv/`, with the 27 `@mcp.tool()`-decorated functions split across `mcpsrv/tools/{papers,taxonomy,bibliography,figures,chunks}.py`. See [MCP_TOOLS.md](MCP_TOOLS.md) for the full tool surface.

## Steering the client session

The server often needs to shape *how* the client uses the corpus — report vs. manuscript structure, cross-validation of synthesized claims, citation discipline — not just answer queries. MCP is **client-driven**: a server can only add to the model's context as a *response* to a client-initiated request (`initialize`, `tools/call`, `prompts/get`, `resources/read`), and can **never unilaterally push a turn** mid-session. Elicitation, sampling, and notifications do *not* inject free-form prompting into the model's context. So all steering rides on three response channels, in increasing specificity:

1. **`instructions`** — `InitializeResult.instructions`, sourced from `<corpuscle>/instructions.md`, injected once at session start. The always-on baseline: tool routing, citation rules, "respect the active output profile." Established.
2. **Tool-result guidance** — short guidance text appended to a tool's response payload, which lands in context the moment that tool is called. This is the idiomatic *just-in-time* nudge channel; it beats front-loading everything into `instructions` (which the model reads once and drifts from). `format_citation`'s provenance warning is the established example. The intended extension is to make this **output-profile-aware** (tracked in [#101](https://github.com/caseywdunn/corpus/issues/101)): e.g. under a `manuscript` profile, figure tools remind the model to verify publishability and emit the attribution string, and synthesis retrievals nudge cross-validation against multiple sources.
3. **MCP Prompts** (`prompts/list` / `prompts/get`, surfaced in Claude Code as `/mcp__corpus__<name>`) — user-invoked structural scaffolds (a manuscript skeleton, a monographic-review recipe, a cross-validation checklist). Zero context cost until invoked. Forward-looking, tracked in #101.

Two rules govern the idiom:

- **Soft steering** (document structure, cross-validation, house style) goes in the nudge channels above, keyed to the client-selected output profile. It is advisory — the model may ignore it.
- **Hard requirements** (figure publishability, citation provenance) are enforced **server-side at the tool boundary**, never delivered as a nudge the model can ignore. This is the #79/#100 lesson: a trust-critical gate must be code, not a prompt.

Tool-result nudges cost tokens on every call and pull against the served-bundle payload-trimming work (#76, #81–86), so keep them terse and conditional rather than emitting on every call.

## Key files

| File / package | Role |
|---|---|
| `pipeline/` | Stage 1 + Pass 3b/3c orchestrator (split per-stage; CLI in `pipeline/main.py`) |
| `process_corpus.py` | Thin CLI shim into `pipeline.main` (kept for backwards compatibility) |
| `update_corpus.py` | Orchestrator that runs pipeline + post-pipeline scripts in dependency order |
| `corpus_status.py` | Single-command rollup of stage completion, quality flags, staleness |
| `embed_chunks.py` | Stage 2 embedding (BGE-M3 → LanceDB) |
| `pipeline/figures.py` | Figure extraction, classification, caption parsing, chunk-figure linking |
| `pipeline/vision.py` | Vision-LLM backends (local Qwen2.5-VL, Claude API) |
| `pipeline/taxa.py` | Taxonomy DB access, taxon-mention extraction, lexicon loaders |
| `pipeline/grobid_client.py` | Grobid TEI parsing for headers, references, in-text citations |
| `pipeline/embeddings.py` | BGE-M3 embedding backends (local, HF, etc.) |
| `pipeline/external.py` | Shared retry + circuit breaker + `--strict-network` mode |
| `pipeline/version.py` | Single-source `__version__` stamped into every artifact |
| `mcpsrv/` | MCP server implementation (FastMCP, eager index, stdio + SSE transports) |
| `mcp_server.py` | Thin CLI shim into `mcpsrv.main` |
| `bib/` | BibTeX parser, importer, exporter (round-trip curation) |
| `config.yaml` | Pipeline configuration (loaded by `pipeline.config.load_config`) |
| `slurm/batch_pipeline.sh` | SLURM orchestrator: chains Grobid, Stage 1 array, cleanup, Pass 3b, Embed |
| `slurm/batch_process_corpus.sh` | SLURM Stage 1 batch script (supports job arrays) |
| `slurm/batch_pass3b.sh` | SLURM Pass 3b batch script (GPU) |
| `slurm/batch_embed.sh` | SLURM embedding batch script (GPU) |
| `slurm/bouchet_paths.sh` | Shared path definitions for all batch scripts |
| `ingest_taxonomy.py` | Build Darwin Core taxonomy SQLite from a DwC file, archive, or the WoRMS API |
| `build_biblio_authority.py` | Bibliographic authority DB (deduplicated works + citation graph) |
| `build_taxon_mentions.py` | Cross-paper taxon mentions SQLite |
| `backfill_intext_citations.py` | TEI body → `intext_citations.json` per paper |
| `reconcile_corpus_to_biblio.py` | Merge ghost cited-references onto corpus papers |
| `package_for_serve.py` | Whitelist + manifest the served bundle for S3 / EC2 deploy |
| `demo/lexicon.yaml` | Example multi-category lexicon (siphonophore anatomy under `anatomy:`) |
