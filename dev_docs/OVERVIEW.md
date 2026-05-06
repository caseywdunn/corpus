# Architecture Overview

The corpus pipeline processes a collection of scientific literature PDFs into a searchable, queryable knowledge base. It handles born-digital papers, scanned historical plates, multilingual text, and 19th-century Fraktur typefaces.

## Repo layout

The top-level entry-point scripts are thin shims; the implementation is grouped into packages:

| Package | Role |
|---|---|
| `pipeline/` | Stage 1 + Pass 3b/3c orchestrator. Split into `scan.py` (OCR), `extract.py` (docling), `metadata.py` (Grobid + bib), `chunking.py`, `annotate.py` (taxa + lexicons), `figure_passes.py`, `runner.py` (per-paper orchestrator), `main.py` (CLI), and supporting `config.py` / `io.py` / `log.py` / `stages.py`. |
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

Figure extraction is a multi-pass process designed for historical taxonomic literature where plates, multi-panel composites, and captions separated from figures are common.

1. **Docling extraction**: `document.pictures` with provenance (page, bbox). Each figure is classified by a transformer model (`DocumentFigureClassifier-v2.5`) into content types.
2. **Caption extraction**: Two-path approach -- first tries docling's structural caption links, then falls back to a proximity heuristic that considers same-page and facing-page text items (important for historical plate monographs where caption pages face plate pages).
3. **Figure classification** (`figures.py:classify_figure`): Categorizes each extracted region as `furniture` (headers, page numbers), `plate_label` (standalone "Plate N" text), `figure` (captioned, numbered, reasonable size), or `unclassified`.
4. **Figure-number parsing**: Multilingual regex (`_FIGURE_PREFIX`) covers Figure, Abb./Abbildung, Plate, Pl., image, illustration, and their abbreviations in English, German, French, Russian, Spanish, and Italian.
5. **PyMuPDF fallback**: When docling finds no savable figures, raw `page.get_images()` extracts embedded images with basic metadata.
6. **Deduplication**: Perceptual hashing (average hash) merges near-duplicate extractions from the two paths.
7. **Pass 3b** (vision LLM): Detects panels and compound figures that regex/heuristic methods miss.
8. **Chunk-figure linking** (`link_chunks_to_figures`): Cross-references figure numbers mentioned in text chunks with extracted figures, enabling "find figures related to this text" queries.

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

Lexicons are inputs you maintain alongside the literature, not part of the tool. Each category's content is fingerprinted independently (SHA-256 over the canonical JSON of just that section) and recorded both inside the `<category>.json` artifact and in the per-paper `pipeline_state.json` completion record. On `--resume`, editing one section re-runs `taxa_anatomy_extraction` against the new fingerprint; sections whose hash didn't change stay cached.

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

## Key files

| File / package | Role |
|---|---|
| `pipeline/` | Stage 1 + Pass 3b/3c orchestrator (split per-stage; CLI in `pipeline/main.py`) |
| `process_corpus.py` | Thin CLI shim into `pipeline.main` (kept for backwards compatibility) |
| `update_corpus.py` | Orchestrator that runs pipeline + post-pipeline scripts in dependency order |
| `corpus_status.py` | Single-command rollup of stage completion, quality flags, staleness |
| `embed_chunks.py` | Stage 2 embedding (BGE-M3 → LanceDB) |
| `figures.py` | Figure extraction, classification, caption parsing, chunk-figure linking |
| `vision.py` | Vision-LLM backends (local Qwen2.5-VL, Claude API) |
| `mcpsrv/` | MCP server implementation (FastMCP, eager index, stdio + SSE transports) |
| `mcp_server.py` | Thin CLI shim into `mcpsrv.main` |
| `bib/` | BibTeX parser, importer, exporter (round-trip curation) |
| `external.py` | Shared retry + circuit breaker + `--strict-network` mode |
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
