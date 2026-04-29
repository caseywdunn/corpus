# Architecture Overview

The corpus pipeline processes a collection of scientific literature PDFs (focused on siphonophore taxonomy spanning the late 18th century to present) into a searchable, queryable knowledge base. It handles born-digital papers, scanned historical plates, multilingual text (English, German, French, Russian, Latin), and Fraktur typefaces.

## Content-addressed storage

Every unique PDF is identified by the first 12 hex characters of its SHA-256 hash. All artifacts for that PDF live under `<output_dir>/documents/<HASH>/`. Identical PDFs found at multiple input paths are processed once; every discovered path is recorded in `summary.json`. The presence of `summary.json` is the completion marker that `--resume` checks.

```
output/
  documents/
    a1b2c3d4e5f6/
      processed.pdf          # OCR'd copy (if scanned)
      docling_doc.json       # full docling parse
      scan_detection.json    # born_digital | scanned | broken_text_layer
      text.json              # extracted text (markdown)
      figures.json           # figure metadata, captions, classifications
      figures/               # extracted figure PNGs
      visualizations/        # per-page QC overlays (red=words, yellow=figures)
      metadata.json          # Grobid bibliographic metadata + references
      chunks.json            # text chunks for embedding
      pipeline.log           # per-document processing log
      summary.json           # provenance + processing status
    ...
  vector_db/
    lancedb/                 # LanceDB vector index
    <HASH>_embedded.done     # per-hash embedding completion markers
```

## Pipeline stages

The pipeline is split into independent stages so that CPU-heavy work, GPU work, and API-cost work can run on different hardware and be resumed independently.

### Stage 1: CPU processing (`process_corpus.py`)

Orchestrated by `run_pdf_processing_pipeline`, five steps per PDF:

| Step | What happens | Output |
|---|---|---|
| **Scan detection** | PyMuPDF text-density heuristic classifies as `born_digital`, `scanned`, or `broken_text_layer` | `scan_detection.json` |
| **PDF preparation** | Copies born-digital PDFs as-is; runs `ocrmypdf` on scanned/broken PDFs with language-appropriate Tesseract models (eng, deu, fra, rus, lat, deu_latf for Fraktur) | `processed.pdf` |
| **Text + figure extraction** | Docling parses the PDF into structured text and figure regions. Figures go through a classification/caption pipeline (see [Figure pipeline](#figure-pipeline) below). Falls back to raw PyMuPDF image extraction when docling finds nothing. | `text.json`, `figures.json`, `figures/*.png`, `visualizations/*.png` |
| **Metadata extraction** | Grobid extracts title, authors, year, DOI, abstract, section structure, and parsed references. Falls back to placeholder when Grobid is unavailable. | `metadata.json` |
| **Chunking** | Splits extracted text into overlapping chunks, annotated with section class, taxon mentions (via the configured DwC taxonomy snapshot), and anatomy-lexicon matches | `chunks.json` |

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

- **Default backend**: local sentence-transformers with [BGE-M3](https://huggingface.co/BAAI/bge-m3) (1024-dim), runs on CUDA/MPS/CPU
- **Alternative**: OpenAI `text-embedding-3-small` (1536-dim, requires API key)
- Embeddings are batched (100 per call). Per-hash completion is marked by `<HASH>_embedded.done`.
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

A user-supplied anatomy lexicon (passed via `--anatomy-lexicon path/to/lexicon.yaml`) similarly tags chunks with anatomical terms specific to the group under study. The lexicon is treated as input data parallel to `--bib`, not part of the tool — see [demo/anatomy_lexicon.yaml](../demo/anatomy_lexicon.yaml) for an example covering siphonophore terms (nectophore, pneumatophore, gastrozooid, etc.).

## MCP server (`mcp_server.py`)

Exposes the processed corpus as an MCP (Model Context Protocol) server that LLM clients can query. The server is a read-only view over per-paper artifacts; it does not store data of its own. See the [tool surface table in README.md](../README.md#tool-surface-14-tools) for the full list of 14+ tools.

## Key files

| File | Role |
|---|---|
| `process_corpus.py` | Stage 1 + Pass 3b/3c orchestrator |
| `embed_chunks.py` | Stage 2 embedding |
| `figures.py` | Figure extraction, classification, caption parsing, chunk-figure linking |
| `vision.py` | Vision-LLM backends (local Qwen, Claude API) |
| `mcp_server.py` | MCP server for LLM-driven corpus queries |
| `config.yaml` | Configuration (partially wired -- some values still hard-coded in scripts) |
| `slurm/batch_pipeline.sh` | SLURM orchestrator: chains Grobid, Stage 1 array, cleanup, Pass 3b, Embed |
| `slurm/batch_process_corpus.sh` | SLURM Stage 1 batch script (supports job arrays) |
| `slurm/batch_pass3b.sh` | SLURM Pass 3b batch script (GPU) |
| `slurm/batch_embed.sh` | SLURM embedding batch script (GPU) |
| `slurm/bouchet_paths.sh` | Shared path definitions for all batch scripts |
| `ingest_taxonomy.py` | Build Darwin Core taxonomy SQLite from a DwC file, archive, or the WoRMS API |
| `demo/anatomy_lexicon.yaml` | Example siphonophore anatomy lexicon (user-supplied input pattern) |
