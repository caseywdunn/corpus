# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

A workflow for interrogating a corpus of scientific literature PDFs, spanning both scanned/old and born-digital papers. Target use cases: searching, citation analysis, and collecting anatomical figures of the same species across different papers (focused on siphonophore literature).

## Environment setup

System dependencies (macOS via Homebrew) are required before Python deps:

```
brew install ghostscript tesseract pngquant jbig2enc
conda create --name corpus python=3.13 && conda activate corpus
pip install -r requirements.txt
```

An OpenAI API key (loaded via `.env` / `python-dotenv`) is required for the embedding step.

## Commands

Main pipeline (recommended path — supersedes Snakemake):

```bash
# Process all PDFs in an input tree; output uses hash-based layout
python process_corpus.py <input_dir> <output_dir>

# Resume, skipping hashes that already have summary.json
python process_corpus.py <input_dir> <output_dir> --resume

# Embed chunks into LanceDB (separate step — expensive, needs OpenAI key)
python embed_chunks.py <output_dir>
python embed_chunks.py <output_dir> --pdf-hash <HASH>   # single doc
python embed_chunks.py <output_dir> --resume
```

Analysis & QC:

```bash
python scripts/analysis/search.py                                     # vector search
python scripts/analysis/chat.py                                       # chat over corpus
python scripts/qc/visualize_docling_parse.py <pdf> --level <char|word|line>
```

Legacy Snakemake workflow (only for the `demo/` directory, kept for reference):

```bash
snakemake --cores 1 --scheduler greedy
# or the fallback wrapper that tries multiple schedulers:
python run_workflow.py
```

## Architecture

### Hash-based, git-like organization

Each unique PDF is identified by the first 8 hex chars of its SHA-256 (`calculate_pdf_hash` in [process_corpus.py](process_corpus.py)). All artifacts for that PDF live under `<output_dir>/documents/<HASH>/`. This gives:

- **Content-addressed dedup**: identical PDFs found at multiple paths are processed once; every discovered path is recorded in that hash's `summary.json` under `relative_paths`.
- **Idempotent resume**: presence of `summary.json` in a hash directory is the completion marker checked by `--resume`.

### Two-stage pipeline split

Heavy CPU work (OCR, docling parsing, figure extraction) is separate from API-cost work (embeddings):

1. `process_corpus.py` — discovery → scan detection → OCR → docling extraction → metadata → chunking. Writes all `*.json` artifacts and figure PNGs. No network required.
2. `embed_chunks.py` — reads `chunks.json` per hash, calls OpenAI `text-embedding-3-small` (1536-dim), writes into `<output_dir>/vector_db/lancedb/` with a per-hash `<HASH>_embedded.done` marker.

This split is deliberate: stage 1 can be run on many machines / subdirectories and the results merged; stage 2 is gated on OpenAI credits and can be run overnight.

### Stage 1 internals ([process_corpus.py](process_corpus.py))

`run_pdf_processing_pipeline` orchestrates five steps per PDF, each writing a dedicated JSON under the hash dir:

| Step | Function | Output |
|---|---|---|
| Scan detection | `detect_scan_type` (PyMuPDF, text density heuristic) | `scan_detection.json` |
| OCR / copy | `prepare_pdf` (shells to `ocrmypdf` if scanned) | `processed.pdf` |
| Text + figure extraction | `extract_docling_content` | `text.json`, `figures.json`, `figures/*.png`, `visualizations/*.png` |
| Metadata | `extract_metadata` | `metadata.json` *(currently a placeholder — real Grobid integration is TODO)* |
| Chunking | `chunk_text` | `chunks.json` *(currently naive 1000-char split — not tiktoken-aware despite config.yaml saying otherwise)* |

Notable details:

- **Figure extraction has a fallback chain**: docling's `document.pictures` is tried first (with picture classification + caption extraction); if nothing savable comes back, it falls back to raw PyMuPDF `page.get_images()`. Only figures that actually land on disk are recorded in `figures.json` — see `append_figure` guard in [process_corpus.py:373](process_corpus.py#L373).
- **Visualizations overlay text-cell boxes (red) and figure bboxes (yellow/orange)** via `create_cell_visualizations`. Coordinates must be Y-flipped from docling's bottom-left origin to PIL's top-left — relevant if you touch that function.
- **Pipeline is resilient, not strict**: most steps swallow exceptions into `processing_summary["errors"]` and continue. Missing docling still produces a placeholder `text.json` so downstream steps don't crash.

### Stage 2 internals ([embed_chunks.py](embed_chunks.py))

Uses LanceDB with two Pydantic models, `ChunkMetadata` and `DocumentChunk` (with a `Vector(1536)` field). Embeddings are batched (100 per call). On API failure, zero-vectors are used as a fallback — watch for this when debugging bad search results.

### Config

[config.yaml](config.yaml) exists and describes chunking/embedding/OCR settings, **but `process_corpus.py` largely does not read it** — values like OCR flags, chunk size, and embedding model are currently hard-coded. If modifying behavior, edit the script directly or plumb config through first.

### Legacy Snakemake workflow

[Snakefile](Snakefile) + [scripts/](scripts/) implement the same pipeline as per-file rules, hard-wired to `DEMO_DIR = "demo"` and `OUTPUT_DIR = "output"` with filename-based (not hash-based) organization. Kept for the demo, but new work should go through `process_corpus.py`. The two systems write to different output layouts and are not interchangeable.
