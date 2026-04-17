# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

A workflow for interrogating a corpus of scientific literature PDFs, spanning both scanned/old and born-digital papers. Target use cases: searching, citation analysis, and collecting anatomical figures of the same species across different papers (focused on siphonophore literature).

## Documentation

- [README.md](README.md) — installation, usage, MCP server setup, examples
- [dev_docs/OVERVIEW.md](dev_docs/OVERVIEW.md) — pipeline architecture, stage internals, figure pipeline, key files
- [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md) — HPC operational runbook (SLURM, Grobid, job arrays)
- [PLAN.md](PLAN.md) — roadmap and design decisions

## Quick reference

```bash
# Full pipeline
python process_corpus.py <input_dir> <output_dir> --resume
python embed_chunks.py <output_dir> --resume --backend local

# Parallel on Bouchet (day partition, 256 PDFs per batch)
NUM_BATCHES=8 bash batch_pipeline.sh

# Analysis
python scripts/analysis/search.py
python scripts/analysis/chat.py
```

## Implementation notes for contributors

- Each unique PDF is identified by the first 12 hex chars of its SHA-256. All artifacts live under `<output_dir>/documents/<HASH>/`. Presence of `summary.json` is the completion marker for `--resume`.
- Stage 1 (`process_corpus.py`) is CPU-only: scan detection, OCR, docling extraction, Grobid metadata, chunking. Stage 2 (`embed_chunks.py`) is GPU: BGE-M3 embeddings into LanceDB.
- Figure extraction has a fallback chain: docling pictures first, then raw PyMuPDF `page.get_images()`. Only figures that land on disk are recorded in `figures.json`.
- Visualizations overlay text-cell boxes (red) and figure bboxes (yellow/orange). Coordinates are Y-flipped from docling's bottom-left origin to PIL's top-left.
- Pipeline is resilient: most steps swallow exceptions into `processing_summary["errors"]` and continue.
- `config.yaml` exists but is only partially wired — many values are still hard-coded in scripts. Edit scripts directly or plumb config through.
- Legacy Snakemake workflow ([Snakefile](Snakefile) + [scripts/](scripts/)) is kept for the `demo/` directory only. New work goes through `process_corpus.py`.
