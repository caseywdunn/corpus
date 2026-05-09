# AGENTS.md

Orientation for AI coding agents working in this repository.

## Purpose

A workflow for interrogating a corpus of scientific literature PDFs, spanning both scanned/old and born-digital papers. Target use cases: searching, citation analysis, and collecting figures of the same species across different papers. Corpus-agnostic by design — the group-specific data (taxonomy snapshot, lexicon, BibTeX, instructions) lives in the per-instance corpuscle directory; the lab's reference deployment is a siphonophore corpus, but every other group plugs into the same machinery.

## Documentation map

- [README.md](README.md) — installation, usage, MCP server setup, examples
- [CONTRIBUTING.md](CONTRIBUTING.md) — branching model, release ritual, version bumps
- [CHANGELOG.md](CHANGELOG.md) — what changed in each release
- [dev_docs/PLAN.md](dev_docs/PLAN.md) — roadmap and design decisions for the active version
- [dev_docs/OVERVIEW.md](dev_docs/OVERVIEW.md) — pipeline architecture, stage internals, figure pipeline, key files
- [dev_docs/MCP_TOOLS.md](dev_docs/MCP_TOOLS.md) — full MCP tool surface (27 tools)
- [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md) — HPC operational runbook (SLURM, Grobid, job arrays)
- [dev_docs/DEPLOY.md](dev_docs/DEPLOY.md) — AWS deploy runbook (S3 bundle + EC2 systemd)
- [dev_docs/TESTING.md](dev_docs/TESTING.md) — quality test suite, ground truth format, evaluation workflow
- [INSTALL.md](INSTALL.md) — optional OCR extras, pip-only fallback, platform notes

## Quick reference

```bash
# Full pipeline + post-pipeline scripts in dependency order, with --resume throughout.
python update_corpus.py <input_dir> <output_dir> --resume

# Or run individual stages — equivalent expanded form:
python process_corpus.py <input_dir> <output_dir> --resume
python embed_chunks.py <output_dir> --resume
python build_biblio_authority.py <output_dir>
python build_taxon_mentions.py <output_dir>
python backfill_intext_citations.py <output_dir>
python reconcile_corpus_to_biblio.py <output_dir>

# Status / staleness rollup
python corpus_status.py <output_dir>

# Parallel Stage 1 on Bouchet (day partition, 256 PDFs per batch)
NUM_BATCHES=8 bash slurm/batch_pipeline.sh

# MCP server
python mcp_server.py <output_dir>
```

## Repo layout

| Path | Role |
| --- | --- |
| `pipeline/` | Stage 1 + Pass 3b/3c orchestrator and shared library modules. CLI in `pipeline/main.py`; per-stage modules: `scan.py`, `extract.py`, `metadata.py`, `chunking.py`, `annotate.py`, `figure_passes.py`, `runner.py`. Library modules: `figures.py`, `taxa.py`, `grobid_client.py`, `embeddings.py`, `vision.py`, `external.py`, `version.py`. |
| `mcpsrv/` | MCP server. FastMCP `app.py` + 27 tools in `tools/{papers,taxonomy,bibliography,figures,chunks}.py`. |
| `bib/` | BibTeX parser, importer, exporter (round-trip curation). |
| `slurm/` | SLURM batch scripts (Bouchet). |
| `deploy/` | CloudFormation, nginx, systemd, sync + update scripts. |
| `tests/` | Ground-truth + corpus-wide consistency tests. |
| `process_corpus.py`, `mcp_server.py` | Thin CLI shims into `pipeline.main` / `mcpsrv.main`. |
| `update_corpus.py` | Dependency-ordered orchestrator over the pipeline + post-pipeline scripts. |
| `corpus_status.py` | Single-command rollup of stage completion, quality flags, staleness. |

## Implementation notes for contributors

- Each unique PDF is identified by the first 12 hex chars of its SHA-256. All artifacts live under `<output_dir>/documents/<HASH>/`. Per-paper `pipeline_state.json` records each stage's `pipeline_version` + `input_fingerprint`; on `--resume` a stage re-runs iff its record is missing or disagrees with the current input. A lexicon edit only re-runs the annotation pass.
- Stage 1 (`pipeline/`) is CPU-only: scan detection, OCR, docling extraction, Grobid metadata, chunking, annotation. Stage 2 (`embed_chunks.py`) is GPU: BGE-M3 embeddings into LanceDB. Pass 3b (vision LLM) is opt-in via `--vision-backend {claude,local}`.
- Figure extraction has a fallback chain: docling pictures first, then raw PyMuPDF `page.get_images()`. Only figures that land on disk are recorded in `figures.json`.
- Visualizations overlay text-cell boxes (red) and figure bboxes (yellow/orange). Coordinates are Y-flipped from docling's bottom-left origin to PIL's top-left.
- Annotation artifacts (`taxa.json` plus one `<category>.json` per category in the `--lexicon` YAML) stamp an `input_fingerprint`; per-stage resume detects mismatch and re-runs the annotation pass automatically.
- Pipeline failures land in structured `summary.json["stage_failures"]` with reason codes (`timeout`, `crash`, `external_unavailable`, `unsupported_format`, `corrupted`, `quality_gate`, `too_large`); silent-failure quality gates emit `quality_flags`. Free-text `errors[]` is being phased out.
- `config.yaml` is loaded by `pipeline.config.load_config`; CLI flags override config values, which override built-in defaults.
- External calls (Grobid, BHL, CrossRef, OpenAlex) flow through `pipeline/external.py` for shared retry + backoff + circuit breaker. `--strict-network` aborts on the first transient failure for release builds.
- `__version__` in `pipeline/version.py` is the single source of truth; CONTRIBUTING.md covers the branching model.
