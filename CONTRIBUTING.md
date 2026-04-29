# Contributing

Thanks for the interest. Everything below assumes the `corpus` conda environment is active — see the [README](README.md) and [dev_docs/INSTALL.md](dev_docs/INSTALL.md).

## Tests

The suite is in [tests/](tests/). Two kinds of check run under the same pytest invocation:

- **Structural / unit tests** — fast, no corpus output required. Cover the bibliography cascade, reconciliation logic, reference parsing, text extraction helpers, figure-extraction helpers.
- **Ground-truth + corpus-wide tests** — parametrized over YAML files in [tests/ground_truth/](tests/ground_truth/) and over every processed paper in the corpus output. Require `CORPUS_OUTPUT_DIR` to point at a processed output tree (or an `output/` dir in the repo root).

Full detail on layouts, fixtures, and how to add a new ground-truth paper: [dev_docs/TESTING.md](dev_docs/TESTING.md).

### Running

```bash
# Structural tests only — fast, no corpus needed
python -m pytest tests/test_biblio_cascade.py tests/test_reconcile.py \
    tests/test_reference_extraction.py tests/test_text_extraction.py \
    tests/test_figure_extraction.py -q

# Ground-truth + corpus-wide (needs a processed output dir)
CORPUS_OUTPUT_DIR=/path/to/output python -m pytest tests/ -q

# A single file, verbose
python -m pytest tests/test_biblio_cascade.py -v
```

### What to run before opening a PR

- The structural test command above must pass (61 tests, sub-second).
- If you touch the pipeline, run at least one end-to-end on the demo papers (`python process_corpus.py demo demo_output`) and verify the affected ground-truth YAMLs still match.
- If you touch the MCP server surface, exercise the relevant tools against `demo_output/`.

## Commit style

- Short, present-tense subject line (<70 chars). Body explains the *why*, not the *what* — the diff covers what.
- Reference GitHub issues as `Closes #N` in the body when a commit resolves one.
- Don't squash unrelated changes — a cleanup and a bug fix belong in separate commits. See recent `main` history for the house style.

## Layout conventions

- Code at the repo root is either a library module imported by the pipeline (`embeddings.py`, `figures.py`, `grobid_client.py`, `taxa.py`, `vision.py`) or a CLI entry point the user runs directly (`process_corpus.py`, `embed_chunks.py`, `mcp_server.py`, `build_biblio_authority.py`, `build_taxon_mentions.py`, `ingest_worms.py`, `reconcile_corpus_to_biblio.py`).
- [slurm/](slurm/) — SLURM batch scripts for Bouchet; documented in [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md).
- [tools/](tools/) — developer helpers not part of the daily pipeline: QC visualizations, the MCP-launcher shell wrapper used by `.mcp.json`.
- [tests/](tests/) — one file per subsystem.
- [resources/](resources/) — regenerable precompiled data (WoRMS SQLite, bibliography authority DB, taxon mentions DB). The SQLite files are `.gitignore`d; the YAML lexicon is checked in.

## Dependencies — two files, on purpose

Two install manifests live at the repo root:

- [environment.yaml](environment.yaml) — conda env used on dev laptops and on Bouchet (YCRC). Pulls native binaries (`tesseract`, `ghostscript`) and the heavier Python deps with bundled C libs (`pymupdf`, `lxml`) from conda-forge in one shot, alongside the rest of the Python deps. This is the path most contributors hit.
- [requirements.txt](requirements.txt) — plain pip requirements used by the AWS deploy target ([dev_docs/DEPLOY.md §5](dev_docs/DEPLOY.md)), which builds a stock `python3.12 -m venv` on Ubuntu and has no conda. Native bins come from `apt`; only Python deps come from this file.

**They must stay in sync.** If you add a Python dependency, add it to *both* — to `environment.yaml` (in the conda block when there's a good conda-forge build, otherwise in the `pip:` section) and to `requirements.txt`. The only entries that legitimately appear in `environment.yaml` alone are native binaries that come from `apt` on the AWS host (`tesseract`, `ghostscript`).

## Reporting bugs + proposing changes

Open a GitHub issue with enough context to reproduce: paper hashes, exact commands, the output tree you're running against. For design changes, link or quote the relevant [PLAN.md](PLAN.md) section.
