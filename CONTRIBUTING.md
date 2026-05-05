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

## Branching model

### Branches

- `main` — tagged releases only. Never commit directly to `main`.
- `dev` — active development for the next unreleased version. All feature
  work targets this branch.
- `vN` (e.g. `v0`, `v1`) — release maintenance branches, created from the
  release tag when a major version ships. Used only for patches to released
  versions.
- Issue branches — created from `dev` (or from `vN` for patches), named after
  the issue they address (e.g. `issue-17`).

### Normal workflow

1. Create an issue branch from `dev` (`git checkout -b issue-NNN dev`)
2. Do the work; ensure the structural test command and any affected
   ground-truth tests pass (see [What to run before opening a PR](#what-to-run-before-opening-a-pr)).
3. Merge to `dev` and push (`git checkout dev && git merge issue-NNN && git push`)
4. Delete the issue branch
5. When ready to release: merge `dev` to `main`, tag the release, create
   a `vN` branch from the tag

Always merge and push completed issue branches to `dev` before starting the
next issue. This keeps `dev` up to date and avoids dependency tangles when
later issue branches need earlier work.

### Patching a released version

If a bug is found in v0.1 while v0.2 is in development on `dev`:

1. Create an issue branch from `v0`
2. Fix the bug
3. Merge to `v0`, tag `v0.1.1`
4. Merge `v0` to `main`
5. Cherry-pick the fix into `dev` so the next version gets it too

### Releasing

1. On `dev`, set `__version__` in [version.py](version.py) to the release
   string (e.g. `"0.2.0"`) and confirm `CHANGELOG.md` has a dated entry
   for it. Commit.
2. Merge `dev` → `main`, push.
3. Tag: `git tag -a vX.Y.Z -m "vX.Y.Z" && git push origin vX.Y.Z`. The tag
   name and `__version__` must match (modulo the `v` prefix).
4. Create the GitHub release from the tag (CHANGELOG entry as release notes).
5. Create the `vN` maintenance branch from the tag if this is a new major.
6. On `dev`, bump `__version__` to the next pre-release (e.g.
   `"0.3.0.dev0"` — PEP 440 suffix) and open a new `## [Unreleased]`
   section in `CHANGELOG.md`.

### Versioning

[version.py](version.py) is the single source of truth for the code
version. It is imported by `mcp_server.py` (surfaced via the
`bundle_info` tool) and `package_for_serve.py` (default for
`--version`, stamped into `bundle_manifest.json`), and gets paired with
the git HEAD SHA on every artifact. So any output tree carries the
exact code that produced it.

Between releases, `dev` carries a PEP 440 pre-release suffix
(`0.3.0.dev0`, `0.3.0a1`, `0.3.0rc1`). The release commit drops the
suffix; the next commit on `dev` reintroduces one for the *next*
target version.

Schema versions for on-disk artifacts (sqlite DBs, LanceDB index) are
tracked separately from code SemVer — a code bug-fix doesn't break old
DBs, but a schema change does. Stamp schema versions inside the file
itself if you change one.

## Layout conventions

- Code at the repo root is either a library module imported by the pipeline (`embeddings.py`, `external.py`, `figures.py`, `grobid_client.py`, `taxa.py`, `vision.py`) or a CLI entry point the user runs directly (`process_corpus.py`, `embed_chunks.py`, `mcp_server.py`, `build_biblio_authority.py`, `build_taxon_mentions.py`, `ingest_taxonomy.py`, `reconcile_corpus_to_biblio.py`, `corpus_status.py`, `update_corpus.py`, `backfill_intext_citations.py`, `package_for_serve.py`, `bib_export.py`, `bib_import.py`).
- [bib/](bib/) — bibliographic round-trip package: `parser.py` (BibTeX parser, `BibIndex`), `export.py` (DB → BibTeX), `importer.py` (BibTeX → DB). Exposed as a single namespace via `from bib import …`. The root `bib_export.py` / `bib_import.py` are thin CLI shims; `bib_metadata.py` at root is a backwards-compat re-export shim slated for removal.
- [slurm/](slurm/) — SLURM batch scripts for Bouchet; documented in [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md).
- [tools/](tools/) — developer helpers not part of the daily pipeline: QC visualizations, the MCP-launcher shell wrapper used by `.mcp.json`.
- [tests/](tests/) — one file per subsystem.
- Per-instance data (SQLites, embeddings, per-paper artifacts) lives inside the user's *corpuscle* directory — passed as the first positional arg to every CLI — not under the repo root. See the [corpuscle layout](README.md#corpuscle-layout) in README.md. The repo no longer ships a `resources/` directory.
- [demo/](demo/) — small bundle for smoke-testing the pipeline: 11 siphonophore PDFs, a matching `siphonophores.bib`, and an example `anatomy_lexicon.yaml`. The lexicon is treated as user input, parallel to `--bib` — not part of the tool.

## Dependencies — two files, on purpose

Two install manifests live at the repo root:

- [environment.yaml](environment.yaml) — conda env used on dev laptops and on Bouchet (YCRC). Pulls native binaries (`tesseract`, `ghostscript`) and the heavier Python deps with bundled C libs (`pymupdf`, `lxml`) from conda-forge in one shot, alongside the rest of the Python deps. This is the path most contributors hit.
- [requirements.txt](requirements.txt) — plain pip requirements used by the AWS deploy target ([dev_docs/DEPLOY.md §5](dev_docs/DEPLOY.md)), which builds a stock `python3.12 -m venv` on Ubuntu and has no conda. Native bins come from `apt`; only Python deps come from this file.

**They must stay in sync.** If you add a Python dependency, add it to *both* — to `environment.yaml` (in the conda block when there's a good conda-forge build, otherwise in the `pip:` section) and to `requirements.txt`. The only entries that legitimately appear in `environment.yaml` alone are native binaries that come from `apt` on the AWS host (`tesseract`, `ghostscript`).

## Reporting bugs + proposing changes

Open a GitHub issue with enough context to reproduce: paper hashes, exact commands, the output tree you're running against. For design changes, link or quote the relevant [PLAN.md](PLAN.md) section.
