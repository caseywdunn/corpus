# Contributing

Thanks for the interest. Everything below assumes the `corpus` conda environment is active — see the [README](README.md) and [INSTALL.md](INSTALL.md).

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
- If you touch the pipeline, run at least one end-to-end on the demo papers (`cd demo && corpus run`) and verify the affected ground-truth YAMLs still match.
- If you touch the MCP server surface, exercise the relevant tools against `demo/output/` (`cd demo && corpus serve`).

### Continuous integration

Four test tiers (#75); the first two run automatically in GitHub Actions, the last two are manual release-time checks.

| Tier | Trigger | Where | What it catches |
|---|---|---|---|
| **T0 — lint + unit** | every push, every branch | [`.github/workflows/lint.yml`](.github/workflows/lint.yml) | pyflakes (NameError-class bugs) + ~314 unit tests with no corpus dependency |
| **T1 — demo build + serve, Linux** | every push, every branch + every PR | [`.github/workflows/integration.yml`](.github/workflows/integration.yml) | `corpus run` on the 4-paper demo against real Grobid + LanceDB, bundle-manifest shape, audit-clean, SSE round-trip, all `corpus_required` parametrized tests, then the 4 + 1 implicit-resume scenario (copy [`tests/fixtures/round2_paper/Siebert_etal2011.pdf`](tests/fixtures/round2_paper/) into `demo/`, re-run, assert `skipped=4, embedded≥1, failed=0` — regression check for #71) |
| **T2 — demo build + serve, macOS arm64** | every push, every branch + every PR | [`.github/workflows/integration.yml`](.github/workflows/integration.yml) | same as T1 on `macos-15`, with `grobid.disable: true` (Docker Desktop isn't on GHA macOS runners) — catches macOS-specific regressions including darwin-specific LanceDB resume behavior |
| **T3 — clean-room EC2** | manual, pre-release | [`dev_docs/ec2_smoke.sh`](dev_docs/ec2_smoke.sh) | full install from absolutely nothing on Ubuntu — catches `environment.yaml` drift that warm GHA caches hide |
| **T4 — operator walkthrough** | manual, when CLI changes | [`dev_docs/clean_install_walkthrough.sh`](dev_docs/clean_install_walkthrough.sh) | every operator verb interactively (`completion`, `--cite`, `status --report`, full `bib export/import` round-trip) |

Local equivalents:

```bash
pytest -m "not corpus_required and not resume_scenario"   # T0
cd demo && corpus run --no-vision                          # T1/T2 round-1 corpus build
pytest -m corpus_required                                  #   ground-truth assertions
cp tests/fixtures/round2_paper/Siebert_etal2011.pdf demo/  # T1/T2 round-2 setup
cd demo && corpus run --no-vision                          #   resume scenario (needs Grobid + ~2 min)
```

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

1. On `dev`, set `__version__` in [pipeline/version.py](pipeline/version.py) to the release
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
7. Clean up [dev_docs/PLAN.md](dev_docs/PLAN.md): tick off items that
   shipped in this release (or strike them through), prune anything
   that's no longer the plan, and open a fresh section for the next
   version's roadmap. The CHANGELOG records what *happened*; PLAN
   records what's *next* — they shouldn't drift.

Worked example — releasing `vX.Y.Z` end-to-end from the shell. Run
each block separately and read the output; don't paste end-to-end.

```bash
# Pre-flight: tests green and the smoke walkthrough still works.
python -m pytest tests/ -q
corpus --version           # confirm `pip install -e .` is current
# (optionally) re-run dev_docs/clean_install_walkthrough.sh against a
# scratch corpus to catch CLI-affecting regressions.

# 1. On dev, bump version + date the CHANGELOG entry.
git checkout dev && git pull --ff-only
$EDITOR pipeline/version.py    # "X.Y.Z.devN"  → "X.Y.Z"
$EDITOR CHANGELOG.md           # "## [Unreleased]" → "## [X.Y.Z] - YYYY-MM-DD"
git add pipeline/version.py CHANGELOG.md
git commit -m "release: vX.Y.Z"
git push origin dev

# 2. Merge dev → main (no fast-forward so the merge commit marks the release).
git checkout main && git pull --ff-only
git merge --no-ff dev -m "Release vX.Y.Z"
git push origin main

# 3. Tag from main; the tag name MUST match __version__ modulo the leading `v`.
git tag -a vX.Y.Z -m "vX.Y.Z"
git push origin vX.Y.Z

# 4. Publish the GitHub release. --notes-from-tag would dump the tag
#    annotation; we want the CHANGELOG section instead. Use --notes-file
#    pointed at a temp file you extracted from CHANGELOG.md.
awk '/^## \[X\.Y\.Z\]/{flag=1; next} /^## \[/{flag=0} flag' CHANGELOG.md > /tmp/notes.md
gh release create vX.Y.Z --title "vX.Y.Z" --notes-file /tmp/notes.md

# 5. New major only — create the v(N-1) maintenance branch from the
#    previous major's tag so post-release fixes have a home (#48).
#    Skip for minor/patch bumps within an existing major.
# git checkout v0.last && git checkout -b v0 && git push origin v0

# 6. Back on dev, open the next pre-release.
git checkout dev
$EDITOR pipeline/version.py   # bump to "X.Y+1.0.dev0" (or major+1)
$EDITOR CHANGELOG.md          # insert a fresh "## [Unreleased]" block
git add pipeline/version.py CHANGELOG.md
git commit -m "post-release: bump dev to X.Y+1.0.dev0"
git push origin dev

# 7. Clean up the roadmap: prune shipped items and open the next
#    version's section.
$EDITOR dev_docs/PLAN.md
git add dev_docs/PLAN.md
git commit -m "plan: clear shipped X.Y.Z items, open X.Y+1.0 planning"
git push origin dev
```

### Versioning

[pipeline/version.py](pipeline/version.py) is the single source of truth for the code
version. It is read dynamically by `pyproject.toml` (so `pip show
corpus`, `corpus --version`, and the bundle manifest never drift),
imported by `mcpsrv` (surfaced via the `bundle_info` tool) and
`mcpsrv.bundle` (default for `--version`, stamped into
`bundle_manifest.json`), and gets paired with the git HEAD SHA on every
artifact. So any output tree carries the exact code that produced it.

Between releases, `dev` carries a PEP 440 pre-release suffix
(`0.3.0.dev0`, `0.3.0a1`, `0.3.0rc1`). The release commit drops the
suffix; the next commit on `dev` reintroduces one for the *next*
target version.

Schema versions for on-disk artifacts (sqlite DBs, LanceDB index) are
tracked separately from code SemVer — a code bug-fix doesn't break old
DBs, but a schema change does. Stamp schema versions inside the file
itself if you change one.

## Layout conventions

- The repo root has **zero Python files** (post-v0.3 / #60). The unified `corpus` binary is the only operator-facing entry point; it lands on PATH via `[project.scripts]` after `pip install -e .` from `pyproject.toml`. The CLI router lives at [pipeline/cli.py](pipeline/cli.py); subcommands dispatch via `python -m <module>` to private modules under `pipeline/`, `bib/`, or `mcpsrv/`.
- [pipeline/](pipeline/) — top-level CLI router (`cli.py` → the `corpus` binary), Stage 1 + Pass 3b/3c orchestrator (`main.py`, `orchestrator.py`), post-pipeline modules (`embed.py`, `taxon_mentions.py`, `intext_citations.py`, `taxonomy_ingest.py`, `status.py`), pydantic config schema + bundled template (`config_schema.py`, `config.template.yaml`), shared rich console layer (`console.py`), and supporting library modules: `scan.py` (OCR), `extract.py` (docling), `metadata.py` (Grobid + bib), `chunking.py`, `annotate.py` (taxa + lexicons), `figure_passes.py`, `runner.py` (per-paper orchestrator), `figures.py`, `taxa.py`, `grobid_client.py`, `embeddings.py`, `vision.py`, `external.py` (shared retry / circuit breaker), `version.py` (single-source `__version__`), `config.py` / `io.py` / `log.py` / `stages.py`.
- [bib/](bib/) — bibliographic round-trip + biblio authority + reconcile: `parser.py` (BibTeX parser + `BibIndex`), `export.py` (DB → BibTeX), `importer.py` (BibTeX → DB; sidecar staging for #67), `authority.py` (the build), `reconcile.py` (ghost merge). Exposed as a single namespace via `from bib import …`.
- [mcpsrv/](mcpsrv/) — MCP server package: `app.py` (FastMCP instance + index accessor), `indexes.py` (`CorpusIndex` / `TaxonMentionDB` / `BiblioAuthority`), `tools/` (papers / taxonomy / figures / chunks / bibliography), `transport.py` (bearer auth + SSE), `main.py` (argparse + dispatch), `bundle.py` (built-bundle → served-bundle distillation). Adding a new MCP tool: extend the appropriate `mcpsrv/tools/<concern>.py` module and import it from `mcpsrv/tools/__init__.py` so the `@mcp.tool()` decorator registers at startup.
- [slurm/](slurm/) — SLURM batch scripts for Bouchet; documented in [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md).
- [tools/](tools/) — developer helpers not part of the daily pipeline: QC visualizations, the MCP-launcher shell wrapper used by `.mcp.json`, one-off maintenance scripts (`dedup_ghost_works.py`, `unify_doi_corpus_key.py`).
- [templates/](templates/) — copy-and-customize starters that operators use, not pipeline inputs. Currently just the optional corpuscle-specific `instructions.md` scaffold; see [Editing client-side instructions](#editing-client-side-instructions) below.
- [tests/](tests/) — one file per subsystem.
- Per-instance data (SQLites, embeddings, per-paper artifacts) lives inside the user's *corpuscle* directory — passed as the first positional arg to every CLI — not under the repo root. See the [corpuscle layout](README.md#corpuscle-layout) in README.md. The repo no longer ships a `resources/` directory.
- [demo/](demo/) — small bundle for smoke-testing the pipeline: 11 siphonophore PDFs, a matching `siphonophores.bib`, and an example multi-category `lexicon.yaml`. The lexicon is treated as user input, parallel to `--bib` — not part of the tool.

## Editing client-side instructions

The MCP server returns a markdown blob to clients in `InitializeResult.instructions`, which well-behaved clients (Claude Desktop, Claude Code) inject into the LLM context at the start of every chat session. There are two layers, joined at server startup:

- **Defaults** — [mcpsrv/default_instructions.md](mcpsrv/default_instructions.md). Ships with the server code and is always served. Edit here for guidance that applies to *any* corpus (defer to the corpus taxonomy / bibliography, preserve historical terminology, etc.). No operator action needed for this to reach clients.
- **Per-corpuscle** — `<corpuscle>/instructions.md` (optional). Edit here for nudges specific to one corpus: a one-paragraph description of what it covers, common misconceptions to correct, domain-specific terminology cautions. Operators copy [templates/instructions.md](templates/instructions.md) as a starting scaffold. The server prepends this file to the defaults so corpus-specific guidance lands first.

Both files are paid for in tokens on every client session — keep additions terse and high-leverage. Test by restarting the MCP server (`corpus serve --output-dir <output_dir>`) and inspecting the `InitializeResult.instructions` payload, or by starting a session in any client.

## Dependencies — two files, on purpose

Two install manifests live at the repo root:

- [environment.yaml](environment.yaml) — conda env used on dev laptops and on Bouchet (YCRC). Pulls native binaries (`tesseract`, `ghostscript`) and the heavier Python deps with bundled C libs (`pymupdf`, `lxml`) from conda-forge in one shot, alongside the rest of the Python deps. This is the path most contributors hit.
- [requirements.txt](requirements.txt) — plain pip requirements used by the AWS deploy target ([DEPLOY.md §5](DEPLOY.md)), which builds a stock `python3.12 -m venv` on Ubuntu and has no conda. Native bins come from `apt`; only Python deps come from this file.

**They must stay in sync.** If you add a Python dependency, add it to *both* — to `environment.yaml` (in the conda block when there's a good conda-forge build, otherwise in the `pip:` section) and to `requirements.txt`. The only entries that legitimately appear in `environment.yaml` alone are native binaries that come from `apt` on the AWS host (`tesseract`, `ghostscript`).

## Reporting bugs + proposing changes

Open a GitHub issue with enough context to reproduce: paper hashes, exact commands, the output tree you're running against. For design changes, link or quote the relevant [dev_docs/PLAN.md](dev_docs/PLAN.md) section.
