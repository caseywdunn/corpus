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

[![T0 lint + unit (dev)](https://github.com/caseywdunn/corpus/actions/workflows/lint.yml/badge.svg?branch=dev)](https://github.com/caseywdunn/corpus/actions/workflows/lint.yml?query=branch%3Adev)
[![T1/T2 integration (dev)](https://github.com/caseywdunn/corpus/actions/workflows/integration.yml/badge.svg?branch=dev)](https://github.com/caseywdunn/corpus/actions/workflows/integration.yml?query=branch%3Adev)

Badges above reflect the `dev` branch — the active development line — so contributors see immediately whether the next release is green. The README's badges are scoped to `main` for the release-state signal users care about.

Five test tiers (#75); the first four run automatically in GitHub Actions — T0/T1/T2 on every push, T3 on a weekly schedule — and the last two are manual release-time checks.

**On caching, because it has misled us once:** `setup-miniconda` does *not* cache the solved conda environment, and this repo adds no `actions/cache` for it. Every T0/T1/T2 run re-resolves `environment.yaml` from scratch; the only caches are integration.yml's HuggingFace model-hub caches. So CI is a genuine dependency-drift detector — its weakness is *timeliness*, not blindness. `mcp` 2.0.0 broke every clean install and went unnoticed for 24 days purely because nobody pushed to `dev` between 2026-07-06 and 2026-07-30 (#156). T3's `schedule:` trigger is what closes that gap.

| Tier | Trigger | Where | What it catches |
|---|---|---|---|
| **T0 — lint + unit** | every push, every branch | [`.github/workflows/lint.yml`](.github/workflows/lint.yml) | pyflakes (NameError-class bugs) + ~314 unit tests with no corpus dependency |
| **T1 — demo build + serve, Linux** | every push, every branch + every PR | [`.github/workflows/integration.yml`](.github/workflows/integration.yml) | `corpus run` on the 4-paper demo against real Grobid + LanceDB, bundle-manifest shape, audit-clean, SSE round-trip, all `corpus_required` parametrized tests, then the 4 + 1 implicit-resume scenario (copy [`tests/fixtures/round2_paper/Siebert_etal2011.pdf`](tests/fixtures/round2_paper/) into `demo/`, re-run, assert `skipped=4, embedded≥1, failed=0` — regression check for #71) |
| **T2 — demo build + serve, macOS arm64** | every push, every branch + every PR | [`.github/workflows/integration.yml`](.github/workflows/integration.yml) | same as T1 on `macos-15`, with `grobid.disable: true` (Docker Desktop isn't on GHA macOS runners) — catches macOS-specific regressions including darwin-specific LanceDB resume behavior |
| **T3 — clean-room install** | **weekly (Mon 06:00 UTC)**, `workflow_dispatch`, **any PR targeting `main`**, or a push touching the lane itself | [`.github/workflows/clean-room.yml`](.github/workflows/clean-room.yml) | the same install → demo → serve path as T1 but on a *clock* rather than a push, and with the HuggingFace cache deliberately disabled so a genuine first-run model download is exercised (the path that 429'd in #140). Also drives the real `docker-compose.yml`. This is the [PLAN.md §0](dev_docs/PLAN.md) release gate. **Note:** `schedule` and `workflow_dispatch` only work for workflows on the default branch, so on a feature branch this lane is unregistered (`gh workflow run` → 404) — that is why it also triggers on PRs to `main`, which is where the release gate actually needs it |
| **T3-bare — clean-room EC2** | manual, pre-release | [`dev_docs/ec2_smoke.sh`](dev_docs/ec2_smoke.sh) | what T3 can't cover: the bare-host bootstrap (apt, miniforge install) from absolutely nothing on a real Ubuntu EC2 instance. Criteria in [`dev_docs/PLATFORM_SMOKE.md`](dev_docs/PLATFORM_SMOKE.md) |
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

1. **Confirm CI on `dev` is green.** Every per-push lane must be ✓ on
   the latest `dev` commit (`gh run list --branch dev --limit 4`, or
   the [Actions tab](https://github.com/caseywdunn/corpus/actions)).
   A red `dev` means the release merge into `main` will be red too —
   fix the cause before proceeding. The README badges are scoped to
   `main`, so anything red here gets stamped on the released version.

   **Confirm the Zenodo webhook is armed**, while there is still time to
   fix it:

   ```bash
   gh api repos/caseywdunn/corpus/hooks \
     --jq '.[] | select(.config.url | contains("zenodo")) |
           "active=\(.active) events=\(.events|join(","))"'
   ```

   Expect `active=true events=release`. No output at all means the repo
   is toggled **off** in Zenodo and the release will not be archived.
   Zenodo silently missed v0.3.0 through v1.0.0 — six releases — because
   nothing checked this and a missing archive looks exactly like a
   successful release from GitHub's side. Re-arm by toggling the repo
   off and on at https://zenodo.org/account/settings/github/, which
   recreates the hook; GitHub then sends a `ping` that Zenodo should
   answer `202`.
2. On `dev`, set `__version__` in [pipeline/version.py](pipeline/version.py) to the release
   string (e.g. `"0.2.0"`) and confirm `CHANGELOG.md` has a dated entry
   for it. Bump the release-pinned install line in
   [INSTALL.md](INSTALL.md) (`pip install git+…@vX.Y.Z`) to the same tag —
   it names a tag that only exists once you reach step 4, so nothing
   catches it going stale, and it has now drifted twice
   (`grep -rn '\.git@v' INSTALL.md`). Commit.
3. **Open a pull request `dev` → `main`** and wait for every check,
   then merge. Do **not** push `main` directly: T3's `push` trigger is
   path-filtered to the clean-room workflow file, so a direct merge runs
   T0/T1/T2 and silently skips **T3 — clean-room install**, which is the
   [PLAN.md §0](dev_docs/PLAN.md) release gate. `main` goes green either
   way, so the omission is invisible. A PR targeting `main` is the only
   trigger that runs T3 on a release that does not happen to edit that
   workflow — v1.0.0 fired it on push only because that merge *added*
   the lane.

   After merging, wait for CI on `main` to go green as well
   (`gh run watch`, or poll `gh run list --branch main`) before
   tagging — once you tag, the version-stamped artifact is fixed; you
   don't want to be the person who released vX.Y.Z and then immediately
   tagged vX.Y.Z+1 to fix a CI failure that the merge surfaced.
4. Tag: `git tag -a vX.Y.Z -m "vX.Y.Z" && git push origin vX.Y.Z`. The tag
   name and `__version__` must match (modulo the `v` prefix).
5. Create the GitHub release from the tag (CHANGELOG entry as release notes).
   **This — not the tag — is what triggers Zenodo.** Verify within a
   minute or two, while the payload is still redeliverable:

   ```bash
   gh api repos/caseywdunn/corpus/hooks/<ID>/deliveries \
     --jq '.[] | "\(.delivered_at) \(.event)/\(.action) code=\(.status_code)"' | head -5
   ```

   One release produces three deliveries — `created`, `published`,
   `released`. Zenodo accepts the first with **202** and answers the
   other two **409**; that pattern is success, not failure. What you are
   checking for is a 202 on at least one of them. A 4xx across all three
   means the hook is stale — fix it and redeliver rather than cutting
   another release:

   ```bash
   gh api -X POST repos/caseywdunn/corpus/hooks/<ID>/deliveries/<DELIVERY_ID>/attempts
   ```

   The record itself takes a few minutes to appear; confirm with
   `curl -sL -H 'Accept: application/json' https://zenodo.org/api/records/19964909`
   and check that `metadata.version` is the tag you just cut. The
   concept DOI in `CITATION.cff` always resolves to the newest archived
   version, so if Zenodo misses a release every citation silently points
   at an older one.
6. Create the `vN` maintenance branch from the tag if this is a new major.
7. On `dev`, bump `__version__` to the next pre-release (e.g.
   `"0.3.0.dev0"` — PEP 440 suffix) and open a new `## [Unreleased]`
   section in `CHANGELOG.md`.
8. Clean up [dev_docs/PLAN.md](dev_docs/PLAN.md): tick off items that
   shipped in this release (or strike them through), prune anything
   that's no longer the plan, and open a fresh section for the next
   version's roadmap. The CHANGELOG records what *happened*; PLAN
   records what's *next* — they shouldn't drift.

Worked example — releasing `vX.Y.Z` end-to-end from the shell. Run
each block separately and read the output; don't paste end-to-end.

```bash
# Pre-flight A: CI on dev is green (both T0 and T1/T2).
gh run list --branch dev --limit 4   # expect ✓ on the top two rows

# Pre-flight B: local tests green and the smoke walkthrough still works.
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

# 2. Open a PR dev → main. Required: a PR is the ONLY trigger that runs
#    T3 (clean-room install) on a release that doesn't edit that
#    workflow. Merging main directly skips the release gate silently.
gh pr create --base main --head dev --title "Release vX.Y.Z" --body "..."
gh pr checks <N> --watch          # all lanes, T3 included
gh pr merge <N> --merge           # merge commit; do NOT --delete-branch (dev is permanent)

# 2a. WAIT FOR CI ON MAIN to be green before tagging — the tag pins
#     the version-stamped artifact in the bundle manifest, and a red
#     main badge on the README would persist until the next release.
gh run watch   # or: gh run list --branch main --limit 4

# 3. Tag from main; the tag name MUST match __version__ modulo the leading `v`.
git tag -a vX.Y.Z -m "vX.Y.Z"
git push origin vX.Y.Z

# 4. Publish the GitHub release. --notes-from-tag would dump the tag
#    annotation; we want the CHANGELOG section instead. Use --notes-file
#    pointed at a temp file you extracted from CHANGELOG.md.
awk '/^## \[X\.Y\.Z\]/{flag=1; next} /^## \[/{flag=0} flag' CHANGELOG.md > /tmp/notes.md
gh release create vX.Y.Z --title "vX.Y.Z" --notes-file /tmp/notes.md --verify-tag

# 4a. VERIFY ZENODO TOOK IT — the release, not the tag, is the trigger.
#     Three deliveries per release: 202 on one and 409 on the others is
#     success (409 = duplicate suppression). All-4xx means a stale hook.
gh api repos/caseywdunn/corpus/hooks/<ID>/deliveries \
  --jq '.[] | "\(.delivered_at) \(.event)/\(.action) code=\(.status_code)"' | head -5

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
version. It is read dynamically by `pyproject.toml` (so `corpus
--version` and the bundle manifest never drift),
imported by `mcpsrv` (surfaced via the `bundle_info` tool) and
`mcpsrv.bundle` (default for `--version`, stamped into
`bundle_manifest.json`), and gets paired with the git HEAD SHA on every
artifact. So any output tree carries the exact code that produced it.

Between releases, `dev` carries a PEP 440 pre-release suffix
(`0.3.0.dev0`, `0.3.0a1`, `0.3.0rc1`). The release commit drops the
suffix; the next commit on `dev` reintroduces one for the *next*
target version. Use PEP 440 spelling — `1.0.0rc1`, not `1.0.0-rc1`,
which setuptools would silently normalize, making this file disagree
with the installed metadata.

**After bumping, re-run `pip install -e .`.** An editable install's pip
metadata is a snapshot taken at install time, so `pip show corpus` keeps
reporting the *old* version until you reinstall — while `corpus
--version`, the import, and every stamped artifact immediately report the
new one. Only `pip show` drifts, and only until reinstall; artifacts are
never wrong.

**A version bump forces a full reprocess.** `_stage_recorded_complete`
(`pipeline/stages.py`) treats a stage record whose `pipeline_version`
differs from the current one as stale, so `--resume` re-runs every stage
on every paper. That is usually what you want when the code changed, but
bump deliberately before a large incremental run rather than during it.

Schema versions for on-disk artifacts (sqlite DBs, LanceDB index) are
tracked separately from code SemVer — a code bug-fix doesn't break old
DBs, but a schema change does. Stamp schema versions inside the file
itself if you change one.

## Layout conventions

- The repo root has **zero Python files** (post-v0.3 / #60). The unified `corpus` binary is the only operator-facing entry point; it lands on PATH via `[project.scripts]` after `pip install -e .` from `pyproject.toml`. The CLI router lives at [pipeline/cli.py](pipeline/cli.py); subcommands dispatch via `python -m <module>` to private modules under `pipeline/`, `bib/`, or `mcpsrv/`.
- [pipeline/](pipeline/) — top-level CLI router (`cli.py` → the `corpus` binary), Stage 1 + Pass 3b/3c orchestrator (`main.py`, `orchestrator.py`), post-pipeline modules (`embed.py`, `taxon_mentions.py`, `intext_citations.py`, `taxonomy_ingest.py`, `status.py`), pydantic config schema + bundled template (`config_schema.py`, `config.template.yaml`), shared rich console layer (`console.py`), and supporting library modules: `scan.py` (OCR), `extract.py` (docling), `metadata.py` (Grobid + bib), `chunking.py`, `annotate.py` (taxa + lexicons), `figure_passes.py`, `runner.py` (per-paper orchestrator), `figures.py`, `taxa.py`, `grobid_client.py`, `embeddings.py`, `vision.py`, `external.py` (shared retry / circuit breaker), `version.py` (single-source `__version__`), `config.py` / `io.py` / `log.py` / `stages.py`.
- [bib/](bib/) — bibliographic round-trip + biblio authority + reconcile: `parser.py` (BibTeX parser + `BibIndex`), `export.py` (DB → BibTeX), `importer.py` (BibTeX → DB; sidecar staging for #67), `authority.py` (the build), `reconcile.py` (ghost merge). Exposed as a single namespace via `from bib import …`.
- [mcpsrv/](mcpsrv/) — MCP server package: `app.py` (MCPServer instance + index accessor), `indexes.py` (`CorpusIndex` / `TaxonMentionDB` / `BiblioAuthority`), `tools/` (papers / taxonomy / figures / chunks / bibliography), `transport.py` (bearer auth + SSE), `main.py` (argparse + dispatch), `bundle.py` (built-bundle → served-bundle distillation). Adding a new MCP tool: extend the appropriate `mcpsrv/tools/<concern>.py` module and import it from `mcpsrv/tools/__init__.py` so the `@mcp.tool()` decorator registers at startup.
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
