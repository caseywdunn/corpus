# AGENTS.md

Orientation for AI coding agents working in this repository.

## Purpose

A workflow for interrogating a corpus of scientific literature PDFs, spanning both scanned/old and born-digital papers. Target use cases: searching, citation analysis, and collecting figures of the same species across different papers. Corpus-agnostic by design — the group-specific data (taxonomy snapshot, lexicon, BibTeX, instructions) lives in the per-instance corpuscle directory; the lab's reference deployment is a siphonophore corpus, but every other group plugs into the same machinery.

## Documentation scope

Corpus is group-agnostic. The general documentation applies to any corpus and any
cluster. The siphonophore / Bouchet files are one worked example — **do not treat
them as requirements or as the only supported workflow**.

| File | Scope | Notes |
|---|---|---|
| `README.md`, `INSTALL.md` | General, **public-facing** | Any corpus, any platform. Written for someone on their own machine. |
| `DEPLOY.md` | General | AWS deploy runbook; any corpus |
| `dev_docs/OVERVIEW.md`, `dev_docs/PLAN.md`, `dev_docs/MCP_TOOLS.md` | General | Architecture, roadmap, tool reference |
| `dev_docs/TESTING.md` | General | Quality test suite (ground truth, eval workflow) |
| `dev_docs/OCR_LANGUAGES.md` | General | What the gold set says about choosing Tesseract packs: when a union helps, when it hurts, why the native pack matters (character coverage, not dictionaries), and pinning vs detection. Stable — cite this from code comments rather than PLAN.md. |
| `tools/smoke_test_sse.py` | General | Programmatic MCP smoke test; works on any corpuscle. Requires a compute node (not login) for full Layer 3 coverage — Layer 3 loads the ~600 MB BGE-M3 embedder for semantic search. |
| `dev_docs/BOUCHET.md` | **Siphonophore + Yale Bouchet (example)** | Runbook for the Dunn-lab siphonophore corpus on Yale's Bouchet HPC. The SLURM scripts, paths, and partition names are Bouchet-specific; the pattern (config-driven `corpus run`, job arrays, pre-download models) is general. Use as a template, not a literal guide. |
| `dev_docs/ACCEPTANCE_PROMPTS.md` | **Siphonophore (example)** | Manual acceptance prompts for the siphonophore corpuscle. Taxon names are siphonophore-specific; adapt for other groups. Referenced from BOUCHET.md. |

When authoring documentation for a new corpus or cluster, add analogues to
`BOUCHET.md` / `ACCEPTANCE_PROMPTS.md` in your own runbook file — don't modify the
general docs to embed group-specific details.

### The siphonophore library's size is dynamic — don't hard-code it

The reference siphonophore library grows as papers are added; any exact count
written into docs or code is wrong by the next `git lfs pull`. Say **`~2,000
papers`**, "the current corpus", or "as the corpus grows" instead of a precise
figure, and phrase anything derived from the size — walltimes, batch counts,
disk sizing — as a rule that scales rather than a number measured once.

This applies to the corpus size only. Genuinely fixed counts stay exact: the
4 + 1-paper demo corpus, the MCP tool count, the 64-PDF `BATCH_SIZE`. Note
`NUM_BATCHES` is *not* fixed — it is `ceil(unique PDFs / BATCH_SIZE)`, so give
the formula alongside any literal. And verify a "fixed" count before repeating
it: every number in this paragraph was wrong when it was checked in Aug 2026 —
the demo was slimmed 11 → 4 + 1 in v0.4, and 256 is `PASS3B_BATCH_SIZE` (the
vision pass), not `BATCH_SIZE`, which defaults to 64.

`CHANGELOG.md` is exempt too: it records what was true at a release, so leave
historical figures alone.

### Don't cite `dev_docs/PLAN.md` from code comments

PLAN.md is ephemeral. It is pruned and renumbered at every release
(CONTRIBUTING.md's release ritual, step 8), so a comment reading
`see dev_docs/PLAN.md §10` points at a section that no longer exists — and says
nothing about what the code does. This has already happened twice: 16 such
comments across `mcpsrv/`, `pipeline/`, `tests/` and BOUCHET.md were left
pointing at a v0.1-era numbering.

Cite a **stable** doc instead — `dev_docs/OVERVIEW.md` for architecture and
stage internals, `DEPLOY.md` for the served surface, `dev_docs/MCP_TOOLS.md` for
the tool catalog — or drop the pointer and say the thing. Link by *section
title*, not number, since those move too. Issue numbers (`#176`) are fine and
preferred: they are permanent and carry the reasoning.

### Keep Yale specifics out of the public docs

`README.md` and `INSTALL.md` are read by people who have never heard of Yale and
are running corpus on a laptop. Keep them crisp and general:

- **No Bouchet-specific content** — no partition names, `module load` lines,
  `/nfs/roberts/...` paths, NetIDs, QoS limits, or cluster-policy claims ("login
  nodes block X"). Cluster policy varies per site and goes stale; asserting one
  cluster's behavior as a general fact is how these docs acquire bugs.
- **One pointer is the budget.** README gets a single mention — *"we use Yale's
  Bouchet, runbook in `dev_docs/BOUCHET.md`"* — plus its entry in the docs index.
  INSTALL may link `BOUCHET.md` as the worked example for a recipe. That's it.
- **Generic HPC comments are fine** in moderation: "on a SLURM cluster", "on an
  HPC cluster, check `module avail` first". Say the general thing, then link the
  worked example rather than inlining it.

Everything site-specific belongs in `dev_docs/BOUCHET.md` and `slurm/`.

## Documentation map

- [README.md](README.md) — installation, usage, MCP server setup, examples
- [CONTRIBUTING.md](CONTRIBUTING.md) — branching model, release ritual, version bumps
- [CHANGELOG.md](CHANGELOG.md) — what changed in each release
- [dev_docs/PLAN.md](dev_docs/PLAN.md) — roadmap and design decisions for the active version. Turns over every release — never cite it from code comments (see above)
- [dev_docs/OVERVIEW.md](dev_docs/OVERVIEW.md) — pipeline architecture, stage internals, figure pipeline, key files
- [dev_docs/MCP_TOOLS.md](dev_docs/MCP_TOOLS.md) — full MCP tool surface + count
- [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md) — HPC runbook for siphonophore corpus on Yale Bouchet (example; adapt for other clusters/corpora)
- [dev_docs/ACCEPTANCE_PROMPTS.md](dev_docs/ACCEPTANCE_PROMPTS.md) — manual acceptance prompts for the siphonophore corpuscle (example; adapt taxon names)
- [DEPLOY.md](DEPLOY.md) — AWS deploy runbook (S3 bundle + EC2 systemd)
- [dev_docs/TESTING.md](dev_docs/TESTING.md) — quality test suite, ground truth format, evaluation workflow
- [INSTALL.md](INSTALL.md) — optional OCR extras, pip-only fallback, platform notes
- [.github/workflows/](.github/workflows/) — CI tiers T0 (lint + unit) and T1/T2 (demo build + serve on Linux + macOS, each followed by the 4 + 1 implicit-resume scenario); all fire on every push and every PR. See CONTRIBUTING.md's "Continuous integration" section for the full tier table.

## Quick reference

```bash
# Scaffold + run a corpuscle (v0.3+; #58–#68).
cd <corpuscle>          # contains config.yaml; see `corpus init` for a starter
corpus run              # full pipeline + post-pipeline + bundle, implicit resume
corpus status           # rollup + report
corpus serve            # MCP server (reads bundle from config.yaml's output_dir)

# Module-level entry points for ad-hoc / debug invocations:
python -m pipeline.main          <input_dir> <output_dir>          # Stage 1
python -m pipeline.embed         <output_dir>                       # Stage 2
python -m bib.authority          <output_dir>                       # biblio DB
python -m pipeline.taxon_mentions <output_dir>                       # taxon DB
python -m pipeline.intext_citations <output_dir>                     # intext refs
python -m bib.reconcile          <output_dir>                       # ghost merge
python -m pipeline.taxonomy_ingest <output_dir> --source <dwc|dwca|worms> ...
python -m mcpsrv.bundle          <output_dir> <serve_dir> --version vX.Y.Z

# Parallel Stage 1 on Bouchet (day partition, 256 PDFs per batch)
NUM_BATCHES=8 bash slurm/batch_pipeline.sh
```

## Repo layout

| Path | Role |
| --- | --- |
| `pipeline/` | Top-level CLI router (`pipeline/cli.py` → the `corpus` binary), Stage 1 + Pass 3b/3c orchestrator, post-pipeline modules (`embed`, `taxon_mentions`, `intext_citations`, `taxonomy_ingest`, `status`, `orchestrator`), and shared library modules. Pydantic config schema in `config_schema.py`; bundled `config.template.yaml`. |
| `mcpsrv/` | MCP server. MCPServer `app.py` + the `@mcp.tool()` surface in `tools/{papers,taxonomy,bibliography,figures,chunks,lexicon}.py` (catalog + count in [dev_docs/MCP_TOOLS.md](dev_docs/MCP_TOOLS.md)); `mcpsrv.bundle` distills a build into a served bundle. |
| `bib/` | BibTeX round-trip + biblio authority + reconcile (`bib.parser`, `bib.export`, `bib.importer`, `bib.authority`, `bib.reconcile`). |
| `slurm/` | SLURM batch scripts (Bouchet). |
| `deploy/` | CloudFormation, nginx, systemd, sync + update scripts. |
| `tools/` | Operator and measurement tooling. Reads corpuscle artifacts; **never imported by `pipeline/`**. `tools/qc/` scores a build against the gold transcription set. |
| `tests/` | Ground-truth + corpus-wide consistency tests. |

### Where library-facing work goes

A *library* is the upstream collection of PDFs and a `.bib`; a *corpuscle* is
what corpus builds from one. Code that helps curate a library — resolving a
language tag to OCR packs, inspecting pages to decide a page range — has an
obvious pull toward living in the library repo that needs it. It should not,
because two libraries then carry two copies of the same knowledge and they
drift. The siphonophore library already carries a temporary local copy of the
BCP-47 → Tesseract table for exactly this reason.

Three tiers, and the question that sorts them is **does this encode knowledge
about PDFs and OCR, or knowledge about one particular collection?**

| kind | example | home |
| --- | --- | --- |
| Generic PDF/OCR knowledge | `bcp47_to_tesseract` (#215) | `pipeline/` — public function, no wall needed |
| Read-only inspection over a library | per-page evidence for `keeppages` (#217) | `tools/` |
| A judgment about *this* collection | "this document is `de-Latf`"; "pages 1–2 are a wrapper" | `skills/` (#178), or the library repo |

**`pipeline/` may never import from `tools/` or `skills/`.** That is the rule
that keeps the tiers real, and `tests/test_import_direction.py` enforces it.
Dependencies point one way: a skill may import from `pipeline/`, never the
reverse.

The trap this avoids is not tidiness. A mapping table duplicated into a
library repo has to agree with what `scan.py` actually does, and nothing
checks that it still does — the same shape as `consolidateCitations` being
settable and never read, or `tesseract_packs` being recorded and then
recomputed by a different route.

## Implementation notes for contributors

- Each unique PDF is identified by the first 12 hex chars of its SHA-256. All artifacts live under `<output_dir>/documents/<HASH>/`. Per-paper `pipeline_state.json` records each stage's `pipeline_version` + `input_fingerprint`; `corpus run` re-runs whichever stages' record is missing or disagrees with the current input (implicit resume per #60; the v0.2 `--resume` flag is gone). A lexicon edit only re-runs the annotation pass.
- Stage 1 (`pipeline.main`) is CPU-only: scan detection, OCR, docling extraction, Grobid metadata, chunking, annotation. Stage 2 (`pipeline.embed`) is GPU: BGE-M3 embeddings into LanceDB. Pass 3b (vision LLM) is opt-in via `vision.backend: {claude, local}` in `config.yaml`.
- Figure extraction has a fallback chain: docling pictures first, then raw PyMuPDF `page.get_images()`. Only figures that land on disk are recorded in `figures.json`.
- Visualizations overlay text-cell boxes (red) and figure bboxes (yellow/orange). Coordinates are Y-flipped from docling's bottom-left origin to PIL's top-left.
- Annotation artifacts (`taxa.json` plus one `<category>.json` per category in the `--lexicon` YAML) stamp an `input_fingerprint`; per-stage resume detects mismatch and re-runs the annotation pass automatically.
- Pipeline failures land in structured `summary.json["stage_failures"]` with reason codes (`timeout`, `crash`, `external_unavailable`, `unsupported_format`, `corrupted`, `quality_gate`, `too_large`); silent-failure quality gates emit `quality_flags`. Free-text `errors[]` is being phased out.
- `config.yaml` is loaded by `pipeline.config.load_config`; CLI flags override config values, which override built-in defaults.
- External calls (Grobid, BHL, CrossRef, OpenAlex) flow through `pipeline/external.py` for shared retry + backoff + circuit breaker. `--strict-network` aborts on the first transient failure for release builds.
- Steering the client session: the MCP server is client-driven — it can only add to the model's context in *response* to a request, never push unsolicited. Soft guidance (report/manuscript structure, cross-validation nudges, house style) rides on `instructions` (always-on, from `<corpuscle>/instructions.md`), tool-result text (just-in-time, keyed to the active output profile per #101; `format_citation`'s provenance warning is the established example), and user-invoked MCP prompts (`/mcp__corpus__<name>`). Keep tool-result nudges terse — they cost tokens on every call (#76, #81–86). Hard requirements (figure publishability, citation provenance) are enforced server-side at the tool boundary, never via a nudge the model can ignore. See OVERVIEW.md "Steering the client session".
- `__version__` in `pipeline/version.py` is the single source of truth; CONTRIBUTING.md covers the branching model.
