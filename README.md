# corpus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19964909.svg)](https://doi.org/10.5281/zenodo.19964909)

## Citing corpus

We are preparing a manuscript that presents corpus. In the mean time, please cite corpus as:

> Dunn, C. W., Manko, M. K., Zapata, F., & Church, S. H. (2026). Corpus: Model Context Protocol (MCP) server for taxonomic literature. <https://doi.org/10.5281/zenodo.19964909>

## What corpus does

Think of corpus as an interface between a body of scientific literature and a Large Language Model (LLM) like Claude (Anthropic). LLMs absorb broad knowledge during training but lack granular detail on specific scientific topics. Corpus supplies that detail from a literature collection of your choice, so the model can answer questions grounded in the papers themselves rather than its training data alone.

The workflow turns a folder of PDFs — born-digital articles alongside centuries-old scans in multiple languages — into a queryable knowledge base, exposed over the [Model Context Protocol (MCP)](https://modelcontextprotocol.io/). MCP is the open standard that LLM clients (Claude Desktop, Claude Code, claude.ai web, Cursor) use to reach external data. The primary target is taxonomic literature.

Corpus itself is general-purpose — not tied to any group of organisms. You aim it at your PDFs and build a **corpuscle**: a self-contained MCP server for one group (siphonophores, drosophila, ferns, …). Run a corpuscle locally for your own work, or deploy it on a server so colleagues can connect. We build corpuscles on Yale's HPC cluster and host public servers on Amazon Web Services (AWS), but those are our choices; any machine with Python and adequate disk will do.

### Implementation

Corpus reads a folder of PDFs and produces a knowledge base built around: **taxa, figures, bibliographies, and the cross-paper relationships between them**. It handles several tasks:

- **OCR for old scans** — including 19th-century German Fraktur and other historical typefaces, with per-paper language detection.
- **Figure extraction** — every plate, photograph, and line drawing pulled out with its caption, then vision-LLM-tagged for taxonomy and anatomy.
- **Bibliography parsing and reconciliation** — references inside each paper are deduplicated across the whole corpus, so, for example, "Totton 1965" is one entity even when twenty papers cite it differently.
- **Taxonomy linking** — every species mention is matched against a Darwin Core taxonomy with full synonymy, so a question about *Apolemia uvaria* finds papers that only ever called it *Stephanomia uvaria*.

The output is a per-paper artifact tree plus cross-paper databases. You query it through an **MCP server** (`corpus serve`, backed by [mcpsrv/](mcpsrv/)) that any [MCP](https://modelcontextprotocol.io/) client (Claude Desktop, Claude Code, claude.ai web) can connect to — so "show me every figure of *Nanomia bijuga* nectophores in the corpus" becomes something you ask in chat instead of grep across a hard drive. The server is a read-only view over the per-paper artifacts, so re-running the pipeline is the only way new data reaches the tools.

### Example uses

- *"List every valid* Apolemia *species and for each list the papers that discuss them."*
- *"Show the bibliography of Pugh 1997 and mark which cited works are also in the corpus."*
- *"Give me every figure showing nectophore morphology in* Nanomia bijuga *."*
- *"What are the most-cited papers that aren't currently in the corpus?"*
- *"Summarize what the corpus says about pneumatophore structure."*
- *"Translate the diagnosis of* Forskalia edwardsii *from Haeckel 1888 (German) into English."*

The full tool surface is in [dev_docs/MCP_TOOLS.md](dev_docs/MCP_TOOLS.md).

## What you bring

### A body of pdfs (required)

A directory of PDFs is the only required input. They can be any mix of born-digital papers and scanned older literature — the pipeline handles both. Two practical points:

- **Quality matters.** Clean high-resolution scans produce dramatically better OCR and figure extraction than re-scans of photocopies. For born-digital papers, the original publisher PDF beats any downstream export.
- **Filename convention.** Name files `<Surname><Year>.pdf` (e.g., `Totton1965a.pdf`). The bibliography reconciler uses this filename as a fallback when Grobid mis-parses the paper reference.

### Bibliographic data for pdfs (optional)

Grobid extracts metadata (authors, title, year, journal) from each PDF's header, but it gets confused on older scans and unusual layouts. If you have a curated BibTeX file, pass it with `--bib` and it overrides Grobid for any paper that an entry references via its `file = {...}` field:

```bibtex
@article{Totton1965a,
  author  = {Totton, A. K.},
  title   = {A Synopsis of the Siphonophora},
  year    = {1965},
  journal = {British Museum (Natural History)},
  file    = {Totton1965a.pdf},
}
```

```yaml
# in your corpuscle's config.yaml
bib: ./references.bib
```

Matching is on the basename (case-insensitive), so the `file` value can be a bare filename or a full path. This is the cleanest way to get accurate per-paper metadata where you've already curated references for your own writing.

### External taxonomic data (optional)

The taxonomy database (`taxonomy.sqlite`) is a [Darwin Core](https://dwc.tdwg.org/) snapshot that drives synonymy resolution — the layer that lets a question about *Apolemia uvaria* find papers that only ever wrote *Stephanomia uvaria*. It's built once on the first `corpus run` from whichever source you point at; subsequent runs reuse the cached SQLite unless you pass `--force-rebuild-taxonomy`.

The `taxonomy:` block in `config.yaml` picks the source. The bundled template ships it commented out — opt in by uncommenting and choosing one of:

```yaml
# WoRMS — marine groups. Fetches an AphiaID subtree from the World
# Register of Marine Species. Look up your group at marinespecies.org.
taxonomy:
  source: worms
  root_id: 1371       # Siphonophorae order (example)

# Darwin Core Archive — non-marine groups. Pull a DwC-A .zip from
# GBIF / ITIS / Catalogue of Life and point at it locally.
taxonomy:
  source: dwca
  path: ./dwca.zip

# Darwin Core directory — rare, an already-unzipped DwC tree.
taxonomy:
  source: dwc
  path: ./dwc/
```

| Source | What it is | When to use |
| --- | --- | --- |
| `worms` | [WoRMS](https://www.marinespecies.org/) subtree fetched by AphiaID | Marine taxa (siphonophores, copepods, fish, ...) |
| `dwca` | [Darwin Core Archive](https://dwc.tdwg.org/text/) `.zip` from [GBIF](https://www.gbif.org/), [ITIS](https://www.itis.gov/), [Catalogue of Life](https://www.catalogueoflife.org/) | Anything else; the broad default |
| `dwc` | Directory of Darwin Core `.tsv` files | An already-unzipped DwC tree |

Or invoke the ingester directly without `corpus run`: `corpus taxonomy ingest --source dwca --input path/to/dwca.zip` (equivalent to `python -m pipeline.taxonomy_ingest <output_dir>`).

The inverse — dumping a corpus's built taxonomy back out as a DwC-A — is `corpus taxonomy export -o taxonomy.zip`. Use it to share a taxonomy snapshot without forcing the recipient to walk WoRMS again, or to commit a small fixture into a downstream repo so CI exercises the `dwca` ingest path without network calls. The round-trip property: `corpus taxonomy ingest --source dwca --input <export.zip>` recovers the same `taxa` row set as the source SQLite.

**Without a `taxonomy:` block** the pipeline still extracts taxon mentions from text — you only lose the synonymy graph that links historical names to current valid names. The default template ships with the block commented out, so leaving it alone is the no-taxonomy path.

### Lexicons (optional)

A lexicon is a small YAML file listing terms you want tagged in figure captions and chunk text. Top-level keys are categories (e.g. `anatomy:`, `biogeography:`, `methods:`); each value is a flat term-to-metadata map. See [demo/lexicon.yaml](demo/lexicon.yaml) for a worked example with siphonophore anatomy (*nectophore*, *pneumatophore*, *gastrozooid* …).

```yaml
anatomy:
  nectophore:
    synonyms: [nectophores, swimming bell]
    translations: {de: [Schwimmglocke]}
biogeography:
  pelagic:
    synonyms: [open water]
```

Set `lexicon: ./lexicon.yaml` in your corpuscle's `config.yaml`. Each category emits its own `<hash>/<category>.json` artifact. Per-category content is fingerprinted independently, so editing a single section only invalidates that category's annotations on the next `corpus run` — the rest stay cached (implicit resume; #60 dropped the `--resume` flag).

Without a `lexicon:` entry, lexicon extraction is skipped entirely. Like `bib:`, the lexicon is an input you maintain alongside your literature, not something the tool ships.

### Instructions for the LLM (optional)

A markdown file at `<corpuscle>/instructions.md` whose contents land in every chat session against the corpus. The MCP server returns it in `InitializeResult.instructions`, and well-behaved clients (Claude Desktop, Claude Code) inject it into the LLM's context at session start.

Use it for per-corpus nudges that no taxonomy or lexicon entry can express — see [demo/instructions.md](demo/instructions.md) for a worked example, which (among other things) tells the model that *Velella* and *Porpita* are not siphonophores even when older literature lumps them in.

Override the default location with `--instructions <path>` when starting the MCP server. If the file is absent, no instructions are sent.

## Computational requirements

- **Disk.** Plan on several times the size of the original PDFs.
- **CPU vs. GPU.** Stage 1 (OCR, layout, Grobid, chunking, annotation) is CPU-bound and parallelizes well across PDFs. Stage 2 (BGE-M3 embeddings) and the optional vision pass (Qwen2.5-VL-7B figure tagging) are GPU-accelerated; embeddings still run on CPU but slowly, and the local vision backend needs a GPU to be usable. A Claude API vision backend is also available — costs API credits but runs anywhere.
- **Scale.** A few dozen PDFs run on a laptop in a couple of hours. A few thousand benefit from an HPC cluster — we use Yale's Bouchet, runbook in [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md).
- **MCP client.** The query interface is MCP, so you'll need a client that speaks it (Claude Desktop, Claude Code, claude.ai web with custom connectors, Cursor, Continue). Most require an Anthropic subscription.
- **Remote deployment.** Serving the corpus to others over the network requires a server. The reference deploy is AWS (EC2 instances behind a shared Application Load Balancer), but any host with Python and an open port works. See [Deploying MCP server remotely](#deploying-mcp-server-remotely) below.

## Corpuscle layout

A *corpuscle* is the on-disk container for one corpus instance — siphonophores, drosophila, or whatever group you're working on. It's a single directory:

```text
<corpuscle>/
├── documents/<HASH>/         # per-paper artifacts (text, chunks, figures, taxa, anatomy, …)
├── vector_db/lancedb/        # embeddings index
├── taxonomy.sqlite           # Darwin Core snapshot, built by `pipeline.taxonomy_ingest`
├── biblio_authority.sqlite   # deduplicated works + citation graph
├── taxon_mentions.sqlite     # cross-paper taxon index
└── _serve/                   # distilled served bundle (#60); contains
                              # bundle_manifest.json + a whitelisted subset
                              # of the above. Pass to `corpus serve --output-dir`
                              # for ship-to-host workflows. Skip with
                              # `corpus run --no-bundle`.
```

Three cross-paper layers are rebuildable independently of the per-paper artifacts (no re-running OCR or extraction):

| Layer | File | Built by |
| --- | --- | --- |
| **Taxonomy** — Darwin Core snapshot with synonymy | `taxonomy.sqlite` | `pipeline.taxonomy_ingest` (sources: `dwc` / `dwca` / `worms`) |
| **Bibliography** — deduplicated works + citation graph | `biblio_authority.sqlite` | `bib.authority` + `bib.reconcile` |
| **Taxon mentions** — cross-paper taxon-to-paper index | `taxon_mentions.sqlite` | `pipeline.taxon_mentions` |

Every CLI takes the corpuscle root as its first positional argument and resolves all per-instance files from there. Run two corpora side-by-side by giving each its own corpuscle directory; they don't share state.

## Installation

```bash
conda env create -f environment.yaml
conda activate corpus
pip install -e .
```

`pip install -e .` puts the `corpus` binary on PATH (via the
`[project.scripts]` entry point in `pyproject.toml`); the package
metadata version stays in sync with `pipeline/version.py` so
`pip show corpus`, `corpus --version`, and the bundle manifest
never drift.

### Supported platforms

| Target | Status |
| --- | --- |
| **linux-x86_64** | Supported. The reference HPC + deploy target (Bouchet, AWS EC2). |
| **macOS arm64** (Apple Silicon) | Supported, with two extra constraints — see below. |
| macOS x86_64 (Intel Mac, or Rosetta on Apple Silicon) | **Not supported.** Apple dropped Intel-mac PyTorch wheels after 2.2; `docling` and `transformers ≥ 5.0` require torch ≥ 2.4, which has no macOS Intel build. |
| linux-aarch64 (Graviton, Ampere) | Not currently supported. Wheels exist for most dependencies but we don't test it; add as a third target if a real Graviton use case appears. |

**On Apple Silicon, two things have to be right:**

1. **Use an arm64-native conda** ([miniforge](https://github.com/conda-forge/miniforge)) to create the corpus env. The default download from anaconda.com is the Intel x86_64 build and keeps running under Rosetta even on an arm64 host — that traps the install in the unsupported macOS Intel matrix above, and `transformers` then silently disables PyTorch. After `conda env create -f environment.yaml`, verify with:

   ```bash
   ~/miniforge3/envs/corpus/bin/python -c "import platform; print(platform.machine())"
   # expect: arm64   (NOT x86_64)
   ```

2. **Grobid still runs under Rosetta.** The official Grobid Docker image (`grobid/grobid:0.8.1`) is `linux/amd64` only. On Apple Silicon, Docker Desktop runs it under x86_64 emulation. Enable **Settings → General → "Use Rosetta for x86_64/amd64 emulation"** in Docker Desktop for the fastest path; without Rosetta, QEMU emulation works but is noticeably slower. Memory budget is the bigger concern — keep the `JAVA_OPTS=-Xmx8g` heap in `docker-compose.yml` only if Docker Desktop has at least 12 GB allocated under **Settings → Resources**, otherwise drop to `-Xmx4g`. Grobid is only used at pipeline build time, not at MCP serve time, so this overhead is bounded to `corpus run`.

See [INSTALL.md](INSTALL.md#supported-platforms) for the optional OCR helper install (pngquant, jbig2enc) and the pip-only fallback.

Grobid runs as a separate service that must be up *before* you call `corpus run`. `docker compose up -d` runs it in the background; leave it running while you work and stop it with `docker compose stop grobid` when you're done. `corpus run` won't try to launch it for you — auto-launching cross-platform is awkward (docker on a laptop, Singularity on Bouchet, neither on a stripped-down host). `corpus check` confirms reachability before you commit to a long pipeline run.

**Docker is a prerequisite.** Install it from <https://docs.docker.com/engine/install/> (or `apt install docker.io` on Debian/Ubuntu, `brew install --cask docker` on macOS) before the commands below. On HPC hosts without Docker, [Apptainer](https://apptainer.org/) can pull the same Grobid image — see [INSTALL.md](INSTALL.md#grobid-on-bouchet).

```bash
docker compose up -d grobid              # start in background; persists across runs
curl http://localhost:8070/api/isalive   # should print "true"
corpus check                             # confirms grobid + GPU + config + disk
```

After creating the env, run `bash tools/install_tessdata.sh` to download the Tesseract language packs (the conda-forge `tesseract` package ships only English data — see [INSTALL.md](INSTALL.md#ocr-language-packs)). Optional `jbig2enc` compression is covered in the same doc.

## Language support

OCR is the only stage where language matters. Each scanned PDF gets its language detected automatically (`langdetect` on the first few pages of extracted text), and the result is routed to a matching [Tesseract](https://github.com/tesseract-ocr/tesseract) language pack — a `<code>.traineddata` file that teaches Tesseract to recognize a particular language or script. Born-digital PDFs with a clean text layer skip OCR entirely, so packs only matter for scans.

**What ships by default.** The conda-forge `tesseract` package ships only the English LSTM model. The default fallback set in `config.yaml` (`ocr.ocr_languages_default`: `eng`, `deu`, `fra`, `rus`, `lat`, `spa`, `por`, `chi_sim`, `chi_tra`, `jpn`, `ell`, `kor`) plus `grc` (Ancient Greek, opt-in) and `deu_latf` (19th-c. German Fraktur) are installed by running [`tools/install_tessdata.sh`](tools/install_tessdata.sh) after `conda env create` — it drops the matching `<code>.traineddata` files into `$CONDA_PREFIX/share/tessdata/` directly from the official [`tesseract-ocr/tessdata_best`](https://github.com/tesseract-ocr/tessdata_best) repo. None of the per-language packs are on conda-forge.

**Without `deu_latf`,** scanned 19th-c. German papers (Goldfuss 1820, Brandt 1837, Pagenstecher 1869, Dönitz 1871, …) OCR to whitespace. The default install bundles it, so no separate step is needed unless you trimmed the language list.

**How to tell if you need more packs.** After a Stage 1 run, `<output_dir>/documents/<HASH>/scan_detection.json` records the detected language for each scanned PDF. Skim those to see what languages your corpus actually contains. If a language is detected but no matching pack is installed, OCR silently falls back to the union of installed packs (English is always appended last) — extraction still succeeds, but accuracy on that paper drops. Adding more languages than you need is safe but slows OCR, so trim `ocr.ocr_languages_default` for a known-narrow corpus.

**Adding a language outside the default set.** The full ISO-to-Tesseract map covers ~40 languages including CJK, Cyrillic, Greek, Arabic, and Indic scripts — see `_ISO_TO_TESSERACT` in [pipeline/scan.py](pipeline/scan.py). Pass the codes to the installer to grab the matching packs:

```bash
bash tools/install_tessdata.sh ara hin tha
```

To make a new language part of the fallback set tried when detection is uncertain, add its code to `ocr.ocr_languages_default` in `config.yaml`.

**Everything past OCR is language-agnostic.** The [BGE-M3](https://huggingface.co/BAAI/bge-m3) embedding model is multilingual (trained on 100+ languages) and built for cross-lingual retrieval — a question in English finds passages in German, Latin, French, and so on. Grobid metadata extraction works across European languages, and the vision LLM that tags figures reads labels in whatever script is on the plate. So once OCR succeeds, no further language configuration is needed.

## Try it on the demo corpus

The repo ships [demo/](demo/) — a regular corpuscle (11 siphonophore PDFs, `siphonophores.bib`, `lexicon.yaml`, `instructions.md`, a pre-built Siphonophorae taxonomy as `taxonomy.zip`, and a `config.yaml` already pointing at all of it). The bundled DwC-A means the first `corpus run` doesn't walk the WoRMS REST API — it ingests the full taxonomy from a local file in seconds. One command runs the full pipeline + cross-paper builds + bundle:

```bash
cd demo && corpus run
```

Output lands in `demo/output/` (gitignored). Re-run anytime — `corpus run` is always idempotent.

## First time run

```bash
mkdir my-corpus && cd my-corpus
corpus init                # scaffold config.yaml from the bundled template
$EDITOR config.yaml        # set input_pdfs, taxonomy, optional bib + lexicon
corpus check               # pre-flight: validates config + probes Grobid/GPU/disk
corpus run                 # full pipeline + post-pipeline + bundle, hands-off
corpus status --report     # stage pass-rates + cross-paper artifact ✓/✗
corpus serve               # local MCP server against the freshly built bundle
```

The work divides into two stages, best run on different hardware: Stage 1 is CPU-heavy and embarrassingly parallel; Stage 2 wants a GPU. Each unique PDF is identified by the first 12 hex chars of its SHA-256, and all artifacts live under `<output_dir>/documents/<HASH>/`. `corpus run` is always idempotent — re-runs only do work whose inputs have changed (per-stage state + content fingerprints from v0.2). `--force-rebuild` covers the rare clean-rebuild case.

For large corpora on Bouchet, the SLURM chain (Stage 1 + embed + finalize) is one line:

```bash
NUM_BATCHES=8 bash slurm/batch_pipeline.sh
```

For long laptop runs that outlast a single SSH session, detach with either tool:

```bash
nohup corpus run > run.log 2>&1 &     # background; tail run.log to follow
# or
tmux new-session -s corpus 'corpus run; bash'   # detach with C-b d
```

## Vision pass (optional)

Set `vision.backend: local` (needs CUDA/MPS) or `vision.backend: claude` (needs `ANTHROPIC_API_KEY`) in `config.yaml` and `corpus run` picks it up automatically. If the backend isn't usable on the host, the pass skips gracefully with a one-line nudge — Pass 3b never hard-fails deep into the run because of a missing GPU or missing key. `--no-vision` is the explicit-skip flag for a config that requests it.

On Bouchet, the SLURM chain (`slurm/batch_pipeline.sh`) runs Pass 3b on `gpu_h200` automatically — see [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md).

## Adding and updating documents

Drop new PDFs into the configured `input_pdfs:` directory and re-run `corpus run`. Implicit resume only re-processes what changed; orphan PDFs (removed from the input set) get their `documents/<HASH>/` and matching LanceDB rows pruned by default, with a 25%-safety-rail (override with `--force-prune`; opt out entirely with `--no-prune`).

Lexicon / taxonomy edits need no extra flag — fingerprints land in `<hash>/<category>.json` and the affected category re-annotates on the next run.

## Curating bibliographic metadata

Grobid mis-parses some references. Round-trip via BibTeX:

```bash
corpus bib export -o my_corpus.bib   # current biblio_authority → BibTeX
$EDITOR my_corpus.bib                # fix titles, years, authors
corpus bib import my_corpus.bib      # apply edits back
```

`corpus bib import` works *before* the first `corpus run` too — it stages the edits to a sidecar that's applied automatically when biblio_authority.sqlite first gets built.

Each entry carries a stable `corpus_hash` field that bib_import uses to match edits back to the right `works` row, with DOI as a fallback. The `.bib` is the source of truth for user edits — preserve it across pipeline rebuilds and re-run `bib.importer` afterward to re-apply your curation. Hand-edited entries override Grobid-derived metadata in every MCP tool that reads from the authority DB.

`--bib` (the input flag described in [Bibliographic data for pdfs](#bibliographic-data-for-pdfs-optional)) and the round-trip above are independent: the input flag overrides Grobid's per-paper header from a curated BibTeX *before* the authority DB is built; the round-trip edits the authority DB *after*. Use whichever fits where you already maintain the metadata.

## Distilling a served bundle

The corpuscle the pipeline emits is the **build bundle**: everything `pipeline.main` produces, including `processed.pdf`, raw docling dumps, per-page QC visualizations, and per-paper logs. For ~2,000 papers that's ~10 GB. The MCP server doesn't need most of it — only the JSON artifacts it reads at startup, the figure PNGs, and the precompiled indices. [`mcpsrv.bundle`](mcpsrv/bundle.py) distills the build bundle down to a **served bundle** (~3 GB for the same corpus) by copying only whitelisted files, scrubbing absolute paths from JSON values, and writing a versioned `bundle_manifest.json` that the MCP `bundle_info` tool surfaces:

```bash
python -m mcpsrv.bundle <output_dir> <serve_bundle_dir> --version v1.0.0
```

You only need this when **shipping a corpus to a different host** — a remote MCP server, a colleague's machine, an S3 bucket. For local serving against your own pipeline output (the [Deploying MCP server locally](#deploying-mcp-server-locally) section below), point the server straight at `<output_dir>` and skip distillation entirely.

The full file whitelist and path-scrubbing contract are documented in the [`mcpsrv.bundle`](mcpsrv/bundle.py) module docstring.

## Deploying MCP server locally

```bash
corpus serve   # reads bundle path from config.yaml's output_dir
```

Point your MCP client at this. A project-scoped [.mcp.json](.mcp.json) ships in this repo for Claude Code; for Claude Desktop, edit `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "corpus": {
      "command": "/opt/anaconda3/envs/corpus/bin/corpus",
      "args": ["serve", "--output-dir", "/path/to/output"]
    }
  }
}
```

Restart the client and the corpus tools appear. From there, the [example queries](#example-uses) above become chat messages.

## Deploying MCP server remotely

For sharing a corpus with colleagues or hosting on AWS — bearer-token SSE startup, smoke-test, three-way client-config matrix (Claude Code, Custom Connectors UI, mcp-remote bridge), full EC2 + ALB + S3 runbook — see [DEPLOY.md](DEPLOY.md). One-liner version: `corpus serve --output-dir <bundle> -- --transport sse --host 127.0.0.1 --port 8080 --auth-token-file <token-file>`.

## Additional documentation and resources

- [AGENTS.md](AGENTS.md) — orientation for AI coding agents working in the repo
- [CONTRIBUTING.md](CONTRIBUTING.md) — branching model, release ritual, version bumps
- [CHANGELOG.md](CHANGELOG.md) — what changed in each release
- [dev_docs/OVERVIEW.md](dev_docs/OVERVIEW.md) — pipeline architecture, stage internals, figure pipeline, key files
- [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md) — HPC operational runbook (SLURM, Grobid, job arrays)
- [DEPLOY.md](DEPLOY.md) — AWS deploy runbook (S3 bundle + EC2 systemd)
- [dev_docs/TESTING.md](dev_docs/TESTING.md) — quality test suite, ground truth format, evaluation workflow
- [INSTALL.md](INSTALL.md) — optional OCR extras, pip-only fallback, platform notes
- [dev_docs/MCP_TOOLS.md](dev_docs/MCP_TOOLS.md) — full MCP tool surface
- [dev_docs/PLAN.md](dev_docs/PLAN.md) — roadmap and design decisions
- [dev_docs/clean_install_walkthrough.sh](dev_docs/clean_install_walkthrough.sh) — copy-paste UX walkthrough: fresh env → build → serve, exercising every operator-facing verb at least once
- [dev_docs/PLATFORM_SMOKE.md](dev_docs/PLATFORM_SMOKE.md) — pre-release platform-portability smoke runbook (macOS arm64 + linux-x86_64 Bouchet, clean env recreate); references [dev_docs/ec2_smoke.sh](dev_docs/ec2_smoke.sh) for the highest-signal clean-room linux validation

External:

- [MCP (Model Context Protocol)](https://modelcontextprotocol.io/) — the standard the query interface speaks
- [Docling tutorial](https://youtu.be/9lBTS5dM27c?si=d8cS4wbY6eGXaHH1), based on [this repo](https://github.com/daveebbelaar/ai-cookbook/tree/main/knowledge/docling)
- [Docling figure export docs](https://docling-project.github.io/docling/examples/export_figures/)
- [docling-parse](https://github.com/docling-project/docling-parse)
