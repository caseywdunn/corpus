# corpus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19964909.svg)](https://doi.org/10.5281/zenodo.19964909)

A workflow for turning a collection of scientific-literature PDFs — spanning born-digital papers and centuries-old scans in multiple languages — into a queryable knowledge base exposed over the [Model Context Protocol (MCP)](https://modelcontextprotocol.io/), so AI agents can interrogate and synthesize the literature directly. The primary target is taxonomic literature.

## Citing corpus

We are preparing a manuscript that presents corpus. In the mean time, please cite corpus as:

> Dunn, C. W., Zapata, F., & Church, S. H. (2026). Corpus: Model Context Protocol (MCP) server for taxonomic literature. <https://doi.org/10.5281/zenodo.19964909>

## What corpus does

corpus reads a folder of PDFs and produces a knowledge base built around: **taxa, figures, bibliographies, and the cross-paper relationships between them**. It handles several tasks:

- **OCR for old scans** — including 19th-century German Fraktur and other historical typefaces, with per-paper language detection.
- **Figure extraction** — every plate, photograph, and line drawing pulled out with its caption, then vision-LLM-tagged for taxonomy and anatomy.
- **Bibliography parsing and reconciliation** — references inside each paper are deduplicated across the whole corpus, so, for example, "Totton 1965" is one entity even when twenty papers cite it differently.
- **Taxonomy linking** — every species mention is matched against a Darwin Core taxonomy with full synonymy, so a question about *Apolemia uvaria* finds papers that only ever called it *Stephanomia uvaria*.

The output is a per-paper artifact tree plus cross-paper databases. You query it through an **MCP server** ([mcp_server.py](mcp_server.py)) that any [MCP](https://modelcontextprotocol.io/) client (Claude Desktop, Claude Code, claude.ai web) can connect to — so "show me every figure of *Nanomia bijuga* nectophores in the corpus" becomes something you ask in chat instead of grep across a hard drive. The server is a read-only view over the per-paper artifacts, so re-running the pipeline is the only way new data reaches the tools.

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

```bash
python process_corpus.py <input_dir> <output_dir> --bib references.bib --resume
```

Matching is on the basename (case-insensitive), so the `file` value can be a bare filename or a full path. This is the cleanest way to get accurate per-paper metadata where you've already curated references for your own writing.

### External taxonomic data (optional)

The taxonomy database (`taxonomy.sqlite`) is a [Darwin Core](https://dwc.tdwg.org/) snapshot that drives synonymy resolution. By default it builds from a Darwin Core Archive of your choice — for siphonophores we use the [WoRMS](https://www.marinespecies.org/) export. Other groups: pull a DwC-A from the relevant authority ([GBIF](https://www.gbif.org/), [ITIS](https://www.itis.gov/), [Catalogue of Life](https://www.catalogueoflife.org/)) and ingest into your corpuscle:

```bash
python ingest_taxonomy.py <output_dir> --source dwca --input path/to/dwca.zip
```

Without external taxonomy the pipeline still extracts taxon mentions from text — you lose the synonymy graph that links historical names to current valid names.

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

Pass the file with `--lexicon path/to/lexicon.yaml`. Each category emits its own `<hash>/<category>.json` artifact. Per-category content is fingerprinted independently, so editing a single section only invalidates that category's annotations on `--resume` — the rest stay cached.

Without `--lexicon`, lexicon extraction is skipped entirely. Like `--bib`, the lexicon is an input you maintain alongside your literature, not something the tool ships.

### Instructions for the LLM (optional)

A markdown file at `<corpuscle>/instructions.md` whose contents land in every chat session against the corpus. The MCP server returns it in `InitializeResult.instructions`, and well-behaved clients (Claude Desktop, Claude Code) inject it into the LLM's context at session start.

Use it for per-corpus nudges that no taxonomy or lexicon entry can express — see [demo/instructions.md](demo/instructions.md) for a worked example, which (among other things) tells the model that *Velella* and *Porpita* are not siphonophores even when older literature lumps them in.

Override the default location with `--instructions <path>` when starting the MCP server. If the file is absent, no instructions are sent.

## Computational requirements

- **Disk.** Plan on several times the size of the original PDFs.
- **CPU vs. GPU.** Stage 1 (OCR, layout, Grobid, chunking, annotation) is CPU-bound and parallelizes well across PDFs. Stage 2 (BGE-M3 embeddings) and the optional vision pass (Qwen2.5-VL-7B figure tagging) are GPU-accelerated; embeddings still run on CPU but slowly, and the local vision backend needs a GPU to be usable. A Claude API vision backend is also available — costs API credits but runs anywhere.
- **Scale.** A few dozen PDFs run on a laptop in a couple of hours. A few thousand benefit from an HPC cluster — we use Yale's Bouchet, runbook in [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md).
- **MCP client.** The query interface is MCP, so you'll need a client that speaks it (Claude Desktop, Claude Code, claude.ai web with custom connectors, Cursor, Continue). Most require an Anthropic subscription.
- **Remote deployment.** Serving the corpus to others over the network requires a server. The reference deploy is AWS (EC2 + nginx + Let's Encrypt), but any host with Python and an open port works. See [Deploying MCP server remotely](#deploying-mcp-server-remotely) below.

## Corpuscle layout

A *corpuscle* is the on-disk container for one corpus instance — siphonophores, drosophila, or whatever group you're working on. It's a single directory:

```text
<corpuscle>/
├── documents/<HASH>/         # per-paper artifacts (text, chunks, figures, taxa, anatomy, …)
├── vector_db/lancedb/        # embeddings index
├── taxonomy.sqlite           # Darwin Core snapshot, built by ingest_taxonomy.py
├── biblio_authority.sqlite   # deduplicated works + citation graph
└── taxon_mentions.sqlite     # cross-paper taxon index
```

Two cross-paper layers are built from the per-paper artifacts and are independently rebuildable without re-running OCR or extraction:

| Layer | File | Built by |
| --- | --- | --- |
| **Taxonomy** — Darwin Core snapshot with synonymy | `taxonomy.sqlite` | `ingest_taxonomy.py` (sources: `dwc` / `dwca` / `worms`) |
| **Bibliography** — deduplicated works + citation graph | `biblio_authority.sqlite` | `build_biblio_authority.py` + `reconcile_corpus_to_biblio.py` |

Every CLI takes the corpuscle root as its first positional argument and resolves all per-instance files from there. Run two corpora side-by-side by giving each its own corpuscle directory; they don't share state.

## Installation

```bash
conda env create -f environment.yaml
conda activate corpus
```

Grobid runs as a separate service:

```bash
docker compose up -d grobid
curl http://localhost:8070/api/isalive   # should print "true"
```

Platform-specific OCR extras (Fraktur for 19th-century German, additional `pngquant`/`jbig2enc` compression) are covered in [dev_docs/INSTALL.md](dev_docs/INSTALL.md).

## Language support

OCR is the only stage where language matters. Each scanned PDF gets its language detected automatically (`langdetect` on the first few pages of extracted text), and the result is routed to a matching [Tesseract](https://github.com/tesseract-ocr/tesseract) language pack — a `<code>.traineddata` file that teaches Tesseract to recognize a particular language or script. Born-digital PDFs with a clean text layer skip OCR entirely, so packs only matter for scans.

**What ships by default.** `environment.yaml` installs `tesseract` plus the language packs that back the default fallback set in `config.yaml` (`ocr.ocr_languages_default`): `eng`, `deu`, `fra`, `rus`, `lat`, `spa`, `por`, `chi_sim`, `chi_tra`, `jpn`, `ell`, `kor`. The conda env also pulls down `grc` (Ancient Greek), but it's intentionally left out of the default fallback union — opt in via `config.yaml` if you have classical sources.

**One manual step required.** 19th-century German Fraktur (`deu_latf`) is part of the default fallback set but isn't packaged on conda-forge, so it has to be downloaded into the conda env's `tessdata` directory by hand after the env is built. Without it, scanned 19th-c. German papers (Goldfuss 1820, Brandt 1837, Pagenstecher 1869, Dönitz 1871, …) OCR to whitespace. Steps in [dev_docs/INSTALL.md](dev_docs/INSTALL.md#additional-ocr-language-packs).

**How to tell if you need more packs.** After a Stage 1 run, `<output_dir>/documents/<HASH>/scan_detection.json` records the detected language for each scanned PDF. Skim those to see what languages your corpus actually contains. If a language is detected but no matching pack is installed, OCR silently falls back to the union of installed packs (English is always appended last) — extraction still succeeds, but accuracy on that paper drops. Adding more languages than you need is safe but slows OCR, so trim `ocr.ocr_languages_default` for a known-narrow corpus.

**Adding a language outside the default set.** The full ISO-to-Tesseract map covers ~40 languages including CJK, Cyrillic, Greek, Arabic, and Indic scripts — see `_ISO_TO_TESSERACT` in [pipeline/scan.py](pipeline/scan.py). Install the matching `tesseract-data-<code>` pack from conda-forge and the detector will pick it up on the next run. To make it part of the fallback set tried when detection is uncertain, add the code to `ocr.ocr_languages_default` in `config.yaml`.

**Everything past OCR is language-agnostic.** The [BGE-M3](https://huggingface.co/BAAI/bge-m3) embedding model is multilingual (trained on 100+ languages) and built for cross-lingual retrieval — a question in English finds passages in German, Latin, French, and so on. Grobid metadata extraction works across European languages, and the vision LLM that tags figures reads labels in whatever script is on the plate. So once OCR succeeds, no further language configuration is needed.

## Try it on the demo corpus

The repo ships with 11 siphonophore PDFs, a matching BibTeX, and an example anatomy lexicon under [demo/](demo/). Running the whole pipeline against them takes a few minutes on a laptop, with two uses:

- **First-run smoke test** — confirm your install works before processing your own corpus.
- **Regression check on top of a real corpus** — `output_demo/` is its own corpuscle, so it lives alongside (and never clobbers) your production corpuscle.

```bash
# Build a small WoRMS taxonomy snapshot into the demo corpuscle.  Skip
# this step if you only want to smoke-test stages 1 + 2.
python ingest_taxonomy.py output_demo --source worms --root-id 1371

# update_corpus.py runs the pipeline + post-pipeline scripts in
# dependency order — equivalent to chaining process_corpus.py,
# embed_chunks.py, build_biblio_authority.py,
# build_taxon_mentions.py, backfill_intext_citations.py, and
# reconcile_corpus_to_biblio.py with --resume throughout.
python update_corpus.py demo output_demo \
    --bib demo/siphonophores.bib \
    --lexicon demo/lexicon.yaml \
    --resume
```

Every cross-paper SQLite lands inside `output_demo/` automatically — no override flags needed.

## First time run

The work divides into two stages, best run on different hardware. Stage 1 is CPU-heavy and embarrassingly parallel; Stage 2 wants a GPU. Each unique PDF is identified by the first 12 hex chars of its SHA-256, and all artifacts live under `<output_dir>/documents/<HASH>/`.

The fastest path is `update_corpus.py`, which runs the pipeline and the four post-pipeline scripts in dependency order with `--resume` throughout — see [update_corpus.py](update_corpus.py) for the full list of steps:

```bash
python update_corpus.py <input_dir> <output_dir> --resume
```

Run-by-run, that command is equivalent to:

```bash
# Stage 1 — OCR, docling, Grobid, chunking, annotation (CPU, parallelizable)
python process_corpus.py <input_dir> <output_dir> --resume

# Stage 2 — embeddings (GPU helpful, not required)
python embed_chunks.py <output_dir> --resume

# Cross-paper databases (fast — minutes at corpus scale)
python build_biblio_authority.py <output_dir>
python build_taxon_mentions.py <output_dir>
python backfill_intext_citations.py <output_dir>
python reconcile_corpus_to_biblio.py <output_dir>
```

`--resume` skips work whose artifact is already on disk — per-stage inside Stage 1 (so a lexicon edit only re-runs the annotation pass, not Docling) and per-paper across Stage 2. For large corpora, run stage 1 on a cluster — the SLURM one-liner on Bouchet is:

```bash
NUM_BATCHES=8 bash slurm/batch_pipeline.sh
```

To enrich pre-DOI references against the Biodiversity Heritage Library, add `--enrich-bhl` to the bibliography step (needs `BHL_API_KEY`, rate-limited at hours not minutes):

```bash
python build_biblio_authority.py <output_dir> --enrich-bhl
```

After a run, [`corpus_status.py`](corpus_status.py) reports stage completion, quality flags, and which papers have stale annotations relative to the current lexicon / taxonomy snapshot:

```bash
python corpus_status.py <output_dir>
```

## Vision pass (optional)

The default pipeline extracts figures with their captions and bounding boxes but doesn't open the images. **Pass 3b** runs a vision-language model over each extracted figure to detect multi-panel structure, recover compound-figure splits, and tag panels with ROI bounding boxes. It's a separate, opt-in pass because it's expensive — GPU minutes per paper for a local model, or API tokens for the remote backend.

```bash
# Local backend — needs a CUDA or MPS GPU; uses Qwen2.5-VL-7B-Instruct
python process_corpus.py <input_dir> <output_dir> --resume --vision-backend local

# Remote backend — needs ANTHROPIC_API_KEY; runs anywhere
python process_corpus.py <input_dir> <output_dir> --resume --vision-backend claude
```

`--refresh-vision` re-runs only Pass 3b on hashes that already have a `summary.json`, so you can layer vision over an existing corpus without re-doing OCR / Docling / Grobid. On Bouchet, [slurm/batch_pass3b.sh](slurm/batch_pass3b.sh) does the same job as a SLURM array on `gpu_h200` — see [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md).

## Adding and updating documents

Drop new PDFs into `<input_dir>` and re-run `update_corpus.py`. `--resume` recognizes existing artifacts by SHA-256 hash, so only new or changed PDFs are reprocessed:

```bash
python update_corpus.py <input_dir> <output_dir> --resume
```

When the lexicon or taxonomy snapshot changes (not the PDFs), just re-run with `--resume`. Each per-paper `pipeline_state.json` records the input fingerprint of every stage; a mismatch forces only the affected stage to re-run, so editing one lexicon section regenerates that category's annotations and leaves everything else cached:

```bash
python update_corpus.py <input_dir> <output_dir> --resume \
    --lexicon path/to/lexicon.yaml
```

To remove papers, delete the PDFs from `<input_dir>` and run `python process_corpus.py <input_dir> <output_dir> --audit-orphans` for a read-only listing of `documents/<HASH>/` directories whose source PDF is gone (deletion stays manual).

The cross-paper databases (bibliography, taxon mentions) rebuild from scratch each time — fast at corpus scale, and the only way they pick up new cross-references.

## Curating bibliographic metadata

Grobid will mis-parse some references — wrong year on a 19th-c. monograph, scrambled author list, dropped subtitle. Rather than hand-editing per-paper `metadata.json` files (which get clobbered on the next pipeline run), corpus supports a **BibTeX round-trip** against `biblio_authority.sqlite`:

```bash
# Export the authority database as BibTeX
python bib_export.py <output_dir> -o my_corpus.bib

# Hand-edit my_corpus.bib (fix titles, years, authors, …)

# Apply the edits back into the authority DB
python bib_import.py <output_dir> my_corpus.bib
```

Each entry carries a stable `corpus_hash` field that bib_import uses to match edits back to the right `works` row, with DOI as a fallback. The `.bib` is the source of truth for user edits — preserve it across pipeline rebuilds and re-run `bib_import.py` afterward to re-apply your curation. Hand-edited entries override Grobid-derived metadata in every MCP tool that reads from the authority DB.

`--bib` (the input flag described in [Bibliographic data for pdfs](#bibliographic-data-for-pdfs-optional)) and the round-trip above are independent: the input flag overrides Grobid's per-paper header from a curated BibTeX *before* the authority DB is built; the round-trip edits the authority DB *after*. Use whichever fits where you already maintain the metadata.

## Distilling a served bundle

The corpuscle the pipeline emits is the **build bundle**: everything `process_corpus.py` produces, including `processed.pdf`, raw docling dumps, per-page QC visualizations, and per-paper logs. For ~2,000 papers that's ~10 GB. The MCP server doesn't need most of it — only the JSON artifacts it reads at startup, the figure PNGs, and the precompiled indices. [`package_for_serve.py`](package_for_serve.py) distills the build bundle down to a **served bundle** (~3 GB for the same corpus) by copying only whitelisted files, scrubbing absolute paths from JSON values, and writing a versioned `bundle_manifest.json` that the MCP `bundle_info` tool surfaces:

```bash
python package_for_serve.py <output_dir> <serve_bundle_dir> --version v1.0.0
```

You only need this when **shipping a corpus to a different host** — a remote MCP server, a colleague's machine, an S3 bucket. For local serving against your own pipeline output (the [Deploying MCP server locally](#deploying-mcp-server-locally) section below), point the server straight at `<output_dir>` and skip distillation entirely.

The full file whitelist and path-scrubbing contract are documented in the [`package_for_serve.py`](package_for_serve.py) module docstring.

## Deploying MCP server locally

Point your MCP client at a local server. A project-scoped [.mcp.json](.mcp.json) ships in this repo for Claude Code; for Claude Desktop, edit `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "corpus": {
      "command": "/opt/anaconda3/envs/corpus/bin/python",
      "args": ["/path/to/corpus/mcp_server.py", "/path/to/output"]
    }
  }
}
```

Restart the client and the corpus tools appear. From there, the [example queries](#example-uses) above become chat messages.

## Deploying MCP server remotely

The same `mcp_server.py` speaks SSE over HTTP for remote clients — useful for sharing a corpus with colleagues or running it from a tablet. The remote host serves a [distilled bundle](#distilling-a-served-bundle), not the full build output, so plan to run `package_for_serve.py` on the build host and ship the resulting `<serve_bundle_dir>` (typically via `aws s3 sync`) before pointing the remote server at it. The full AWS runbook is in [dev_docs/DEPLOY.md](dev_docs/DEPLOY.md); the steps below cover the minimum to bring up an SSE server on any host you already have.

**1. Generate a bearer token.** Use `secrets.token_urlsafe` for a strong URL-safe random string, save to a mode-600 file, never pass it on the CLI (would leak via `ps`):

```bash
python -c 'import secrets; print(secrets.token_urlsafe(32))' > ~/corpus-mcp.token
chmod 600 ~/corpus-mcp.token
```

**2. Start the server over SSE.** `output` is your pipeline output directory — replace with an absolute path if running from elsewhere:

```bash
python mcp_server.py output \
    --transport sse --host 127.0.0.1 --port 8080 \
    --auth-token-file ~/corpus-mcp.token
```

Clients send `Authorization: Bearer <token>` on every request. Without `--auth-token-file` or `CORPUS_MCP_TOKEN` the server runs open and logs a loud warning — fine for localhost experiments, never safe for public-facing deploys.

**3. Smoke-test before exposing to anyone else.** [tools/smoke_test_sse.py](tools/smoke_test_sse.py) launches its own server on a free port, generates its own token, and drives the full stack — 401 without auth, 200 with auth, MCP `initialize` → `list_tools` → `bundle_info` → `list_papers`:

```bash
python tools/smoke_test_sse.py output
```

All seven checks should pass locally before you move to a real server. The deployment pattern (EC2 + ALB + CloudFront, bundle pulled from S3) is in [dev_docs/DEPLOY.md](dev_docs/DEPLOY.md) and [dev_docs/PLAN.md §10](dev_docs/PLAN.md).

**4. Connect a client.** Remote-MCP support exists natively in Claude Desktop, claude.ai web, and Claude Code. The paths differ in friction and in which transport + auth they target. Recommended in order:

- **Claude Code (CLI)** — tested and working against our server:

  ```bash
  claude mcp add corpus-remote http://127.0.0.1:18080/sse \
      --transport sse \
      --scope user \
      --header "Authorization: Bearer $(cat ~/corpus-mcp.token)"
  ```

  Use `--scope project` to tie to one repo, `claude mcp remove corpus-remote --scope user` to undo. `/mcp` in any session lists connected servers and their tools.

- **Claude Desktop / claude.ai web (Custom Connectors UI)** — Settings → Connectors → "Add custom connector", paste the URL. Available on Free / Pro / Max / Team / Enterprise plans ([Anthropic docs](https://support.claude.com/en/articles/11175166-get-started-with-custom-connectors-using-remote-mcp)). **Caveat:** Custom Connectors target the **Streamable HTTP** transport with OAuth-style auth, while our server speaks **SSE** with a static bearer token. Whether the UI accepts our URL depends on how strict the client is — worth trying. If it rejects the SSE endpoint, fall back to the bridge below.

- **Claude Desktop fallback — mcp-remote bridge.** Edit `~/Library/Application Support/Claude/claude_desktop_config.json` to launch [`mcp-remote`](https://www.npmjs.com/package/mcp-remote) as a stdio subprocess that bridges to the SSE endpoint:

  ```json
  {
    "mcpServers": {
      "corpus-remote": {
        "command": "/opt/homebrew/bin/npx",
        "args": ["-y", "mcp-remote",
                 "http://127.0.0.1:18080/sse",
                 "--header", "Authorization: Bearer <token>"]
      }
    }
  }
  ```

A cleaner Custom Connectors experience would add a **Streamable HTTP** transport (`mcp.streamable_http_app()` in FastMCP) with OAuth discovery, but this is deferred indefinitely — SSE + bearer token works for the ~20-collaborator deploy target.

## Additional documentation and resources

- [AGENTS.md](AGENTS.md) — orientation for AI coding agents working in the repo
- [CONTRIBUTING.md](CONTRIBUTING.md) — branching model, release ritual, version bumps
- [CHANGELOG.md](CHANGELOG.md) — what changed in each release
- [dev_docs/OVERVIEW.md](dev_docs/OVERVIEW.md) — pipeline architecture, stage internals, figure pipeline, key files
- [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md) — HPC operational runbook (SLURM, Grobid, job arrays)
- [dev_docs/DEPLOY.md](dev_docs/DEPLOY.md) — AWS deploy runbook (S3 bundle + EC2 systemd)
- [dev_docs/TESTING.md](dev_docs/TESTING.md) — quality test suite, ground truth format, evaluation workflow
- [dev_docs/INSTALL.md](dev_docs/INSTALL.md) — optional OCR extras, pip-only fallback, platform notes
- [dev_docs/MCP_TOOLS.md](dev_docs/MCP_TOOLS.md) — full MCP tool surface
- [dev_docs/PLAN.md](dev_docs/PLAN.md) — roadmap and design decisions

External:

- [MCP (Model Context Protocol)](https://modelcontextprotocol.io/) — the standard the query interface speaks
- [Docling tutorial](https://youtu.be/9lBTS5dM27c?si=d8cS4wbY6eGXaHH1), based on [this repo](https://github.com/daveebbelaar/ai-cookbook/tree/main/knowledge/docling)
- [Docling figure export docs](https://docling-project.github.io/docling/examples/export_figures/)
- [docling-parse](https://github.com/docling-project/docling-parse)
