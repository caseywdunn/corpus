# corpus

A workflow for turning a collection of scientific-literature PDFs — spanning born-digital papers and centuries-old scans in multiple languages — into a queryable knowledge base exposed over MCP. The primary target is taxonomic literature.

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
- *"What are the most-cited siphonophore papers that aren't currently in the corpus?"*
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

### Anatomy lexicon (optional)

A small YAML file that lists the anatomical terms you want tagged in figure captions and chunk text — *nectophore*, *pneumatophore*, *gastrozooid* for siphonophores, or whatever vocabulary fits your group. Pass it with `--anatomy-lexicon`; without it, anatomy extraction is skipped.

The format is a flat term-to-metadata map; see [demo/anatomy_lexicon.yaml](demo/anatomy_lexicon.yaml) for a worked example. Like `--bib`, the lexicon is something you maintain alongside your literature, not something the tool ships.

## Computational requirements

- **Disk.** Plan on several times the size of the original PDFs.
- **CPU vs. GPU.** Stage 1 (OCR, layout, Grobid, chunking) is CPU-bound and parallelizes well across PDFs. Stage 2 (vision-LLM figure analysis and BGE-M3 embeddings) is GPU-accelerated; without a GPU these steps still run but are much slower.
- **Scale.** A few dozen PDFs run on a laptop in a couple of hours. A few thousand benefit from an HPC cluster — we use Yale's Bouchet, runbook in [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md).
- **MCP client.** The query interface is MCP, so you'll need a client that speaks it (Claude Desktop, Claude Code, claude.ai web with custom connectors, Cursor, Continue). Most require an Anthropic subscription.
- **Remote deployment.** Serving the corpus to others over the network requires a server. The reference deploy is AWS (EC2 + ALB + CloudFront), but any host with Python and an open port works. See [Deploying remotely](#deploying-remotely) below.

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

## Try it on the demo corpus

The repo ships with 11 siphonophore PDFs, a matching BibTeX, and an example anatomy lexicon under [demo/](demo/). Running the whole pipeline against them takes a few minutes on a laptop, with two uses:

- **First-run smoke test** — confirm your install works before processing your own corpus.
- **Regression check on top of a real corpus** — `output_demo/` is its own corpuscle, so it lives alongside (and never clobbers) your production corpuscle.

```bash
# Build a small WoRMS taxonomy snapshot into the demo corpuscle.  Skip
# this step if you only want to smoke-test stages 1 + 2.
python ingest_taxonomy.py output_demo --source worms --root-id 1371

python process_corpus.py demo output_demo \
    --bib demo/siphonophores.bib \
    --anatomy-lexicon demo/anatomy_lexicon.yaml \
    --resume
python embed_chunks.py output_demo --resume

python build_biblio_authority.py output_demo
python reconcile_corpus_to_biblio.py output_demo
python build_taxon_mentions.py output_demo
```

Every cross-paper SQLite lands inside `output_demo/` automatically — no override flags needed.

## First time run

The work divides into two stages, best run on different hardware. Stage 1 is CPU-heavy and embarrassingly parallel; Stage 2 wants a GPU. Each unique PDF is identified by the first 12 hex chars of its SHA-256, and all artifacts live under `<output_dir>/documents/<HASH>/`.

```bash
# Stage 1 — OCR, docling, Grobid, chunking (CPU, parallelizable)
python process_corpus.py <input_dir> <output_dir> --resume

# Stage 2 — embeddings (GPU helpful, not required)
python embed_chunks.py <output_dir> --resume
```

`--resume` skips PDFs that already have a `summary.json`, so you can interrupt and restart without losing work. For large corpora, run stage 1 on a cluster — the SLURM one-liner on Bouchet is:

```bash
NUM_BATCHES=8 bash slurm/batch_pipeline.sh
```

Once Stage 1 and 2 are done, build the cross-paper databases. These are fast (minutes against a few thousand papers) and run comfortably on a laptop:

```bash
python build_biblio_authority.py <output_dir>
python reconcile_corpus_to_biblio.py <output_dir>
python build_taxon_mentions.py <output_dir>
```

To enrich references that lack DOIs against the Biodiversity Heritage Library, add `--enrich-bhl`:

```bash
python build_biblio_authority.py <output_dir> --enrich-bhl
```

This needs `BHL_API_KEY` and is rate-limited (hours, not minutes). Skip it for a fast local build against in-corpus references only.

## Adding and updating documents

Drop new PDFs into `<input_dir>` and re-run the same commands. `--resume` recognizes existing artifacts by SHA-256 hash, so only new or changed PDFs are processed:

```bash
python process_corpus.py <input_dir> <output_dir> --resume
python embed_chunks.py <output_dir> --resume
python build_biblio_authority.py <output_dir>
python reconcile_corpus_to_biblio.py <output_dir>
python build_taxon_mentions.py <output_dir>
```

The cross-paper databases (bibliography, taxon mentions) rebuild from scratch each time — fast at corpus scale, and the only way they pick up new cross-references.

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

The same `mcp_server.py` speaks SSE over HTTP for remote clients — useful for sharing a corpus with colleagues or running it from a tablet.

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

All seven checks should pass locally before you move to a real server. The deployment pattern (EC2 + ALB + CloudFront, bundle pulled from S3) is in [dev_docs/DEPLOY.md](dev_docs/DEPLOY.md) and [PLAN.md §10](PLAN.md).

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

For a cleaner Custom Connectors experience in the future, the server would add a **Streamable HTTP** transport (`mcp.streamable_http_app()` in FastMCP) with OAuth discovery. That's a post-deploy item, tracked against the AWS rollout.

## Gotchas



## Additional documentation and resources

- [dev_docs/OVERVIEW.md](dev_docs/OVERVIEW.md) — pipeline architecture, stage internals, figure pipeline, key files
- [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md) — HPC operational runbook (SLURM, Grobid, job arrays)
- [dev_docs/DEPLOY.md](dev_docs/DEPLOY.md) — AWS deploy runbook (S3 bundle + EC2 systemd)
- [dev_docs/TESTING.md](dev_docs/TESTING.md) — quality test suite, ground truth format, evaluation workflow
- [dev_docs/INSTALL.md](dev_docs/INSTALL.md) — optional OCR extras, pip-only fallback, platform notes
- [dev_docs/MCP_TOOLS.md](dev_docs/MCP_TOOLS.md) — full MCP tool surface
- [PLAN.md](PLAN.md) — roadmap and design decisions

External:

- [MCP (Model Context Protocol)](https://modelcontextprotocol.io/) — the standard the query interface speaks
- [Docling tutorial](https://youtu.be/9lBTS5dM27c?si=d8cS4wbY6eGXaHH1), based on [this repo](https://github.com/daveebbelaar/ai-cookbook/tree/main/knowledge/docling)
- [Docling figure export docs](https://docling-project.github.io/docling/examples/export_figures/)
- [docling-parse](https://github.com/docling-project/docling-parse)
