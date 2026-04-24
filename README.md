# corpus

A workflow for turning a collection of scientific-literature PDFs — spanning born-digital papers and centuries-old scans in multiple languages — into a queryable knowledge base exposed over MCP. The primary target is the siphonophore literature, but the pipeline is discipline-agnostic as long as a taxonomic / entity layer relevant to the domain can be supplied.

## What this gives you

Once the pipeline has run, the corpus is available through the **MCP server** ([mcp_server.py](mcp_server.py)). [MCP (Model Context Protocol)](https://modelcontextprotocol.io/) is a standard way to hand LLM clients a set of tools they can call; Claude Desktop, Claude Code, Cursor, and Continue all speak it. The server is a read-only view over the per-paper artifacts — no separate store, no sync issues — so re-running the pipeline is the only way new data reaches the tools.

Representative queries you can drive from a chat:

- *"List every valid* Apolemia *species and for each list the papers that discuss them."*
- *"Show the bibliography of Pugh 1997 and mark which cited works are also in the corpus."*
- *"Give me every figure showing nectophore morphology in* Nanomia bijuga *."*
- *"What are the most-cited siphonophore papers that aren't currently in the corpus?"*
- *"Summarize what the corpus says about pneumatophore structure."*

The full tool surface lives in [dev_docs/MCP_TOOLS.md](dev_docs/MCP_TOOLS.md).

## Cross-paper sources of truth

Three precompiled SQLite databases under `resources/` link entity mentions across every paper in the corpus:

| Layer | File | Built by | Status |
| --- | --- | --- | --- |
| **Taxonomy** — WoRMS snapshot with synonymy | `worms.sqlite` | `ingest_worms.py` | operational |
| **Bibliography** — deduplicated works + citation graph | `biblio_authority.sqlite` | `build_biblio_authority.py` + `reconcile_corpus_to_biblio.py` | operational |
| **Geography** — location mentions with coordinates | `locations.sqlite` (planned) | TBD | [PLAN.md §12](PLAN.md) |

Each is independently rebuildable and does not require re-running OCR or extraction. See [PLAN.md §12](PLAN.md) for the design.

## Workflow

The work divides naturally into four stages, best run on different hardware:

### 1. Collect the corpus

PDFs live outside this repo. For the siphonophore corpus I use [`dunnlab/siphonophores`](https://github.com/dunnlab/siphonophores) (Git LFS). Quality matters: **clean high-resolution scans** produce dramatically better OCR and figure extraction than re-scans of photocopies. For born-digital papers, the original publisher PDF beats any downstream export. Filenames should carry `<Surname><Year>` (e.g., `Totton1965a.pdf`) — the reconciler uses this to link corpus papers to their citations when Grobid mis-parses the header.

### 2. Process PDFs on a cluster (the big lift)

Stage 1 (OCR + docling + Grobid + chunking) is CPU-heavy; Pass 3b (vision-LLM figure analysis) and embedding are GPU. One pass over ~2,000 PDFs takes tens of hours wall-clock, even with job-array parallelization. Run it on an HPC cluster.

On Yale's Bouchet (YCRC): [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md) is the operational runbook — SLURM partitions, storage paths, dependency chain, known traps. The one-liner:

```bash
NUM_BATCHES=8 bash slurm/batch_pipeline.sh
```

### 3. Post-process locally

The fast, low-resource steps — bibliography authority DB, taxon mentions, reconciliation, QC — run comfortably on a laptop against the output from stage 2. Copy `output/` down from the cluster, then:

```bash
# Bibliography authority DB.  --enrich-bhl turns on Biodiversity
# Heritage Library lookups for references without a DOI; needs
# BHL_API_KEY and is rate-limited, so long-running (hours).  Skip
# the flag for a ~minute-scale local build against corpus refs only.
python build_biblio_authority.py output --enrich-bhl
python reconcile_corpus_to_biblio.py output

# Taxon mentions DB
python build_taxon_mentions.py output
```

### 4. Deploy the MCP server

**Personal use** — point your MCP client at a local server. A project-scoped [.mcp.json](.mcp.json) ships in the repo for Claude Code; for Claude Desktop, edit `~/Library/Application Support/Claude/claude_desktop_config.json`:

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

**Shared access** — the same `mcp_server.py` speaks SSE over HTTP for remote clients. Three steps:

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

All seven checks should pass locally before you move to AWS. The deployment pattern (EC2 + ALB + CloudFront, bundle pulled from S3) is spelled out in [PLAN.md §10](PLAN.md).

**4. Connect a client.** Remote-MCP support now exists natively in Claude Desktop, claude.ai web, and Claude Code. The paths differ in friction and in which transport + auth they target. Recommended in order:

- **Claude Code (CLI)** — tested and working against our server. One-liner:

  ```bash
  claude mcp add corpus-remote http://127.0.0.1:18080/sse \
      --transport sse \
      --scope user \
      --header "Authorization: Bearer $(cat ~/corpus-mcp.token)"
  ```

  Use `--scope project` to tie to one repo, `claude mcp remove corpus-remote --scope user` to undo. `/mcp` in any Claude Code session lists connected servers and their tools. The CLI accepts `--transport sse` with custom headers, which maps cleanly to our current server.

- **Claude Desktop / claude.ai web (Custom Connectors UI)** — Settings → Connectors → "Add custom connector", paste the URL. Available on Free / Pro / Max / Team / Enterprise plans ([Anthropic docs](https://support.claude.com/en/articles/11175166-get-started-with-custom-connectors-using-remote-mcp)). **Caveat for this repo's server:** Custom Connectors target the **Streamable HTTP** transport with OAuth-style auth, while our server speaks **SSE** with a static bearer token. Whether the UI accepts our URL depends on how strict the client is about the handshake — worth trying (Desktop → ⌘Q restart → add custom connector → paste URL), since it's reversible. If it rejects the SSE endpoint, fall back to the `mcp-remote` bridge below.

- **Claude Desktop fallback — mcp-remote bridge.** If the Custom Connectors UI doesn't connect to our SSE + bearer server, edit `~/Library/Application Support/Claude/claude_desktop_config.json` to launch [`mcp-remote`](https://www.npmjs.com/package/mcp-remote) as a stdio subprocess that bridges to our SSE endpoint:

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

For a cleaner Custom Connectors experience in the future, the server would add a **Streamable HTTP** transport (`mcp.streamable_http_app()` in FastMCP) with OAuth discovery. That's a post-deploy item; tracked against the AWS rollout, not this initial release.

The repo's project-scoped [.mcp.json](.mcp.json) stays useful for local development — it's stdio, no running server or token needed. Use it for coding against the corpus; use one of the above for testing or consuming the deployed pattern.

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

## Quick reference

```bash
# Full pipeline on a small local collection
python process_corpus.py <input_dir> <output_dir> --resume
python embed_chunks.py <output_dir> --resume

# MCP server
python mcp_server.py <output_dir>

# Rebuild the bibliography authority DB (after new papers land)
python build_biblio_authority.py <output_dir>
python reconcile_corpus_to_biblio.py <output_dir>
```

## Documentation

- [dev_docs/OVERVIEW.md](dev_docs/OVERVIEW.md) — pipeline architecture, stage internals, figure pipeline, key files
- [dev_docs/BOUCHET.md](dev_docs/BOUCHET.md) — HPC operational runbook (SLURM, Grobid, job arrays)
- [dev_docs/TESTING.md](dev_docs/TESTING.md) — quality test suite, ground truth format, evaluation workflow
- [dev_docs/INSTALL.md](dev_docs/INSTALL.md) — optional OCR extras, pip-only fallback, platform notes
- [dev_docs/MCP_TOOLS.md](dev_docs/MCP_TOOLS.md) — full MCP tool surface
- [PLAN.md](PLAN.md) — roadmap and design decisions

## Resources

- [Docling tutorial](https://youtu.be/9lBTS5dM27c?si=d8cS4wbY6eGXaHH1), based on [this repo](https://github.com/daveebbelaar/ai-cookbook/tree/main/knowledge/docling)
- [Docling figure export docs](https://docling-project.github.io/docling/examples/export_figures/)
- [docling-parse](https://github.com/docling-project/docling-parse)
