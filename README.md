# corpus

A workflow for turning a collection of scientific-literature PDFs — spanning born-digital papers and centuries-old scans in multiple languages — into a queryable knowledge base exposed over MCP. The primary target is the siphonophore literature, but the pipeline is discipline-agnostic as long as a taxonomic / entity layer relevant to the domain can be supplied.

## What this gives you

Once the pipeline has run, the corpus is available through the **MCP server** ([mcp_server.py](mcp_server.py)). [MCP (Model Context Protocol)](https://modelcontextprotocol.io/) is a standard way to hand LLM clients a set of tools they can call; Claude Desktop, Claude Code, Cursor, and Continue all speak it. The server is a read-only view over the per-paper artifacts — no separate store, no sync issues — so re-running the pipeline is the only way new data reaches the tools.

Representative queries you can drive from a chat:

- *"List every valid* Apolemia *species and flag which have papers in the corpus."*
- *"Show the bibliography of Dunn 2005 and mark which cited works are also in the corpus."*
- *"Give me every figure showing nectophore morphology in* Nanomia bijuga *papers."*
- *"What are the most-cited siphonophore papers that aren't currently in the corpus?"*
- *"Summarize what the corpus says about pneumatophore structure, with paper hashes."*

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

**Shared access** — the server uses stdio transport, so exposing it to others means fronting it with a small HTTP bridge and deploying to a cloud VM. An AWS deployment pattern (EC2 + nginx + auth) is on the roadmap; see [PLAN.md §10](PLAN.md) for the auth/access control plan.

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
