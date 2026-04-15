# corpus

A workflow for interrogating a corpus of scientific literature pdfs, including old (scanned, potentially spanning centuries of formatting conventions) and new (born-digital) papers. My own goal is analysis of the siphonophore literature, such as:

- Searching
- Citation analysis (including discovery of relevant manuscripts not in the corpus)
- Collecting anatomical figures of the same species across different papers

## System Overview

The corpus processing system uses **hash-based organization** similar to Git, where each unique PDF (identified by SHA-256 hash) gets processed once regardless of how many copies exist in the input directory structure.

### Key Features

- **Smart Deduplication**: Identical PDFs are only processed once
- **Complete Provenance**: Tracks all locations where each PDF was found
- **Flexible Input**: Works with any directory structure
- **Resumable Processing**: Can resume interrupted workflows
- **Git-like Organization**: Hash-based file naming for unique identification

### Processing Pipeline

1. **PDF Discovery**: Recursively finds all PDFs in input directory
2. **Hash-based Deduplication**: Groups identical PDFs by content hash
3. **Scan Detection**: Determines if PDF needs OCR using text extraction analysis
4. **PDF Preparation**: Runs OCR with `ocrmypdf` if needed
5. **Content Extraction**: Uses `docling` to extract text and figures
6. **Metadata Extraction**: Extracts bibliographic metadata with `Grobid`
7. **Text Chunking**: Splits text into searchable chunks
8. **Vector Embedding**: Creates embeddings and stores in vector database


## Installation

Create the conda environment — this installs Python, tesseract, ghostscript, pymupdf, docling, lancedb, and the rest of the Python dependencies:

```bash
conda env create -f environment.yaml
conda activate corpus
```

Updating an existing environment after `environment.yaml` changes:

```bash
conda env update -f environment.yaml --prune
```

If a previous attempt left a partial `corpus` environment behind, remove it first with `conda env remove -n corpus`.

### Optional: pngquant + jbig2enc for higher OCR optimization

`ocrmypdf` is installed via pip from the environment above. It works without any extra native tools, but enabling its highest optimization level requires `pngquant` and `jbig2enc`, which aren't on conda-forge for all platforms (notably `osx-arm64`). Install them at the system level if you want smaller output PDFs:

```bash
# macOS
brew install pngquant jbig2enc

# Debian/Ubuntu
sudo apt-get install pngquant jbig2enc
```

On Bouchet, `module avail pngquant` will tell you whether a module is available; if not, it's safe to skip — `ocrmypdf` falls back gracefully.

### Additional OCR language packs

The base `tesseract` install includes only English. For the multilingual / historical portion of the corpus, add the language data you need:

```bash
# Modern European languages + Latin
conda install -c conda-forge tesseract-data-deu tesseract-data-fra tesseract-data-rus tesseract-data-lat
```

**19th-century German Fraktur** (e.g., Haeckel, Schneider) is not on conda-forge. Download `deu_latf.traineddata` from the Tesseract project and drop it into the env's `tessdata` directory:

```bash
TESSDATA="$CONDA_PREFIX/share/tessdata"
curl -L -o "$TESSDATA/deu_latf.traineddata" \
  https://github.com/tesseract-ocr/tessdata_best/raw/main/deu_latf.traineddata
```

### OpenAI API key (transitional)

The current `embed_chunks.py` uses OpenAI; a local open-weights backend run on a GPU is in progress (see [PLAN.md §7](PLAN.md)). Until that lands, put your key in a `.env` file at the repo root:

```
OPENAI_API_KEY=sk-...
```

### On Bouchet (YCRC HPC)

Same file, same commands, after loading the module system:

```bash
module load miniconda
conda env create -f environment.yaml
conda activate corpus
```

See [PLAN.md §7](PLAN.md) for GPU partition guidance when running the embedding stage on Bouchet.

### Running Grobid

The pipeline uses [Grobid](https://grobid.readthedocs.io/) to extract bibliographic metadata and references from each PDF. Start it with Docker Compose:

```bash
docker compose up -d grobid
# Wait ~60 s for first-run startup; check status with:
docker compose ps
curl http://localhost:8070/api/isalive   # should print "true"
```

Stop it with `docker compose down` when you're done.

If Grobid isn't reachable when `process_corpus.py` runs, the pipeline logs a warning and falls back to placeholder metadata — the rest of the steps (text, figures, chunks) still produce useful output. To explicitly skip Grobid (for rapid iteration on non-metadata stages):

```bash
python process_corpus.py <input> <output> --no-grobid
```

On Bouchet, [Singularity](https://docs.sylabs.io/) can pull the same image:

```bash
singularity build grobid.sif docker://lfoppiano/grobid:0.8.1
singularity run --bind $HOME grobid.sif &
# then point the pipeline at the service's URL
python process_corpus.py <input> <output> --grobid-url http://localhost:8070
```

### Pip-only fallback

If you can't use conda, `requirements.txt` still works but you'll need to install the system tools yourself (`brew install ghostscript tesseract pngquant jbig2enc` on macOS):

```bash
pip install -r requirements.txt
```

## Usage

### Basic Processing

Process all PDFs in an input directory:

```bash
python process_corpus.py /path/to/input/directory /path/to/output/directory
```

**Input Directory**: Can contain PDFs in any nested structure. Examples:
- `/Users/researcher/papers/` - Directory with PDFs and subdirectories
- `/Volumes/PDFs/siphonophore_literature/` - External drive with paper collection
- `./my_pdfs/` - Relative path to local PDF collection

**Output Directory**: Will be created with this structure:
```
output_directory/
├── documents/           # Processed documents organized by hash
│   ├── A1B2C3D4/       # Hash-based directory for unique PDF
│   │   ├── processed.pdf        # OCR'd version (if needed)
│   │   ├── figures/             # Extracted figure images
│   │   ├── visualizations/      # Page visualization PNGs
│   │   │   ├── page_1_visualization.png  # Red boxes for words, yellow/orange for figures
│   │   │   └── page_2_visualization.png
│   │   ├── text.json           # Extracted text content
│   │   ├── figures.json        # Figure metadata and captions
│   │   ├── metadata.json       # Document metadata (title, authors, etc.)
│   │   ├── chunks.json         # Text chunks for vector search
│   │   ├── scan_detection.json # OCR requirement analysis
│   │   └── summary.json        # Processing summary and provenance
│   └── E5F6G7H8/       # Another unique PDF
└── vector_db/           # Vector database files
    ├── A1B2C3D4_embedded.done  # Embedding completion markers
    └── lancedb/                 # Actual vector database
```

### Querying the corpus from an LLM (MCP)

After processing PDFs, expose the output as an MCP server that any MCP-capable client (Claude Desktop, Claude Code, Cursor, Continue, …) can drive. The server is a thin view over the per-paper artifacts; it does not store data of its own. See [PLAN.md §3](PLAN.md) for design context and [PLAN.md §8](PLAN.md) for the queries it's designed to support.

#### Quickstart: Claude Code against the demo output

A project-scoped [.mcp.json](.mcp.json) ships in the repo root and points at `demo_output/` with the local conda env's Python. From a Claude Code session started in this repo, run `/mcp` to list / enable / reload MCP servers. After approving `corpus`, the 15 tools are available in-session.

This works as long as:

- You've run `python process_corpus.py demo demo_output` (and ideally `python embed_chunks.py demo_output` for the `get_chunks_for_topic` semantic-search tool).
- Your conda env lives at `/opt/anaconda3/envs/corpus` — edit `.mcp.json` if yours is somewhere else.
- You've pointed at a different corpus than `demo_output/`, in which case change the second arg in `.mcp.json`.

Restart the MCP server (via `/mcp` → reload) after re-running `process_corpus.py` or `embed_chunks.py` so the index picks up new papers / new vectors.

#### Claude Desktop

Edit `~/Library/Application Support/Claude/claude_desktop_config.json` (absolute paths — Claude Desktop doesn't know about your repo CWD):

```json
{
  "mcpServers": {
    "corpus": {
      "command": "/opt/anaconda3/envs/corpus/bin/python",
      "args": [
        "/Users/you/repos/corpus/mcp_server.py",
        "/Users/you/repos/corpus/demo_output"
      ]
    }
  }
}
```

Restart Claude Desktop after editing.

#### Running the server by hand

Mostly useful for testing the server's stdio handshake. Clients normally launch it themselves through the config snippets above.

```bash
python mcp_server.py demo_output
# or, with an alternate WoRMS snapshot or embedding backend:
python mcp_server.py demo_output \
    --worms-sqlite /path/to/worms.sqlite \
    --embedding-backend local --embedding-model BAAI/bge-m3
```

#### Tool surface (14 tools)

| Tool | What it returns |
|---|---|
| `list_papers(year_from, year_to, limit)` | Every paper with its metadata + annotation counts |
| `get_paper(paper_hash)` | One paper's full metadata + top-10 taxa / anatomy rollups |
| `search_taxon(name)` | WoRMS lookup, resolving through synonymy |
| `get_papers_for_taxon(taxon_name)` | Papers mentioning a taxon, ordered by mention count |
| `get_chunks_for_taxon(taxon_name, paper_hash?, limit)` | Every chunk that mentions the taxon (exhaustive, synonymy-aware) |
| `get_chunks_for_topic(query, k, paper_hash?)` | Semantic search via the LanceDB vector index (needs `embed_chunks.py` to have been run) |
| `translate_chunk(paper_hash, chunk_id, target_language?)` | On-demand translation via Claude API; caches the result (needs `ANTHROPIC_API_KEY`) |
| `get_figures_for_taxon(taxon_name, paper_hash?)` | Figures from mentioning papers, ranked by caption overlap |
| `get_figures_for_anatomy(anatomy_term, paper_hash?)` | Figures whose captions mention an anatomy lexicon term |
| `get_chunks_by_section(paper_hash, section_class)` | Section-filtered chunks (introduction / description / references / …) |
| `get_bibliography(paper_hash)` | Parsed Grobid references |
| `get_papers_by_author(surname)` | Papers by a given author surname |
| `list_valid_species_under(parent_taxon)` | WoRMS-tree filter with per-species corpus coverage |
| `get_figure(paper_hash, figure_id)` | One figure's full record incl. image path |
| `get_chunk(paper_hash, chunk_id)` | One chunk's full record |

#### Example prompts

In a Claude Code session with the `corpus` server connected:

- *"What papers in the corpus mention Bargmannia? List them with year + first author."*
- *"List every currently-valid Apolemia species, and flag which ones already have papers in the corpus."*
- *"Give me the figures in Pugh 2001 that show nectophore morphology — include image paths so I can open them."*
- *"Summarize what the corpus says about the structure and function of pneumatophores, citing paper hashes."*
- *"Read the bibliography of Dunn 2005 and tell me which cited works are also in the corpus as their own papers."*

### Advanced Usage

**Resume interrupted processing** (skip already processed PDFs):
```bash
python process_corpus.py /path/to/input /path/to/output --resume
```

**Process embeddings separately** (useful for large corpora):
```bash
# First, run main processing
python process_corpus.py /path/to/input /path/to/output

# Then embed chunks. Default backend is local sentence-transformers
# (BGE-M3, 1024-dim) running on the best available accelerator —
# CUDA on Bouchet, MPS on Apple Silicon, CPU otherwise. The first
# run downloads the model (~2 GB) and caches it under
# ~/.cache/huggingface/.
python embed_chunks.py /path/to/output

# Pick a different model
python embed_chunks.py /path/to/output --model BAAI/bge-large-en-v1.5

# Force a device
python embed_chunks.py /path/to/output --device cpu

# Use OpenAI text-embedding-3-small (transitional, requires
# OPENAI_API_KEY in .env)
python embed_chunks.py /path/to/output --backend openai

# Embed one specific document by hash
python embed_chunks.py /path/to/output --pdf-hash a1b2c3d4e5f6

# Resume (skip docs already embedded with this backend+model)
python embed_chunks.py /path/to/output --resume

# Drop the existing table and re-embed (use when switching models)
python embed_chunks.py /path/to/output --rebuild
```

Once embedded, the MCP server's `get_chunks_for_topic` tool becomes
available for semantic search over the corpus.

### Example Workflows

**Small local collection:**
```bash
python process_corpus.py ./my_papers ./corpus_output
```

**Large research archive:**
```bash
# Process documents first
python process_corpus.py /Volumes/Research/PDFs ./large_corpus_output --resume

# Run embeddings overnight (requires OpenAI credits)
python embed_chunks.py ./large_corpus_output
```

**Distributed processing:**
```bash
# Process different subdirectories on different machines
python process_corpus.py /shared/pdfs/biology /shared/output/biology
python process_corpus.py /shared/pdfs/physics /shared/output/physics

# Then combine and embed
python embed_chunks.py /shared/output/biology
python embed_chunks.py /shared/output/physics
```

### Legacy Usage (Snakemake - Deprecated)

The old Snakemake workflow is still available for the demo directory:

```bash
snakemake --cores 1 --scheduler greedy
```

**Note**: The new `process_corpus.py` system is recommended as it provides better deduplication, provenance tracking, and handles arbitrary input directory structures.


### Quality Control

Tools to spot check the preprocessing steps:
1. **OCR Inspection**: Review processed PDFs for OCR quality
2. **Docling Visualization**: View parsing results with bounding boxes
3. **Metadata Validation**: Check extracted bibliographic information
4. **Chunk Review**: Assess text chunking quality

**Built-in Visualizations**: Each processed PDF automatically generates page visualization PNGs in the `visualizations/` directory showing:
- **Red boxes**: Word-level text detection
- **Yellow/Orange regions**: Detected figures and images

**Manual Visualization**: To create custom visualizations on any PDF:

```bash
python scripts/qc/visualize_docling_parse.py <path_to_pdf> --level <char|word|line>
```

### Analysis

The processed corpus supports:

- **Vector Search**: Semantic search across document chunks
- **Citation Analysis**: Discovery of relevant manuscripts
- **Figure Collection**: Anatomical figures of same species across papers

Analysis tools:
```bash
python scripts/analysis/search.py    # Vector search interface
python scripts/analysis/chat.py     # Chat with your corpus
```

## Output File Descriptions

Each processed document creates these files in `output_dir/documents/HASH/`:

- **`processed.pdf`**: OCR'd version of original (if OCR was needed)
- **`summary.json`**: Processing metadata and all input file paths
- **`scan_detection.json`**: Analysis of whether PDF needed OCR
- **`text.json`**: Full extracted text content in markdown format  
- **`figures.json`**: Metadata for all extracted figures
- **`figures/`**: Directory containing extracted figure images
- **`visualizations/`**: PNG files showing word boxes (red) and figure regions (yellow/orange)
- **`metadata.json`**: Bibliographic metadata (title, authors, DOI, etc.)
- **`chunks.json`**: Text split into searchable chunks with metadata

## Deduplication Example

If you have duplicate PDFs in your input directory:

```
input/
├── paper1.pdf
├── research/
│   ├── paper1.pdf          # Same content as above
│   └── different_paper.pdf
└── archive/
    └── paper1_copy.pdf     # Also same content
```

The system will:
1. Process `paper1.pdf` content only once (saves time and storage)
2. Store it in `output/documents/A1B2C3D4/` (example hash)
3. Record all three file locations in `summary.json`:
   ```json
   {
     "pdf_hash": "A1B2C3D4",
     "relative_paths": [
       "paper1.pdf",
       "research/paper1.pdf", 
       "archive/paper1_copy.pdf"
     ],
     "total_copies_found": 3
   }
   ```

## Resources

- [This tutorial](https://youtu.be/9lBTS5dM27c?si=d8cS4wbY6eGXaHH1), based on [this repo](https://github.com/daveebbelaar/ai-cookbook/tree/main/knowledge/docling)
- https://docling-project.github.io/docling/examples/export_figures/
- https://github.com/docling-project/docling-parse