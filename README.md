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

```
brew install ghostscript
brew install tesseract
brew install pngquant
brew install jbig2enc
conda create --name corpus python=3.13
conda activate corpus
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

### Advanced Usage

**Resume interrupted processing** (skip already processed PDFs):
```bash
python process_corpus.py /path/to/input /path/to/output --resume
```

**Process embeddings separately** (useful for large corpora):
```bash
# First, run main processing without embeddings taking too long
python process_corpus.py /path/to/input /path/to/output

# Then run embeddings separately (requires OpenAI API key)
python embed_chunks.py /path/to/output

# Or embed specific documents by hash
python embed_chunks.py /path/to/output --pdf-hash A1B2C3D4

# Resume embedding (skip already embedded documents)
python embed_chunks.py /path/to/output --resume
```

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