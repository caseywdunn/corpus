# Generalized Corpus Processing

The corpus processing system has been refactored to be more general and use hash-based file organization similar to Git.

## New Architecture

### Input/Output Structure

```
Input Directory (any structure):
├── paper1.pdf
├── subfolder/
│   ├── paper2.pdf
│   └── paper1.pdf  # duplicate (same content)
└── another_folder/
    └── paper3.pdf

Output Directory:
├── documents/
│   ├── A1B2C3D4/  # hash of paper1.pdf content
│   │   ├── processed.pdf
│   │   ├── figures/
│   │   │   ├── figure_1.png
│   │   │   └── figure_2.png
│   │   ├── text.json
│   │   ├── figures.json
│   │   ├── metadata.json
│   │   ├── chunks.json
│   │   ├── scan_detection.json
│   │   └── summary.json
│   ├── E5F6G7H8/  # hash of paper2.pdf
│   └── I9J0K1L2/  # hash of paper3.pdf
└── vector_db/
    ├── A1B2C3D4_embedded.done
    ├── E5F6G7H8_embedded.done
    └── lancedb/  # actual vector database
```

### Key Features

1. **Hash-based Deduplication**: PDFs with identical content (detected via SHA-256 hash) are only processed once, regardless of how many copies exist in the input directory.

2. **Preserves Source Information**: The `summary.json` file in each document directory tracks all input paths where the PDF was found.

3. **Git-like Organization**: Each unique PDF gets a hash-based directory (first 8 characters of SHA-256, like Git).

4. **Complete Processing Pipeline**: Each document directory contains all processed artifacts:
   - `processed.pdf`: OCR'd version if needed
   - `figures/`: Extracted figure images
   - `text.json`: Extracted text content
   - `figures.json`: Figure metadata and captions
   - `metadata.json`: Document metadata (title, authors, etc.)
   - `chunks.json`: Text chunks for vector search
   - `scan_detection.json`: Whether document needed OCR
   - `summary.json`: Input paths and processing metadata

## Usage

### Basic Processing

```bash
python process_corpus.py /path/to/input/pdfs /path/to/output
```

### Resume Processing (skip already processed)

```bash
python process_corpus.py /path/to/input/pdfs /path/to/output --resume
```

### Embed Chunks (can be run separately)

```bash
# Embed all documents
python embed_chunks.py /path/to/output

# Embed specific document by hash
python embed_chunks.py /path/to/output --pdf-hash A1B2C3D4

# Resume embedding (skip already embedded)
python embed_chunks.py /path/to/output --resume
```

## Example Output Structure

For a PDF with hash `A1B2C3D4`, the `summary.json` would look like:

```json
{
  "pdf_hash": "A1B2C3D4",
  "input_dir": "/path/to/input",
  "relative_paths": [
    "paper1.pdf",
    "subfolder/paper1.pdf"
  ],
  "total_copies_found": 2,
  "processing_summary": {
    "status": "success",
    "processing_steps": [
      "scan_detection",
      "pdf_preparation", 
      "docling_extraction",
      "metadata_extraction",
      "text_chunking"
    ],
    "files_created": [...]
  },
  "output_directory": "/path/to/output/documents/A1B2C3D4"
}
```

## Benefits

1. **Efficient**: Duplicate PDFs are only processed once
2. **Traceable**: Full provenance of where each PDF came from
3. **Resumable**: Can resume interrupted processing
4. **Modular**: Embedding can be run separately from main processing
5. **Scalable**: Works with arbitrarily large and nested input directories
6. **Git-like**: Familiar hash-based organization for developers

## Migration from Old System

The old system used the `demo/` directory with Snakemake. The new system:

- Can process any input directory structure
- Uses direct Python scripts instead of Snakemake
- Provides better deduplication and organization
- Maintains compatibility with the same underlying processing tools

To migrate, simply run the new `process_corpus.py` script on your existing PDF collection.