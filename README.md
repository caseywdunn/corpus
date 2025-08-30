# corpus

A workflow for interrogating a corpus of scientific literature pdfs, including old (scanned, potentially spanning centuries of formatting conventions) and new (born-digital) papers. My own goal is analysis of the siphonophore literature, such as:

- Searching
- Citation analysis (including discovery of relevant manuscripts not in the corpus)
- Collecting anatomical figures of the same species across different papers

## Overview

### Preprocessing

A snakemake worklflow that:

1. Does OCR if needed with `ocrmypdf`.
2. Converts pdf to docling format, extracting text and figures.
3. Extracts metadata (title, authors, year, journal, doi) with `Grobid`.
4. Chunks text per [this](https://github.com/daveebbelaar/ai-cookbook/blob/main/knowledge/docling/2-chunking.py)
5. Tokenize and vector embed text, and put it into a vector database. Following [this](https://github.com/daveebbelaar/ai-cookbook/blob/main/knowledge/docling/3-embedding.py)

### Quality Control

Tools to spot check the preprocessing steps:
1. Inspection of OCR output.
2. Visualization of docling parse with `docling_parse/visualize.py`
3. Validation of extracted metadata.
4. Review of chunking and tokenization quality.

This relies on a set of demo documents that have been manually assessed.

### Aanalysis

Includes:

- search, per [this](https://github.com/daveebbelaar/ai-cookbook/blob/main/knowledge/docling/4-search.py)
- chat, per [this](https://github.com/daveebbelaar/ai-cookbook/blob/main/knowledge/docling/5-chat.py)

And more, to be explored later.


## Installation

```
conda create --name corpus python=3.13
pip install -r requirements.txt
```

## Usage

```
snakemake --cores 1

```


## Resources

- [This tutorial](https://youtu.be/9lBTS5dM27c?si=d8cS4wbY6eGXaHH1), based on [this repo](https://github.com/daveebbelaar/ai-cookbook/tree/main/knowledge/docling)
- https://docling-project.github.io/docling/examples/export_figures/
- https://github.com/docling-project/docling-parse