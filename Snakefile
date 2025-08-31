
import os
from pathlib import Path

# Configuration
DEMO_DIR = "demo"
OUTPUT_DIR = "output"
PDFS = glob_wildcards(f"{DEMO_DIR}/{{pdf}}.pdf").pdf

rule all:
    input:
        expand(f"{OUTPUT_DIR}/vector_db/{{pdf}}_embedded.done", pdf=PDFS)

rule detect_scan:
    input:
        f"{DEMO_DIR}/{{pdf}}.pdf"
    output:
        f"{OUTPUT_DIR}/scan_detection/{{pdf}}_detection.json"
    script:
        "scripts/detect_scan.py"

rule prepare_pdf:
    input:
        pdf=f"{DEMO_DIR}/{{pdf}}.pdf",
        detection=f"{OUTPUT_DIR}/scan_detection/{{pdf}}_detection.json"
    output:
        f"{OUTPUT_DIR}/processed/{{pdf}}_processed.pdf"
    script:
        "scripts/prepare_pdf.py"

rule extract_docling:
    input:
        f"{OUTPUT_DIR}/processed/{{pdf}}_processed.pdf"
    output:
        text=f"{OUTPUT_DIR}/docling/{{pdf}}_text.json",
        figures=f"{OUTPUT_DIR}/docling/{{pdf}}_figures.json"
    script:
        "scripts/extract_docling.py"

rule extract_metadata:
    input:
        f"{OUTPUT_DIR}/processed/{{pdf}}_processed.pdf"
    output:
        f"{OUTPUT_DIR}/metadata/{{pdf}}_metadata.json"
    script:
        "scripts/extract_metadata.py"

rule chunk_text:
    input:
        text=f"{OUTPUT_DIR}/docling/{{pdf}}_text.json",
        metadata=f"{OUTPUT_DIR}/metadata/{{pdf}}_metadata.json"
    output:
        f"{OUTPUT_DIR}/chunks/{{pdf}}_chunks.json"
    script:
        "scripts/chunk_text.py"

rule embed_and_store:
    input:
        f"{OUTPUT_DIR}/chunks/{{pdf}}_chunks.json"
    output:
        f"{OUTPUT_DIR}/vector_db/{{pdf}}_embedded.done"
    script:
        "scripts/embed_and_store.py"
