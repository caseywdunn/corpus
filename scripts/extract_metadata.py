#!/usr/bin/env python3

import os
import json
import subprocess
import tempfile
from pathlib import Path

def extract_metadata_with_grobid(pdf_path):
    """Extract metadata using Grobid (placeholder - requires Grobid service)"""
    # This is a placeholder implementation
    # In practice, you'd need to have Grobid running as a service
    # and make HTTP requests to it
    
    filename = Path(pdf_path).stem
    
    # Basic metadata extraction from filename patterns
    metadata = {
        "title": "",
        "authors": [],
        "year": None,
        "journal": "",
        "doi": "",
        "filename": filename
    }
    
    # Try to extract year from filename
    import re
    year_match = re.search(r'(\d{4})', filename)
    if year_match:
        metadata["year"] = int(year_match.group(1))
    
    # Extract potential author names from filename
    # This is a simple heuristic - real Grobid would do much better
    parts = filename.replace('_', ' ').replace('-', ' ').split()
    potential_authors = []
    for part in parts:
        if part.istitle() and len(part) > 2 and not part.isdigit():
            potential_authors.append(part)
    
    if potential_authors:
        metadata["authors"] = potential_authors[:3]  # Take first 3 as authors
    
    return metadata

def main():
    input_pdf = snakemake.input[0]
    output_metadata = snakemake.output[0]
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_metadata), exist_ok=True)
    
    try:
        # Extract metadata
        metadata = extract_metadata_with_grobid(input_pdf)
        
        # Save metadata
        with open(output_metadata, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)
        
        print(f"Extracted metadata from {input_pdf}")
        
    except Exception as e:
        print(f"Error extracting metadata from {input_pdf}: {e}")
        # Create empty metadata file to not break the workflow
        empty_metadata = {
            "title": "",
            "authors": [],
            "year": None,
            "journal": "",
            "doi": "",
            "filename": Path(input_pdf).stem,
            "error": str(e)
        }
        
        with open(output_metadata, 'w', encoding='utf-8') as f:
            json.dump(empty_metadata, f, indent=2)

if __name__ == "__main__":
    main()