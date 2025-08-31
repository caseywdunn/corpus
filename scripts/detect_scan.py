#!/usr/bin/env python3

import sys
import subprocess
import json
from pathlib import Path

def has_meaningful_text(pdf_path):
    """Check if PDF has meaningful extractable text content"""
    try:
        import fitz  # PyMuPDF
        doc = fitz.open(pdf_path)
        total_text = ""
        total_chars = 0
        
        # Sample first few pages to determine text quality
        pages_to_check = min(3, len(doc))
        
        for page_num in range(pages_to_check):
            page = doc[page_num]
            page_text = page.get_text()
            total_text += page_text
            total_chars += len(page_text.strip())
        
        doc.close()
        
        # Remove whitespace and count meaningful characters
        meaningful_text = ''.join(total_text.split())
        
        # Calculate text density (characters per page)
        avg_chars_per_page = total_chars / pages_to_check if pages_to_check > 0 else 0
        
        # Consider it meaningful if:
        # 1. Has substantial text (>500 chars total meaningful text)
        # 2. Has reasonable density (>100 chars per page on average)
        has_substantial_text = len(meaningful_text) > 500
        has_good_density = avg_chars_per_page > 100
        
        return has_substantial_text and has_good_density
        
    except Exception as e:
        print(f"Error checking text content: {e}")
        # If we can't determine, assume it needs OCR to be safe
        return False

def main():
    # Get inputs/outputs from Snakemake
    pdf_path = snakemake.input[0]
    output_json = snakemake.output[0]
    
    # Create output directory
    Path(output_json).parent.mkdir(parents=True, exist_ok=True)
    
    # Check if file has meaningful text
    has_text = has_meaningful_text(pdf_path)
    
    result = {
        "filename": Path(pdf_path).name,
        "has_text": has_text,
        "needs_ocr": not has_text,
        "file_type": "born_digital" if has_text else "scanned"
    }
    
    # Save result
    with open(output_json, 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"{Path(pdf_path).name}: {result['file_type']}")

if __name__ == "__main__":
    main()