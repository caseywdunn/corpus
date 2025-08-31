#!/usr/bin/env python3

import json
import shutil
import subprocess
import sys
from pathlib import Path

def main():
    try:
        pdf_file = snakemake.input.pdf
        detection_file = snakemake.input.detection
        output_file = snakemake.output[0]
        
        print(f"Processing: {pdf_file}")
        print(f"Detection file: {detection_file}")
        print(f"Output file: {output_file}")
        
        # Create output directory
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        # Load detection results
        with open(detection_file, 'r') as f:
            detection = json.load(f)
        
        print(f"Detection result: {detection}")
        
        needs_ocr = detection.get('needs_ocr', False)
        
        if needs_ocr:
            # Run OCR on scanned PDF
            print(f"Running OCR on {Path(pdf_file).name} (detected as scanned)")
            cmd = [
                'ocrmypdf', 
                '--force-ocr',  # Force OCR even if some text exists
                '--optimize', '2',
                '--color-conversion-strategy', 'RGB',
                '--output-type', 'pdf',  # Skip PDF/A to avoid color space issues
                pdf_file, 
                output_file
            ]
            print(f"OCR command: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"OCR failed with return code {result.returncode}")
                print(f"STDOUT: {result.stdout}")
                print(f"STDERR: {result.stderr}")
                sys.exit(1)
        else:
            # Just copy born-digital PDF
            print(f"Copying {Path(pdf_file).name} (detected as born-digital)")
            shutil.copy2(pdf_file, output_file)
        
        print(f"Successfully processed: {Path(pdf_file).name} -> {Path(output_file).name}")
        
    except Exception as e:
        print(f"Error in prepare_pdf.py: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()