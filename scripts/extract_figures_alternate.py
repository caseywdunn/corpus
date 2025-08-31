#!/usr/bin/env python3

import fitz  # PyMuPDF
import os
from pathlib import Path
import json

def extract_images_with_pymupdf(pdf_path, output_dir):
    """Extract images directly from PDF using PyMuPDF"""
    
    figures_data = []
    pdf_name = Path(pdf_path).stem
    
    try:
        doc = fitz.open(pdf_path)
        
        for page_num in range(len(doc)):
            page = doc.load_page(page_num)
            image_list = page.get_images()
            
            print(f"Page {page_num + 1}: Found {len(image_list)} images")
            
            for img_index, img in enumerate(image_list):
                # Get image data
                xref = img[0]  # xref number
                pix = fitz.Pixmap(doc, xref)
                
                if pix.n - pix.alpha < 4:  # Skip if not RGB/GRAY
                    figure_filename = f"{pdf_name}_page{page_num+1}_img{img_index+1}.png"
                    figure_path = Path(output_dir) / figure_filename
                    
                    pix.save(str(figure_path))
                    
                    # Get image dimensions and position
                    img_rect = page.get_image_rects(img)[0] if page.get_image_rects(img) else None
                    
                    figure_info = {
                        "figure_id": f"pymupdf_p{page_num+1}_i{img_index+1}",
                        "filename": figure_filename,
                        "file_path": str(figure_path),
                        "page": page_num + 1,
                        "width": pix.width,
                        "height": pix.height,
                        "bbox": str(img_rect) if img_rect else None,
                        "extraction_method": "pymupdf"
                    }
                    figures_data.append(figure_info)
                    
                pix = None  # Free memory
        
        doc.close()
        
    except Exception as e:
        print(f"Error extracting images with PyMuPDF: {e}")
    
    return figures_data

def main():
    # This can be called from the main extract_docling.py or standalone
    pdf_path = "demo/Dunn-etal2005_Marrus.pdf"  # Example path
    output_dir = "output/test_figures"
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    figures = extract_images_with_pymupdf(pdf_path, output_dir)
    
    print(f"Extracted {len(figures)} figures")
    for fig in figures:
        print(f"  - {fig['filename']}: {fig['width']}x{fig['height']}")

if __name__ == "__main__":
    main()