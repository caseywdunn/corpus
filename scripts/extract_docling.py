#!/usr/bin/env python3

import os
import json
from pathlib import Path
from docling.document_converter import DocumentConverter, PdfFormatOption
from docling.datamodel.base_models import InputFormat
from docling.datamodel.pipeline_options import PdfPipelineOptions

def main():
    input_pdf = snakemake.input[0]
    text_output = snakemake.output.text
    figures_output = snakemake.output.figures
    
    # Ensure output directories exist
    os.makedirs(os.path.dirname(text_output), exist_ok=True)
    os.makedirs(os.path.dirname(figures_output), exist_ok=True)
    
    # Create figures output directory
    pdf_name = Path(input_pdf).stem
    figures_dir = Path(figures_output).parent / f"{pdf_name}_figures"
    figures_dir.mkdir(exist_ok=True)
    
    # Check if this is a scanned document based on detection results  
    pdf_detection_file = str(input_pdf).replace('processed/', 'scan_detection/').replace('_processed.pdf', '_detection.json')
    is_scanned = False
    
    try:
        if Path(pdf_detection_file).exists():
            with open(pdf_detection_file, 'r') as f:
                detection = json.load(f)
                is_scanned = detection.get('needs_ocr', False)
    except Exception as e:
        print(f"Could not read detection file: {e}")
    
    print(f"Document detected as: {'scanned' if is_scanned else 'born-digital'}")
    
    # Configure pipeline based on document type
    if is_scanned:
        # For scanned documents: enable OCR and proper figure segmentation
        pipeline_options = PdfPipelineOptions(
            do_ocr=True,  # Enable OCR for scanned content
            do_table_structure=True,
            generate_picture_images=True,  # Extract figure images instead of full pages
            generate_page_images=False,   # Don't extract full page as images
            do_picture_classification=True,  # Classify image content
        )
        pipeline_options.table_structure_options.do_cell_matching = True
        
        # Use PdfFormatOption to properly configure the format options
        pdf_format_option = PdfFormatOption(pipeline_options=pipeline_options)
        converter = DocumentConverter(
            format_options={InputFormat.PDF: pdf_format_option}
        )
        print("Configured docling for scanned document with OCR and figure segmentation")
    else:
        # For born-digital documents: standard processing
        converter = DocumentConverter()
        print("Using standard docling processing for born-digital document")
    
    # Convert the PDF
    result = converter.convert(input_pdf)
    document = result.document
    
    # Debug: Print document structure
    print(f"Document attributes: {dir(document)}")
    if hasattr(document, 'body'):
        print(f"Body attributes: {dir(document.body)}")
        if hasattr(document.body, 'elements'):
            print(f"Found {len(document.body.elements)} elements")
            for i, elem in enumerate(document.body.elements[:5]):  # Show first 5
                print(f"Element {i}: type={getattr(elem, 'type', 'unknown')}, text preview={str(elem)[:50]}...")
    
    # Extract text content
    text_content = {
        "title": document.name,
        "text": document.export_to_markdown(),
        "pages": len(document.pages) if hasattr(document, 'pages') else None
    }
    
    # Save text content
    with open(text_output, 'w', encoding='utf-8') as f:
        json.dump(text_content, f, indent=2, ensure_ascii=False)
    
    # Extract and save figures
    figures_data = []
    try:
        # Method 0: Check if document has pictures directly
        if hasattr(document, 'pictures') and len(document.pictures) > 0:
            print(f"Found {len(document.pictures)} pictures in document.pictures")
            for pic_idx, picture in enumerate(document.pictures):
                print(f"Picture {pic_idx}: {type(picture).__name__}")
                
                figure_filename = f"{pdf_name}_docling_{len(figures_data)+1}.png"
                figure_path = figures_dir / figure_filename
                
                figure_info = {
                    "figure_id": f"docling_{len(figures_data)+1}",
                    "filename": figure_filename,
                    "file_path": str(figure_path),
                    "type": type(picture).__name__,
                    "caption": getattr(picture, 'text', ''),
                    "bbox": str(getattr(picture, 'bbox', None)),
                    "extraction_method": "docling_pictures"
                }
                figures_data.append(figure_info)
                
                # Try to save the picture
                try:
                    if hasattr(picture, 'get_image'):
                        image = picture.get_image(document) if callable(picture.get_image) else picture.get_image
                        if image:
                            image.save(str(figure_path))
                            print(f"Saved docling picture: {figure_filename}")
                    elif hasattr(picture, 'image'):
                        if picture.image:
                            picture.image.save(str(figure_path))
                            print(f"Saved docling picture: {figure_filename}")
                except Exception as e:
                    print(f"Could not save docling picture {pic_idx}: {e}")
        
        # Method 1: Look through document elements for images/figures
        elif hasattr(document, 'body') and hasattr(document.body, 'children'):
            print(f"Checking {len(document.body.children)} document body children for figures")
            
            for elem_idx, element in enumerate(document.body.children):
                element_type = str(type(element).__name__)
                if elem_idx < 10:  # Only show first 10 to avoid spam
                    print(f"Element {elem_idx}: {element_type} - {str(element)[:100]}")
                
                # Check for picture/image elements (broader check)
                if ('picture' in element_type.lower() or 'image' in element_type.lower() or
                    'Picture' in element_type or 'Figure' in element_type):
                    figure_filename = f"{pdf_name}_docling_{len(figures_data)+1}.png"
                    figure_path = figures_dir / figure_filename
                    
                    figure_info = {
                        "figure_id": f"docling_{len(figures_data)+1}",
                        "filename": figure_filename,
                        "file_path": str(figure_path),
                        "type": element_type,
                        "caption": getattr(element, 'text', ''),
                        "bbox": str(getattr(element, 'bbox', None)),
                        "extraction_method": "docling"
                    }
                    figures_data.append(figure_info)
                    
                    # Try to extract and save the image
                    try:
                        if hasattr(element, 'get_image'):
                            image = element.get_image(document) if callable(element.get_image) else element.get_image
                            if image:
                                image.save(str(figure_path))
                                print(f"Saved docling image: {figure_filename}")
                    except Exception as e:
                        print(f"Could not save docling image {elem_idx}: {e}")
        else:
            print("No document.pictures or document.body.children found")
        
        # Method 2: Use PyMuPDF as fallback to extract images directly from PDF
        print("Trying PyMuPDF extraction as fallback...")
        try:
            import fitz
            doc = fitz.open(input_pdf)
            
            for page_num in range(len(doc)):
                page = doc.load_page(page_num)
                image_list = page.get_images()
                
                if image_list:
                    print(f"Page {page_num + 1}: Found {len(image_list)} images")
                
                for img_index, img in enumerate(image_list):
                    try:
                        # Get image data
                        xref = img[0]  # xref number
                        pix = fitz.Pixmap(doc, xref)
                        
                        if pix.n - pix.alpha < 4:  # RGB or GRAY
                            figure_filename = f"{pdf_name}_page{page_num+1}_img{img_index+1}.png"
                            figure_path = figures_dir / figure_filename
                            
                            pix.save(str(figure_path))
                            
                            # Get image position if available
                            try:
                                img_rects = page.get_image_rects(img)
                                img_rect = img_rects[0] if img_rects else None
                            except:
                                img_rect = None
                            
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
                            print(f"Saved PyMuPDF image: {figure_filename} ({pix.width}x{pix.height})")
                        
                        pix = None  # Free memory
                        
                    except Exception as e:
                        print(f"Error saving image {img_index} from page {page_num + 1}: {e}")
            
            doc.close()
            
        except ImportError:
            print("PyMuPDF not available for fallback extraction")
        except Exception as e:
            print(f"Error with PyMuPDF extraction: {e}")
            
    except Exception as e:
        print(f"Warning: Could not extract figures: {e}")
        import traceback
        traceback.print_exc()
        figures_data = []
    
    print(f"Extracted and saved {len(figures_data)} figures to {figures_dir}")
    
    # Add figures directory info to the JSON
    figures_info = {
        "figures": figures_data,
        "figures_directory": str(figures_dir),
        "total_figures": len(figures_data)
    }
    
    # Save figures data
    with open(figures_output, 'w', encoding='utf-8') as f:
        json.dump(figures_info, f, indent=2)
    
    print(f"Extracted text and figures from {input_pdf}")

if __name__ == "__main__":
    main()