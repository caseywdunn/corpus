#!/usr/bin/env python3
"""
Generalized corpus processing script that takes input_dir and output_dir arguments.
Recursively finds PDFs in input_dir, processes them with hash-based organization.
"""

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Set
import tempfile
import os


def create_cell_visualizations(input_pdf: Path, output_dir: Path, pdf_name: str, document):
    """Create cell visualization PNGs using docling-parse with figure regions highlighted"""
    try:
        from docling_core.types.doc.page import TextCellUnit
        from docling_parse.pdf_parser import DoclingPdfParser, PdfDocument
        from PIL import ImageDraw
        
        print(f"    Creating cell visualizations...")
        
        pdf_parser = DoclingPdfParser()
        pdf_doc: PdfDocument = pdf_parser.load(path_or_stream=str(input_pdf))
        
        # Get figure bounding boxes from the docling document
        figure_bboxes_by_page = {}
        if hasattr(document, 'pictures') and document.pictures:
            for picture in document.pictures:
                if hasattr(picture, 'prov') and picture.prov:
                    for prov_item in picture.prov:
                        if hasattr(prov_item, 'page_no') and hasattr(prov_item, 'bbox'):
                            page_no = prov_item.page_no
                            if page_no not in figure_bboxes_by_page:
                                figure_bboxes_by_page[page_no] = []
                            figure_bboxes_by_page[page_no].append(prov_item.bbox)
        
        # Create visualization for word-level cells
        cell_unit = TextCellUnit.WORD
        
        for page_no, pred_page in pdf_doc.iterate_pages():
            img = pred_page.render_as_image(cell_unit=cell_unit)
            draw = ImageDraw.Draw(img)
            img_height = img.height
            
            # First, highlight figure regions in yellow
            if page_no in figure_bboxes_by_page:
                for bbox in figure_bboxes_by_page[page_no]:
                    # Convert bbox coordinates (assumes bottom-left origin) to PIL coordinates (top-left origin)
                    x0 = bbox.l
                    x1 = bbox.r
                    y0_raw = bbox.t  # top in document coordinates
                    y1_raw = bbox.b  # bottom in document coordinates
                    # Flip Y coordinates for PIL (top-left origin)
                    y0 = img_height - y0_raw
                    y1 = img_height - y1_raw
                    
                    # Ensure coordinates are in correct order
                    x0, x1 = min(x0, x1), max(x0, x1)
                    y0, y1 = min(y0, y1), max(y0, y1)
                    
                    rect = (x0, y0, x1, y1)
                    draw.rectangle(rect, fill="yellow", outline="orange", width=3)
            
            # Then, draw red rectangles around text cells (these will show over non-figure areas)
            for cell in pred_page.iterate_cells(unit_type=cell_unit):
                # Convert cell.rect to (x0, y0, x1, y1) and flip y for PIL
                x0 = min(getattr(cell.rect, "r_x0", 0), getattr(cell.rect, "r_x2", 0))
                x1 = max(getattr(cell.rect, "r_x0", 0), getattr(cell.rect, "r_x2", 0))
                y0_raw = min(getattr(cell.rect, "r_y0", 0), getattr(cell.rect, "r_y2", 0))
                y1_raw = max(getattr(cell.rect, "r_y0", 0), getattr(cell.rect, "r_y2", 0))
                y0 = img_height - y1_raw
                y1 = img_height - y0_raw
                rect = (x0, y0, x1, y1)
                
                # Check if this cell overlaps with any figure region
                is_in_figure = False
                if page_no in figure_bboxes_by_page:
                    for bbox in figure_bboxes_by_page[page_no]:
                        fig_x0, fig_x1 = bbox.l, bbox.r
                        fig_y0_raw, fig_y1_raw = bbox.t, bbox.b
                        fig_y0 = img_height - fig_y0_raw
                        fig_y1 = img_height - fig_y1_raw
                        fig_x0, fig_x1 = min(fig_x0, fig_x1), max(fig_x0, fig_x1)
                        fig_y0, fig_y1 = min(fig_y0, fig_y1), max(fig_y0, fig_y1)
                        
                        # Check for overlap
                        if (x0 < fig_x1 and x1 > fig_x0 and y0 < fig_y1 and y1 > fig_y0):
                            is_in_figure = True
                            break
                
                # Only draw cell boundaries for text regions (not in figures)
                if not is_in_figure:
                    draw.rectangle(rect, outline="red", width=2)
            
            # Save the visualization
            viz_filename = f"page_{page_no}_visualization.png"
            viz_path = output_dir / viz_filename
            img.save(str(viz_path))
            
        print(f"    Saved {len(list(pdf_doc.iterate_pages()))} visualization PNGs")
            
    except ImportError as e:
        print(f"    Warning: Could not create visualizations - missing dependencies: {e}")
    except Exception as e:
        print(f"    Warning: Could not create visualizations: {e}")


def calculate_pdf_hash(pdf_path: Path) -> str:
    """Calculate SHA-256 hash of PDF file for unique identification."""
    hasher = hashlib.sha256()
    with open(pdf_path, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hasher.update(chunk)
    return hasher.hexdigest()[:8].upper()  # Use first 8 chars like git


def find_all_pdfs(input_dir: Path) -> Dict[str, List[Path]]:
    """
    Recursively find all PDFs in input directory.
    Returns dict mapping PDF hash to list of paths (for duplicate detection).
    """
    pdf_map = {}
    
    for pdf_path in input_dir.rglob("*.pdf"):
        if pdf_path.is_file():
            try:
                pdf_hash = calculate_pdf_hash(pdf_path)
                if pdf_hash not in pdf_map:
                    pdf_map[pdf_hash] = []
                pdf_map[pdf_hash].append(pdf_path)
            except Exception as e:
                print(f"Warning: Could not hash {pdf_path}: {e}")
    
    return pdf_map


def create_output_structure(output_dir: Path):
    """Create the output directory structure."""
    documents_dir = output_dir / "documents"
    vector_db_dir = output_dir / "vector_db"
    
    documents_dir.mkdir(parents=True, exist_ok=True)
    vector_db_dir.mkdir(parents=True, exist_ok=True)
    
    return documents_dir, vector_db_dir


def get_relative_paths(pdf_paths: List[Path], input_dir: Path) -> List[str]:
    """Get relative paths of PDFs from input directory."""
    return [str(path.relative_to(input_dir)) for path in pdf_paths]


def run_pdf_processing_pipeline(pdf_path: Path, hash_dir: Path, temp_dir: Path) -> Dict:
    """
    Run the PDF processing pipeline for a single PDF.
    Returns summary information about the processing.
    """
    pdf_name = pdf_path.stem
    temp_pdf = temp_dir / f"{pdf_name}.pdf"
    
    # Copy PDF to temp directory for processing
    shutil.copy2(pdf_path, temp_pdf)
    
    # Create subdirectories in hash directory
    figures_dir = hash_dir / "figures"
    figures_dir.mkdir(exist_ok=True)
    visualizations_dir = hash_dir / "visualizations"
    visualizations_dir.mkdir(exist_ok=True)
    
    processing_summary = {
        "original_pdf": str(pdf_path),
        "processing_steps": [],
        "files_created": [],
        "errors": []
    }
    
    try:
        # Step 1: Detect scan
        print(f"  Detecting scan type...")
        detection_result = detect_scan_type(temp_pdf)
        detection_file = hash_dir / "scan_detection.json"
        with open(detection_file, 'w') as f:
            json.dump(detection_result, f, indent=2)
        processing_summary["files_created"].append(str(detection_file))
        processing_summary["processing_steps"].append("scan_detection")
        
        # Step 2: Prepare PDF (OCR if needed)
        print(f"  Preparing PDF...")
        processed_pdf = hash_dir / f"processed.pdf"
        prepare_pdf_result = prepare_pdf(temp_pdf, detection_result, processed_pdf)
        processing_summary["files_created"].append(str(processed_pdf))
        processing_summary["processing_steps"].append("pdf_preparation")
        
        # Step 3: Extract text and figures with docling
        print(f"  Extracting text and figures...")
        text_file = hash_dir / "text.json"
        figures_file = hash_dir / "figures.json"
        extract_docling_content(processed_pdf, text_file, figures_file, figures_dir, visualizations_dir)
        processing_summary["files_created"].extend([str(text_file), str(figures_file)])
        processing_summary["processing_steps"].append("docling_extraction")
        
        # Step 4: Extract metadata
        print(f"  Extracting metadata...")
        metadata_file = hash_dir / "metadata.json"
        extract_metadata(processed_pdf, metadata_file)
        processing_summary["files_created"].append(str(metadata_file))
        processing_summary["processing_steps"].append("metadata_extraction")
        
        # Step 5: Chunk text
        print(f"  Chunking text...")
        chunks_file = hash_dir / "chunks.json"
        chunk_text(text_file, metadata_file, chunks_file)
        processing_summary["files_created"].append(str(chunks_file))
        processing_summary["processing_steps"].append("text_chunking")
        
        processing_summary["status"] = "success"
        
    except Exception as e:
        processing_summary["status"] = "error"
        processing_summary["errors"].append(str(e))
        print(f"  Error processing PDF: {e}")
    
    return processing_summary


def detect_scan_type(pdf_path: Path) -> Dict:
    """Detect if PDF is scanned or born-digital."""
    try:
        import fitz  # PyMuPDF
        doc = fitz.open(pdf_path)
        total_text = ""
        total_chars = 0
        
        pages_to_check = min(3, len(doc))
        
        for page_num in range(pages_to_check):
            page = doc[page_num]
            page_text = page.get_text()
            total_text += page_text
            total_chars += len(page_text.strip())
        
        doc.close()
        
        meaningful_text = ''.join(total_text.split())
        avg_chars_per_page = total_chars / pages_to_check if pages_to_check > 0 else 0
        
        has_substantial_text = len(meaningful_text) > 500
        has_good_density = avg_chars_per_page > 100
        has_text = has_substantial_text and has_good_density
        
        return {
            "filename": pdf_path.name,
            "has_text": has_text,
            "needs_ocr": not has_text,
            "file_type": "born_digital" if has_text else "scanned"
        }
        
    except ImportError:
        # If PyMuPDF isn't available, avoid expensive OCR; proceed as born-digital.
        print(f"    Warning: PyMuPDF not available; treating as born-digital (no OCR)")
        return {
            "filename": pdf_path.name,
            "has_text": True,
            "needs_ocr": False,
            "file_type": "born_digital"
        }
    except Exception as e:
        print(f"    Warning: Error checking text content: {e}")
        return {
            "filename": pdf_path.name,
            "has_text": False,
            "needs_ocr": True,
            "file_type": "scanned"
        }


def prepare_pdf(input_pdf: Path, detection_result: Dict, output_pdf: Path):
    """Prepare PDF by running OCR if needed."""
    needs_ocr = detection_result.get('needs_ocr', False)
    
    if needs_ocr:
        # Check if ocrmypdf is available
        if shutil.which('ocrmypdf') is None:
            print(f"    Warning: ocrmypdf not found, copying original PDF")
            shutil.copy2(input_pdf, output_pdf)
            return
            
        print(f"    Running OCR on {input_pdf.name} (detected as scanned)")
        cmd = [
            'ocrmypdf', 
            '--force-ocr',
            '--optimize', '2',
            '--color-conversion-strategy', 'RGB',
            '--output-type', 'pdf',
            str(input_pdf), 
            str(output_pdf)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"    Warning: OCR failed, copying original PDF")
            print(f"    OCR Error: {result.stderr}")
            shutil.copy2(input_pdf, output_pdf)
        else:
            print(f"    OCR completed successfully")
    else:
        print(f"    Copying {input_pdf.name} (detected as born-digital)")
        shutil.copy2(input_pdf, output_pdf)


def extract_docling_content(pdf_path: Path, text_output: Path, figures_output: Path, figures_dir: Path, visualizations_dir: Path):
    """Extract text and figures using docling, with PyMuPDF fallback for figures.

    Guarantees that only successfully saved figures are listed in figures.json.
    """
    document = None
    text_content = None

    # Ensure figure directory exists
    figures_dir.mkdir(parents=True, exist_ok=True)

    try:
        from docling.document_converter import DocumentConverter, PdfFormatOption
        from docling.datamodel.base_models import InputFormat
        from docling.datamodel.pipeline_options import PdfPipelineOptions
        # Configure converter to generate picture images explicitly
        pipeline_options = PdfPipelineOptions(
            do_ocr=False,
            do_table_structure=True,
            generate_picture_images=True,
            generate_page_images=False,
            do_picture_classification=True,
        )
        pdf_format_option = PdfFormatOption(pipeline_options=pipeline_options)
        converter = DocumentConverter(format_options={InputFormat.PDF: pdf_format_option})
        result = converter.convert(str(pdf_path))
        document = result.document

        # Extract text from docling if available
        text_content = {
            "title": document.name,
            "text": document.export_to_markdown(),
            "pages": len(document.pages) if hasattr(document, "pages") else None,
        }
    except ImportError as e:
        print(f"    Warning: Docling not available ({e}), will fallback for figures")
    except Exception as e:
        print(f"    Warning: Docling extraction failed ({e}), will fallback for figures")

    # If we couldn't get text via docling, write a placeholder text file
    if text_content is None:
        text_content = {
            "title": pdf_path.stem,
            "text": f"# {pdf_path.stem}\n\n[Text extraction unavailable — see logs]",
            "pages": None,
        }

    # Save text content (from docling or placeholder)
    with open(text_output, "w", encoding="utf-8") as f:
        json.dump(text_content, f, indent=2, ensure_ascii=False)

    # Accumulate figure metadata only when file is saved
    figures_data = []

    # Helper to append a figure entry only if file exists
    def append_figure(figure_id: str, figure_path: Path, caption: str, meta: dict):
        if figure_path.exists() and figure_path.stat().st_size > 0:
            entry = {
                "figure_id": figure_id,
                "filename": figure_path.name,
                "file_path": str(figure_path),
                "caption": caption or "",
            }
            entry.update(meta or {})
            figures_data.append(entry)
        else:
            print(f"    Warning: Skipping missing/empty figure file: {figure_path}")

    # First try: use docling pictures if document is available
    if document is not None and hasattr(document, "pictures"):
        saved_count = 0
        for idx, picture in enumerate(document.pictures or []):
            figure_path = figures_dir / f"figure_{idx + 1}.png"

            # Caption (best-effort)
            caption = ""
            try:
                if hasattr(picture, "caption_text"):
                    caption = picture.caption_text(document) if callable(picture.caption_text) else (
                        str(picture.caption_text) if picture.caption_text else ""
                    )
            except Exception as e:
                print(f"    Note: caption extraction failed for figure {idx+1}: {e}")

            # Save image if available
            try:
                image = None
                if hasattr(picture, "get_image"):
                    image = picture.get_image(document) if callable(picture.get_image) else picture.get_image
                if image is not None and hasattr(image, "save"):
                    image.save(str(figure_path))
                    saved_count += 1
                    append_figure(
                        figure_id=f"docling_{idx + 1}",
                        figure_path=figure_path,
                        caption=caption,
                        meta={"extraction_method": "docling"},
                    )
                else:
                    print(f"    Warning: Docling did not return a savable image for figure {idx+1}")
            except Exception as e:
                print(f"    Warning: Could not save docling figure {idx+1}: {e}")

        if saved_count == 0:
            print("    No figures saved via docling; will try PyMuPDF fallback")

    # Fallback: use PyMuPDF to extract embedded images from the PDF
    if len(figures_data) == 0:
        try:
            import fitz  # PyMuPDF
            doc = fitz.open(str(pdf_path))
            for page_num in range(len(doc)):
                page = doc.load_page(page_num)
                for img_index, img in enumerate(page.get_images()):
                    try:
                        xref = img[0]
                        pix = fitz.Pixmap(doc, xref)
                        if pix.n - pix.alpha < 4:  # RGB or GRAY
                            figure_path = figures_dir / f"page{page_num + 1}_img{img_index + 1}.png"
                            pix.save(str(figure_path))
                            append_figure(
                                figure_id=f"pymupdf_p{page_num + 1}_i{img_index + 1}",
                                figure_path=figure_path,
                                caption="",
                                meta={
                                    "extraction_method": "pymupdf",
                                    "page": page_num + 1,
                                    "width": pix.width,
                                    "height": pix.height,
                                },
                            )
                        pix = None
                    except Exception as e:
                        print(f"    Warning: Failed to save PyMuPDF image p{page_num+1} i{img_index+1}: {e}")
            doc.close()
        except ImportError:
            print("    Warning: PyMuPDF not available; cannot run fallback image extraction")
        except Exception as e:
            print(f"    Warning: PyMuPDF fallback failed: {e}")

    # Write figures.json
    figures_info = {
        "figures": figures_data,
        "figures_directory": str(figures_dir),
        "total_figures": len(figures_data),
    }
    with open(figures_output, "w", encoding="utf-8") as f:
        json.dump(figures_info, f, indent=2)

    # Create cell visualizations only if we have a document from docling
    if document is not None:
        pdf_name = pdf_path.stem
        try:
            create_cell_visualizations(pdf_path, visualizations_dir, pdf_name, document)
        except Exception as e:
            print(f"    Warning: creating visualizations failed: {e}")


def extract_metadata(pdf_path: Path, metadata_output: Path):
    """Extract metadata using Grobid."""
    # Placeholder - would need to integrate actual Grobid extraction
    # For now, create a basic metadata structure
    metadata = {
        "filename": pdf_path.name,
        "title": "",
        "authors": [],
        "year": "",
        "journal": "",
        "doi": "",
        "extraction_method": "placeholder"
    }
    
    with open(metadata_output, 'w') as f:
        json.dump(metadata, f, indent=2)


def chunk_text(text_file: Path, metadata_file: Path, chunks_output: Path):
    """Chunk the extracted text."""
    # Load text and metadata
    with open(text_file, 'r', encoding='utf-8') as f:
        text_data = json.load(f)
    
    with open(metadata_file, 'r') as f:
        metadata = json.load(f)
    
    # Simple chunking (placeholder - could integrate more sophisticated chunking)
    text = text_data.get('text', '')
    
    # Split into chunks of roughly 1000 characters
    chunk_size = 1000
    chunks = []
    
    for i in range(0, len(text), chunk_size):
        chunk_text = text[i:i + chunk_size]
        chunks.append({
            "chunk_id": f"chunk_{len(chunks)}",
            "text": chunk_text,
            "start_char": i,
            "end_char": min(i + chunk_size, len(text))
        })
    
    chunks_data = {
        "metadata": metadata,
        "total_chunks": len(chunks),
        "chunks": chunks
    }
    
    with open(chunks_output, 'w', encoding='utf-8') as f:
        json.dump(chunks_data, f, indent=2, ensure_ascii=False)


def ingest_to_vector_db(chunks_file: Path, vector_db_dir: Path, pdf_hash: str):
    """Ingest chunks into vector database."""
    # Placeholder for vector database ingestion
    # This would typically involve embedding the text and storing in a vector database
    
    # Create a simple marker file for now
    ingestion_marker = vector_db_dir / f"{pdf_hash}_embedded.done"
    
    with open(ingestion_marker, 'w') as f:
        json.dump({
            "pdf_hash": pdf_hash,
            "chunks_file": str(chunks_file),
            "ingestion_timestamp": str(Path(chunks_file).stat().st_mtime),
            "status": "completed"
        }, f, indent=2)


def create_summary_json(pdf_hash: str, pdf_paths: List[Path], input_dir: Path, 
                       hash_dir: Path, processing_summary: Dict):
    """Create summary JSON with all metadata."""
    relative_paths = get_relative_paths(pdf_paths, input_dir)
    
    summary = {
        "pdf_hash": pdf_hash,
        "input_dir": str(input_dir),
        "relative_paths": relative_paths,
        "total_copies_found": len(pdf_paths),
        "processing_summary": processing_summary,
        "output_directory": str(hash_dir)
    }
    
    summary_file = hash_dir / "summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    return summary_file


def main():
    parser = argparse.ArgumentParser(description="Process a corpus of PDFs with hash-based organization")
    parser.add_argument("input_dir", type=Path, help="Input directory containing PDFs")
    parser.add_argument("output_dir", type=Path, help="Output directory for processed files")
    parser.add_argument("--resume", action="store_true", help="Resume processing, skip already processed PDFs")
    
    args = parser.parse_args()
    
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    
    if not input_dir.exists():
        print(f"Error: Input directory {input_dir} does not exist")
        sys.exit(1)
    
    print(f"Processing PDFs from: {input_dir}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory structure
    documents_dir, vector_db_dir = create_output_structure(output_dir)
    
    # Find all PDFs and group by hash
    print("Discovering PDFs...")
    pdf_map = find_all_pdfs(input_dir)
    
    print(f"Found {sum(len(paths) for paths in pdf_map.values())} PDF files")
    print(f"Unique PDFs (by hash): {len(pdf_map)}")
    
    # Process each unique PDF
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = Path(temp_dir)
        
        for pdf_hash, pdf_paths in pdf_map.items():
            print(f"\nProcessing PDF hash {pdf_hash} ({len(pdf_paths)} copies)")
            
            # Create hash-based directory
            hash_dir = documents_dir / pdf_hash
            
            # Skip if already processed and resume flag is set
            if args.resume and (hash_dir / "summary.json").exists():
                print(f"  Skipping {pdf_hash} (already processed)")
                continue
            
            hash_dir.mkdir(exist_ok=True)
            
            # Use the first copy for processing (they're all identical by hash)
            primary_pdf = pdf_paths[0]
            
            print(f"  Primary file: {primary_pdf.relative_to(input_dir)}")
            if len(pdf_paths) > 1:
                print(f"  Additional copies:")
                for path in pdf_paths[1:]:
                    print(f"    {path.relative_to(input_dir)}")
            
            # Run the processing pipeline
            processing_summary = run_pdf_processing_pipeline(primary_pdf, hash_dir, temp_dir)
            
            # Create summary JSON
            summary_file = create_summary_json(pdf_hash, pdf_paths, input_dir, hash_dir, processing_summary)
            print(f"  Created summary: {summary_file}")
            
            # Ingest into vector database if processing was successful
            if processing_summary.get("status") == "success":
                chunks_file = hash_dir / "chunks.json"
                if chunks_file.exists():
                    print(f"  Ingesting into vector database...")
                    ingest_to_vector_db(chunks_file, vector_db_dir, pdf_hash)
            
    print(f"\nProcessing complete! Results saved to: {output_dir}")
    print(f"  Documents: {documents_dir}")
    print(f"  Vector DB: {vector_db_dir}")


if __name__ == "__main__":
    main()
