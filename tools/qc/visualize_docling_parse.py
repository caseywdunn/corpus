# from https://github.com/docling-project/docling-parse
import argparse
from docling_core.types.doc.page import TextCellUnit
from docling_parse.pdf_parser import DoclingPdfParser, PdfDocument

def main(args):
    pdf_parser = DoclingPdfParser()
    pdf_doc: PdfDocument = pdf_parser.load(path_or_stream=args.filename)

    # Map string to TextCellUnit
    level_map = {
        "char": TextCellUnit.CHAR,
        "word": TextCellUnit.WORD,
        "line": TextCellUnit.LINE
    }
    cell_unit = level_map[args.level]

    from PIL import ImageDraw
    for page_no, pred_page in pdf_doc.iterate_pages():
        print(f"Page {page_no}")
        img = pred_page.render_as_image(cell_unit=cell_unit)
        draw = ImageDraw.Draw(img)
        img_height = img.height
        for cell in pred_page.iterate_cells(unit_type=cell_unit):
            print(cell.rect, ": ", cell.text)
            # Convert cell.rect to (x0, y0, x1, y1) and flip y for PIL
            x0 = min(getattr(cell.rect, "r_x0", 0), getattr(cell.rect, "r_x2", 0))
            x1 = max(getattr(cell.rect, "r_x0", 0), getattr(cell.rect, "r_x2", 0))
            y0_raw = min(getattr(cell.rect, "r_y0", 0), getattr(cell.rect, "r_y2", 0))
            y1_raw = max(getattr(cell.rect, "r_y0", 0), getattr(cell.rect, "r_y2", 0))
            y0 = img_height - y1_raw
            y1 = img_height - y0_raw
            rect = (x0, y0, x1, y1)
            draw.rectangle(rect, outline="red", width=2)
        img.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize Docling parse of a PDF.")
    parser.add_argument("filename", type=str, help="Path to the PDF file.")
    parser.add_argument("--level", type=str, choices=["char", "word", "line"], default="word", help="Cell unit level to visualize.")
    args = parser.parse_args()
    main(args)