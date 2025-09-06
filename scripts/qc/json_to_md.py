import json
from pathlib import Path
from docling_core.datamodel.document import DoclingDocument
import argparse

def main(input_path, output_path):
    input_path = Path(input_path)
    output_path = Path(output_path)
    print(f"Loading Docling Document from: {input_path}")
    with input_path.open("r") as fp:
        doc_dict = json.load(fp)
    docling_doc = DoclingDocument.model_validate(doc_dict)
    print("Exporting to Markdown...")
    markdown_text = docling_doc.export_to_markdown()
    with output_path.open("w", encoding="utf-8") as fp:
        fp.write(markdown_text)
    print(f"Conversion successful. Markdown saved to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Docling JSON to Markdown.")
    parser.add_argument("input_path", help="Path to the Docling JSON file.")
    parser.add_argument("output_path", help="Path to save the Markdown file.")
    args = parser.parse_args()
    main(args.input_path, args.output_path)