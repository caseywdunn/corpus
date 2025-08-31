#!/usr/bin/env python3

import json
import argparse
from pathlib import Path

def visualize_docling_parse(text_file, figures_file):
    """Visualize the results of docling parsing"""
    
    # Load text data
    with open(text_file, 'r', encoding='utf-8') as f:
        text_data = json.load(f)
    
    # Load figures data
    with open(figures_file, 'r', encoding='utf-8') as f:
        figures_data = json.load(f)
    
    print(f"=== Document: {text_data.get('title', 'Unknown')} ===")
    print(f"Pages: {text_data.get('pages', 'Unknown')}")
    print(f"Figures found: {len(figures_data)}")
    print()
    
    print("=== TEXT PREVIEW ===")
    text = text_data.get('text', '')
    preview = text[:500] + "..." if len(text) > 500 else text
    print(preview)
    print()
    
    print("=== FIGURES ===")
    for figure in figures_data:
        print(f"Figure {figure.get('figure_id', 'unknown')} (Page {figure.get('page', '?')})")
        caption = figure.get('caption', '')
        if caption:
            print(f"  Caption: {caption[:100]}...")
        print(f"  BBox: {figure.get('bbox', 'N/A')}")
        print()

def main():
    parser = argparse.ArgumentParser(description='Visualize docling parsing results')
    parser.add_argument('--text', required=True, help='Path to text JSON file')
    parser.add_argument('--figures', required=True, help='Path to figures JSON file')
    
    args = parser.parse_args()
    
    visualize_docling_parse(args.text, args.figures)

if __name__ == "__main__":
    main()