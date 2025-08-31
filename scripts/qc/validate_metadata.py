#!/usr/bin/env python3

import json
import argparse
from pathlib import Path
import os

def validate_metadata(metadata_dir):
    """Validate extracted metadata files"""
    
    metadata_files = list(Path(metadata_dir).glob("*_metadata.json"))
    
    print(f"=== METADATA VALIDATION REPORT ===")
    print(f"Found {len(metadata_files)} metadata files")
    print()
    
    total_files = 0
    files_with_title = 0
    files_with_authors = 0
    files_with_year = 0
    files_with_doi = 0
    files_with_errors = 0
    
    for metadata_file in metadata_files:
        total_files += 1
        
        with open(metadata_file, 'r', encoding='utf-8') as f:
            metadata = json.load(f)
        
        filename = metadata.get('filename', metadata_file.stem)
        
        print(f"--- {filename} ---")
        
        # Check for errors
        if 'error' in metadata:
            files_with_errors += 1
            print(f"  ERROR: {metadata['error']}")
        
        # Check fields
        title = metadata.get('title', '')
        if title and title.strip():
            files_with_title += 1
            print(f"  Title: {title[:80]}...")
        else:
            print(f"  Title: [MISSING]")
        
        authors = metadata.get('authors', [])
        if authors:
            files_with_authors += 1
            print(f"  Authors: {', '.join(authors)}")
        else:
            print(f"  Authors: [MISSING]")
        
        year = metadata.get('year')
        if year:
            files_with_year += 1
            print(f"  Year: {year}")
        else:
            print(f"  Year: [MISSING]")
        
        doi = metadata.get('doi', '')
        if doi and doi.strip():
            files_with_doi += 1
            print(f"  DOI: {doi}")
        else:
            print(f"  DOI: [MISSING]")
        
        print()
    
    # Summary statistics
    print("=== SUMMARY ===")
    print(f"Total files: {total_files}")
    print(f"Files with errors: {files_with_errors} ({files_with_errors/total_files*100:.1f}%)")
    print(f"Files with title: {files_with_title} ({files_with_title/total_files*100:.1f}%)")
    print(f"Files with authors: {files_with_authors} ({files_with_authors/total_files*100:.1f}%)")
    print(f"Files with year: {files_with_year} ({files_with_year/total_files*100:.1f}%)")
    print(f"Files with DOI: {files_with_doi} ({files_with_doi/total_files*100:.1f}%)")

def main():
    parser = argparse.ArgumentParser(description='Validate extracted metadata')
    parser.add_argument('--metadata-dir', default='output/metadata', 
                       help='Directory containing metadata JSON files')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.metadata_dir):
        print(f"Metadata directory {args.metadata_dir} does not exist")
        return
    
    validate_metadata(args.metadata_dir)

if __name__ == "__main__":
    main()