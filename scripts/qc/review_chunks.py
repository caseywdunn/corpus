#!/usr/bin/env python3

import json
import argparse
from pathlib import Path
import os

def review_chunks(chunks_dir):
    """Review chunking and tokenization quality"""
    
    chunks_files = list(Path(chunks_dir).glob("*_chunks.json"))
    
    print(f"=== CHUNKING QUALITY REVIEW ===")
    print(f"Found {len(chunks_files)} chunk files")
    print()
    
    total_documents = 0
    total_chunks = 0
    total_tokens = 0
    
    for chunks_file in chunks_files:
        total_documents += 1
        
        with open(chunks_file, 'r', encoding='utf-8') as f:
            chunks_data = json.load(f)
        
        chunks = chunks_data.get('chunks', [])
        doc_name = chunks_data.get('source_metadata', {}).get('filename', chunks_file.stem)
        
        print(f"--- {doc_name} ---")
        
        if 'error' in chunks_data:
            print(f"  ERROR: {chunks_data['error']}")
            continue
        
        doc_chunks = len(chunks)
        total_chunks += doc_chunks
        
        if doc_chunks == 0:
            print(f"  No chunks created")
            continue
        
        # Calculate statistics
        chunk_sizes = [chunk.get('tokens', len(chunk.get('text', '').split())) for chunk in chunks]
        avg_tokens = sum(chunk_sizes) / len(chunk_sizes) if chunk_sizes else 0
        total_tokens += sum(chunk_sizes)
        
        print(f"  Chunks: {doc_chunks}")
        print(f"  Avg tokens/chunk: {avg_tokens:.1f}")
        print(f"  Min tokens: {min(chunk_sizes) if chunk_sizes else 0}")
        print(f"  Max tokens: {max(chunk_sizes) if chunk_sizes else 0}")
        
        # Show sample chunks
        print(f"  Sample chunks:")
        for i, chunk in enumerate(chunks[:2]):  # Show first 2 chunks
            text_preview = chunk.get('text', '')[:100] + "..."
            tokens = chunk.get('tokens', 'unknown')
            print(f"    Chunk {i+1} ({tokens} tokens): {text_preview}")
        
        print()
    
    # Summary statistics
    print("=== SUMMARY ===")
    print(f"Total documents: {total_documents}")
    print(f"Total chunks: {total_chunks}")
    print(f"Avg chunks per document: {total_chunks/total_documents:.1f}" if total_documents > 0 else "N/A")
    print(f"Total tokens: {total_tokens:,}")
    print(f"Avg tokens per chunk: {total_tokens/total_chunks:.1f}" if total_chunks > 0 else "N/A")

def main():
    parser = argparse.ArgumentParser(description='Review chunking quality')
    parser.add_argument('--chunks-dir', default='output/chunks', 
                       help='Directory containing chunk JSON files')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.chunks_dir):
        print(f"Chunks directory {args.chunks_dir} does not exist")
        return
    
    review_chunks(args.chunks_dir)

if __name__ == "__main__":
    main()