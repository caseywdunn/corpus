#!/usr/bin/env python3

import os
import json
import tiktoken
from pathlib import Path
from docling.chunking import HybridChunker
from docling.document_converter import DocumentConverter

class CustomTokenizer:
    def __init__(self, model_name="text-embedding-3-large"):
        self.tokenizer = tiktoken.encoding_for_model("gpt-4")
    
    def count_tokens(self, text):
        return len(self.tokenizer.encode(text))

def main():
    text_file = snakemake.input.text
    metadata_file = snakemake.input.metadata
    output_chunks = snakemake.output[0]
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_chunks), exist_ok=True)
    
    # Load text and metadata
    with open(text_file, 'r', encoding='utf-8') as f:
        text_data = json.load(f)
    
    with open(metadata_file, 'r', encoding='utf-8') as f:
        metadata = json.load(f)
    
    # Configuration for chunking
    MAX_TOKENS = 8191
    
    try:
        # Initialize tokenizer and chunker
        tokenizer = CustomTokenizer()
        
        # For now, implement simple text splitting since HybridChunker needs a docling document
        text = text_data.get('text', '')
        
        # Simple chunking approach - split by paragraphs and combine up to max tokens
        paragraphs = text.split('\n\n')
        chunks = []
        current_chunk = ""
        current_tokens = 0
        
        for paragraph in paragraphs:
            paragraph = paragraph.strip()
            if not paragraph:
                continue
                
            para_tokens = tokenizer.count_tokens(paragraph)
            
            # If adding this paragraph would exceed max tokens, save current chunk
            if current_tokens + para_tokens > MAX_TOKENS and current_chunk:
                chunks.append({
                    "text": current_chunk.strip(),
                    "tokens": current_tokens,
                    "chunk_id": len(chunks)
                })
                current_chunk = paragraph
                current_tokens = para_tokens
            else:
                if current_chunk:
                    current_chunk += "\n\n" + paragraph
                    current_tokens += para_tokens
                else:
                    current_chunk = paragraph
                    current_tokens = para_tokens
        
        # Add the last chunk
        if current_chunk:
            chunks.append({
                "text": current_chunk.strip(),
                "tokens": current_tokens,
                "chunk_id": len(chunks)
            })
        
        # Add metadata to each chunk
        for chunk in chunks:
            chunk["metadata"] = {
                "filename": metadata.get("filename", ""),
                "title": metadata.get("title", text_data.get("title", "")),
                "authors": metadata.get("authors", []),
                "year": metadata.get("year"),
                "journal": metadata.get("journal", ""),
                "doi": metadata.get("doi", ""),
                "total_pages": text_data.get("pages")
            }
        
        # Save chunks
        output_data = {
            "chunks": chunks,
            "total_chunks": len(chunks),
            "source_metadata": metadata
        }
        
        with open(output_chunks, 'w', encoding='utf-8') as f:
            json.dump(output_data, f, indent=2, ensure_ascii=False)
        
        print(f"Created {len(chunks)} chunks from {text_file}")
        
    except Exception as e:
        print(f"Error chunking text from {text_file}: {e}")
        # Create empty chunks file to not break the workflow
        empty_chunks = {
            "chunks": [],
            "total_chunks": 0,
            "source_metadata": metadata,
            "error": str(e)
        }
        
        with open(output_chunks, 'w', encoding='utf-8') as f:
            json.dump(empty_chunks, f, indent=2)

if __name__ == "__main__":
    main()