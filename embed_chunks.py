#!/usr/bin/env python3
"""
Standalone script to embed and store chunks from the new hash-based structure.
Can be used as part of the main processing pipeline or run separately.
"""

import os
import json
from pathlib import Path
from typing import List, Dict
import lancedb
from lancedb.pydantic import LanceModel, Vector
import numpy as np

# Load environment variables
from dotenv import load_dotenv
load_dotenv()

class ChunkMetadata(LanceModel):
    pdf_hash: str
    filename: str = ""
    title: str = ""
    authors: List[str] = []
    year: int = None
    journal: str = ""
    doi: str = ""
    total_pages: int = None
    chunk_id: str = ""
    relative_paths: List[str] = []

class DocumentChunk(LanceModel):
    text: str
    vector: Vector(1536)  # text-embedding-3-small has 1536 dimensions
    metadata: ChunkMetadata

def get_embeddings(texts: List[str]) -> List[List[float]]:
    """Get embeddings for a list of texts using OpenAI API"""
    try:
        from openai import OpenAI
        client = OpenAI()
        
        # Process in batches to avoid API limits
        batch_size = 100
        all_embeddings = []
        
        for i in range(0, len(texts), batch_size):
            batch = texts[i:i + batch_size]
            print(f"  Getting embeddings for batch {i//batch_size + 1}/{(len(texts) + batch_size - 1)//batch_size}")
            
            response = client.embeddings.create(
                model="text-embedding-3-small",
                input=batch
            )
            
            batch_embeddings = [embedding.embedding for embedding in response.data]
            all_embeddings.extend(batch_embeddings)
        
        return all_embeddings
        
    except Exception as e:
        print(f"Error getting embeddings: {e}")
        # Return zero vectors as fallback
        return [np.zeros(1536).tolist() for _ in texts]

def load_document_data(hash_dir: Path) -> Dict:
    """Load all data for a document from its hash directory."""
    data = {}
    
    # Load summary
    summary_file = hash_dir / "summary.json"
    if summary_file.exists():
        with open(summary_file, 'r') as f:
            data['summary'] = json.load(f)
    
    # Load chunks
    chunks_file = hash_dir / "chunks.json"
    if chunks_file.exists():
        with open(chunks_file, 'r', encoding='utf-8') as f:
            data['chunks'] = json.load(f)
    
    # Load metadata
    metadata_file = hash_dir / "metadata.json"
    if metadata_file.exists():
        with open(metadata_file, 'r') as f:
            data['metadata'] = json.load(f)
    
    return data

def embed_document(hash_dir: Path, vector_db_dir: Path) -> bool:
    """Embed a single document and store in vector database."""
    pdf_hash = hash_dir.name
    print(f"Embedding document {pdf_hash}...")
    
    # Load document data
    doc_data = load_document_data(hash_dir)
    
    if 'chunks' not in doc_data:
        print(f"  No chunks found for {pdf_hash}")
        return False
    
    chunks_data = doc_data['chunks']
    chunks = chunks_data.get('chunks', [])
    
    if not chunks:
        print(f"  No chunks to embed for {pdf_hash}")
        return False
    
    # Extract metadata
    summary = doc_data.get('summary', {})
    metadata = doc_data.get('metadata', {})
    
    # Prepare texts for embedding
    texts = [chunk['text'] for chunk in chunks]
    
    print(f"  Getting embeddings for {len(texts)} chunks...")
    embeddings = get_embeddings(texts)
    
    # Connect to LanceDB
    db = lancedb.connect(str(vector_db_dir))
    table_name = "document_chunks"
    
    # Prepare records
    records = []
    for chunk, embedding in zip(chunks, embeddings):
        chunk_metadata = ChunkMetadata(
            pdf_hash=pdf_hash,
            filename=summary.get('relative_paths', [''])[0] if summary.get('relative_paths') else '',
            title=metadata.get('title', ''),
            authors=metadata.get('authors', []),
            year=metadata.get('year'),
            journal=metadata.get('journal', ''),
            doi=metadata.get('doi', ''),
            total_pages=chunks_data.get('metadata', {}).get('total_pages'),
            chunk_id=chunk.get('chunk_id', ''),
            relative_paths=summary.get('relative_paths', [])
        )
        
        record = DocumentChunk(
            text=chunk['text'],
            vector=embedding,
            metadata=chunk_metadata
        )
        records.append(record)
    
    # Store in database
    if table_name in db.table_names():
        table = db.open_table(table_name)
        table.add(records)
    else:
        table = db.create_table(table_name, records)
    
    print(f"  Stored {len(records)} chunks in vector database")
    
    # Create completion marker
    marker_file = vector_db_dir / f"{pdf_hash}_embedded.done"
    with open(marker_file, 'w') as f:
        json.dump({
            "pdf_hash": pdf_hash,
            "chunks_count": len(records),
            "relative_paths": summary.get('relative_paths', []),
            "status": "completed"
        }, f, indent=2)
    
    return True

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Embed document chunks into vector database")
    parser.add_argument("output_dir", type=Path, help="Output directory containing documents/ and vector_db/")
    parser.add_argument("--pdf-hash", help="Process only specific PDF hash")
    parser.add_argument("--resume", action="store_true", help="Skip already embedded documents")
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    documents_dir = output_dir / "documents"
    vector_db_dir = output_dir / "vector_db"
    
    if not documents_dir.exists():
        print(f"Error: Documents directory {documents_dir} does not exist")
        return 1
    
    vector_db_dir.mkdir(exist_ok=True)
    
    # Find documents to process
    if args.pdf_hash:
        hash_dirs = [documents_dir / args.pdf_hash]
        if not hash_dirs[0].exists():
            print(f"Error: PDF hash directory {args.pdf_hash} not found")
            return 1
    else:
        hash_dirs = [d for d in documents_dir.iterdir() if d.is_dir()]
    
    print(f"Found {len(hash_dirs)} document(s) to process")
    
    success_count = 0
    for hash_dir in hash_dirs:
        pdf_hash = hash_dir.name
        marker_file = vector_db_dir / f"{pdf_hash}_embedded.done"
        
        # Skip if already processed and resume flag is set
        if args.resume and marker_file.exists():
            print(f"Skipping {pdf_hash} (already embedded)")
            continue
        
        try:
            if embed_document(hash_dir, vector_db_dir):
                success_count += 1
        except Exception as e:
            print(f"Error embedding document {pdf_hash}: {e}")
    
    print(f"\nEmbedding complete! Successfully embedded {success_count}/{len(hash_dirs)} documents")
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())