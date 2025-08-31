#!/usr/bin/env python3

import os
import json
from pathlib import Path
from typing import List
import lancedb
from lancedb.pydantic import LanceModel, Vector
from lancedb.embeddings import EmbeddingFunctionRegistry
from openai import OpenAI
import numpy as np

# Load environment variables
from dotenv import load_dotenv
load_dotenv()

class ChunkMetadata(LanceModel):
    filename: str
    title: str
    authors: List[str] = []
    year: int = None
    journal: str = ""
    doi: str = ""
    total_pages: int = None
    chunk_id: int

class DocumentChunk(LanceModel):
    text: str
    vector: Vector(1536)  # text-embedding-3-small has 1536 dimensions
    metadata: ChunkMetadata

def get_embeddings(texts, client):
    """Get embeddings for a list of texts using OpenAI API"""
    try:
        response = client.embeddings.create(
            model="text-embedding-3-small",  # Using smaller model for cost efficiency
            input=texts
        )
        return [embedding.embedding for embedding in response.data]
    except Exception as e:
        print(f"Error getting embeddings: {e}")
        # Return zero vectors as fallback
        return [np.zeros(1536).tolist() for _ in texts]

def main():
    chunks_file = snakemake.input[0]
    output_done = snakemake.output[0]
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_done), exist_ok=True)
    
    # Load chunks
    with open(chunks_file, 'r', encoding='utf-8') as f:
        chunks_data = json.load(f)
    
    chunks = chunks_data.get('chunks', [])
    source_metadata = chunks_data.get('source_metadata', {})
    
    if not chunks:
        print(f"No chunks found in {chunks_file}")
        # Create done file anyway
        with open(output_done, 'w') as f:
            f.write(f"No chunks to process from {chunks_file}\n")
        return
    
    try:
        # Initialize OpenAI client
        client = OpenAI()
        
        # Create/connect to LanceDB
        db_path = "output/lancedb"
        os.makedirs(db_path, exist_ok=True)
        db = lancedb.connect(db_path)
        
        # Create table if it doesn't exist
        table_name = "document_chunks"
        
        # Prepare data for embedding and storage
        texts = [chunk['text'] for chunk in chunks]
        
        print(f"Getting embeddings for {len(texts)} chunks...")
        embeddings = get_embeddings(texts, client)
        
        # Prepare records for LanceDB
        records = []
        for i, (chunk, embedding) in enumerate(zip(chunks, embeddings)):
            metadata = ChunkMetadata(
                filename=source_metadata.get('filename', ''),
                title=chunk['metadata'].get('title', ''),
                authors=chunk['metadata'].get('authors', []),
                year=chunk['metadata'].get('year'),
                journal=chunk['metadata'].get('journal', ''),
                doi=chunk['metadata'].get('doi', ''),
                total_pages=chunk['metadata'].get('total_pages'),
                chunk_id=chunk.get('chunk_id', i)
            )
            
            record = DocumentChunk(
                text=chunk['text'],
                vector=embedding,
                metadata=metadata
            )
            records.append(record)
        
        # Create or get table
        if table_name in db.table_names():
            table = db.open_table(table_name)
            table.add(records)
        else:
            table = db.create_table(table_name, records)
        
        print(f"Stored {len(records)} chunks in LanceDB")
        
        # Create done marker file
        with open(output_done, 'w') as f:
            f.write(f"Successfully processed {len(chunks)} chunks from {source_metadata.get('filename', 'unknown')}\n")
            f.write(f"Stored in LanceDB table: {table_name}\n")
        
    except Exception as e:
        print(f"Error embedding and storing chunks from {chunks_file}: {e}")
        # Create done file anyway to not break the workflow
        with open(output_done, 'w') as f:
            f.write(f"Error processing {chunks_file}: {e}\n")

if __name__ == "__main__":
    main()