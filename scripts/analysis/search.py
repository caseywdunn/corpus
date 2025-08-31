#!/usr/bin/env python3

import argparse
import lancedb
import numpy as np
from openai import OpenAI
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

def search_corpus(query, top_k=5, db_path="output/lancedb"):
    """Search the document corpus using semantic similarity"""
    
    # Initialize OpenAI client
    client = OpenAI()
    
    # Get embedding for the query
    response = client.embeddings.create(
        model="text-embedding-3-small",
        input=[query]
    )
    query_embedding = response.data[0].embedding
    
    # Connect to LanceDB
    try:
        db = lancedb.connect(db_path)
        table = db.open_table("document_chunks")
    except Exception as e:
        print(f"Error connecting to database: {e}")
        return []
    
    # Perform similarity search
    try:
        results = table.search(query_embedding).limit(top_k).to_list()
        return results
    except Exception as e:
        print(f"Error performing search: {e}")
        return []

def format_results(results, query):
    """Format search results for display"""
    
    print(f"=== SEARCH RESULTS for: '{query}' ===")
    print(f"Found {len(results)} results")
    print()
    
    for i, result in enumerate(results, 1):
        metadata = result.get('metadata', {})
        text = result.get('text', '')
        score = result.get('_distance', 'N/A')
        
        print(f"--- Result {i} (Score: {score}) ---")
        print(f"Document: {metadata.get('filename', 'Unknown')}")
        print(f"Title: {metadata.get('title', 'Unknown')}")
        
        authors = metadata.get('authors', [])
        if authors:
            print(f"Authors: {', '.join(authors)}")
        
        year = metadata.get('year')
        if year:
            print(f"Year: {year}")
        
        # Show text preview
        preview = text[:300] + "..." if len(text) > 300 else text
        print(f"Text: {preview}")
        print()

def main():
    parser = argparse.ArgumentParser(description='Search the document corpus')
    parser.add_argument('query', help='Search query')
    parser.add_argument('--top-k', type=int, default=5, 
                       help='Number of results to return (default: 5)')
    parser.add_argument('--db-path', default='output/lancedb',
                       help='Path to LanceDB database')
    
    args = parser.parse_args()
    
    # Perform search
    results = search_corpus(args.query, args.top_k, args.db_path)
    
    # Display results
    format_results(results, args.query)

if __name__ == "__main__":
    main()