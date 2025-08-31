#!/usr/bin/env python3

import argparse
import lancedb
from openai import OpenAI
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

def search_relevant_chunks(query, top_k=3, db_path="output/lancedb"):
    """Search for relevant chunks to use as context"""
    
    client = OpenAI()
    
    # Get embedding for the query
    response = client.embeddings.create(
        model="text-embedding-3-small",
        input=[query]
    )
    query_embedding = response.data[0].embedding
    
    # Connect to LanceDB and search
    try:
        db = lancedb.connect(db_path)
        table = db.open_table("document_chunks")
        results = table.search(query_embedding).limit(top_k).to_list()
        return results
    except Exception as e:
        print(f"Error searching database: {e}")
        return []

def chat_with_corpus(question, db_path="output/lancedb"):
    """Chat with the corpus using RAG (Retrieval-Augmented Generation)"""
    
    # Search for relevant context
    relevant_chunks = search_relevant_chunks(question, top_k=3, db_path=db_path)
    
    if not relevant_chunks:
        print("No relevant documents found in the corpus.")
        return
    
    # Build context from relevant chunks
    context_parts = []
    sources = []
    
    for chunk in relevant_chunks:
        text = chunk.get('text', '')
        metadata = chunk.get('metadata', {})
        
        context_parts.append(text)
        source = f"{metadata.get('filename', 'Unknown')} - {metadata.get('title', 'Unknown')}"
        if source not in sources:
            sources.append(source)
    
    context = "\n\n".join(context_parts)
    
    # Create prompt with context
    system_prompt = """You are a helpful research assistant that answers questions based on the provided scientific literature context. 
    Always base your answers on the information in the context. If the context doesn't contain enough information to answer the question, say so.
    When possible, reference specific details from the papers."""
    
    user_prompt = f"""Context from scientific papers:
    
{context}

Question: {question}

Please provide a comprehensive answer based on the context above."""
    
    # Get response from OpenAI
    client = OpenAI()
    
    try:
        response = client.chat.completions.create(
            model="gpt-4",
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt}
            ],
            temperature=0.1
        )
        
        answer = response.choices[0].message.content
        
        # Display results
        print(f"=== QUESTION ===")
        print(question)
        print()
        
        print(f"=== ANSWER ===")
        print(answer)
        print()
        
        print(f"=== SOURCES ===")
        for source in sources:
            print(f"• {source}")
        
    except Exception as e:
        print(f"Error getting response from OpenAI: {e}")

def interactive_chat(db_path="output/lancedb"):
    """Start an interactive chat session"""
    
    print("=== INTERACTIVE CORPUS CHAT ===")
    print("Ask questions about your document corpus. Type 'quit' to exit.")
    print()
    
    while True:
        try:
            question = input("Question: ").strip()
            
            if question.lower() in ['quit', 'exit', 'q']:
                print("Goodbye!")
                break
            
            if not question:
                continue
            
            print()
            chat_with_corpus(question, db_path)
            print("\n" + "="*50 + "\n")
            
        except KeyboardInterrupt:
            print("\nGoodbye!")
            break
        except Exception as e:
            print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(description='Chat with your document corpus')
    parser.add_argument('--question', help='Single question to ask')
    parser.add_argument('--interactive', action='store_true', 
                       help='Start interactive chat session')
    parser.add_argument('--db-path', default='output/lancedb',
                       help='Path to LanceDB database')
    
    args = parser.parse_args()
    
    if args.question:
        chat_with_corpus(args.question, args.db_path)
    elif args.interactive:
        interactive_chat(args.db_path)
    else:
        # Default to interactive if no question provided
        interactive_chat(args.db_path)

if __name__ == "__main__":
    main()