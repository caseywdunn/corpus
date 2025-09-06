#!/usr/bin/env python3

import os
import json
import subprocess
import tempfile
import re
import time
from pathlib import Path
import urllib.parse
import urllib.request
from urllib.error import HTTPError, URLError

def extract_pdf_metadata(pdf_path):
    """Extract metadata using PyPDF2 or similar PDF parsing"""
    try:
        import PyPDF2
        with open(pdf_path, 'rb') as file:
            reader = PyPDF2.PdfReader(file)
            pdf_metadata = reader.metadata or {}
            
            # Extract basic PDF metadata
            extracted = {}
            if pdf_metadata.get('/Title'):
                extracted['title'] = str(pdf_metadata.get('/Title'))
            if pdf_metadata.get('/Author'):
                extracted['author'] = str(pdf_metadata.get('/Author'))
            if pdf_metadata.get('/Subject'):
                extracted['subject'] = str(pdf_metadata.get('/Subject'))
            if pdf_metadata.get('/Creator'):
                extracted['creator'] = str(pdf_metadata.get('/Creator'))
            if pdf_metadata.get('/Producer'):
                extracted['producer'] = str(pdf_metadata.get('/Producer'))
            if pdf_metadata.get('/CreationDate'):
                extracted['creation_date'] = str(pdf_metadata.get('/CreationDate'))
            
            return extracted
    except ImportError:
        print("PyPDF2 not available, trying alternative methods")
        return {}
    except Exception as e:
        print(f"Error extracting PDF metadata: {e}")
        return {}

def extract_metadata_from_text(text_content):
    """Extract metadata from the text content using pattern matching"""
    metadata = {}
    
    if not text_content:
        return metadata
    
    # Split into lines for easier processing
    lines = text_content.split('\n')
    first_50_lines = lines[:50]  # Look in first 50 lines for metadata
    text_sample = '\n'.join(first_50_lines)
    
    # Extract title - look for patterns that suggest a title
    # Title is often the first substantial line, or after author names
    title_candidates = []
    for i, line in enumerate(first_50_lines):
        line = line.strip()
        if len(line) > 20 and len(line) < 200:  # Reasonable title length
            # Skip lines that look like headers, page numbers, etc.
            if not re.match(r'^(page \d+|abstract|introduction|keywords|references)', line.lower()):
                # Skip lines that are all caps or mostly punctuation
                if not line.isupper() and len(re.sub(r'[^a-zA-Z0-9\s]', '', line)) > 10:
                    title_candidates.append((line, i))
    
    # Pick the best title candidate
    if title_candidates:
        # Prefer candidates from lines 0-10 that aren't author lines
        for title, line_num in title_candidates[:5]:
            if not re.search(r'\b(by|author|et al\.)', title.lower()):
                metadata['title'] = title
                break
        
        # Fallback to first candidate
        if 'title' not in metadata:
            metadata['title'] = title_candidates[0][0]
    
    # Extract authors - look for common author patterns
    author_patterns = [
        r'^([A-Z]{1,2}\.[A-Z]{1,2}\.?\s+[A-Z][A-Z]+)$',  # P.R. PUGH pattern
        r'^([A-Z][a-z]+\s+[A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)$',  # John Smith or John A Smith
        r'(?:by|author[s]?:?\s+)([A-Z][a-z]+(?:\s+[A-Z][a-z]*)*(?:\s+and\s+[A-Z][a-z]+(?:\s+[A-Z][a-z]*)*)*)',
        r'^([A-Z][a-z]+(?:\s+[A-Z]\.?\s*[A-Z][a-z]+)*(?:\s*,\s*[A-Z][a-z]+(?:\s+[A-Z]\.?\s*[A-Z][a-z]+)*)*)',
        r'([A-Z][a-z]+\s+et\s+al\.)',
    ]
    
    authors = []
    for pattern in author_patterns:
        for line in first_50_lines[:10]:  # Check first 10 lines for authors
            matches = re.findall(pattern, line.strip())
            for match in matches:
                if len(match) > 3 and len(match) < 100:  # Reasonable author name length
                    authors.extend(re.split(r',\s*|\s+and\s+', match))
    
    if authors:
        # Clean up author names
        cleaned_authors = []
        for author in authors[:5]:  # Max 5 authors
            author = author.strip().replace('et al.', '').strip()
            if len(author) > 2 and not author.lower() in ['by', 'author', 'authors']:
                cleaned_authors.append(author)
        metadata['authors'] = cleaned_authors
    
    # Extract year
    year_matches = re.findall(r'\b(19|20)\d{2}\b', text_sample)
    if year_matches:
        # Prefer years between 1980-2030
        valid_years = [int(match) for match in year_matches if 1980 <= int(match) <= 2030]
        if valid_years:
            metadata['year'] = valid_years[0]
    
    # Extract DOI
    doi_match = re.search(r'doi:?\s*(10\.\d+/[^\s]+)', text_sample, re.IGNORECASE)
    if doi_match:
        metadata['doi'] = doi_match.group(1)
    
    # Extract journal name - look for common journal patterns
    journal_patterns = [
        r'(?:published in|journal:|in)\s+([A-Z][a-zA-Z\s&]+?)(?:\s+vol|\s+\d+|\s*,|\s*$)',
        r'^([A-Z][a-zA-Z\s&]+Journal[a-zA-Z\s]*)',
    ]
    
    for pattern in journal_patterns:
        for line in first_50_lines[:20]:
            match = re.search(pattern, line.strip())
            if match:
                journal = match.group(1).strip()
                if 5 < len(journal) < 100:  # Reasonable journal name length
                    metadata['journal'] = journal
                    break
        if 'journal' in metadata:
            break
    
    return metadata

def query_crossref_api(title, authors=None, year=None):
    """Query Crossref API for paper metadata"""
    try:
        # Build query string - use bibliographic query for best results
        query_parts = []
        if title:
            query_parts.append(title.strip())
        if authors:
            for author in authors[:2]:  # Use first 2 authors
                query_parts.append(author.strip())
        if year:
            query_parts.append(str(year))
        
        query = ' '.join(query_parts)
        if not query:
            return None
        
        # URL encode the query
        encoded_query = urllib.parse.quote(query)
        
        # Build Crossref API URL
        url = f"https://api.crossref.org/works?query.bibliographic={encoded_query}&rows=3&mailto=researcher@example.com"
        
        # Make request
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Academic Metadata Extractor (mailto:researcher@example.com)')
        
        with urllib.request.urlopen(req, timeout=10) as response:
            data = json.loads(response.read().decode())
            
            if data.get('status') == 'ok' and data.get('message', {}).get('items'):
                items = data['message']['items']
                
                # Find the best match
                for item in items:
                    # Score based on title similarity and year match
                    item_title = item.get('title', [''])[0].lower() if item.get('title') else ''
                    title_similarity = 0
                    if title and item_title:
                        # Simple word overlap scoring
                        title_words = set(title.lower().split())
                        item_words = set(item_title.split())
                        if title_words and item_words:
                            title_similarity = len(title_words.intersection(item_words)) / len(title_words.union(item_words))
                    
                    year_match = True
                    if year and item.get('published-print', {}).get('date-parts'):
                        item_year = item['published-print']['date-parts'][0][0]
                        year_match = abs(int(year) - item_year) <= 1
                    elif year and item.get('published-online', {}).get('date-parts'):
                        item_year = item['published-online']['date-parts'][0][0]
                        year_match = abs(int(year) - item_year) <= 1
                    
                    # Accept if title similarity > 0.3 or if year matches closely
                    if title_similarity > 0.3 or (year_match and title_similarity > 0.1):
                        return format_crossref_metadata(item)
                
                # Fallback: return first result if no good match
                return format_crossref_metadata(items[0])
        
        return None
        
    except (HTTPError, URLError, json.JSONDecodeError, Exception) as e:
        print(f"Crossref API error: {e}")
        return None

def query_bhl_api(title, authors=None, year=None):
    """Query Biodiversity Heritage Library API for paper metadata"""
    try:
        # BHL API requires an API key - for now, skip if no key available
        # In production, you'd get a key from https://www.biodiversitylibrary.org/getapikey.aspx
        api_key = os.environ.get('BHL_API_KEY')
        if not api_key:
            return None
        
        # Build search parameters
        params = {}
        if title:
            params['title'] = title.strip()
        if authors:
            params['authorname'] = authors[0].strip()  # Use first author
        if year:
            params['year'] = str(year)
        
        if not params:
            return None
        
        # Build URL
        params['apikey'] = api_key
        params['format'] = 'json'
        query_string = urllib.parse.urlencode(params)
        url = f"https://www.biodiversitylibrary.org/api3?op=PublicationSearchAdvanced&{query_string}"
        
        # Make request
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Academic Metadata Extractor')
        
        with urllib.request.urlopen(req, timeout=10) as response:
            data = json.loads(response.read().decode())
            
            if data.get('Status') == 'ok' and data.get('Result'):
                results = data['Result']
                if results:
                    return format_bhl_metadata(results[0])
        
        return None
        
    except (HTTPError, URLError, json.JSONDecodeError, Exception) as e:
        print(f"BHL API error: {e}")
        return None

def format_crossref_metadata(item):
    """Convert Crossref API response to our metadata format"""
    metadata = {}
    
    # Title
    if item.get('title'):
        metadata['title'] = item['title'][0]
    
    # Authors
    authors = []
    if item.get('author'):
        for author in item['author'][:5]:  # Max 5 authors
            given = author.get('given', '')
            family = author.get('family', '')
            if given and family:
                authors.append(f"{given} {family}")
            elif family:
                authors.append(family)
    metadata['authors'] = authors
    
    # Year
    if item.get('published-print', {}).get('date-parts'):
        metadata['year'] = item['published-print']['date-parts'][0][0]
    elif item.get('published-online', {}).get('date-parts'):
        metadata['year'] = item['published-online']['date-parts'][0][0]
    
    # Journal
    if item.get('container-title'):
        metadata['journal'] = item['container-title'][0]
    
    # DOI
    if item.get('DOI'):
        metadata['doi'] = item['DOI']
    
    # Volume, Issue, Pages
    if item.get('volume'):
        metadata['volume'] = item['volume']
    if item.get('issue'):
        metadata['issue'] = item['issue']
    if item.get('page'):
        metadata['pages'] = item['page']
    
    # URL
    if item.get('URL'):
        metadata['url'] = item['URL']
    
    return metadata

def format_bhl_metadata(item):
    """Convert BHL API response to our metadata format"""
    metadata = {}
    
    # Extract relevant fields from BHL response
    if item.get('FullTitle'):
        metadata['title'] = item['FullTitle']
    if item.get('Authors'):
        metadata['authors'] = [author.get('Name', '') for author in item['Authors'] if author.get('Name')]
    if item.get('Date'):
        try:
            metadata['year'] = int(item['Date'])
        except (ValueError, TypeError):
            pass
    if item.get('ContainerTitle'):
        metadata['journal'] = item['ContainerTitle']
    if item.get('Url'):
        metadata['url'] = item['Url']
    
    return metadata

def generate_bibtex(metadata, entry_type='article'):
    """Generate BibTeX entry from metadata"""
    if not metadata.get('title'):
        return None
    
    # Generate citation key
    year = metadata.get('year', 'NoYear')
    first_author = 'Anonymous'
    if metadata.get('authors'):
        first_author = metadata['authors'][0].split()[-1]  # Last name
        first_author = re.sub(r'[^a-zA-Z]', '', first_author)
    
    citekey = f"{first_author}{year}"
    
    # Start BibTeX entry
    bibtex_lines = [f"@{entry_type}{{{citekey},"]
    
    # Add fields
    if metadata.get('title'):
        bibtex_lines.append(f'  title = {{{metadata["title"]}}},')
    
    if metadata.get('authors'):
        authors_str = ' and '.join(metadata['authors'])
        bibtex_lines.append(f'  author = {{{authors_str}}},')
    
    if metadata.get('journal'):
        bibtex_lines.append(f'  journal = {{{metadata["journal"]}}},')
    
    if metadata.get('year'):
        bibtex_lines.append(f'  year = {{{metadata["year"]}}},')
    
    if metadata.get('volume'):
        bibtex_lines.append(f'  volume = {{{metadata["volume"]}}},')
    
    if metadata.get('issue'):
        bibtex_lines.append(f'  number = {{{metadata["issue"]}}},')
    
    if metadata.get('pages'):
        bibtex_lines.append(f'  pages = {{{metadata["pages"]}}},')
    
    if metadata.get('doi'):
        bibtex_lines.append(f'  doi = {{{metadata["doi"]}}},')
    
    if metadata.get('url'):
        bibtex_lines.append(f'  url = {{{metadata["url"]}}},')
    
    # Close entry
    bibtex_lines.append("}")
    
    return '\n'.join(bibtex_lines)

def extract_metadata_comprehensive(pdf_path):
    """Extract metadata using multiple methods and combine results"""
    filename = Path(pdf_path).stem
    
    # Initialize with filename-based metadata
    metadata = {
        "title": "",
        "authors": [],
        "year": None,
        "journal": "",
        "doi": "",
        "filename": filename,
        "extraction_methods": []
    }
    
    # Method 1: Extract from PDF metadata
    pdf_meta = extract_pdf_metadata(pdf_path)
    if pdf_meta:
        metadata["extraction_methods"].append("pdf_metadata")
        if pdf_meta.get('title'):
            metadata['title'] = pdf_meta['title']
        if pdf_meta.get('author'):
            metadata['authors'] = [pdf_meta['author']]
    
    # Method 2: Extract from docling text content if available
    text_json_path = str(pdf_path).replace('/processed/', '/docling/').replace('_processed.pdf', '_text.json')
    try:
        if os.path.exists(text_json_path):
            with open(text_json_path, 'r', encoding='utf-8') as f:
                text_data = json.load(f)
                text_content = text_data.get('text', '') or text_data.get('content', '')
                
                if text_content:
                    metadata["extraction_methods"].append("text_analysis")
                    text_meta = extract_metadata_from_text(text_content)
                    
                    # Use text-based metadata if PDF metadata is empty
                    if not metadata.get('title') and text_meta.get('title'):
                        metadata['title'] = text_meta['title']
                    if not metadata.get('authors') and text_meta.get('authors'):
                        metadata['authors'] = text_meta['authors']
                    if not metadata.get('year') and text_meta.get('year'):
                        metadata['year'] = text_meta['year']
                    if not metadata.get('journal') and text_meta.get('journal'):
                        metadata['journal'] = text_meta['journal']
                    if not metadata.get('doi') and text_meta.get('doi'):
                        metadata['doi'] = text_meta['doi']
    except Exception as e:
        pass  # Silently continue if text content is not available
    
    # Method 3: External API lookups for enhanced metadata
    if metadata.get('title') or metadata.get('authors'):
        # Try Crossref first (free, comprehensive)
        crossref_meta = query_crossref_api(
            title=metadata.get('title'), 
            authors=metadata.get('authors'), 
            year=metadata.get('year')
        )
        
        if crossref_meta:
            metadata["extraction_methods"].append("crossref_api")
            # Use Crossref data to fill gaps or improve existing data
            if not metadata.get('title') and crossref_meta.get('title'):
                metadata['title'] = crossref_meta['title']
            if not metadata.get('authors') and crossref_meta.get('authors'):
                metadata['authors'] = crossref_meta['authors']
            if not metadata.get('year') and crossref_meta.get('year'):
                metadata['year'] = crossref_meta['year']
            if not metadata.get('journal') and crossref_meta.get('journal'):
                metadata['journal'] = crossref_meta['journal']
            if not metadata.get('doi') and crossref_meta.get('doi'):
                metadata['doi'] = crossref_meta['doi']
            
            # Add additional fields from Crossref
            if crossref_meta.get('volume'):
                metadata['volume'] = crossref_meta['volume']
            if crossref_meta.get('issue'):
                metadata['issue'] = crossref_meta['issue']
            if crossref_meta.get('pages'):
                metadata['pages'] = crossref_meta['pages']
            if crossref_meta.get('url'):
                metadata['url'] = crossref_meta['url']
            
            # Rate limiting - be polite to Crossref
            time.sleep(0.1)
        
        # Try BHL for biodiversity papers (if API key available)
        if not crossref_meta and (
            'biodiversity' in metadata.get('title', '').lower() or
            'species' in metadata.get('title', '').lower() or
            'taxonomy' in metadata.get('title', '').lower() or
            'genus' in metadata.get('title', '').lower()
        ):
            bhl_meta = query_bhl_api(
                title=metadata.get('title'),
                authors=metadata.get('authors'),
                year=metadata.get('year')
            )
            
            if bhl_meta:
                metadata["extraction_methods"].append("bhl_api")
                # Use BHL data to fill gaps
                if not metadata.get('title') and bhl_meta.get('title'):
                    metadata['title'] = bhl_meta['title']
                if not metadata.get('authors') and bhl_meta.get('authors'):
                    metadata['authors'] = bhl_meta['authors']
                if not metadata.get('year') and bhl_meta.get('year'):
                    metadata['year'] = bhl_meta['year']
                if not metadata.get('journal') and bhl_meta.get('journal'):
                    metadata['journal'] = bhl_meta['journal']
                if bhl_meta.get('url'):
                    metadata['url'] = bhl_meta['url']
    
    # Method 4: Filename-based extraction as fallback
    if not metadata.get('year'):
        year_match = re.search(r'(\d{4})', filename)
        if year_match:
            metadata["year"] = int(year_match.group(1))
            if "filename_parsing" not in metadata["extraction_methods"]:
                metadata["extraction_methods"].append("filename_parsing")
    
    if not metadata.get('authors'):
        parts = filename.replace('_', ' ').replace('-', ' ').split()
        potential_authors = []
        for part in parts:
            if part.istitle() and len(part) > 2 and not part.isdigit():
                potential_authors.append(part)
        
        if potential_authors:
            metadata["authors"] = potential_authors[:3]
            if "filename_parsing" not in metadata["extraction_methods"]:
                metadata["extraction_methods"].append("filename_parsing")
    
    # Generate BibTeX record
    bibtex_entry = generate_bibtex(metadata)
    if bibtex_entry:
        metadata['bibtex'] = bibtex_entry
    
    return metadata

def main():
    input_pdf = snakemake.input[0]
    output_metadata = snakemake.output[0]
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_metadata), exist_ok=True)
    
    try:
        # Extract metadata using comprehensive approach
        metadata = extract_metadata_comprehensive(input_pdf)
        
        # Save metadata
        with open(output_metadata, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)
        
        print(f"Extracted metadata from {input_pdf}")
        if metadata.get("extraction_methods"):
            print(f"  Methods used: {', '.join(metadata['extraction_methods'])}")
        if metadata.get("title"):
            print(f"  Title: {metadata['title'][:100]}...")
        if metadata.get("authors"):
            print(f"  Authors: {', '.join(metadata['authors'][:3])}")
        
    except Exception as e:
        print(f"Error extracting metadata from {input_pdf}: {e}")
        import traceback
        traceback.print_exc()
        
        # Create empty metadata file to not break the workflow
        empty_metadata = {
            "title": "",
            "authors": [],
            "year": None,
            "journal": "",
            "doi": "",
            "filename": Path(input_pdf).stem,
            "extraction_methods": ["error"],
            "error": str(e)
        }
        
        with open(output_metadata, 'w', encoding='utf-8') as f:
            json.dump(empty_metadata, f, indent=2)

if __name__ == "__main__":
    main()