# MCP tool surface

The MCP server ([mcp_server.py](../mcp_server.py)) currently exposes 23 tools. This table is generated from the docstrings in the source; when the server definition changes, regenerate with:

```bash
python3 -c "
import re
src = open('mcp_server.py').read()
for m in re.finditer(r'@mcp\.tool\(\).*?def (\w+)\([^)]*\).*?:\s*\"\"\"([^.]+)\.', src, re.DOTALL):
    name, summary = m.group(1), re.sub(r'\s+', ' ', m.group(2).strip())
    print(f'| \`{name}\` | {summary}. |')
"
```

## Discovery

| Tool | Returns |
| --- | --- |
| `list_papers` | Every paper in the corpus with bibliographic + annotation counts. |
| `get_paper` | Full metadata for one paper: title, authors, year, abstract, DOI, plus top-10 taxa and top anatomy terms. |
| `get_chunk` | One chunk's full record: text, headings, section_class, figure_refs. |
| `get_chunks_by_section` | Chunks of a paper filtered by section class. |

## Taxonomy

| Tool | Returns |
| --- | --- |
| `search_taxon` | Resolve a taxon name against the WoRMS snapshot. |
| `get_papers_for_taxon` | Papers mentioning a taxon, resolved through synonymy. |
| `get_chunks_for_taxon` | Every chunk that mentions the taxon (resolved through synonymy). |
| `get_taxon_mentions` | All text-span mentions of a taxon across the corpus. |
| `list_valid_species_under` | All currently-valid species that descend from the given taxon in the WoRMS snapshot. |
| `get_original_description` | Find the original-description paper for a taxon. |

## Bibliography + citation graph

| Tool | Returns |
| --- | --- |
| `get_bibliography` | Parsed references for one paper (from Grobid TEI). |
| `get_citation_graph` | Citation graph around a work. |
| `resolve_reference` | Resolve a free-text bibliographic reference to a work in the authority database. |
| `get_missing_references` | Works cited by corpus papers that are NOT in the corpus. |
| `get_works_by_author` | All works by an author across the full bibliographic authority database. |
| `get_papers_by_author` | Papers authored by the given surname (case-insensitive). |

## Figures

| Tool | Returns |
| --- | --- |
| `get_figures_for_taxon` | Figures from papers that mention the taxon, ranked by caption relevance. |
| `get_figures_for_anatomy` | Figures whose captions mention an anatomy term. |
| `get_figure` | One figure's full record: caption, page, bbox, image path, cross-references. |
| `list_figure_rois` | Per-panel / per-subfigure ROIs annotated on a figure. |
| `get_figure_roi_image` | Crop a panel ROI out of a figure image and return the crop's path. |

## Semantic search + translation

Requires `embed_chunks.py` to have been run (for `get_chunks_for_topic`) and `ANTHROPIC_API_KEY` (for `translate_chunk`).

| Tool | Returns |
| --- | --- |
| `get_chunks_for_topic` | Semantic search over chunks via the LanceDB vector index. |
| `translate_chunk` | Translate one chunk to the target language (default English), via the Anthropic Claude API. |
