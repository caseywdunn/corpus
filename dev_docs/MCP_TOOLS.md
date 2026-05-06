# MCP tool surface

The MCP server exposes 27 `@mcp.tool()`-decorated functions, split across `mcpsrv/tools/{papers,taxonomy,bibliography,figures,chunks}.py`. The top-level [mcp_server.py](../mcp_server.py) is a thin shim into `mcpsrv.main`.

This table is generated from the docstrings in the source; when the server definition changes, regenerate with:

```bash
python3 -c "
import re, pathlib
for f in sorted(pathlib.Path('mcpsrv/tools').glob('*.py')):
    if f.name.startswith('_'): continue
    src = f.read_text()
    for m in re.finditer(r'@mcp\\.tool\\(\\)\\s+(?:async\\s+)?def (\\w+)\\([^)]*\\)[^:]*:\\s*\"\"\"(.*?)\"\"\"', src, re.DOTALL):
        name, doc = m.group(1), m.group(2).strip()
        first = re.sub(r'\\s+', ' ', doc.split('.')[0].strip())
        print(f'| \`{name}\` | {first}. |')
"
```

## Discovery + bundle metadata

| Tool | Returns |
| --- | --- |
| `bundle_info` | Bundle version, server name + version, paper / chunk / figure counts, embedding model, pipeline git SHA. Lets clients detect stale endpoints and cite a corpus version in downstream work. |
| `list_papers` | Every paper in the corpus with bibliographic + annotation counts. Optional `year_from` / `year_to` filters. |
| `get_paper` | Full metadata for one paper: title, authors, year, abstract, DOI, plus top taxa and anatomy terms. |
| `get_chunk` | One chunk's full record: text, headings, section_class, figure_refs. |
| `get_chunks_by_section` | Chunks of a paper filtered by section class. |

## Taxonomy

| Tool | Returns |
| --- | --- |
| `search_taxon` | Resolve a taxon name against the configured Darwin Core taxonomy snapshot. |
| `get_papers_for_taxon` | Papers mentioning a taxon, resolved through synonymy. |
| `get_chunks_for_taxon` | Every chunk that mentions the taxon (resolved through synonymy). |
| `get_taxon_mentions` | All text-span mentions of a taxon across the corpus, with surrounding context. |
| `list_valid_species_under` | All currently-valid species descending from the given taxon in the configured taxonomy snapshot. |
| `get_papers_by_author` | Papers authored by the given surname (case-insensitive). |

## Bibliography + citation graph

| Tool | Returns |
| --- | --- |
| `get_bibliography` | Parsed references for one paper (from Grobid TEI). |
| `get_intext_citations` | In-text `<ref type="bibr">` markers for one paper, with deduplicated paragraph excerpts and section context. |
| `get_excerpts_citing` | Cross-corpus: every passage citing a given work — surface text, section, and the surrounding paragraph. |
| `get_citation_graph` | Citation graph around a work or paper (in / out / both). |
| `resolve_reference` | Resolve a free-text bibliographic reference to a work in the authority database. |
| `get_missing_references` | Works cited by corpus papers that are NOT in the corpus. |
| `get_works_by_author` | All works by an author across the full bibliographic authority database (corpus papers + cited references + taxonomic-authority stubs). |
| `get_original_description` | Find the original-description paper for a taxon. |

## Figures

| Tool | Returns |
| --- | --- |
| `get_figures_for_taxon` | Figures from papers that mention the taxon, ranked by caption relevance. |
| `get_figures_for_anatomy` | Figures whose captions mention an anatomy term. |
| `get_figure` | One figure's full record: caption, page, bbox, image path, cross-references. |
| `get_figure_image` | A figure (or panel crop) returned as inline PNG bytes. |
| `list_figure_rois` | Per-panel / per-subfigure ROIs annotated on a figure. |
| `get_figure_roi_image` | Crop a panel ROI out of a figure image and return the crop's path. |

## Semantic search + translation

Requires `embed_chunks.py` to have been run (for `get_chunks_for_topic`) and `ANTHROPIC_API_KEY` (for `translate_chunk`).

| Tool | Returns |
| --- | --- |
| `get_chunks_for_topic` | Semantic search over chunks via the LanceDB vector index. |
| `translate_chunk` | Translate one chunk to the target language (default English), via the Anthropic Claude API. |
