# MCP tool surface

The MCP server exposes 40 `@mcp.tool()`-decorated functions, split across `mcpsrv/tools/{papers,taxonomy,bibliography,figures,chunks,lexicon}.py`. The top-level [mcp_server.py](../mcp_server.py) is a thin shim into `mcpsrv.main`.

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
| `corpus_summary` | One-call orientation: paper counts by decade, lexicon coverage per category (top terms), top taxa, figure totals, bundle identity. Server-side join over the in-memory indexes — fixed-shape payload (~2–5 k tokens) regardless of corpus size. Caps on `top_taxa` + `top_terms_per_category`. |
| `list_papers` | Every paper in the corpus with bibliographic + annotation counts. Optional `year_from` / `year_to` filters. |
| `get_paper` | Full metadata for one paper: title, authors, year, abstract, DOI, plus top taxa and anatomy terms. |
| `get_chunk` | One chunk's full record: text, headings, section_class, figure_refs. |
| `get_papers` | Batched `get_paper` for many hashes at once with optional field whitelist. Output is in input order. Synthesizes `first_author` so the caller doesn't walk `authors[0].surname` themselves. |
| `get_chunks` | Batched chunk fetch for one paper — drill-down pair to every `*_dossier` tool. Pass `chunk_ids=[...]` for a subset; `with_text=False` emits metadata-only (~80 chars/chunk). |
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
| `get_taxon_dossier` | One-call comprehensive view of a taxon across the corpus: metadata, papers (sorted by mention count), chunk_index (IDs only — pair with `get_chunks`), figure_index, top lexicon terms per category, cooccurring taxa. Supersedes the `search_taxon` + `get_papers_for_taxon` + N× `get_paper` + N× `get_chunks_for_taxon` + `get_figures_for_taxon` chain (~45 round-trips → 1). `include=[...]` trims sections. |
| `get_taxon_lexicon_slice` | Lexicon coverage for one taxon under one category, joined at the chunk level. Same-chunk co-occurrence — tighter than `corpus_summary` / dossier rollups: a term only counts when it appears in a chunk where the taxon is also mentioned. Returns `{term, n_chunks, n_papers, paper_examples}` per term. Category-agnostic; unknown category returns the available list. |
| `get_taxon_subtree_dossier` | Walk a clade's accepted species/subspecies via the DwC `parent_name_usage_id` tree; for each species with corpus coverage, return a capsule (paper count, mention count, authorship). Plus a deduplicated aggregate paper list across the subtree with `n_species_covered` per paper. Supersedes the p07 monographic pattern of `list_valid_species_under` + N× `get_papers_for_taxon`. |

## Bibliography + citation graph

| Tool | Returns |
| --- | --- |
| `get_bibliography` | Parsed references for one paper (from Grobid TEI). |
| `get_intext_citations` | In-text `<ref type="bibr">` markers for one paper, with deduplicated paragraph excerpts and section context. |
| `get_excerpts_citing` | Cross-corpus: every passage citing a given work — surface text, section, and the surrounding paragraph. |
| `get_citation_graph` | Citation graph around a work or paper (in / out / both). |
| `resolve_reference` | Resolve a free-text bibliographic reference to a work in the authority database. |
| `format_citation` | Fully-assembled citation string for a work in the authority DB, plus provenance tier (`bib` / `grobid_reconciled` / `unresolved`) and verbatim warning footnote. The route for every citation an LLM client emits — never recombine fields client-side. |
| `format_citations` | Batched `format_citation` (#88): pass one of `queries` / `work_ids` / `paper_hashes` (a list); returns `citations[]` in input order, each the `format_citation` payload or a per-item error. Prefer this over N single calls when emitting a reference list. |
| `get_missing_references` | Works cited by corpus papers that are NOT in the corpus. |
| `get_works_by_author` | All works by an author across the full bibliographic authority database (corpus papers + cited references + taxonomic-authority stubs). |
| `get_original_description` | Find the original-description paper for a taxon. |

## Figures

| Tool | Returns |
| --- | --- |
| `get_figures_for_taxon` | Figures from papers that mention the taxon, ranked by caption relevance. |
| `get_figures_for_lexicon_term` | Figures whose captions mention a term from one lexicon category (anatomy, biogeography, …). |
| `get_figure_dossier_for_taxon` | Figures linked to a taxon, each with `linked_chunks` (chunk IDs that reference the figure via `chunks.json:figure_refs`) + summarized ROIs. Single call replaces `get_figures_for_taxon` + per-figure `list_figure_rois` + cross-ref against `get_chunks_for_taxon`. |
| `get_figure_dossier_for_term` | Same shape, for figures whose captions match a lexicon term. Category-agnostic. |
| `get_figure` | One figure's full record: caption, page, bbox, image path, cross-references. |
| `get_figure_image` | A figure (or panel crop) returned as inline PNG bytes. |
| `get_figure_url` | A bearer-gated HTTP URL (plus `auth_header` + license fields) the caller can `curl -o` to land the figure PNG on disk without loading its bytes into the model's context — for pandoc / LaTeX / PDF assembly. |
| `list_figure_rois` | Per-panel / per-subfigure ROIs annotated on a figure. |
| `get_figure_roi_image` | Crop a panel ROI out of a figure image and return the crop's path. |

## Semantic search + translation

Requires `embed_chunks.py` to have been run (for `get_chunks_for_topic`) and `ANTHROPIC_API_KEY` (for `translate_chunk`).

| Tool | Returns |
| --- | --- |
| `get_chunks_for_topic` | Semantic search over chunks via the LanceDB vector index. Pass `with_text=False` for a metadata-only scan (#82): ~80 chars/row vs ~600 with full text, then drill down with `get_chunks(paper_hash, chunk_ids=[...])`. Pass `with_cites=True` (#88) to attach `cited_work_ids` (the chunk's parent paper's in-text citation targets) — feed to `format_citations`. |
| `translate_chunk` | Translate one chunk to the target language (default English), via the Anthropic Claude API. |

## Lexicon

Category-agnostic — every tool takes `category` as an argument and returns `{"error": "unknown_category", "available": [...]}` if the bundle doesn't declare it.

| Tool | Returns |
| --- | --- |
| `lexicon_matrix` | Lexicon-coverage view for one category. **Default `detail=False`** returns compact per-term totals (`term_totals[]` with `total_mentions` + `papers_with_mentions`) over the selected papers. **`detail=True`** returns the full paper × term mention-count grid (`rows[]`), which is O(papers × terms) and was a multi-MB runaway, hence opt-in (#88). Caller-controlled columns (`terms=`) or top-N by total mention count; caller-controlled paper set (`paper_hashes=`) or all papers, optionally year-filtered. |
| `get_lexicon_term_dossier` | Per-term cross-corpus view: rollup counts, top papers by mention count, chunk examples (IDs only — pair with `get_chunks`), term description. |

## Prompt-cache integration (Anthropic API clients)

The tool catalog + system prompt this server emits are static across a corpus build — #76 principle 4 (deterministic ordering and shape) means cache breakpoints stay valid until a new bundle ships. A client consuming the surface via Anthropic's Messages API directly (the `tests/test_prompt_quality.py` harness in #79 is the reference implementation; Claude Desktop / Claude Code handle caching themselves) cuts cached-token pricing ~10× by placing `cache_control` breakpoints after the static prefix:

```python
client.messages.create(
    model="claude-sonnet-4-6",
    system=[
        {"type": "text", "text": default_instructions_md,
         "cache_control": {"type": "ephemeral"}},   # ← end of system prompt
    ],
    tools=[
        *tool_definitions[:-1],
        {**tool_definitions[-1],
         "cache_control": {"type": "ephemeral"}},   # ← end of tool catalog
    ],
    messages=[...],
)
```

Two breakpoints land below the ~5 k-token cache-eligibility floor on typical bundles:

- **System prompt** (`mcpsrv/default_instructions.md` concatenated with any corpuscle-specific `instructions.md`): ~500–1500 tokens.
- **Tool catalog** (40 tools × ~100 tokens after the #81 docstring trim): ~4 k tokens.

Together ~5 k tokens get cached across all turns of a session — a non-trivial saving on conversations that fan out into many tool-use rounds. Cache lives ~5 minutes by default; subsequent sessions against the same build hit the same cache lines.

The catalog is stable across builds *unless* you bump `pipeline/version.py` or a docstring changes. PR diffs that touch tool docstrings should expect a cache miss on the next session — a routine cost.
