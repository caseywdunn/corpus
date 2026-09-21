# MCP tool surface

The MCP server exposes 38 `@mcp.tool()`-decorated functions, split across `mcpsrv/tools/{papers,taxonomy,bibliography,figures,chunks,lexicon,profiles}.py`. `corpus serve` is the user-facing entry point; it wraps `mcpsrv.main`.

This surface is frozen as of 1.0. What that commits us to — additive vs. breaking, and how anything gets removed — is [API_STABILITY.md](API_STABILITY.md).

All explicit list parameters accept at most 500 items, 4,096 characters per
item and 65,536 characters in total per list. Larger requests return
`invalid_argument` in the tool's normal error shape before corpus lookup;
split them into batches. Lists are never silently truncated. Omitted and empty
lists retain each tool's documented meaning (including whole-paper selection
by `get_chunks(chunk_ids=None)`); this is an input budget, not new pagination.

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
| `get_papers` | Full metadata for one or many papers (title, authors, year, abstract, DOI, top taxa and lexicon terms) with optional field whitelist. Output is in input order. Synthesizes `first_author` so the caller doesn't walk `authors[0].surname` themselves. Pass a single-element `hashes` list for the one-paper case. |
| `get_chunks` | Chunk fetch for one paper — drill-down pair to every `*_dossier` tool. Pass `chunk_ids=[...]` for a subset; `with_text=False` omits chunk prose. Carries materialized treatment/source context, with optional exact `treatment_name` and `section_type` filters (#319). |
| `get_chunks_by_section` | Chunks of a paper filtered by existing `section_class`, optional materialized `section_type` (such as `diagnosis`), and exact resolved `treatment_name`; filters intersect before `limit` is applied (#319). |


### Bounded diagnosis retrieval

Use `get_chunks_by_section(paper_hash="…", section_type="diagnosis", limit=20,
with_text=False)` to inspect diagnosis IDs and stored treatment names. Add
`treatment_name="Genus species"` using the exact returned resolved name to
restrict the result, then fetch the selected IDs with
`get_chunks(paper_hash="…", chunk_ids=["…"])`. `treatment_name` is an exact,
case-sensitive comparison of the build's resolved treatment name; it does not
resolve taxonomy synonyms or search for names in prose. The existing
`section_class="description"` remains valid and includes materialized diagnoses.
An explicit genus or subgenus treatment carries `rank="genus"` or
`rank="subgenus"`; its name is not an inferred species assignment. A later
evidenced species heading replaces that context. Missing species evidence does
not justify choosing a species from nearby prose.

Both tools expose `treatment_context`, `section_type` and `source_items` (source
item, page, box and character span). Bounded producer observations appear as
`text_integrity`, `tables` and `key_branches`, including partial-table and
unverified-spelling flags. Source prose fields (including original/replacement
text, headings and table-cell text) use `{preview, char_count,
preview_charspan, truncated}` objects: at most 32 characters with
`with_text=False`, or 96 with text. Preview offsets refer to that stored source
field; repair/source spans retain their original item-relative coordinates.
Full evidence remains in build artifacts. The server does not recompute it.

Newly chunked text keys also expose `key_branches[].chunk_scope`:
`coverage` is `complete`, `partial` or `unknown` when exact source-text
alignment is unavailable; unknown coverage leaves completeness flags null.
`complete_branch` says whether this chunk contains the entire associated lead
and destination; `continuation` marks partial context. `destination_in_text`
distinguishes an endpoint visible in this chunk from one carried only by its
source association. Zero-based `fragment_index`, `fragment_count` and adjacent
`previous_chunk_id`/`next_chunk_id` support bounded continuation: fetch one
adjacent ID with `get_chunks` and follow the same `item_ref`. The source
association's status, including unverified spelling, remains separate from
fragment completeness. Legacy chunks lacking this scope have unavailable
split context until rechunked; the server does not infer it. These fields use
the same context projection limits as other producer evidence.

`context_projection` reports record counts (`available`, `returned`, `truncated`)
for each evidence array and a combined truncation flag. At most six records per
array and six entries per nested list are projected, with nested list counts
reported as `<field>_scope`. New context fields together are capped at 8 KiB per
row, measured as compact UTF-8 JSON. Nonempty optional evidence arrays share a
64 KiB budget across the response. These limits apply to the new context fields,
not existing chunk prose, headings, or the whole response: required context
status/count notices still scale with the number of requested rows. Requested
rows are retained when evidence is omitted, with explicit counts.

`status="resolved"` carries the enclosing treatment name and its heading
evidence; `status="unknown"` means the build examined context but could not
assign it. Captions or comparative mentions do not establish a treatment.
Literal-name routes such as `get_chunks_for_taxon` retain their mention semantics;
a diagnosis may be retrieved through its treatment even when it never repeats
the species name.

Older artifacts remain readable through ordinary chunk and `section_class`
retrieval, reporting `treatment_context.status="unavailable"`. Requests using
`treatment_name` or `section_type` on artifacts without the materialized context
policy return `code="rebuild_required"` with regeneration guidance. Rebuild
extraction/chunks and their annotation, embedding and bundle descendants to
provide the new context. The tools do not repair an old bundle while serving.

## Taxonomy

| Tool | Returns |
| --- | --- |
| `search_taxon` | Resolve a taxon name against the configured Darwin Core taxonomy snapshot. Pass `parent_chain=True` (#88) to add the accepted taxon's ancestry (immediate parent → root, each `{taxon_id, scientific_name, rank}`). |
| `get_papers_for_taxon` | Papers mentioning a taxon, resolved through synonymy. |
| `get_chunks_for_taxon` | Every chunk that mentions the taxon (resolved through synonymy). |
| `get_taxon_mentions` | All text-span mentions of a taxon across the corpus, with surrounding context. |
| `list_valid_species_under` | Accepted species/subspecies descending from a taxon, resolving synonyms in the configured snapshot. Optional `limit`/`offset` pagination preserves the species-list payload and reports continuation in MCP `_meta.pagination`; oversized unpaged calls fail explicitly. |
| `get_papers_by_author` | Papers authored by the given surname (case-insensitive). |
| `get_taxon_dossier` | One-call comprehensive view of a taxon across the corpus: metadata, papers (sorted by mention count), chunk_index (IDs only — pair with `get_chunks`), figure_index, top lexicon terms per category, cooccurring taxa. Supersedes the `search_taxon` + `get_papers_for_taxon` + `get_papers` + N× `get_chunks_for_taxon` + `get_figures_for_taxon` chain (~45 round-trips → 1). `include=[...]` trims sections. Lexicon/cooccurrence totals cover the capped paper selection; `aggregate_scope` reports selected/available paper counts and the deterministic selection policy, even when the paper list is omitted (#318). |
| `get_taxon_lexicon_slice` | Lexicon coverage for one taxon under one category, joined at the chunk level. Same-chunk co-occurrence — tighter than `corpus_summary` / dossier rollups: a term only counts when it appears in a chunk where the taxon is also mentioned. Returns `{term, n_chunks, n_papers, paper_examples}` per term. Category-agnostic; unknown category returns the available list. |
| `get_taxon_subtree_dossier` | Walk a clade's accepted species/subspecies via the DwC `parent_name_usage_id` tree; for each species with corpus coverage, return a capsule (paper count, mention count, authorship). Plus a deduplicated aggregate paper list across the subtree with `n_species_covered` per paper. Supersedes the p07 monographic pattern of `list_valid_species_under` + N× `get_papers_for_taxon`. |

### Species-list pagination

`list_valid_species_under(parent_taxon_name, limit=None, offset=0)` preserves
the existing list of species objects. With no limit, a complete list is
returned if its MCP result fits 256 KiB. **Safety change:** a larger unpaged
request now returns a bounded error row (`code="invalid_argument"`,
`reason="pagination_required"`, `suggested_limit=100`, `next_offset=0`),
rather than emitting an oversized event or silently returning a prefix.

Pass a positive `limit` to opt into pagination. The shared cap is 500;
`limit=0` is invalid, and a nonzero `offset` requires an explicit limit.
Ordering is by scientific name, then taxon ID, within the immutable snapshot.
Synonyms still resolve to their accepted parent; only accepted species and
subspecies below that parent are returned, with existing fields unchanged.

The MCP result retains its per-row text blocks and its
`structuredContent.result` list. One additive metadata object,
`_meta.pagination`, contains `offset`, effective `limit`, `returned`,
`total_available`, `next_offset`, `truncated`, `truncated_reason`,
`result_bytes` and `max_result_bytes`. For example, a client requests
`limit=3, offset=0`, reads the same species list, then requests the metadata's
`next_offset` until it is null. Read metadata even when fewer than `limit`
rows arrive: the byte budget can shorten a page (`response_bytes` reason).
Empty successful pages also carry metadata. The budget measures UTF-8 JSON
for the entire MCP result, including both content representations and
pagination metadata; the small JSON-RPC/SSE framing is additional.

Rows and taxon names are never clipped. An individual row exceeding the
budget returns an explicit `unavailable` error with
`reason="row_exceeds_response_budget"` and `blocked_offset`; pagination
reports zero returned and no next offset. This is a failed page, not the end
of a complete enumeration. That row requires inspection of the source
snapshot; the server does not silently skip it or repeat a nonadvancing
continuation. Ordinary configuration/argument failures retain error rows.

## Bibliography + citation graph

| Tool | Returns |
| --- | --- |
| `get_bibliography` | Parsed references for one paper (from Grobid TEI). |
| `get_intext_citations` | In-text `<ref type="bibr">` markers for one paper, with deduplicated paragraph excerpts and section context. |
| `get_excerpts_citing` | Cross-corpus citation markers with surface text, section, surrounding paragraph and source marker/paragraph indices. Rebuilt records expose author/year spans, target validation, and whether citation text came from PDF coordinates or TEI; `get_intext_citations` also retains raw TEI marker observations (#309, #317). `limit` (default 50) + `offset` paginate in paper-hash/source-marker order; follow `next_offset`. `excerpts_available` / `excerpts_returned` count **markers**, so multiple rows can share a paragraph. Reports `truncated`, `truncated_reason`, and exact `response_bytes`; JSON is bounded by `CORPUS_EXCERPTS_MAX_BYTES` (default 128 KiB). A single oversized marker retains its identity with explicitly marked text previews (`truncated_fields` records original character counts); retrieve its source via `get_intext_citations`. `limit=0` returns counts only (#325). |
| `get_citation_graph` | Citation graph around a work or paper (`citing` / `cited_by` / `both`; other values return `invalid_argument`). Bounded breadth (#87): `max_edges_per_node` (per-node fan-out, survivors ranked by `cited_by_count`) + a shared `max_total_edges` cap across both directions (`citing` first; includes budgets 0 and 1). Reports `edges_available` / `edges_returned` / `response_bytes` alongside `truncated`, because that flag says only whether *this tool* cut edges — it reads `false` on a 145 kB payload a client may not deliver (#166). `response_bytes` is also a real ceiling (`CORPUS_CITATION_GRAPH_MAX_BYTES`, default 256 kB). Generous defaults. |
| `resolve_reference` | Resolve an author/year reference to a work in the authority database. Shares query parsing with `format_citations`, including `et al.`, `and`, `&`, surname particles and diacritics (#310). A failed query does not establish that a publication is absent; refine the query or use a known identifier. Title-only queries remain unsupported. |
| `format_citations` | Fully-assembled citation strings for works in the authority DB — the route for every citation an LLM client emits; never recombine fields client-side (#88). Pass one of `queries` / `work_ids` / `paper_hashes` (a list); returns `{style, count, citations[]}` in input order, each entry a citation payload (`work_id`, `formatted`, `inline`, provenance tier `bib` / `grobid_reconciled` / `unresolved`, verbatim warning footnote) or a per-item error. Batch a whole reference list into one call. |
| `get_missing_references` | Candidate works cited by corpus papers that are not mapped to an in-corpus work, ranked by citation count. This is an acquisition lead, not proof of absence: inspect build-time evidence with `tools/qc/reference_reconciliation.py` before curating the library (#155). Best-effort by construction (#155): a work can still appear as several rows where its citation strings did not reconcile, so verify against `resolve_reference` before treating a row as a gap. Rows with neither title nor year are withheld — parse debris, not leads (477 of 6,953 on the reference corpus). |
| `get_works_by_author` | All works by an author across the full bibliographic authority database (corpus papers + cited references + taxonomic-authority stubs). |
| `get_original_description` | Find the original-description paper for a taxon. |

## Figures

Every response that carries a figure caption or caption-derived ROI also
carries `caption_status` (`bound` / `uncertain` / `unbound`),
`caption_confidence` (`high` / `medium` / `low` / null), and
`caption_page_distance`, plus `caption_kind` (`prose_caption` / `bare_label` /
`unlabelled_caption`; `unknown` on a legacy record). These are build-time
association facts. The server
only normalizes the same facts for bundles made before the fields existed; it
does not recompute caption binding at query time.

| Tool | Returns |
| --- | --- |
| `get_figures_for_taxon` | Caption-name matches rank before paper-only associations; paper mention counts order each group, with paper hash and figure ID breaking ties. `caption_has_taxon: false` identifies paper-only evidence. The legacy numeric `score` can be higher for a nonmatch and is not the primary rank; use `caption_only=True` or `caption_has_taxon` to select caption matches. Neither evidence type verifies what the image depicts. `caption_text` is a preview by default (#85); `full_caption=True` for the verbatim caption. |
| `get_figures_for_lexicon_term` | Figures whose captions mention a lexicon term. The query is resolved to its canonical term and matched against **every surface form observed in the corpus** (#143), so `wing` finds captions saying `ala` and vice versa — including across languages, since a Russian paper's `нектофор` maps to `nectophore`. Rows carry `matched_surfaces` + `canonical`; an unrecognized term degrades to a literal search with `resolved: false`. `caption_text` previewed by default (#85). |
| `get_figure_dossier_for_taxon` | Figures linked to a taxon, using the same caption-first rank and `caption_has_taxon` evidence distinction as `get_figures_for_taxon`; paper association does not establish depiction. Each has `linked_chunks` (chunk IDs that reference the figure via `chunks.json:figure_refs`) + summarized ROIs. Single call replaces `get_figures_for_taxon` + per-figure `list_figure_rois` + cross-ref against `get_chunks_for_taxon`. |
| `get_figure_dossier_for_term` | Same shape, for figures whose captions match a lexicon term. Synonym-aware on the same basis as `get_figures_for_lexicon_term` (#143); the response reports `canonical`, `resolved`, and `surfaces_searched`. |
| `get_figure` | One figure's full record: caption, page, bbox, image path, cross-references, plus `license` / `license_url` / `attribution` from the parent work. The publication-clearance *determination* is included only under a strict `profile=`, or on `include_licensing=True` (#154) — under the permissive default the server has already authorized the figure, so it isn't shipped. **`license: null` is normal and does not mean "no licensing information"** — most historical works carry no explicit license string, and their clearance is *derived from publication year* against `licensing.pd_cutoff_years`. The derived answer lives in `publishable` / `license_source` (`age_based_pd`) / `clearance_state`; ask for `include_licensing=True` to see it. A pre-cutoff work reads `license: null` while being fully cleared as public domain. |
| `get_figure_image` | A figure (or panel crop) returned as inline PNG bytes. Figure-licensing gate keyed to the per-call `profile=` (#101): a strict profile (manuscript/presentation) refuses a figure whose publication clearance can't be established; the default `report` allows it. Refusals retain MCP `isError: true` and provide structured `error`/`code` fields; licensing refusals also identify `profile`, `publication_clearance` and `license_source`, matching URL delivery (#327). |
| `get_figure_url` | A five-minute signed HTTP download URL, scoped to this figure, optional panel and resolved output profile. Fetch directly with `curl -o`; the retained `auth_header` field is `null`, never the shared MCP bearer token. Licensing is checked again at download. Links also expire on restart. Reverse proxies must forward `/figures/` and configure a client-reachable public base (see [DEPLOY.md](../DEPLOY.md)). License/attribution fields retain their existing shape. |
| `list_figure_rois` | Per-panel / per-figure ROIs annotated on an image. Caption-derived A/B/C panels and numeric figures on a shared historical plate are separate fields and carry `pass3_target_kind`; plate targets also report the independent pre-expansion `missing_figures` cross-check. |
| `get_figure_roi_image` | Crop a lettered-panel or numbered-figure ROI and return its logical cache path. Crops live outside the immutable bundle; retrieve bytes with `get_figure_image` or `get_figure_url`, not by joining this path to the bundle directory. Honors the per-call `profile=` licensing gate, including whole-figure fallback when no pixel ROI was detected. |

## Semantic search

Requires the embedding stage to have been run (`corpus run --only embed`), for `get_chunks_for_topic`.

| Tool | Returns |
| --- | --- |
| `get_chunks_for_topic` | Semantic search over chunks via the LanceDB vector index. Pass `with_text=False` for a metadata-only scan (#82): ~80 chars/row vs ~600 with full text, then drill down with `get_chunks(paper_hash, chunk_ids=[...])`. Pass `with_cites=True` (#88) to attach `cited_work_ids` (the chunk's parent paper's in-text citation targets) — feed to `format_citations`. |

## Output profiles

Output type is a per-call client/session property (#101): pass `profile=` (`report` / `manuscript` / `presentation`) to the gated figure tools to select the figure-licensing policy. These tools expose the vocabulary.

| Tool | Returns |
| --- | --- |
| `list_output_profiles` | The built-in output profiles and their policies (`figure_licensing`, `require_attribution`, `citation_provenance`, `excerpt_max_words`) plus the server's fallback `profile`. |
| `get_active_profile` | The server's fallback output profile — the one applied to calls that omit `profile=`. Informational; the authoritative selection is per-call. |

## Lexicon

Category-agnostic — every tool takes `category` as an argument and returns `{"error": "unknown_category", "available": [...]}` if the bundle doesn't declare it.

| Tool | Returns |
| --- | --- |
| `lexicon_matrix` | Lexicon-coverage view for one category. **Default `detail=False`** returns compact per-term totals (`term_totals[]` with `total_mentions` + `papers_with_mentions`) over the selected papers. **`detail=True`** returns the full paper × term mention-count grid (`rows[]`), which is O(papers × terms) and was a multi-MB runaway, hence opt-in (#88). Caller-controlled columns (`terms=`, preserving caller order) or top-N by total mention count **within the selected paper set**, with alphabetical ties; caller-controlled paper set (`paper_hashes=`) or all papers, optionally year-filtered. The grid is bounded by `CORPUS_LEXICON_MATRIX_MAX_BYTES` (default 128 kB) and reports `rows_available` / `rows_returned` / `response_bytes` / `truncated`; #88 made it opt-in but never bounded it, and it runs 382 kB over 1,775 rows. Narrow with `paper_hashes` or `terms` (#83). |
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
- **Tool catalog** (38 tools × ~100 tokens after the #81 docstring trim): ~4 k tokens.

Together ~5 k tokens get cached across all turns of a session — a non-trivial saving on conversations that fan out into many tool-use rounds. Cache lives ~5 minutes by default; subsequent sessions against the same build hit the same cache lines.

The catalog is stable across builds *unless* you bump `pipeline/version.py` or a docstring changes. PR diffs that touch tool docstrings should expect a cache miss on the next session — a routine cost.
