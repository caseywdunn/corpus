# PLAN.md — Corpus pipeline

A roadmap toward a production-quality corpus pipeline. The pipeline is taxon-agnostic: domain specifics (taxonomy snapshot, anatomy lexicon, target queries, instructions to MCP clients) live in the per-corpus *corpuscle* directory, and only the pipeline code lives in this repo. The first deployment targets a corpus of ~2000 papers spanning late-18th-century to modern scientific literature in multiple languages (English, German, French, Russian, and others).

## 1. Assessment of current repo

The initial-prototype bug-fix punch list (eleven items: stub `extract_metadata`, naive chunker, unloaded config, silent zero-vector embeddings, unlinked figure captions, serial processing, 8-char hash collisions, English-only OCR, dual Snakemake/Python pipelines, missing PyMuPDF figure metadata, missing per-paper logs) has been fully delivered. See [`CHANGELOG.md`](CHANGELOG.md) for the per-release history. The "what works well" items below remain the foundation:

- **Hash-based content addressing** — SHA-256 prefix directories give idempotent resume, natural dedup across messy input trees, and survive rename/reorg of source collections.
- **Two-stage split** (CPU extraction vs. GPU embeddings) — embeddings are the resource-bound step that benefits from rerunnability and parallelism.
- **Figure QC visualizations** (red word boxes + yellow/orange figure bboxes per page) — the intermediate artifact that catches extraction problems early.
- **PyMuPDF fallback for figure extraction** — docling doesn't always deliver; the fallback now classifies by size threshold so MCP-visible types are right.

## 2. Is this the right toolchain?

A candid reassessment for the stated goals.

### Docling
Works well on modern born-digital papers. Struggles on: plate-heavy 19th-century monographs (figures separate from captions), scanned Fraktur German, two-column taxonomy journals with irregular layouts. **Keep for modern papers; don't assume it will carry the historical tail.**

### ocrmypdf + Tesseract
Usable if given the right language packs. Needs: `tesseract-lang` (deu, fra, rus, lat), and specifically `deu_latf` / Fraktur models for 19th-c. German. For very old or damaged scans, Tesseract is weak — **Kraken** (historical documents, fraktur, Cyrillic) or a vision LLM may do better.

### Grobid
Still the best general-purpose scientific-PDF metadata/reference extractor. Runs as a Docker (or Singularity on Bouchet) service, called per PDF. Extracts references and section structure, foundational for the citation-analysis goal.

### Alternatives worth piloting
- **Marker** (surya + texify) — often better than docling on mixed-quality PDFs, native markdown output with structure.
- **Nougat** (Meta) — good on modern arxiv-like PDFs but trained on that distribution; poor on historical.
- **Claude / GPT-4o vision** on page images — expensive per page, but dramatically better on hard pages (fraktur, plates, hand-lettered captions). Reasonable as a **fallback for pages that fail structural parsing**, not as the default.
- **BHL (Biodiversity Heritage Library) API** — a sizable fraction of historical zoological literature is already in BHL with OCR and bibliographic metadata. **Query BHL before OCRing a scan yourself.** This alone could save hundreds of papers' worth of OCR effort and provide cleaner ground truth.
- **CrossRef / OpenAlex** for bibliographic metadata by DOI. Faster and more accurate than Grobid when a DOI is available.

### Recommended toolchain
- **Modern born-digital** → docling (primary) + Grobid (metadata + references) + OpenAlex/CrossRef (metadata override when DOI known).
- **Modern scans** → ocrmypdf (multi-lang) → docling + Grobid.
- **Historical (pre-1950)** → BHL lookup first; if not available, Kraken OCR (or vision LLM for the hardest) + custom layout handling for plate/caption association.
- **Embeddings** → switch to a local open-weights model run on Bouchet GPUs. See §7.

## 3. Re-architecture

Current shape is fine; the key moves are to (a) separate extraction concerns more crisply, (b) make figure+caption a first-class joint object, (c) layer an annotation stage on top for the domain-specific goals, and (d) expose the whole corpus as an MCP server over open-format indices so any LLM client — and any non-LLM tool — can interrogate it.

### Proposed pipeline stages

```
Stage 0  Acquisition & identification
         └─ hash, DOI lookup (CrossRef/OpenAlex), BHL lookup, language detect

Stage 1  Normalization
         └─ scan classification → OCR (lang-aware) → processed.pdf

Stage 2  Structural extraction
         ├─ text + sections + references        (docling or Grobid)
         ├─ figures + bboxes + captions joined  (docling pictures + caption refs)
         └─ QC artifacts (per-page visualizations, figure contact sheets)

Stage 3  Chunking
         └─ HybridChunker (respects headings/sections), tiktoken-aware

Stage 4  Annotation (new)
         ├─ taxonomic names        (gnfinder + DwC verification, e.g. WoRMS/GBIF)
         ├─ geographic mentions    (spaCy / GLiNER + GeoNames)
         ├─ anatomy terms          (corpus-specific domain lexicon)
         └─ figure classification  (what anatomical structure is depicted? — vision LLM)

Stage 5  Indexing (all outputs in open, file-based formats)
         ├─ vector index (LanceDB)               — semantic / topic retrieval
         ├─ taxon mention index (parquet)        — accepted-name resolved, chunk offsets
         ├─ DwC taxonomy backbone (sqlite)       — genera / species / synonymy / authorities
         ├─ bibliographic index (parquet)        — authors, year, title, journal, refs (Grobid + CSL-JSON export)
         ├─ reference graph (parquet)            — citations between papers
         ├─ section-type index (parquet)         — chunk → section class (from Grobid TEI)
         ├─ anatomy-term index (parquet)         — domain anatomy NER offsets
         ├─ geographic index (parquet)           — location mentions + coordinates
         ├─ figure index (parquet + PNGs)        — captions, taxon tags, anatomy tags
         └─ trait / character matrix (parquet)   — species × character (deferred; see §8 Q3)

Stage 6  MCP tools layer
         └─ MCP server exposing tools over the indices (see "Endpoint" below)
```

### Corpus-level inputs (shared across all papers)

Not every input comes from the PDF directory. The pipeline also depends on a small set of **corpus-level reference data** that must be materialized alongside the per-paper artifacts:

- **DwC taxonomy snapshot** — the authoritative backbone of accepted names, synonymies, and authority strings for the focal taxonomic group. Sourced from any Darwin Core archive (e.g. WoRMS for marine taxa, GBIF for general biodiversity) and ingested as `taxonomy.sqlite` once, refreshed periodically. Consumed by every taxon-keyed query (§8) and required by the taxon-mention resolver.
- **Domain anatomy lexicon** — curated terminology specific to the corpus's focal group, used by the anatomy-term NER pass and the figure-anatomy classifier. Maintained as a YAML file in the corpuscle directory.
- **(optional) BHL bibliographic index** — pre-queried BHL metadata/OCR for the historical tail; see §4 BHL lookup.

These are distinct from per-paper stages because they apply across the whole corpus and are not content-addressed by PDF hash.

### Figure+caption as a first-class object

Replace the current loose `figures.json` with:

```json
{
  "figure_id": "docling_3",
  "page": 14,
  "bbox": [x0, y0, x1, y1],
  "bbox_coord_system": "pdf_pts_bottom_left",
  "image_file": "figures/figure_3.png",
  "caption_text": "Fig. 3. <example caption text> ...",
  "caption_page": 14,
  "caption_bbox": [...],
  "caption_source": "docling_caption_link" | "heuristic_proximity" | "manual",
  "figure_number_parsed": "3",
  "referenced_in_chunks": ["chunk_7", "chunk_12"]
}
```

This requires: (a) preserving bbox in docling extraction (already have it for visualizations — just emit it), (b) pulling caption text + bbox from docling's `TextItem` children where available, (c) a heuristic fallback (nearest caption-style text block on same or facing page) for cases where docling doesn't link them, and (d) a linker pass that resolves "Fig. 3" mentions in body text to figure objects.

### Endpoint: MCP server over open-format indices

The end-user interface is **not** a chat UI backed by vector search. It is an **MCP (Model Context Protocol) server** that exposes a curated set of tools over the indices above, plus raw-file access for deeper inspection. Any MCP-capable client (Claude Desktop, Claude Code, Cursor, Continue, and a growing list of others) connects and drives synthesis with these tools instead of relying on a single top-K retrieval.

**Why MCP:**

- It is the emerging standard specifically for tool/resource exposure to LLMs — one server, many clients.
- It decouples the corpus from any single front-end. Swap Claude Desktop for a custom agent, or for a future model, without touching the corpus.
- Tools can return **exhaustive** result sets (e.g., every chunk mentioning a given taxon) — the LLM does the synthesis, not the retrieval layer. Top-K truncation is a per-tool opt-in, not a default.
- Outside-lab tools can use the exact same corpus: either via MCP, or — for non-LLM use cases — by reading the underlying parquet/SQLite/LanceDB files directly.

**Initial tool surface:**

| Tool | Purpose | Backs queries |
|---|---|---|
| `search_taxon(name)` | resolve a string against the DwC taxonomy snapshot to accepted name + synonyms + rank | all taxon-keyed |
| `get_chunks_for_taxon(accepted_name, include_synonyms=True)` | exhaustive chunk list | Q1, Q2, Q5 |
| `get_chunks_for_topic(query, k=None)` | semantic vector search; `k=None` returns all above a threshold | Q6, Q8 |
| `get_figures_for_taxon(accepted_name)` | figure records incl. image paths, captions, bboxes | Q2, Q4 |
| `get_figures_for_anatomy(term)` | figures tagged with an anatomy term | Q8, Q4 |
| `get_papers_by_author(name)` | author-filtered paper list | Q5 |
| `get_chunks_by_section(paper_hash, section_class)` | chunks of a given section type ("description", "embryology", …) | Q2, Q6 |
| `get_bibliography(paper_hash)` | refs + citations from Grobid TEI | Q2 |
| `aggregate_taxonomic_acts(rank, group_by="decade")` | analytic rollup over authority strings | Q7 |
| `list_valid_species_under(taxon)` | DwC-snapshot children filter | Q4 |

The server is thin: each tool is a small function over the parquet/SQLite/LanceDB indices. The storage layer is intentionally format-neutral so non-MCP consumers (a Jupyter notebook, a Snakemake rule, another lab's pipeline) can query the same corpus directly via DuckDB / SQLite / LanceDB without going through the MCP layer at all.

## 4. Concrete next steps, ordered

The Phase A–F roadmap that lived here — Grobid wiring, HybridChunker, language-aware OCR, BHL lookup, DwC taxonomy ingest, anatomy lexicon, taxon mention resolution, reference graph, section-type index, MCP server, the golden test set, and the rest — has been fully delivered. See [`CHANGELOG.md`](CHANGELOG.md) for the per-feature record and §6 for current status. Open items have been moved to the GitHub issue tracker.

## 5. Scale considerations for ~2000 papers

- Rough cost envelope: docling is ~30s–3min/paper on CPU, OCR adds 1–5min/paper for scans. At 2000 papers, plan on **2–5 days of wallclock on a single machine** for stage 1–2 serial; hours with parallelism.
- Embeddings: running a local open-weights model on Bouchet GPUs (see §7) eliminates the OpenAI cost and rate-limit path entirely. Full corpus (~300M tokens) embeds in well under an hour on a single GPU; no API bill. The motivation for switching is reproducibility and removing the silent zero-vector failure mode, not cost — OpenAI would have been ~$6.
- Grobid is a Java server; batch it via its concurrency API rather than spawning per-PDF.
- Keep intermediate artifacts on disk, not in memory — the current file-per-step layout is correct, just add per-paper logs.
- Add a corpus-level `manifest.parquet` or SQLite that rolls up one row per hash (status, pages, language, taxa count, figure count, errors). This is what you'll actually query when inspecting corpus health.

## 6. Status snapshot

Implementation status of the items in §4 is reflected by `[x]`/`[ ]` checkboxes in that section. Per-release rollups live in [`CHANGELOG.md`](CHANGELOG.md). This document focuses on the design and roadmap; for "what is currently done?" consult those two surfaces.

### v0.1 release punch list

Tracked work before the v0.1 tag and full corpus rebuild on Bouchet. Issues are tracked in GitHub; checkboxes here are the running status.

**Bouchet (rebuild dependencies):**

- [ ] [#10](https://github.com/caseywdunn/corpus/issues/10) — Diagnose why ocrmypdf produces zero-text PDFs on ~33 docs. Run flag matrix on Stepanjants/MilosMaley/Pagenstecher; install pngquant; apply winning flags to `prepare_pdf`. Likely recovers ~25–30 docs.
- [ ] [#9](https://github.com/caseywdunn/corpus/issues/9) — Install `deu_latf` Fraktur Tesseract pack on Bouchet. Recovers ~6–8 19th-c. German scans (Goldfuss 1820, Pagenstecher 1869, Brandt 1837, Donitz 1871, Stechow 1921, Hoeven 1836, Schmidtlein 1881, Doflein 1906).
- [x] [#8](https://github.com/caseywdunn/corpus/issues/8) — Vendor-watermark detection in `detect_scan_type`. Recovers Karplus 2014, Browne 1905, Fleming 1828. Implemented in commit `403040b`.
- [ ] [#11](https://github.com/caseywdunn/corpus/issues/11) — Queue `slurm/batch_pass3b.sh` so vision Pass 3b + Pass 3c run as part of the rebuild. Resolves the bulk of 6,841 missing-figure records and applies compound-figure splits across the corpus.
- [ ] **Trigger the full rebuild**: `slurm/batch_process_corpus.sh` (Stage 1) → `batch_pass3b.sh` (vision) → `batch_embed.sh` (BGE-M3) → `batch_grobid.sh` if needed → `batch_biblio.sh`. All chainable via `--dependency=afterok`.
- [ ] **Audit the rebuild output** before tagging: zero-text-extraction count, `pass3c_status` coverage, `missing_figures[]` reduction, figure-type coverage on PyMuPDF rescues.

**Local (no Bouchet needed):**

- [ ] [#5](https://github.com/caseywdunn/corpus/issues/5) — Streamable HTTP transport + OAuth for native Custom Connectors UI (defer if not required for v0.1 collaborator access).
- [ ] [#6](https://github.com/caseywdunn/corpus/issues/6) — `deploy/stack.yaml` default-VPC assumption; either fix or document as a known prerequisite.
- [x] **§10 design item: no absolute paths in served JSON.** Implemented in commit `304257a` as a `package_for_serve.py` scrub + audit pass.

**Tagging:**

- [ ] **Date `CHANGELOG.md`** — replace `2026-04-XX` with the tag date once the rebuild lands and the bundle is uploaded.
- [ ] **Distill served bundle**: `python package_for_serve.py "$OUTPUT_DIR" "$BOUCHET_PROJECT/serve_bundle" --version v0.1.0`.
- [ ] **Sync to S3**: `aws s3 sync serve_bundle/ s3://<corpus-name>/v0.1.0/`.
- [ ] **`git tag v0.1.0`** + push tag + create GitHub release pointing at the changelog section.
- [ ] **EC2 reload**: `systemctl reload corpus-mcp` on the served instance to pick up the new bundle.

**Deferred to post-v0.1** (not blocking the tag):

Tracked in the GitHub issue tracker, not in this document:

- [#12](https://github.com/caseywdunn/corpus/issues/12) — Stage 2 parallelism: assess whether further work is needed beyond `batch_embed.sh`.
- [#13](https://github.com/caseywdunn/corpus/issues/13) — Geographic extraction (§12 Layer 3): NER + GeoNames + locality table.
- [#14](https://github.com/caseywdunn/corpus/issues/14) — Trait extraction + identification keys (Q3).
- [#15](https://github.com/caseywdunn/corpus/issues/15) — Refactor `mcp_server.py`: split 2050-line single file into per-concern modules.
- [#7 part 3](https://github.com/caseywdunn/corpus/issues/7) — any remaining work on the in-text citation graph beyond parts 1 and 2 already merged.

## 7. Switch embeddings from OpenAI → local open-weights on Bouchet

The production run will execute on Bouchet (YCRC), which has idle GPU capacity we've already been allocated. Replacing the OpenAI embedding call with a local open-weights model is a small code change with three real benefits and one cost.

### Why switch

- **Reproducibility.** Pinned model weights produce the same vectors forever. OpenAI has already deprecated embedding models in the past (`text-embedding-ada-002`); we don't want a published corpus whose index can't be rebuilt identically.
- **Removes the silent-failure mode.** Current code inserts zero vectors on API failure (bug #4) — harder to hit when there is no API in the loop at all. We should still raise on embed failures regardless.
- **No network / credits dependency.** The whole pipeline becomes runnable from any Bouchet allocation with only `module load` + conda.
- **Throughput.** An H200 embeds the whole corpus in minutes; OpenAI rate limits and 100-chunk batches are slower in practice.

Cost is **not** the motivation. At `text-embedding-3-small` pricing, the full 2000-paper corpus is ~$6 — trivial. Pay the code-change cost only because the first two reasons warrant it.

### Partition choice

Use `--partition=gpu` (RTX 5000 ADA, 32 GB VRAM, 4/node) — enough for any sentence-transformers model up to ~1–2B parameters at fp16. Reserve `gpu_h200` for jobs that actually need H200-class compute; embeddings don't. If we pilot a 7B embedding model (E5-Mistral, Qwen3-Embedding-8B), then `gpu_h200` is justified.

### Model choice

Default: **BAAI/bge-m3** — 568M params, 1024-dim, 8k context, strong multilingual retrieval (covers our German/French/Russian tail), widely used, stable on HF. Dense+sparse+multi-vector modes; start with dense for drop-in compatibility with LanceDB.

Alternatives to keep in mind:

- `nomic-ai/nomic-embed-text-v1.5` — 137M, 768-dim, Matryoshka, good English. Smaller/faster but weaker multilingual.
- `intfloat/multilingual-e5-large` — 560M, 1024-dim, proven multilingual baseline.
- `Qwen/Qwen3-Embedding-8B` — top of MTEB, 4096-dim, needs `gpu_h200` and more care; overkill unless we find BGE-M3 wanting.

### Implementation status

OpenAI removed entirely (commit `202ff04`); BGE-M3 is the only embedding backend (1024-dim). SLURM submission via `slurm/batch_embed.sh` on `--partition=gpu`. Stage 1 stays CPU-only on `day`. Embed failures raise rather than silently zero-vector.

## 8. Target queries & end-user interface

This section defines what the repo is *for*. Every design choice above should be justifiable against the queries below. If a decision doesn't help serve at least one of these, it's speculative.

### Framing: the endpoint is not a RAG app

Top-K vector retrieval is suited to "answer a question using the best few passages." Most of the target queries below are **enumerative** (find everything matching criteria, then synthesize), **analytic** (aggregate over structured metadata), or **trait-structured** (produce a comparative table). Top-K would either truncate the evidence or return irrelevant near-neighbors.

The corpus therefore produces **a set of indices**, and the user-facing layer is a **query-routing + synthesis layer** that:

- picks the right index (or combination) for the query shape,
- pulls **exhaustive** matches from structured indices (not top-K), then
- hands the full evidence set to an LLM (Opus, 1M-token context) for synthesis, or to a plotting/table function for analytic queries.

Natural implementation shape: an agent with tools over the indices (`search_taxon`, `get_chunks_for_taxon`, `get_figures`, `get_bibliography`, `aggregate_by_year`, `run_sql`), not a single `chat.py` backed by vector search. Vector search is one tool among several, used for semantic/topic queries where the match criterion is fuzzy.

### Eight target query patterns

Generic shapes, traced as **Shape → Path → Output → Gaps in current plan**. Concrete instantiations (specific taxa, anatomy terms, authors, topics) live in the corpuscle's `instructions.md`, not here. Placeholders below: `<species>` / `<genus>` for taxa, `<anatomy>` for anatomy terms, `<topic>` for a free-text topic, `<author X>` / `<author Y>` for personal names.

#### Q1. "List all collection locations of `<species>`."

- **Shape:** enumerative, entity-keyed, structured output.
- **Path:** resolve `<species>` to accepted DwC taxon (+ synonyms) → taxon index returns all chunks/captions/tables mentioning the species → geographic NER on those chunks → dedupe + group by paper.
- **Output:** table (location, coordinates if present, paper, page, quote).
- **Gaps:** geographic NER at chunk scope tied to taxon co-occurrence ([#13](https://github.com/caseywdunn/corpus/issues/13)); tabular output recipe (not yet scoped). Taxon↔chunk offset index is in place (must remain exhaustive, not retrieval-ranked).

#### Q2. "Compose a monographic review of `<genus>`."

- **Shape:** enumerative, genus-level, long-form synthesis.
- **Path:** genus → all child species (incl. synonyms) via the DwC snapshot → all chunks mentioning any of them → section-filtered retrieval (descriptions, habitat, distribution, remarks) via Grobid TEI → all figures + captions → all references → Opus synthesis in one or a few passes.
- **Output:** structured long-form document with inline citations, figures, and a reference list.
- **Gaps:** synthesis recipe for monographic output (not yet scoped). The supporting indices — DwC taxonomy backbone, section-type tagging, figure+caption joint object — are all in place.

#### Q3. "Make a key to identify species in `<genus>`."

- **Shape:** trait-structured, comparative.
- **Path:** genus → all valid species → per-species diagnostic traits (meristics, counts, dimensions, colors, structural presence/absence) from species-description sections → pivot into species × character matrix → generate dichotomous (or polyclave) key.
- **Output:** a usable identification key.
- **Gaps:** trait extraction for identification is deferred to post-v0.1 ([#14](https://github.com/caseywdunn/corpus/issues/14)). A key requires (a) a domain-specific character schema, (b) structured LLM extraction against that schema per species description, (c) a trait index as a first-class artifact, (d) key-generation logic. Substantial enough to warrant its own plan section when we return to it; out of scope for the first production run.

#### Q4. "List all currently valid species in the focal group with a one-paragraph summary and diagnostic figures."

- **Shape:** enumerative over a closed taxonomic set; per-entity synthesis; figure selection.
- **Path:** DwC snapshot → list of currently-valid species → per species: all mentioning chunks → Opus one-paragraph summary → figure index: find figures captioned/tagged with the species and a diagnostic anatomy term → assemble.
- **Output:** a catalog document / website.
- **Gaps:** vision-model figure-anatomy classification at corpus scale (Pass 3b on the v0.1 rebuild — see [#11](https://github.com/caseywdunn/corpus/issues/11)); a per-species batch synthesis recipe (not yet scoped); concept of "diagnostic structure" — either user-curated per genus or learned from description sections.

#### Q5. "Summarize all of `<author X>`'s comments about `<author Y>`."

- **Shape:** author-filtered + entity-filtered; meta-commentary extraction.
- **Path:** author index (from Grobid headers) → all papers authored by `<author X>` → chunks mentioning `<author Y>` (as person, not just cited work) → LLM filter: which mentions are substantive commentary vs. routine citation → synthesis.
- **Output:** narrative summary with direct quotes + citations.
- **Gaps:** explicit per-author rollup over Grobid headers; person-name NER that distinguishes cited-author-as-entity from cited-work; "meta-commentary vs. citation" classification (LLM in the loop at query time is fine — don't pre-compute).

#### Q6. "Write a summary of `<topic>` across the corpus."

- **Shape:** topic-semantic, corpus-wide.
- **Path:** (a) section-type index: chunks from sections classified to a relevant section type + (b) vector search for the topic and related concepts across all chunks → union → Opus synthesis with citations.
- **Output:** narrative summary with references.
- **Gaps:** this is the query the vector index was built for. Also benefits from Grobid-provided section labels. No new scope items — just needs what's already planned to work well.

#### Q7. "Plot number of species in the focal group described per decade."

- **Shape:** analytic aggregation over structured metadata. No LLM required.
- **Path:** for every currently-valid species in the DwC snapshot: parse authority string (`"Author, YYYY"` → year) → cross-check against the paper in the corpus where available → aggregate by decade → plot.
- **Output:** a figure.
- **Gaps:** authority-string parsing and taxonomic-act index. Much of the data comes from the DwC snapshot directly — the corpus only needs to supply the original-description PDFs for cross-reference. Consider whether this query even requires the corpus or can be answered from the DwC snapshot alone; if the latter, it's a good sanity-check baseline.

#### Q8. "Summarize what is known about the structure and function of `<anatomy>`."

- **Shape:** topic-semantic, anatomy-keyed.
- **Path:** anatomy index returns all chunks tagged with `<anatomy>` → plus vector search for adjacent concepts (synonyms, related structures from the lexicon) to catch non-lexical matches → figures tagged with `<anatomy>` → Opus synthesis.
- **Output:** narrative summary with citations and figures.
- **Gaps:** none structural — anatomy-term index and vision-based figure anatomy tagging are both in place. Quality of vision tagging at corpus scale lands with the v0.1 rebuild ([#11](https://github.com/caseywdunn/corpus/issues/11)).

### Indices this implies

Vector embeddings are one of ten indices. Only Q6 and Q8 treat vector search as a primary retrieval path; the others use it only as a supplementary recall mechanism, if at all.

| # | Index | Produced by | Used by |
|---|---|---|---|
| 1 | Taxon mention index (chunk-offset) | gnfinder + DwC resolution | Q1, Q2, Q4, Q5 (indirectly), Q7 |
| 2 | DwC taxonomy backbone (genera/species/synonyms/authorities) | DwC snapshot | Q1, Q2, Q3, Q4, Q7 |
| 3 | Bibliographic index (author, year, title, journal) | Grobid headers + CrossRef/OpenAlex | Q2, Q5, Q7 |
| 4 | Reference graph | Grobid bibliography | Q2 |
| 5 | Section-type index (descriptions, distributions, etc.) | Grobid TEI structure | Q2, Q6 |
| 6 | Anatomy-term index | domain lexicon + NER | Q1 (anatomy qualifier), Q3, Q8, Q4 (diagnostic) |
| 7 | Geographic index | NER + GeoNames | Q1 |
| 8 | Figure index (image + caption + taxon + anatomy tags) | docling + vision LLM classifier | Q2, Q4, Q8 |
| 9 | Trait / character matrix (species × character) | structured LLM extraction against schema | Q3 |
| 10 | Vector index (semantic fallback) | local embeddings (§7) | Q6, Q8 (primary); others (recall boost) |

### Coverage of these queries by the current plan

The scope items these queries demand — DwC backbone, anatomy lexicon + NER, authority-string / taxonomic-act parsing, section-type tagging, figure anatomy classification, reference graph, and the MCP tools layer — are all in place. The two genuine gaps remaining are:

- **Trait extraction + key generation (Q3).** Deferred to post-v0.1 ([#14](https://github.com/caseywdunn/corpus/issues/14)). Will get its own plan section when we pick it up.
- **Synthesis recipes (Q2, Q4).** The indices + MCP tools are in place, but the specific Opus-driven prompts and map-reduce patterns for monographic review (Q2) and per-entity catalog (Q4) are not yet scoped. Cheap to iterate on once collaborators are using the served bundle — prototype against the golden test set before committing to any one recipe.

## 9. Three-pass figure extraction with ROI sidecar

The initial figure pipeline (Phase D, §3 "Figure+caption as a first-class object") handles the common case — docling extracts one `Picture` per real figure and we classify + dedupe. Two harder cases surface in the demo set and drive this design:

1. **Adjacent figures merged into one image** — docling sometimes extracts two logically-separate figures as a single bbox when there's insufficient whitespace or heading text between them. Concrete example from the demo set: a paper with adjacent figures (Fig. 3, panels A–B; Fig. 4, panels A–D) comes out as one image whose caption is parsed as `Fig. 3`.
2. **Panels inside one figure image** — docling treats multi-panel figures as a single image, because the panel labels (A, B, C, …) are embedded content, not structural PDF elements. Panel arrangement can be non-rectilinear (L-shapes, insets), so sorting crops by reading order is brittle.

Both have the same underlying need: **content-aware region labels inside an already-extracted figure image**, with labels derived from what's physically in the image (panel letters, embedded "Fig. N" text) rather than guessed from bbox geometry.

### Architecture

Three cleanly-separable passes over the per-paper artifacts. Passes 1–2 are the current pipeline; 2.5 is cheap and always-on; 3 is the new content-aware layer and is opt-in.

| Pass | Input | Output | Cost |
|---|---|---|---|
| 1. Extraction | PDF | raw docling `Picture` objects + captions + bboxes | existing |
| 2. Classification | raw items | `figure_type` + dedupe + semantic filenames + coequal-panel letters | existing |
| 2.5. Caption parsing + gap detection | `figures.json` + `text.json` | `panels_from_caption[]` per figure, `missing_figures[]` at doc level | cheap; always on |
| 3a. Tesseract ROI pass | figure PNGs + captions | `rois[]` sidecar per figure, compound-figure renames | ~1 s / figure; opt-in via `--content-aware-figures` |
| 3b. Vision-LLM ROI pass (future) | Pass 3a's unresolved figures | higher-quality ROIs from Claude over the image | pennies per figure; opt-in |

Each pass is rerunnable independently. Pass 3a can be skipped on first ingest and added later without re-parsing PDFs or re-embedding.

### Sidecar schema (inline in `figures.json`)

A single `rois[]` array per figure object covers both concerns with different `type` values. Panels and sub-figures live in the same list so MCP tools can navigate without branching:

```json
{
  "figure_id": "docling_5",
  "filename": "fig_3-4.png",
  "figure_number": "3-4",
  "image_size_px": [1015, 1305],
  "caption_text": "Fig. 3 A. upper and B. lower views of <subject>. Scale bar 5 mm.",
  "rois": [
    {"type": "figure", "figure_number": "3",
     "roi_px": [0, 0, 500, 1305],
     "caption_text": "Fig. 3 A. upper and B. lower views of <subject>.",
     "source": "ocr:tesseract"},
    {"type": "figure", "figure_number": "4",
     "roi_px": [500, 0, 1015, 1305],
     "caption_text": "Fig. 4 <Genus species>. A. <variant 1> and B., C. <variant 2>…",
     "source": "ocr:tesseract"},
    {"type": "panel", "parent_figure_number": "3", "label": "A",
     "roi_px": [0, 0, 500, 600],
     "description_from_caption": "upper view of <subject>",
     "source": "ocr:tesseract"},
    {"type": "panel", "parent_figure_number": "3", "label": "B",
     "roi_px": [0, 600, 500, 1305],
     "description_from_caption": "lower view of <subject>",
     "source": "ocr:tesseract"}
  ],
  "pass3_status": "completed",
  "previous_filenames": ["fig_3.png"]
}
```

Why inline: keeps the joint object self-contained; MCP tools don't cross-reference a second file; ROI annotations travel with the figure through any pipeline reorg.

### Filename convention

Filenames carry the figure-number(s) so `ls figures/` is self-documenting. Range notation matches bibliographic convention:

| Case | Filename |
|---|---|
| Single figure, no panels | `fig_3.png` |
| Multi-panel coequal (Pass 2 dedupe splits panels) | `fig_3_a.png`, `fig_3_b.png`, … |
| Compound, contiguous range | `fig_3-4.png` |
| Compound, non-contiguous | `fig_3+7.png` |
| Plate | `plate_2.png` (+ `-3`, `_a` variants) |

When Pass 3a renames a file (Pass 2 saved it as `fig_3.png`, Pass 3a detects the embedded "Fig. 4" and renames to `fig_3-4.png`), the previous name is recorded in `previous_filenames[]` so any cached URL can be traced to its current file.

Panels inside a compound image live only in the ROI sidecar — not as separate panel files — since generating six crops (Fig 3 A+B + Fig 4 A+D) for content nobody may query is wasteful. Crop-on-demand via the MCP tool handles that.

### Why not split compound images into separate files

- **Reversible**: a better tool (vision LLM in 3b) can re-segment later without losing the first attempt
- **Preserves provenance**: what docling gave us stays what docling gave us; any split is derived data
- **Crops stay cheap**: PIL crop of a ROI at retrieval time is sub-second; zero cost for panels never queried
- **Honest about uncertainty**: ROI boundaries inferred from a detected panel label are approximate; we don't fabricate crisp image files that imply more precision than we have

### Pass 2.5 — caption parsing + gap detection (cheap, always on)

Two small passes over existing artifacts:

- **Caption-panel parser**: regex-based extraction of `[{label, description}]` from caption text. Handles the common conventions: `A.` / `(A)` / `A:` / `A–C.`, multilingual prefixes (A, A-C, Рис. 3A). Recorded as `panels_from_caption[]` + `panel_count_from_caption` on every figure.
- **Missing-figures scan**: scan `text.json` for `Fig. N` / `Figure N` mentions; any N that's not present as a `figure_number` in `figures.json` gets recorded in a per-paper `missing_figures[]` list with the caption text (parsed out of the running text). Closes the merged-adjacent-figures observability gap even without Pass 3.

Pass 2.5 does **not** change filenames or images — it only annotates.

### Pass 3a — Tesseract ROI pass (opt-in)

For each figure PNG:

1. Load caption's expected panel labels (from Pass 2.5).
2. Run `tesseract image_to_data` on the PNG, get text tokens + bboxes + font-size.
3. Filter for panel-label candidates: single capital letter or letter+period/paren, isolated (whitespace on both sides), font size above median for the image, within ~30 px of an edge of some visual region.
4. Filter for "Fig. N" labels: any `Fig\.?\s*\d+` string inside the image indicates the image contains multiple figures.
5. Cross-check detected panel labels against caption-expected set; emit ROIs when consistent.
6. For compound detection: if two or more `Fig. N` labels are found, split the image region between them (Voronoi by label position works for most layouts), emit `type: "figure"` ROIs per number, and rename the file to range notation.
7. Mark `pass3_status` on the figure: `"not_run"` (skipped), `"completed"`, `"ocr_unresolved"` (ran but couldn't reconcile with caption), `"no_compound_or_panels"` (trivially resolved).

OCR cost on the corpus: ~1 s per figure × ~20 figures/paper × 2000 papers ≈ 11 hours. Tolerable as a single Bouchet run; skippable on dev.

### Pass 3b — vision-LLM fallback (later)

Invoked only for `pass3_status = "ocr_unresolved"` figures (or globally with an explicit flag). One Claude call per figure with the image + caption + OCR candidates, asking for `{label, bbox_px, description}` JSON. Pennies per figure at Haiku; adds up at corpus scale, so gated behind a flag.

### MCP tool surface additions

Two new tools once Pass 3a lands:

```text
list_figure_rois(paper_hash, figure_id)
  → [{type, label/figure_number, roi_px, caption_text, description_from_caption, source}, ...]

get_figure_roi_image(paper_hash, figure_id, roi_index, save_to=None)
  → {image_path: cropped.png, roi_px, caption_text}
    # PIL crop on demand, cached to <hash_dir>/figures/crops/ keyed on
    # (figure_id, roi_index) so repeated asks are free
```

`get_figures_for_taxon` / `get_figures_for_anatomy` gain a `prefer_panels: bool = False` flag so queries asking about a specific anatomy can return panel-level crops when they exist.

## 10. Deployment to AWS (post-Bouchet, for collaborator access)

Production runs on Bouchet, but the served corpus needs a long-lived public-ish endpoint so collaborators at other institutions and lab members can query it without pulling 10 GB of artifacts down to their laptops. AWS is the target; the design needs to be settled now (not after deployment) because the per-paper artifact set is the contract that the distillation, the MCP server, and downstream non-MCP consumers all key off.

### Two-bundle model

Bouchet produces a "build" bundle. AWS holds a "served" bundle. The two are not the same.

| Concern | Build bundle (Bouchet output) | Served bundle (S3 / EC2) |
|---|---|---|
| Audience | the pipeline, golden-test scripts, re-runs | MCP clients, DuckDB users, web tools |
| Per-paper size | ~5 MB | ~1.5 MB (≈70% smaller; ~95% if PDFs excluded) |
| Includes | everything: `processed.pdf`, `docling_doc.json`, `visualizations/`, `grobid.tei.xml`, `pipeline.log` | only the JSONs that MCP tools read, the figure PNGs, the LanceDB index, the DwC taxonomy snapshot, and a manifest |
| Lifecycle | rebuilt each pipeline run | versioned snapshot (`v1.0.0`, …) |

**Files in the served-bundle whitelist** — this is the contract:

- `summary.json`, `metadata.json`, `references.json`
- `text.json`, `chunks.json`
- `figures.json`, `taxa.json`, `anatomy.json`
- `figures/*.png` (all extracted figure images, including Pass 3c-renamed compounds)
- `vector_db/lancedb/` (the embedded chunks)
- `taxonomy.sqlite` (DwC) and any other corpus-level indices
- `bundle_manifest.json` (top-level — see below)

**Excluded** from the served bundle (build-only):

- `processed.pdf` — keep on Bouchet + cold-storage S3 for re-runs; not needed by MCP tools. Optional flag to include (~3 GB extra) if a future tool wants page-image retrieval.
- `docling_doc.json` — raw docling dump; only useful for re-extraction.
- `visualizations/*.png` — pure QC overlays.
- `grobid.tei.xml` — already distilled into `metadata.json` + `references.json`.
- `pipeline.log`, `scan_detection.json`, `figures_report.html` — debug / human-QC.

For 2000 papers: build bundle ~10 GB, served bundle ~3 GB (with PDFs ~6 GB). Small enough to fit comfortably on a $1/mo EBS volume.

### Distillation step

A new `package_for_serve.py` on Bouchet, run after `batch_embed.sh`:

```bash
python package_for_serve.py "$OUTPUT_DIR" "$BOUCHET_PROJECT/serve_bundle" \
    --version v1.0.0 \
    --include-pdfs false
```

It walks `documents/<HASH>/`, copies whitelisted files into `serve_bundle/documents/<HASH>/`, copies the LanceDB and `taxonomy.sqlite`, and writes `serve_bundle/bundle_manifest.json`:

```json
{
  "bundle_version": "v1.0.0",
  "created_at": "2026-04-16T12:00:00Z",
  "pipeline_git_sha": "<commit>",
  "embedding_model": "BAAI/bge-m3",
  "embedding_dim": 1024,
  "vision_backend": "vision:qwen2.5-vl-7b-instruct",
  "taxonomy_snapshot_date": "2026-03-15",
  "paper_count": 2014,
  "figure_count": 38291,
  "chunk_count": 412677,
  "includes_pdfs": false
}
```

The manifest is what collaborators cite ("queried against corpus v1.0.0"); it's also what the MCP server reports via a new `bundle_info` tool so clients can detect stale endpoints.

### AWS layout

Three-tier, deliberately boring:

1. **S3 bucket — `s3://<corpus-name>/`** — source of truth, versioned. Each release lives at `s3://<corpus-name>/v1.0.0/` and is immutable. Bouchet pushes via `aws s3 sync serve_bundle/ s3://…/v1.0.0/`. Lifecycle rule moves unused versions to Glacier after 90 days.
   - Public-read on the parquet/SQLite/LanceDB files lets non-MCP collaborators query directly with DuckDB / pandas — zero infra cost on our end.
2. **MCP server on EC2 — `t3.small` + 10 GB gp3 EBS** — single instance (~$15/mo + $1/mo storage), pulls the bundle from S3 on deploy, runs the MCP server as a systemd service over SSE/HTTP transport. Faster than Fargate (no cold starts, persistent LanceDB in memory). One CloudWatch dashboard, structured JSON logs to S3.
3. **Edge — CloudFront → ALB → EC2**, HTTPS, bearer-token auth (see *Authentication & access control* below). Rate-limit at CloudFront / WAF.

**Updates** = re-run pipeline on Bouchet → distill → `aws s3 sync` to a new version → SSH to EC2 and `systemctl reload corpus-mcp` (which pulls the new version from S3). Old versions stay in S3 for reproducibility.

**Total ongoing cost**: ~$20/mo. Grows slowly with collaborator volume since reads are cheap.

### Authentication & access control

Target audience is ~20 trusted collaborators. The threat model is runaway bills and random traffic, not data protection — the corpus is read-only and the content is already-published literature.

**Shared bearer token first.** A Starlette/FastAPI middleware in `mcp_server.py` checks `Authorization: Bearer <TOKEN>` against a single secret (SSM Parameter Store on EC2, or `/etc/corpus/token` mode 600 as a fallback). Missing or wrong header → 401. Distribute the token over Slack / 1Password; collaborators paste `"headers": {"Authorization": "Bearer xxx"}` into their `.mcp.json` alongside the server URL. Rotation = generate new secret, update SSM, resend. Afternoon of work, all of it one-time.

**Per-user tokens when revocation becomes a real need.** Replace the single-secret check with a `tokens.json` lookup (hashed token → user, issued-at, last-used). Add a small CLI — `corpus-admin {issue,revoke}-token <email>`. ~20 lines of change, still no AWS services.

**API Gateway + usage plans only when quotas or per-key analytics matter.** Adds IAM, Terraform, and cold-start worries; not worth the tax for 20 friendly users.

**Defense in depth, independent of auth.**

- AWS Budgets alarm at $50/mo with email notification — the cheapest way to notice a runaway.
- CloudFront + WAF rate-limit by IP (~100 req/min). Caps a misbehaving client before it hits the EC2.
- If collaborators are at a stable set of institutions, an EC2 security group whitelisting their CIDRs is a second gate that costs nothing.

### Design constraints (now satisfied)

Three constraints had to flow back into the Bouchet pipeline before AWS deployment so the served-bundle contract was stable: (1) a single source-of-truth whitelist of served-bundle files (`PER_PAPER_FILES` in `package_for_serve.py`), (2) no absolute paths in any served JSON (distillation-time scrub + JSON-walk audit, see `tests/test_package_for_serve.py`), and (3) a per-bundle version stamp written into `bundle_manifest.json` (`bundle_version`, `pipeline_git_sha`, `embedding_model`, `taxonomy_snapshot_date`) so collaborators can cite a corpus version. All three are in place.

### What this section deliberately does not pin down

- Whether the MCP server eventually gets a thin HTML/web UI on top (Streamlit, Next.js, …). Out of scope until the MCP-only experience has actual users.
- Multi-region failover, autoscaling, blue/green deploys. Single-instance is fine until it isn't.
- Whether to mirror the served bundle on a Cloudflare R2 / Backblaze B2 alternative for cost. Defer until S3 egress actually shows up on a bill.
- Authentication beyond bearer tokens. If we expose this to a wider audience or need institutional SSO, evaluate Cognito / sign-in-with-ORCID then, not now.

## 11. Bibliographic authority database

The corpus has tens of thousands of parsed references across the full paper set, but without cross-paper deduplication there is no citation graph and no link between taxonomic authority strings and their original-description papers. Typically only a small fraction of historical references carry DOIs — the rest need a fallback identifier. This section defines a single bibliographic authority that unifies corpus papers, cited references, and DwC taxonomic authorities into one queryable graph stored as `<corpuscle>/biblio_authority.sqlite` alongside the per-paper artifacts.

### GUID scheme

Every bibliographic entity gets a single globally-unique identifier, assigned in priority order:

1. **DOI** (normalized lowercase, no URL prefix): `10.1007/s12526-013-0148-x`
2. **BHL Part/Item ID**: `bhl:part/12345` or `bhl:item/31879`
3. **Normalized citation key**: `corpus:smith|1995|introduction to the focal group`

The fallback key is the normalized string `normalize(first_author_surname)|year|normalize(title_first_40_chars)` — lowercase, diacritics stripped (NFD), punctuation removed, whitespace collapsed. Human-readable and reversible. Titles shorter than 40 characters use their full length. This means "Smith 1995", "Smith, J. (1995)", and "SMITH, J., 1995" all produce the same key, which is the main dedup mechanism for references without DOIs.

### Schema

```sql
-- One row per unique bibliographic work
CREATE TABLE works (
    work_id        TEXT PRIMARY KEY,   -- DOI | bhl:* | corpus:*
    guid_type      TEXT NOT NULL,      -- 'doi' | 'bhl' | 'corpus_key'
    title          TEXT,
    year           INTEGER,
    journal        TEXT,
    doi            TEXT,
    bhl_item_id    TEXT,
    bhl_part_id    TEXT,
    openalex_id    TEXT,
    corpus_hash    TEXT,               -- 12-char SHA-256 if in corpus, else NULL
    in_corpus      INTEGER NOT NULL DEFAULT 0,
    source         TEXT NOT NULL,      -- 'corpus_paper' | 'cited_reference' | 'taxon_authority'
    confidence     REAL DEFAULT 1.0,
    created_at     REAL NOT NULL,
    updated_at     REAL NOT NULL
);

-- Ordered author list per work
CREATE TABLE work_authors (
    work_id            TEXT NOT NULL REFERENCES works(work_id),
    position           INTEGER NOT NULL,
    surname            TEXT NOT NULL,
    surname_normalized TEXT NOT NULL,   -- lowercase, no diacritics
    forename           TEXT,
    PRIMARY KEY (work_id, position)
);

-- Directed citation edges
CREATE TABLE citations (
    citing_work_id     TEXT NOT NULL REFERENCES works(work_id),
    cited_work_id      TEXT NOT NULL REFERENCES works(work_id),
    citing_corpus_hash TEXT NOT NULL,   -- which paper's references.json
    grobid_xml_id      TEXT,
    raw_citation       TEXT,
    match_method       TEXT NOT NULL,   -- 'doi_exact' | 'alias_exact' | 'bhl_lookup' | 'title_fuzzy' | 'author_year_only'
    match_score        REAL DEFAULT 1.0,
    PRIMARY KEY (citing_work_id, cited_work_id, citing_corpus_hash)
);

-- Alternate citation forms that resolved to the same work
CREATE TABLE work_aliases (
    alias_key  TEXT NOT NULL,          -- normalized "smith|1995|introduction to the foca"
    work_id    TEXT NOT NULL REFERENCES works(work_id),
    PRIMARY KEY (alias_key, work_id)
);

-- Links Darwin Core taxa (taxonID from the configured taxonomy snapshot)
-- to their original-description works.
CREATE TABLE taxon_work_links (
    taxon_id   TEXT NOT NULL,
    work_id    TEXT NOT NULL REFERENCES works(work_id),
    link_type  TEXT NOT NULL,          -- 'original_description' | 'authority_match'
    confidence REAL DEFAULT 1.0,
    PRIMARY KEY (taxon_id, work_id, link_type)
);

CREATE TABLE build_meta (key TEXT PRIMARY KEY, value TEXT);
```

Plus indexes on: `works(doi)`, `works(corpus_hash)`, `works(year)`, `works(in_corpus)`, `work_authors(surname_normalized)`, `citations(cited_work_id)`, `citations(citing_work_id)`, `taxon_work_links(work_id)`.

### Build pipeline: `build_biblio_authority.py`

Reads existing per-paper artifacts only — no Grobid re-runs needed. Three phases:

**Phase 1 — Seed from corpus papers.** Walk `output/documents/*/metadata.json`. Create a `works` row per paper with `in_corpus=1`. Use DOI as `work_id` if available, else compute `corpus:` hash. Insert authors into `work_authors`, register alias key in `work_aliases`.

**Phase 2 — Ingest cited references + build citation graph.** Walk `output/documents/*/references.json`. For each reference, match to an existing work via a cascade:
1. *DOI exact* — if reference has a DOI, normalize and look up.
2. *Alias key exact* — compute normalized key (first-author-surname + year + title prefix), check `work_aliases`. Handles ~60% of non-DOI matches.
3. *BHL lookup* (optional, gated behind `--enrich-bhl`) — query BHL API by author + year + title. If matched, use `bhl:` as the GUID. Strong coverage of pre-1950 zoological literature, so many historical references that would otherwise fall through to `corpus:` hashes get proper BHL identifiers. Requires network access; skipped in offline builds.
4. *Fuzzy title match* — query candidates by author surname + year, rank by `rapidfuzz.fuzz.token_set_ratio` on title. Threshold 80 = auto-match; 60–80 = lower confidence.
5. *Author+year only* — last resort for references with no parsed title. Match if exactly one candidate; confidence 0.6.

Insert citation edge into `citations` with match method recorded. Expected: ~15,000–25,000 unique works after dedup.

**Phase 3 — Link taxonomic authorities to works.** Parse `scientificNameAuthorship` strings from the configured DwC taxonomy snapshot (e.g., "Eschscholtz, 1829", "(Huxley, 1859)"). Search `works` by author+year. Insert into `taxon_work_links`. Create stub works for unmatched authorities.

Idempotent via `INSERT OR IGNORE` / `INSERT OR REPLACE`. Dependency: `rapidfuzz`.

### New MCP tools this enables

| Tool | Purpose |
|---|---|
| `get_citation_graph(work_id, direction, depth)` | Citing/cited-by works around a given work |
| `resolve_reference(query, author, year)` | Free-text bibliographic lookup against the authority DB |
| `get_missing_references(min_citations)` | Works cited by corpus but not in corpus, ranked by citation count |
| `get_original_description(taxon_name)` | DwC scientificNameAuthorship → work → in_corpus? |
| `get_works_by_author(surname)` | Full bibliography across authority DB, not just corpus papers |

Enhancement to existing `get_bibliography`: add `resolved=True` flag to join each reference to its authority DB work record (work_id, in_corpus, corpus_hash, cited_by_count), making bibliographies navigable.

### Queries this unlocks

- "What information about a focal taxon exists in papers that cite this paper?" → `get_citation_graph` → `get_chunks_for_taxon` on each citing paper
- "Which cited references are missing from the corpus?" → `get_missing_references`
- "Which papers are the original descriptions for species in a focal genus?" → `list_valid_species_under` → `get_original_description` per species
- "Plot species descriptions per decade" → join `taxon_work_links` to `works`, group by `year / 10`

### Integration with §10 (AWS served bundle)

Add `biblio_authority.sqlite` to the served-file whitelist alongside `taxonomy.sqlite`.

## 12. Text-span mention layers: bibliography, taxonomy, geography

The corpus currently stores structured data about papers (metadata, references, figures) but does not track *where in the text* those entities are mentioned. Three parallel annotation layers — bibliographic, taxonomic, and geographic — will deep-link structured entity databases into specific text spans, making every mention queryable and cross-referenceable.

### Design principle: external tables, not inline annotation

Mentions are stored in dedicated database tables that point into text, not as annotations baked into the per-paper JSON artifacts. Reasons:

- **Independent rebuilds.** Re-running the taxonomy linker doesn't require re-running OCR. Updating geographic parsing doesn't touch bibliography. Each layer has its own lifecycle.
- **Immutable paper artifacts.** The pipeline produces per-paper JSONs (`chunks.json`, `references.json`, etc.) that are content-addressed and stable. Layering annotations into those files couples extraction to resolution.
- **Cross-paper queries are database queries.** "All records of a focal species" joins across papers — that's a table scan, not a document search.

The `biblio_authority.sqlite` pattern from §11 is the model. Taxonomy and geography follow the same shape.

### Common mention schema

Each layer has its own mention table, but all share the same spine:

```sql
-- Bibliographic mentions (in-text citations → works)
CREATE TABLE biblio_mentions (
    mention_id     INTEGER PRIMARY KEY,
    work_id        TEXT NOT NULL,          -- FK to works(work_id) in biblio_authority
    corpus_hash    TEXT NOT NULL,
    chunk_index    INTEGER NOT NULL,
    char_start     INTEGER NOT NULL,
    char_end       INTEGER NOT NULL,
    mention_text   TEXT NOT NULL,          -- surface form: "(Dunn et al., 2005)" or "[23]"
    confidence     REAL DEFAULT 1.0,
    method         TEXT NOT NULL           -- 'grobid_tei_ref', 'regex_fallback'
);

-- Taxonomic mentions (species/genus names → Darwin Core taxonID)
CREATE TABLE taxon_mentions (
    mention_id     INTEGER PRIMARY KEY,
    taxon_id       TEXT,                   -- DwC accepted taxonID (NULL if unresolved)
    matched_name   TEXT NOT NULL,          -- the name as resolved (accepted or synonym)
    corpus_hash    TEXT NOT NULL,
    chunk_index    INTEGER NOT NULL,
    char_start     INTEGER NOT NULL,
    char_end       INTEGER NOT NULL,
    mention_text   TEXT NOT NULL,          -- surface form (abbreviated or full binomial)
    taxon_rank     TEXT,                   -- 'species', 'genus', 'family', etc.
    confidence     REAL DEFAULT 1.0,
    method         TEXT NOT NULL           -- 'regex_taxonomy', 'gnfinder', 'manual'
);

-- Geographic mentions (localities → coordinates)
CREATE TABLE geo_mentions (
    mention_id     INTEGER PRIMARY KEY,
    locality_id    TEXT,                   -- FK to a localities table (NULL if unresolved)
    corpus_hash    TEXT NOT NULL,
    chunk_index    INTEGER NOT NULL,
    char_start     INTEGER NOT NULL,
    char_end       INTEGER NOT NULL,
    mention_text   TEXT NOT NULL,          -- surface form: "Villefranche-sur-Mer", "Station 42"
    latitude       REAL,
    longitude      REAL,
    depth_m        REAL,
    confidence     REAL DEFAULT 1.0,
    method         TEXT NOT NULL           -- 'spacy_ner', 'geonames_match', 'regex_coords'
);

-- Normalized locality entities (deduplicated across papers)
CREATE TABLE localities (
    locality_id    TEXT PRIMARY KEY,
    name           TEXT NOT NULL,
    latitude       REAL,
    longitude      REAL,
    country        TEXT,
    geonames_id    INTEGER,
    source         TEXT NOT NULL           -- 'geonames', 'manual', 'inferred'
);
```

### Layer 1: Bibliographic mentions

**Source data.** Grobid TEI XML already tags in-text citations as `<ref type="bibr">` with `target` attributes pointing to bibliography entries. This data is parsed during Stage 2 but currently discarded — only the bibliography entries themselves flow into `references.json`.

**Extraction.** Parse `<ref type="bibr">` elements from `grobid.tei.xml`, map each to a chunk offset in `chunks.json` (by character position in the full text), and resolve the `target` attribute to a `work_id` in `biblio_authority.sqlite`.

**What this enables.**
- "Which papers cite Totton 1965, and what do they say about it?" → `biblio_mentions` WHERE `work_id = <totton_1965>` → return surrounding chunk text
- "Show me the context of every citation of this paper" → all mention spans with their chunk text
- Closes the loop between "this paper cites X" (references.json) and "here is where it cites X" (mention span)

**Feasibility.** Lowest effort of the three layers. Grobid has already done the hard work; this is plumbing.

### Layer 2: Taxonomic mentions

**Source data.** Running text in `chunks.json` contains binomial names (e.g., *Genus species*), abbreviated forms (*G. species*), genus-only references, and higher taxa. The configured Darwin Core taxonomy snapshot (`taxonomy.sqlite`, built by `ingest_taxonomy.py` from any DwC source) provides the authority for resolution.

**Extraction pipeline.**
1. **Detection** — `gnfinder` (Global Names) for raw name detection, strong on historical variants and abbreviations. Supplement with regex for italicized binomials in born-digital PDFs.
2. **Resolution** — match each detected name against accepted names and synonyms in `taxonomy.sqlite`. Record `taxon_id` (DwC `taxonID`) for resolved names, flag unresolved names for manual review.
3. **Abbreviation expansion** — track the most recent genus mention per paper to expand abbreviated forms (e.g., *G. species* → *Genus species* when *Genus* was last mentioned).

**What this enables.**
- "Which papers mention a focal species?" → `taxon_mentions` WHERE `taxon_id = <species>`
- "Make a map of all records of a focal species" → join `taxon_mentions` to `geo_mentions` by (corpus_hash, chunk_index) proximity
- All taxon-keyed queries in §8 become precise span-level lookups rather than text search

**Feasibility.** Medium effort. A focused taxonomic domain (typically a single order or family with hundreds of valid species) bounds the taxonomy. `gnfinder` handles the NER; resolution against the DwC snapshot is a lookup. Main challenge: abbreviated forms and OCR-garbled names in historical text.

### Layer 3: Geographic mentions

**Source data.** Running text mentions sampling localities in several forms: named places ("Villefranche-sur-Mer", "Bay of Naples"), coordinates ("32°15′N, 64°40′W"), station references ("Station 42", "Sta. 1247"), depth ranges ("200–400 m"), and vessel names ("R/V *Dana*", "HMS *Challenger*").

**Extraction pipeline.**
1. **Named location NER** — spaCy `GPE`/`LOC` entities as a starting point. Filter for geographic relevance (not all GPEs are sampling locations).
2. **Geocoding** — resolve named locations to coordinates via GeoNames API or a local GeoNames dump. Record `geonames_id` for traceability.
3. **Coordinate extraction** — regex for degree-minute-second and decimal-degree patterns. Associate with nearby taxon mentions for taxon × locality co-occurrence.
4. **Structured sampling events** (future) — extract (location, depth, date, vessel, taxon) tuples from methods sections. This is harder and may benefit from LLM extraction.

**What this enables.**
- "List all taxa collected at a focal locality" → `geo_mentions` WHERE `locality_id = <locality>` → join to nearby `taxon_mentions`
- "Make a map of all records of a focal species" → `taxon_mentions` for the species → co-occurring `geo_mentions` with coordinates → plot
- GBIF-style occurrence maps built from literature rather than specimen databases

**Feasibility.** Highest effort, lowest initial precision. Named locations are detectable but associating them with the right taxon and distinguishing sampling locations from other geographic references (e.g., "the Mediterranean fauna") requires context. Start with named locations + coordinates only; structured sampling events are a later enrichment.

### Sequencing

1. **Bibliographic mentions** — first, since Grobid TEI already has the data. Closes the citation-context loop.
2. **Taxonomic mentions** — second. Enables all taxon-keyed queries and is prerequisite for taxon × locality co-occurrence in layer 3.
3. **Geographic mentions** — third. Depends on layer 2 for meaningful taxon × locality joins.

### MCP tool surface

| Tool | Layer | Purpose |
|---|---|---|
| `get_citation_context(work_id)` | biblio | All in-text mentions of a cited work, with surrounding chunk text |
| `get_taxon_mentions(name, include_synonyms=True)` | taxon | All text spans mentioning a taxon, resolved to accepted name |
| `get_locality_mentions(locality, radius_km=None)` | geo | All text spans mentioning a locality (or within radius) |
| `get_taxon_localities(name)` | taxon+geo | Co-occurring taxon × locality pairs across the corpus |
| `get_locality_taxa(locality)` | geo+taxon | All taxa mentioned near a given locality |
