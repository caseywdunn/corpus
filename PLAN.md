# PLAN.md — Siphonophore corpus pipeline

A roadmap from the current prototype toward a production-quality corpus of ~2000 siphonophore papers spanning late-18th-century to modern literature in English, German, French, Russian, and others.

## 1. Assessment of current repo

### What works well
- **Hash-based content addressing** is the right foundation for this project. SHA-256 prefix directories give idempotent resume, natural dedup across messy input trees, and survive rename/reorg of source collections.
- **Two-stage split** (local extraction vs. OpenAI embeddings) is correct — embeddings are the expensive, rate-limited, and re-runnable step.
- **Figure QC visualizations** (red word boxes + yellow/orange figure bboxes per page) are exactly the kind of intermediate artifact this corpus needs. Keep and extend these.
- **PyMuPDF fallback for figure extraction** acknowledges docling doesn't always deliver.

### Bugs and gaps to fix before scaling

| # | Issue | Location | Severity |
|---|---|---|---|
| 1 | `extract_metadata` is a stub that writes empty fields | [process_corpus.py:476-491](process_corpus.py#L476-L491) | High — no bibliographic data flows downstream |
| 2 | `chunk_text` is a naive 1000-char character split — ignores `config.yaml`, breaks across sentences, words, headings | [process_corpus.py:494-526](process_corpus.py#L494-L526) | High — bad chunks → bad RAG |
| 3 | `config.yaml` exists but is not loaded anywhere | repo root | Medium — misleading; values are hard-coded |
| 4 | OpenAI failures silently fall back to zero vectors, poisoning the index | [embed_chunks.py:60-63](embed_chunks.py#L60-L63) | High — silent data corruption |
| 5 | Figures and captions are not structurally linked — `figures.json` has a `caption` field but no page/bbox, and captions-as-text blocks can appear in `text.json` disconnected from their figure | [process_corpus.py:318-473](process_corpus.py#L318-L473) | High — defeats the "figures of X species" goal |
| 6 | Serial processing, per-step exceptions swallowed into a summary log; no resumability *within* a PDF, only across PDFs | [process_corpus.py:183-228](process_corpus.py#L183-L228) | High at scale |
| 7 | 8-char SHA-256 prefix: collision probability is low but nonzero at 2000 items and bites forever if it happens. Use 12 chars and compare full hash on reuse | [process_corpus.py:120](process_corpus.py#L120) | Low but cheap to fix |
| 8 | OCR always calls `ocrmypdf` with fixed flags — no language detection, English-only Tesseract by default. Old German/French/Russian text will OCR poorly. | [process_corpus.py:283-315](process_corpus.py#L283-L315) | High for historical corpus |
| 9 | Two parallel implementations (Snakefile + scripts/ vs. process_corpus.py) with different output layouts — confusing, invites drift | repo root + scripts/ | Medium |
| 10 | `figures.json` lacks page number and bbox for fallback (PyMuPDF) figures — can't be cross-referenced with visualizations | [process_corpus.py:438-448](process_corpus.py#L438-L448) | Medium |
| 11 | No per-paper error/log file; `print` only. At 2000 papers you can't diagnose individual failures after the fact. | all | Medium |

## 2. Is this the right toolchain?

A candid reassessment for the stated goals.

### Docling
Works well on modern born-digital papers. Struggles on: plate-heavy 19th-century monographs (figures separate from captions), scanned Fraktur German, two-column taxonomy journals with irregular layouts. **Keep for modern papers; don't assume it will carry the historical tail.**

### ocrmypdf + Tesseract
Usable if given the right language packs. Needs: `tesseract-lang` (deu, fra, rus, lat), and specifically `deu_latf` / Fraktur models for 19th-c. German. For very old or damaged scans, Tesseract is weak — **Kraken** (historical documents, fraktur, Cyrillic) or a vision LLM may do better.

### Grobid (currently stubbed)
Still the best general-purpose scientific-PDF metadata/reference extractor. Run it as a Docker service and call it per PDF. It also extracts references and section structure, which is foundational for the citation-analysis goal. **Re-enable this now, don't ship without it.**

### Alternatives worth piloting
- **Marker** (surya + texify) — often better than docling on mixed-quality PDFs, native markdown output with structure.
- **Nougat** (Meta) — good on modern arxiv-like PDFs but trained on that distribution; poor on historical.
- **Claude / GPT-4o vision** on page images — expensive per page, but dramatically better on hard pages (fraktur, plates, hand-lettered captions). Reasonable as a **fallback for pages that fail structural parsing**, not as the default.
- **BHL (Biodiversity Heritage Library) API** — a sizable fraction of historical siphonophore literature (Haeckel's Challenger report, early Agassiz, Huxley, Claus, Chun, Moser, Totton) is already in BHL with OCR and bibliographic metadata. **Query BHL before OCRing a scan yourself.** This alone could save hundreds of papers' worth of OCR effort and provide cleaner ground truth.
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
         ├─ taxonomic names        (gnfinder + WoRMS / GBIF verification)
         ├─ geographic mentions    (spaCy / GLiNER + GeoNames)
         ├─ anatomy terms          (domain lexicon — nectophore, bract, palpon, ...)
         └─ figure classification  (what anatomical structure is depicted? — vision LLM)

Stage 5  Indexing (all outputs in open, file-based formats)
         ├─ vector index (LanceDB)               — semantic / topic retrieval
         ├─ taxon mention index (parquet)        — accepted-name resolved, chunk offsets
         ├─ WoRMS taxonomic backbone (sqlite)    — genera / species / synonymy / authorities
         ├─ bibliographic index (parquet)        — authors, year, title, journal, refs (Grobid + CSL-JSON export)
         ├─ reference graph (parquet)            — citations between papers
         ├─ section-type index (parquet)         — chunk → section class (from Grobid TEI)
         ├─ anatomy-term index (parquet)         — siphonophore anatomy NER offsets
         ├─ geographic index (parquet)           — location mentions + coordinates
         ├─ figure index (parquet + PNGs)        — captions, taxon tags, anatomy tags
         └─ trait / character matrix (parquet)   — species × character (deferred; see §8 Q3)

Stage 6  MCP tools layer
         └─ MCP server exposing tools over the indices (see "Endpoint" below)
```

### Corpus-level inputs (shared across all papers)

Not every input comes from the PDF directory. The pipeline also depends on a small set of **corpus-level reference data** that must be materialized alongside the per-paper artifacts:

- **WoRMS snapshot** — the authoritative backbone of accepted siphonophore names, synonymies, and authority strings. Ingested as a SQLite/parquet resource once, refreshed periodically. Consumed by every taxon-keyed query (§8 Q1/Q2/Q4/Q7) and required by the taxon-mention resolver.
- **Siphonophore anatomy lexicon** — curated terminology (nectophore, bract, palpon, gonophore, pneumatophore, nectosome, siphosome, cnidoband, tentacle, …) used by the anatomy-term NER pass and the figure-anatomy classifier. Maintained as a YAML file in the repo.
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
  "caption_text": "Fig. 3. Nectophore of Nanomia cara ...",
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
- Tools can return **exhaustive** result sets (e.g., all 180 chunks mentioning *Bargmannia*) — the LLM does the synthesis, not the retrieval layer. Top-K truncation is a per-tool opt-in, not a default.
- Outside-lab tools can use the exact same corpus: either via MCP, or — for non-LLM use cases — by reading the underlying parquet/SQLite/LanceDB files directly.

**Initial tool surface:**

| Tool | Purpose | Backs queries |
|---|---|---|
| `search_taxon(name)` | resolve a string to WoRMS accepted name + synonyms + rank | all taxon-keyed |
| `get_chunks_for_taxon(accepted_name, include_synonyms=True)` | exhaustive chunk list | Q1, Q2, Q5 |
| `get_chunks_for_topic(query, k=None)` | semantic vector search; `k=None` returns all above a threshold | Q6, Q8 |
| `get_figures_for_taxon(accepted_name)` | figure records incl. image paths, captions, bboxes | Q2, Q4 |
| `get_figures_for_anatomy(term)` | figures tagged with an anatomy term | Q8, Q4 |
| `get_papers_by_author(name)` | author-filtered paper list | Q5 |
| `get_chunks_by_section(paper_hash, section_class)` | chunks of a given section type ("description", "embryology", …) | Q2, Q6 |
| `get_bibliography(paper_hash)` | refs + citations from Grobid TEI | Q2 |
| `aggregate_taxonomic_acts(rank, group_by="decade")` | analytic rollup over authority strings | Q7 |
| `list_valid_species_under(taxon)` | WoRMS children filter | Q4 |

The server is thin: each tool is a small function over the parquet/SQLite/LanceDB indices. The storage layer is intentionally format-neutral so non-MCP consumers (a Jupyter notebook, a Snakemake rule, another lab's pipeline) can query the same corpus directly via DuckDB / SQLite / LanceDB without going through the MCP layer at all.

## 4. Concrete next steps, ordered

### Immediate (fix what's there)

- [ ] **Stop writing zero-vector embeddings on failure.** Raise or mark the chunk as failed and skip insertion. [embed_chunks.py:60-63](embed_chunks.py#L60-L63)
- [ ] **Load `config.yaml`** in `process_corpus.py` and actually use its values (OCR flags, chunk tokens, model name). Or delete it. Pick one.
- [ ] **Wire up Grobid.** Add a `docker-compose.yml` with the Grobid image; replace the stub `extract_metadata` with a client call to `/api/processHeaderDocument` (header) and `/api/processReferences` (bibliography). Cache TEI-XML output alongside `metadata.json`. Preserve section structure — downstream stages depend on it.
- [ ] **Replace naive `chunk_text` with docling's `HybridChunker`** (tokenizer-aware, respects structure) — the dependency is already installed.
- [ ] **Language-aware OCR.** Call `ocrmypdf -l eng+deu+fra+rus+lat` as default; detect language per page with langdetect or from existing text density, and narrow the lang list when confident. Install Tesseract lang packs in setup.
- [ ] **Make OCR optional per page, not per document.** Many mixed PDFs have cover pages with text and body pages that are scanned. `ocrmypdf --skip-text` covers this; prefer it over `--force-ocr`.
- [ ] **Bbox + page in `figures.json`** for both docling and PyMuPDF paths. [process_corpus.py:373](process_corpus.py#L373) `append_figure`.
- [ ] **Per-paper log file** at `documents/<HASH>/pipeline.log`. Migrate `print` calls to `logging`.
- [ ] **Use 12-char hash prefix** and verify no collision before reusing a directory.
- [ ] **Switch embeddings to local open-weights on Bouchet GPU** (see §7).

### Near-term (unlock scale & the goals)

- [ ] **Ingest WoRMS taxonomic backbone as a corpus-level input.** Download a WoRMS snapshot (order Siphonophorae + all child taxa, with synonymies and authority strings) and store as SQLite/parquet alongside the per-paper artifacts. This is a **prerequisite**, not an annotation step: Q1/Q2/Q4/Q7 and the taxon-mention resolver all depend on it. Refresh on a schedule; pin version in the corpus manifest.
- [ ] **Curate the siphonophore anatomy lexicon.** A YAML file in the repo listing anatomy terms (nectophore, bract, palpon, gonophore, pneumatophore, nectosome, siphosome, cnidoband, tentacle, stem, …) plus common synonyms and language variants. Used by both the anatomy-term NER and the figure-anatomy classifier.
- [ ] **Figure+caption linking** as described in §3, plus a human-reviewable `figures_report.html` per paper (thumbnails + captions side-by-side). This is your QC surface for the "figure extraction success" goal.
- [ ] **Parallelize stage 2.** Either bring back Snakemake (now targeting the hash-based layout) or use a simple `concurrent.futures.ProcessPoolExecutor` with a file-based lock per hash dir. Snakemake is a good fit at 2000 items — it gives free DAG tracking, retries, cluster submission, and plays naturally with the per-hash file markers you already emit.
- [ ] **BHL lookup stage** before OCR. If a BHL match exists (fuzzy title + year + author), prefer BHL's OCR text as a parallel artifact and flag confidence. Likely to materially help historical papers.
- [ ] **Golden test set.** Pick ~10 papers spanning: (a) modern born-digital, (b) modern scan, (c) 19th-c. German monograph (a Haeckel volume is a good stress test), (d) Russian-language paper, (e) French paper, (f) plate-heavy paper with end-matter figures. Snapshot-test extraction quality. Re-run on every pipeline change.
- [ ] **Deprecate the Snakefile legacy path** once stage 1 reaches parity, or rewrite it against the hash layout. Don't leave two divergent implementations in place.

### Domain-specific (annotation layer + query layer for the stated downstream uses)

- [ ] **Taxonomic mention extraction + resolution.** Integrate `gnfinder` (Global Names — handles historical variants well) for raw detection; resolve every mention against the WoRMS backbone (above) to an accepted name. Emit `taxa.json` per paper with chunk offsets, accepted name, and match confidence. Roll up into a corpus-wide `taxon_mentions.parquet`.
- [ ] **Authority-string parsing + taxonomic-act index.** Parse "Author, Year" authority strings from WoRMS and from in-corpus descriptions. Emit a `taxonomic_acts.parquet` (taxon, author, year, original-description paper hash when resolvable). Backs Q7 directly.
- [ ] **Geographic extraction.** GLiNER or spaCy + GeoNames/Nominatim. Store as `locations.json` per paper; roll up into a corpus-wide `locations.parquet` keyed by taxon co-occurrence. Backs Q1.
- [ ] **Anatomy-term NER.** Apply the curated lexicon to every chunk; emit offsets. Start with exact-match + stemming; extend to a small fine-tuned NER if recall is poor on historical text.
- [ ] **Anatomy-aware figure search.** Start with caption-text keyword + taxon co-occurrence. For deeper coverage, run a vision model (Claude or an open VLM) over figure crops to tag anatomical structures. Highest-leverage place to use an LLM: captions are short, visual classification is cheap per image, directly backs Q4/Q8.
- [ ] **Reference graph.** From Grobid bibliographies, build a `references.parquet` + in-corpus citation graph. Backs Q2 and opens the door to citation-based queries.
- [ ] **Section-type index.** Use Grobid TEI structural labels to tag chunks with a section class (description / distribution / embryology / discussion / remarks / …). Backs Q2 and Q6.
- [x] **MCP server.** Implemented in [mcp_server.py](mcp_server.py) with 13 tools backed by the per-PDF JSON indexes: `list_papers`, `get_paper`, `search_taxon`, `get_papers_for_taxon`, `get_chunks_for_taxon`, `get_figures_for_taxon`, `get_figures_for_anatomy`, `get_chunks_by_section`, `get_bibliography`, `get_papers_by_author`, `list_valid_species_under`, `get_figure`, `get_chunk`. Built on FastMCP, stdio transport, eager in-memory index at startup. `get_chunks_for_topic` (vector-search-backed) is deferred until Phase E swaps in the local embedding backend.
- [ ] **Trait extraction for identification — DEFERRED.** Backs Q3 (keys). Requires a siphonophore-specific character schema and structured LLM extraction per species description. Substantial enough to warrant its own plan section; captured here so it isn't forgotten. Not on the path for the first production run.

## 5. Scale considerations for ~2000 papers

- Rough cost envelope: docling is ~30s–3min/paper on CPU, OCR adds 1–5min/paper for scans. At 2000 papers, plan on **2–5 days of wallclock on a single machine** for stage 1–2 serial; hours with parallelism.
- Embeddings: running a local open-weights model on Bouchet GPUs (see §7) eliminates the OpenAI cost and rate-limit path entirely. Full corpus (~300M tokens) embeds in well under an hour on a single GPU; no API bill. The motivation for switching is reproducibility and removing the silent zero-vector failure mode, not cost — OpenAI would have been ~$6.
- Grobid is a Java server; batch it via its concurrency API rather than spawning per-PDF.
- Keep intermediate artifacts on disk, not in memory — the current file-per-step layout is correct, just add per-paper logs.
- Add a corpus-level `manifest.parquet` or SQLite that rolls up one row per hash (status, pages, language, taxa count, figure count, errors). This is what you'll actually query when inspecting corpus health.

## 6. What to do first

If picking only three items from this plan, start with:

- [ ] **Grobid integration** — unlocks every downstream goal that needs bibliographic, section-structure, or citation data.
- [ ] **WoRMS backbone ingest + figure+caption joint object** — together these back the "nectophore of X species across all papers" goal and every taxon-keyed query in §8.
- [ ] **Golden test set + per-paper logs** — without these, every subsequent change is a leap of faith across 2000 inputs.

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

### Code changes required

[embed_chunks.py](embed_chunks.py) is small; the diff is contained.

1. Replace OpenAI client with `sentence_transformers.SentenceTransformer` (or `FlagEmbedding.BGEM3FlagModel` for BGE-M3-native features). Load once at process start, move to CUDA.
2. Change embedding dim: `Vector(1536)` → `Vector(1024)` (BGE-M3). Any existing LanceDB index has to be rebuilt — it's stage 2 output, cheap to regenerate.
3. Replace per-batch-of-100 API calls with a single `.encode(chunks, batch_size=64, normalize_embeddings=True, show_progress_bar=True)` call per PDF (or across PDFs, if we batch globally).
4. **Fix bug #4 while we're here:** embed failures must raise, not zero-vector. With a local model, failure means OOM or driver issue — both should stop the run loudly.
5. Gate the model/device behind config so we can fall back to OpenAI or swap models without editing code. Finally consume `config.yaml` (bug #3).

### SLURM submission

New `batch_embed.sh` alongside `embed_chunks.py`:

```bash
#!/bin/bash
#SBATCH --job-name=corpus-embed
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=4:00:00
#SBATCH --output=logs/slurm-%j.out

module purge
module load miniconda CUDA
conda activate corpus
python embed_chunks.py <output_dir>
```

Stage 1 (`process_corpus.py`) stays CPU-only and can run on `day` — no reason to hold a GPU during OCR/docling.

### Rollout order

- [ ] Add the local-embedding path behind a `--backend {openai,local}` flag; keep OpenAI working.
- [ ] Run on the existing demo corpus; compare retrieval quality on a handful of query→expected-chunk pairs.
- [ ] If acceptable, flip the default to `local` and drop the OpenAI dependency from `requirements.txt` (or leave as optional extra).
- [ ] Execute the full production run on Bouchet.

## 8. Target queries & end-user interface

This section defines what the repo is *for*. Every design choice above should be justifiable against the queries below. If a decision doesn't help serve at least one of these, it's speculative.

### Framing: the endpoint is not a RAG app

Top-K vector retrieval is suited to "answer a question using the best few passages." Most of the target queries below are **enumerative** (find everything matching criteria, then synthesize), **analytic** (aggregate over structured metadata), or **trait-structured** (produce a comparative table). Top-K would either truncate the evidence or return irrelevant near-neighbors.

The corpus therefore produces **a set of indices**, and the user-facing layer is a **query-routing + synthesis layer** that:

- picks the right index (or combination) for the query shape,
- pulls **exhaustive** matches from structured indices (not top-K), then
- hands the full evidence set to an LLM (Opus, 1M-token context) for synthesis, or to a plotting/table function for analytic queries.

Natural implementation shape: an agent with tools over the indices (`search_taxon`, `get_chunks_for_taxon`, `get_figures`, `get_bibliography`, `aggregate_by_year`, `run_sql`), not a single `chat.py` backed by vector search. Vector search is one tool among several, used for semantic/topic queries where the match criterion is fuzzy.

### The eight target queries

Each is traced as **Shape → Path → Output → Gaps in current plan**.

#### Q1. "List all collection locations of *Agalma elegans*."

- **Shape:** enumerative, entity-keyed, structured output.
- **Path:** resolve "Agalma elegans" to accepted WoRMS taxon (+ synonyms) → taxon index returns all chunks/captions/tables mentioning the species → geographic NER on those chunks → dedupe + group by paper.
- **Output:** table (location, coordinates if present, paper, page, quote).
- **Gaps:** taxon↔chunk offset index (§4 *Taxonomic mention extraction + resolution* — must be exhaustive, not retrieval-ranked); geographic NER at chunk scope tied to taxon co-occurrence (§4 *Geographic extraction*); tabular output recipe (not yet scoped).

#### Q2. "Compose a monographic review of *Apolemia*."

- **Shape:** enumerative, genus-level, long-form synthesis.
- **Path:** genus → all child species (incl. synonyms) via WoRMS → all chunks mentioning any of them → section-filtered retrieval (descriptions, habitat, distribution, remarks) via Grobid TEI → all figures + captions → all references → Opus synthesis in one or a few passes.
- **Output:** structured long-form document with inline citations, figures, and a reference list.
- **Gaps:** WoRMS backbone as a first-class pipeline artifact (§4 *Ingest WoRMS taxonomic backbone*); section-type tagging of chunks from Grobid TEI (§4 *Section-type index*); synthesis recipe for monographic output (not yet scoped); reliable figure+caption joint object (§4 *Figure+caption linking*).

#### Q3. "Make a key to identify *Forskalia* species."

- **Shape:** trait-structured, comparative.
- **Path:** genus → all valid species → per-species diagnostic traits (meristics, counts, dimensions, colors, structural presence/absence) from species-description sections → pivot into species × character matrix → generate dichotomous (or polyclave) key.
- **Output:** a usable identification key.
- **Gaps:** §4 lists *Trait extraction for identification* as **deferred**. A key requires (a) a siphonophore-specific character schema, (b) structured LLM extraction against that schema per species description, (c) a trait index as a first-class artifact, (d) key-generation logic. Substantial enough to warrant its own plan section when we return to it; out of scope for the first production run.

#### Q4. "List all currently valid siphonophore species with a one-paragraph summary and diagnostic figures."

- **Shape:** enumerative over a closed taxonomic set; per-entity synthesis; figure selection.
- **Path:** WoRMS → list of currently-valid siphonophore species (~200) → per species: all mentioning chunks → Opus one-paragraph summary → figure index: find figures captioned/tagged with the species and a diagnostic anatomy term → assemble.
- **Output:** a catalog document / website.
- **Gaps:** figure index with **anatomy tags** (§4 *Anatomy-aware figure search* — needs vision-model figure classification); a per-species batch synthesis recipe (not yet scoped); concept of "diagnostic structure" — either user-curated per genus or learned from description sections.

#### Q5. "Summarize all of Phil Pugh's comments about Haeckel."

- **Shape:** author-filtered + entity-filtered; meta-commentary extraction.
- **Path:** author index (from Grobid headers) → all papers authored by Pugh → chunks mentioning "Haeckel" (as person, not just cited work) → LLM filter: which mentions are substantive commentary vs. routine citation → synthesis.
- **Output:** narrative summary with direct quotes + citations.
- **Gaps:** author-indexed corpus (falls out of §4 *Wire up Grobid* but needs explicit per-author rollup); person-name NER that distinguishes cited-author-as-entity from cited-work; "meta-commentary vs. citation" classification (LLM in the loop at query time is fine — don't pre-compute).

#### Q6. "Write a summary of siphonophore embryology."

- **Shape:** topic-semantic, corpus-wide.
- **Path:** (a) section-type index: chunks from sections titled/classified "development / embryology / larval stages" + (b) vector search for "embryology" and related concepts across all chunks → union → Opus synthesis with citations.
- **Output:** narrative summary with references.
- **Gaps:** this is the query the vector index was built for. Also benefits from Grobid-provided section labels. No new scope items — just needs what's already planned to work well.

#### Q7. "Plot number of siphonophore species described per decade."

- **Shape:** analytic aggregation over structured metadata. No LLM required.
- **Path:** for every currently-valid siphonophore species in WoRMS: parse authority string ("Eschscholtz, 1829" → 1829) → cross-check against the paper in our corpus where available → aggregate by decade → plot.
- **Output:** a figure.
- **Gaps:** authority-string parsing and taxonomic-act index (not in current plan). Much of the data comes from WoRMS directly — the corpus only needs to supply the original description PDFs for cross-reference. Consider whether this query even requires the corpus or can be answered from WoRMS alone; if the latter, it's a good sanity-check baseline.

#### Q8. "Summarize what is known about the structure and function of pneumatophores."

- **Shape:** topic-semantic, anatomy-keyed.
- **Path:** anatomy index returns all chunks tagged with "pneumatophore" → plus vector search for adjacent concepts (pneumatosaccus, gas gland, float, buoyancy) to catch non-lexical matches → figures tagged with pneumatophore → Opus synthesis.
- **Output:** narrative summary with citations and figures.
- **Gaps:** anatomy-term index (§4 *Anatomy-term NER*); vision-based figure anatomy tagging (§4 *Anatomy-aware figure search*).

### Indices this implies

Vector embeddings are one of ten indices. Only Q6 and Q8 treat vector search as a primary retrieval path; the others use it only as a supplementary recall mechanism, if at all.

| # | Index | Produced by | Used by |
|---|---|---|---|
| 1 | Taxon mention index (chunk-offset) | gnfinder + WoRMS resolution | Q1, Q2, Q4, Q5 (indirectly), Q7 |
| 2 | WoRMS taxonomic backbone (genera/species/synonyms/authorities) | WoRMS snapshot | Q1, Q2, Q3, Q4, Q7 |
| 3 | Bibliographic index (author, year, title, journal) | Grobid headers + CrossRef/OpenAlex | Q2, Q5, Q7 |
| 4 | Reference graph | Grobid bibliography | Q2 |
| 5 | Section-type index (descriptions, embryology, etc.) | Grobid TEI structure | Q2, Q6 |
| 6 | Anatomy-term index | domain lexicon + NER | Q1 (anatomy qualifier), Q3, Q8, Q4 (diagnostic) |
| 7 | Geographic index | NER + GeoNames | Q1 |
| 8 | Figure index (image + caption + taxon + anatomy tags) | docling + vision LLM classifier | Q2, Q4, Q8 |
| 9 | Trait / character matrix (species × character) | structured LLM extraction against schema | Q3 |
| 10 | Vector index (semantic fallback) | local embeddings (§7) | Q6, Q8 (primary); others (recall boost) |

### Coverage of these queries by the current plan

The scope items these queries demand — WoRMS backbone as a first-class input, anatomy lexicon + NER, authority-string / taxonomic-act parsing, section-type tagging, figure anatomy classification, reference graph, and the MCP tools layer — are now all reflected as line items in §4. The two genuine gaps remaining are:

- **Trait extraction + key generation (Q3).** Deferred; noted in §4 and called out in Q3's gaps line above. Will get its own plan section when we pick it up.
- **Synthesis recipes (Q2, Q4).** The indices + MCP tools will be in place, but the specific Opus-driven prompts and map-reduce patterns for monographic review (Q2) and per-species catalog (Q4) are not yet scoped. These are cheap to iterate on once the substrate exists — prototype against the golden test set before committing to any one recipe.

## 9. Three-pass figure extraction with ROI sidecar

The initial figure pipeline (Phase D, §3 "Figure+caption as a first-class object") handles the common case — docling extracts one `Picture` per real figure and we classify + dedupe. Two harder cases surface in the demo set and drive this design:

1. **Adjacent figures merged into one image** — docling sometimes extracts two logically-separate figures as a single bbox when there's insufficient whitespace or heading text between them. Pugh 2001's Figures 3 (panels A–B) and 4 (panels A–D) come out as one image whose caption is parsed as `Fig. 3`.
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
  "caption_text": "Fig. 3 A. upper and B. lower views of young nectophore. Scale bar 5 mm.",
  "rois": [
    {"type": "figure", "figure_number": "3",
     "roi_px": [0, 0, 500, 1305],
     "caption_text": "Fig. 3 A. upper and B. lower views of young nectophore.",
     "source": "ocr:tesseract"},
    {"type": "figure", "figure_number": "4",
     "roi_px": [500, 0, 1015, 1305],
     "caption_text": "Fig. 4 Erenna richardi. A. first and B., C. second types of bract…",
     "source": "ocr:tesseract"},
    {"type": "panel", "parent_figure_number": "3", "label": "A",
     "roi_px": [0, 0, 500, 600],
     "description_from_caption": "upper view of young nectophore",
     "source": "ocr:tesseract"},
    {"type": "panel", "parent_figure_number": "3", "label": "B",
     "roi_px": [0, 600, 500, 1305],
     "description_from_caption": "lower view of young nectophore",
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
- **Missing-figures scan**: scan `text.json` for `Fig. N` / `Figure N` mentions; any N that's not present as a `figure_number` in `figures.json` gets recorded in a per-paper `missing_figures[]` list with the caption text (parsed out of the running text). Closes the Pugh 2001 Fig 4 observability gap even without Pass 3.

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

### Rollout

Implemented in three commits:

1. **Pass 2.5** (commit `a794bdc`) — caption parser + missing-figures scan. Cheap, always-on. Adds `panels_from_caption[]` / `panel_count_from_caption` / `missing_figures[]` to artifacts. No image or filename changes.
2. **Pass 3a** (commit `72e0633`) — Tesseract ROI pass, MCP crop-on-demand tools. Gated behind `--content-aware-figures` flag. Recall turned out to be ~20–40% on line-art scientific figures — always-partial, rarely the primary signal.
3. **Pass 3b** (commit `9e98455`) — Claude-vision backend, same result schema as 3a. Flagged via `--vision-backend claude` (+ `--vision-model …` override). On the demo set: 10/11 Pugh completed + 1 `completed_compound` (fig_3 → the Fig 3+4 case), 5/5 Siebert completed. Haiku 4.5 at ~$0.003/fig.

### Status & immediate next steps (as of commit `9e98455`)

**Done**: Phases A–F, MCP server (15 tools), filename-year fallback, figure dedup + position letters, Pass 2.5, Pass 3a, Pass 3b (Claude vision). All work committed; demo/ fully reprocessed.

**Three natural follow-ups**, in rough order of user-facing value:

- [ ] **Compound file rename + sub-figure split** (Phase 3c). When `pass3_status = completed_compound`, attribute the extra panels to a `missing_figures[]` entry, rename the file to `fig_M-N.png`, record `previous_filenames[]`, and emit a second logical figure record so MCP tools can list "Fig 4" separately even though its image is embedded in `fig_3-4.png`. Small — maybe 150 lines.
- [ ] **Local VLM backend** (`vision.LocalVLMBackend` with Qwen2.5-VL-7B on CUDA/MPS). API contract already set; just implement `detect_figure_panels()` via HuggingFace transformers. Needed for the Bouchet production run so no per-figure API cost at the 20k-figure scale.
- [ ] **Bouchet prep** — SLURM batch scripts (`batch_process_corpus.sh`, `batch_pass3b.sh` on `gpu_h200`, `batch_embed.sh` on `gpu`), Grobid via Singularity, model cache pre-download into `~/project/cache/huggingface/`, corpus `git lfs pull` to `~/project`, dry-run on 20-50 papers before the full 2000. Pair with the local-VLM backend above.

**Fourth item deferred indefinitely**: `get_figure_roi_image` panel disambiguation when two A's exist in a compound (minor MCP polish).

Pass 3c handles compound file rename; Phases A–F + 2.5 + 3a + 3b are otherwise considered complete as of this checkpoint.
