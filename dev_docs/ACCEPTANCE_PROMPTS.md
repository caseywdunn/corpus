# Siphonophore Corpus — MCP Acceptance Prompt Checklist

A set of realistic natural-language prompts for manual (or Claude-assisted) acceptance
testing of the siphonophore corpus MCP server.  Run these in order against the live
server after a new corpuscle build to confirm that all databases, indexes, and tool
paths are working correctly.

Each prompt names the **primary tool(s)** it exercises and lists observable **pass
criteria** — things a human tester can confirm without domain expertise.

> **Tip:** connect Claude to the MCP server (e.g. via the Claude desktop app or
> `claude --mcp-config`) and paste each prompt verbatim.  A production corpus built
> from 1 700+ papers should satisfy all pass criteria.

---

## 1  Corpus orientation

**Prompt:**
> How many papers are in this corpus?  How many text chunks?  Give me a one-paragraph
> summary of the bundle.

**Tools:** `bundle_info`, `corpus_summary`

**Pass criteria:**
- Paper count is in the right order of magnitude (≥ 100 for any real build).
- Chunk count is reported and is larger than the paper count.
- Bundle name / version present in the response.

---

## 2  List and paginate papers

**Prompt:**
> List the ten oldest papers in the corpus by year.

**Tools:** `list_papers`

**Pass criteria:**
- Returns exactly 10 papers with titles, authors, and years.
- Earliest year is plausibly old (pre-1900 papers exist in this corpus).

---

## 3  Taxonomy — species lookup

**Prompt:**
> What is the accepted name and taxonomic rank of *Physalia physalis*?  Who described
> it originally, and in what year?

**Tools:** `search_taxon`, `get_original_description`

**Pass criteria:**
- Accepted name returned (likely *Physalia physalis* or a synonym note).
- Rank is "Species".
- `get_original_description` returns an authorship string containing a year.

---

## 4  Taxon dossier — workhorse tool

**Prompt:**
> Give me a full dossier on *Nanomia bijuga*: how many corpus papers mention it, which
> papers have the most mentions, what anatomy terms co-occur most often, and whether
> any figures are associated with it.

**Tools:** `get_taxon_dossier`

**Pass criteria:**
- Response contains `taxon`, `papers`, `lexicon_aggregated`, and `figure_index` sections.
- At least one paper is listed.
- Anatomy terms appear if the lexicon is configured (category `anatomy` present).

---

## 5  Species inventory

**Prompt:**
> List all valid species in the order Siphonophorae that appear in this corpus.

**Tools:** `list_valid_species_under`, `get_papers_for_taxon`

**Pass criteria:**
- Returns a list of species-rank taxa; at least several dozen entries for a full build.
- No Python errors; "not found" is acceptable on a demo corpus.

---

## 6  Subtree dossier

**Prompt:**
> Summarise the literature coverage for the family Physaliidae: how many species, how
> many papers, and which species has the most corpus mentions?

**Tools:** `get_taxon_subtree_dossier`

**Pass criteria:**
- Response contains species list and per-species paper/mention counts.
- Top species is plausibly *Physalia physalis*.

---

## 7  Semantic search

**Prompt:**
> Find the five most relevant text passages about nectophore swimming kinematics.

**Tools:** `get_chunks_for_topic`

**Pass criteria:**
- Returns 5 chunks (or fewer if corpus is small).
- Each chunk contains a `text` field with multiple sentences.
- Passages are thematically related to swimming / nectophores.

---

## 8  Section-based reading

**Prompt:**
> Show me the Methods section of the paper by Dunn et al. on siphonophore phylogenetics.

**Tools:** `list_papers`, `get_chunks_by_section`

**Pass criteria:**
- Correct paper hash is retrieved via `list_papers`.
- `get_chunks_by_section` returns chunks with section headings containing "Method" or
  "Material".
- Text is coherent prose (not garbled OCR).

---

## 9  Citation list for a paper

**Prompt:**
> Give me the full reference list for the Huxley 1859 paper, and tell me which of
> those cited works are also in this corpus.

**Tools:** `list_papers`, `get_bibliography`

**Pass criteria:**
- Reference list is non-empty.
- Each reference has an author string and (usually) a year.
- The `in_corpus` boolean is present on each reference.
- At least some references have `in_corpus: true` in a full build.

---

## 10  In-text citation context

**Prompt:**
> In the Huxley 1859 paper, show me the first ten in-text citation markers and the
> paragraphs they appear in.

**Tools:** `get_intext_citations`

**Pass criteria:**
- Response contains `citations` and `paragraphs` lists.
- Each citation has a `surface` string (the cited text as it appears in the paper).
- `para_index` values correctly index into the returned `paragraphs` list.

---

## 11  Cross-corpus citation retrieval

**Prompt:**
> Find every passage in the corpus that cites Huxley's 1859 paper on siphonophores.
> Show me the citing paper titles and the surrounding text.

**Tools:** `get_works_by_author`, `get_excerpts_citing`

**Pass criteria:**
- `get_works_by_author("Huxley")` returns at least one result.
- `get_excerpts_citing` returns passages with `citing_paper_title` and `paragraph` fields.
- At least one excerpt text mentions "Huxley" or a related concept.

---

## 12  Citation graph

**Prompt:**
> How many times is the Dunn 2005 siphonophore phylogeny paper cited within this corpus?
> Which papers cite it most?

**Tools:** `get_citation_graph`

**Pass criteria:**
- Response includes a `cited_by` count.
- `citing_papers` list is non-empty if the paper is present in the corpus.

---

## 13  Anatomy lexicon co-occurrence

**Prompt:**
> Show me the anatomy term co-occurrence matrix for the top 10 terms, limited to the
> 20 most-relevant papers.  Which terms appear together most often?

**Tools:** `lexicon_matrix`

**Pass criteria:**
- Response contains `terms`, `papers`, and `matrix` (or equivalent) fields.
- Matrix has non-zero entries.
- Term names are recognisable siphonophore anatomy vocabulary (e.g. nectophore,
  pneumatophore, gastrozooid).

---

## 14  Lexicon term dossier

**Prompt:**
> Give me a dossier for the anatomy term "nectophore": which papers use it most, in
> what sections, and are there associated figures?

**Tools:** `get_lexicon_term_dossier`

**Pass criteria:**
- Response contains a papers list and section breakdown.
- Figure list present (may be empty on demo corpus).
- Top papers by mention count seem domain-appropriate.

---

## 15  Figures for a taxon

**Prompt:**
> Show me figures associated with *Physalia* across the corpus.  For the first result,
> describe what the figure caption says.

**Tools:** `get_figures_for_taxon`, `get_figure`

**Pass criteria:**
- `get_figures_for_taxon` returns at least one figure record with `paper_hash`,
  `figure_id`, and a caption excerpt.
- `get_figure` retrieves full metadata for that figure including `caption` text.

---

## 16  Figure dossier (visual survey)

**Prompt:**
> Give me a visual survey of *Nanomia bijuga* figures across the corpus: thumbnails or
> captions grouped by figure type (diagram, photo, etc.).

**Tools:** `get_figure_dossier_for_taxon`

**Pass criteria:**
- Response groups figures by `figure_type`.
- At least one figure type group is non-empty.
- Caption text is legible and mentions the taxon name.

---

## 17  Output profiles

**Prompt:**
> What output profiles are configured for this server?  Which one is currently active?

**Tools:** `list_output_profiles`, `get_active_profile`

**Pass criteria:**
- At least one profile listed.
- Active profile name returned.

---

## 18  End-to-end: researcher workflow

**Prompt:**
> I'm researching the evolution of the nectophore across siphonophore families.
> Point me to the most relevant papers, the key anatomy terms that co-occur with
> "nectophore", and any figures that show nectophore morphology.

**Tools:** `get_chunks_for_topic`, `get_taxon_dossier` (on Siphonophorae),
`lexicon_matrix`, `get_figures_for_lexicon_term`

**Pass criteria:**
- Semantic search returns coherent passages about nectophore biology.
- Anatomy lexicon highlights terms like "nectosac", "velum", "ostium".
- At least one figure caption mentions morphological features.
- Response reads as a useful research starting point, not a list of errors.

---

## Running as a checklist

For each prompt, record:

| # | Tool(s) hit | Pass? | Notes |
|---|-------------|-------|-------|
| 1 | bundle_info, corpus_summary | ☐ | |
| 2 | list_papers | ☐ | |
| 3 | search_taxon, get_original_description | ☐ | |
| 4 | get_taxon_dossier | ☐ | |
| 5 | list_valid_species_under | ☐ | |
| 6 | get_taxon_subtree_dossier | ☐ | |
| 7 | get_chunks_for_topic | ☐ | |
| 8 | get_chunks_by_section | ☐ | |
| 9 | get_bibliography | ☐ | |
| 10 | get_intext_citations | ☐ | |
| 11 | get_works_by_author, get_excerpts_citing | ☐ | |
| 12 | get_citation_graph | ☐ | |
| 13 | lexicon_matrix | ☐ | |
| 14 | get_lexicon_term_dossier | ☐ | |
| 15 | get_figures_for_taxon, get_figure | ☐ | |
| 16 | get_figure_dossier_for_taxon | ☐ | |
| 17 | list_output_profiles, get_active_profile | ☐ | |
| 18 | (multi-tool) | ☐ | |
