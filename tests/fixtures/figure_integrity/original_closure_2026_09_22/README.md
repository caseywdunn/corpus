These receipts complete the named local acceptance checks for #322, #324 and
#332. They distinguish saved-source replay from fresh extraction. No OCR,
Docling conversion, vision model, embedding job or source-library edit ran.

`captions.json` records five real figure records through current caption/panel
helpers, production bundling, direct/dossier tools and MCP image conversion:

- **#322:** the saved actual Sutherland page-5 Docling capture now binds the
  printed Figure 5 side caption to both extracted plots. The Results paragraph
  remains rejected evidence. Direct lookup and dossier agree on binding and
  `caption_completeness=unverified`. The completed Hosia Figures 4/5 source,
  fragment, scale and serving acceptance remains in
  [the earlier receipt](../hosia_gold_2026_09_18/README.md).
- **#324:** unchanged stored Siebert captions were checked against original
  PDF pages 5, 7 and 16. Figures 4/12 yield A–U/A–N; Figure 2 A/B retain their
  distinct species descriptions and shared in-situ context. Requests for U/N
  return explicit whole-image fallback, with returned PNG bytes identical to
  the bundled whole image and no invented coordinates. The source Sutherland
  caption also retains A/B despite inline `N. bijuga` abbreviations.

The source visual review was performed by Codex team agent `/root/fix_figures`,
not a human domain reviewer. PDF hashes, page regions, crop hashes and original
artifact hashes identify the evidence. The small Siebert caption fixture next
to this directory preserves complete original stored text; its whitespace and
source spelling were not corrected. Original source PDFs and crops are not
vendored. The two named Siebert late-letter cases have no stored pixel ROIs;
this is panel-inventory/fallback acceptance, not #305 geometry acceptance.

Two additional defects exposed by these checks are covered by executable tests:
newly recovered expected labels retire an old `completed` status when those
labels lack pixel ROIs; and undeclared dotted initials after `by`/`with` no
longer destroy the explicit Sutherland A/B inventory. No source text or ROI
coordinates are rewritten. A successful subsequent panel pass retires the
stale coverage note. The figure policy invalidates materialization and crossrefs
without rerunning OCR, extraction or chunks. Existing multilingual, bare-label,
cross-page, author-initial, glossary and cross-figure controls remain covered.
The API catalog now distinguishes full stored caption text from verified source
completeness.

`counts.json` and the complete `count_census.tsv` record **#332** acceptance:

| Available snapshot | Documents | Figure records | Stale stored totals | All three API totals |
| --- | ---: | ---: | ---: | ---: |
| Retained v1.2.1, created 2026-09-01 | 1,775 | 22,171 | 0 | 22,171 |
| Completed frozen gold build | 35 | 605 | 0 | 605 |

These fixed snapshot counts are not a permanent library size. The retained
snapshot is **not** the September-9 v1.4 audit snapshot; this does not reproduce
or reclassify its reported 48 mismatches. All four named papers are present in
the retained census, with their actual earlier snapshot totals recorded.

Every figure artifact was copied into a scratch build and passed through real
`_pass25_annotate_figures`, reading its original text without modifying it.
Production `package` validated totals and made count-only bundles. Real
`CorpusIndex`, `bundle_info`, `corpus_summary`, and every page of `list_papers`
then agreed. The TSV retains each original and rematerialized artifact hash,
stored total and array length. Original figure files stayed byte-identical.
Counts include furniture and records sharing an image: the retained snapshot
has 21,521 unique named image references versus 22,171 records; gold has 392
versus 605. Count-only bundles deliberately omit images/chunks/embeddings and
are not deployable replacement corpuscles.

The count replay started before the supplemental inline-abbreviation change;
its exact intermediate producer identities are retained. The counting code and
count rule are unchanged. Its coverage-status changes are producer observations,
not visually adjudicated missing-panel prevalence. `captions.json` records the
final caption producer separately. Fresh final-producer gold/deployment checks
remain overall release gates.

To reproduce, use the source roots and artifact hashes in the receipts in a
new scratch directory. For Sutherland, load `/tmp/sutherland-p5-docling.json`
as `DoclingDocument`, match the two page-5 pictures to original image records
by source bbox, and run `extract_caption_info`/`classify_figure`. For Siebert,
use the unchanged stored captions and images. Run `_pass25_annotate_figures`,
copy actual metadata/taxa and taxonomy, and call production `package`. Query
the five recorded IDs, compare direct/dossier caption evidence, and request
U/N through `get_figure_roi_image` and MCP `get_figure_image`. No synthetic
caption, taxonomy or pixel geometry should be substituted.

For counts, copy each source `figures.json` and `metadata.json`, run Pass 2.5
with the original read-only `text.json`, validate with
`_count_figures_and_chunks`, package, and sum every `list_papers` page alongside
the two aggregate tools. Temporary scripts and raw outputs were retained under
`/tmp/corpus-original-figure-closure/`; those are historical provenance, not
portable dependencies. A new run must report its own source/producer identities
instead of replacing these dated observations.

Final focused checks: **164 passed** across figure integrity, plate admission,
caption association/fragment recovery, configuration updates, compound-caption
evidence, caption binding/evidence and caption preview tests. These include
append/removal, shared-plate totals, stale-bundle rejection and deterministic
legacy counting. Ruff and whitespace checks pass. The bibliography agent's
independent read-only code review found no blocker; it did not run models or
independently grade the complete source population.
