# v1.5 source acceptance — siphonophore worked example

This records the September 2026 audit against the existing siphonophore gold
corpuscle. It is a reference deployment's validation inventory, not a required
corpus or a new fixture collection. General fidelity methodology remains in
[TESTING.md](TESTING.md).

The existing 35-document set remains the rebuild/scoring fixture. Its manifest
is `transcriptions/sources.json` in the siphonophore library. Presence in that
manifest proves that a source is available for gold scoring; it does not prove
that a new citation, geometry or rights assertion already exists. Hydractinia
is not added as a fixture corpus.

Inventory checked 2026-09-18 at library revision `5ad0164e6840ea16adb6e54a3a5712e20b234feb`; manifest SHA-256
`6efcf5b9f2a5c537e1bebcedcbab7bc730854912a128c14af99f4c59f303e378`.

## Explicit source hashes in the issue reports

This table maps the PDF hashes explicitly named in the in-scope issue bodies.
Issue bookkeeping hashes and HTML audit markers are excluded. An empty cell
means no explicitly named hash in that category, not that the behavior is
untested. Name-only examples and candidate libraries still need source review.
Rows are an inventory, not completion claims.

| Issue | Named source already in gold | Named source outside gold |
| --- | --- | --- |
| #296 | — | Manko_Pugh2018.pdf (`cd108374456b`) |
| #299 | Chun1898b | — |
| #300 | — | delle Chiaje1841Volume5.pdf (`4146f5feabcb`); delle Chiaje1841Volume1.pdf (`437c17ceb0d9`); Moore_etal1953.pdf (`663f82651371`); Moore1953.pdf (`d0b8cb7b011d`) |
| #301 | — | Churchetal2015.pdf (`202b537d0a02`); Siebert_etal2013.pdf (`c9e7e8ae50a2`) |
| #302 | Hosiaetal2024 | — |
| #303 | — | Kidwai_Amjad2000.pdf (`392857a972b9`); Hauss_et2016.pdf (`71c01d252e97`); Pakhomov_etal2000.pdf (`8dd620bc32ce`); Sutherland_etal2019b.pdf (`e893ef271e78`); Haddock_etal2005.pdf (`ffc7ebff40b4`) |
| #304 | — | Church_etal2025.pdf (`a0acf5935174`); Haddock_etal2005.pdf (`ffc7ebff40b4`) |
| #305 | — | Deevey_Brooks1971.pdf (`225ea525a91f`); Siebert_etal2013.pdf (`c9e7e8ae50a2`); Hays_etal2018.pdf (`d67cb07d2f0b`) |
| #306 | — | Hung2002.pdf (`146fce7eee02`); Yu2006.pdf (`28d1eb27e929`); Tung2003.pdf (`db7338ea1867`) |
| #307 | — | Daniel1974.pdf (`b4425fd6ae7a`); Pugh_Haddock2016.pdf (`da370f1e0434`) |
| #308 | — | Daniel1985.pdf (`b4146be29447`) |
| #309 | — | Oderberg2020.pdf (`c59691fa845a`) |
| #310 | Totton1965a | — |
| #311 | — | Siebert_etal2013.pdf (`c9e7e8ae50a2`) |
| #312 | — | Boysen-Ennen1987.pdf (`88f1c8dfbd5c`) |
| #313 | Totton1965a | — |
| #314 | — | Pugh1990.pdf (`a16337443af7`); Pugh1974.pdf (`b5a7af6140ca`); Alvarinoetal1990.pdf (`b8ae22a47c6a`); Pugh1992a_Nectopyramidinae.pdf (`d0894c24715f`); Tung2003.pdf (`db7338ea1867`) |
| #315 | — | Mapstone2009.pdf (`4efbb4af134b`); Alvarinoetal1990.pdf (`b8ae22a47c6a`) |
| #316 | — | Fernandezetal2007.pdf (`07f414383157`) |
| #317 | — | Pugh1974.pdf (`b5a7af6140ca`) |
| #319 | — | Mapstone2009.pdf (`4efbb4af134b`); Siebert_etal2013.pdf (`c9e7e8ae50a2`) |
| #320 | Hosiaetal2024 | Church_etal2025.pdf (`a0acf5935174`); Siebert_etal2013.pdf (`c9e7e8ae50a2`); Pugh_Haddock2016.pdf (`da370f1e0434`); Tung2003.pdf (`db7338ea1867`) |
| #321 | — | Yu_thesis2006.pdf (`5545486bbef8`) |
| #322 | Hosiaetal2024 | Sutherland_etal2019a.pdf (`a65c78792a7e`) |
| #323 | Hosiaetal2024 | — |
| #324 | — | Siebert_etal2013.pdf (`c9e7e8ae50a2`) |
| #326 | — | Siebert_etal2013.pdf (`c9e7e8ae50a2`) |
| #327 | — | Siebert_etal2013.pdf (`c9e7e8ae50a2`) |
| #329 | — | Pugh_Haddock2016.pdf (`da370f1e0434`) |
| #332 | — | Nagata_etal2014.pdf (`117921bf43e4`); Mackie1960.pdf (`9be678a89090`); Siebert_etal2013.pdf (`c9e7e8ae50a2`); Mayer1906.pdf (`dc72b945a53f`) |
| #334 | — | Hissmann2005.pdf (`25eb8dd8d21d`); Mapstone2009.pdf (`4efbb4af134b`); DuClos_etal2022.pdf (`5408f4a1fb6f`) |

## How to cover the demonstrated gaps

- Bibliographic identity/locators use the small source-derived BibTeX selection
  in `tests/fixtures/bibliographic_integrity`, with the original library revision,
  source hashes and a compact new-species excerpt. No extra PDFs are vendored.
- Figure tests use captured caption geometry and source identifiers in
  `tests/fixtures/figure_integrity`. Hosia is already in gold. Sutherland's
  side-caption binding and the Erenna edge crop need separate source checks;
  `CORPUS_LIBRARY_DIR=/path/to/siphonophores/library` enables the latter's actual
  PDF replay in `tests/test_figure_integrity_sources.py`.
- Citation-span tests use compact TEI fragments and source-coordinate evidence
  from the specifically failing passages in Pugh, Oderberg and Mapstone. Their
  provenance and measured limits are documented in
  `tests/fixtures/citation_spans/README.md`.
- Encoding/accent assertions and compact native-glyph evidence are in
  `tests/fixtures/text_encoding`; scientific notation and treatment/order
  expectations are in `tests/fixtures/text_integrity`; logical tables, key
  branches and word boundaries are in `tests/fixtures/table_structure`.
  Their READMEs distinguish fresh reproductions, controlled historical forms,
  source-printed anomalies and unresolved readings. Full PDFs remain external.

## Integrated source pilots

The [September 22 normal-pipeline receipt](examples/source_pages_2026_09_22/README.md)
adds fresh OCR/extraction, real BGE-M3 embeddings, annotation and stdio MCP
acceptance for nine documents/13 selected physical pages. Chinese and German
figure captions and image bytes pass; surname evidence and the named spacing
relations survive materialization and serving. Unchanged extraction preserves
99 artifact hashes and all nine documents skip re-embedding. Original Russian
exponent evidence survives full-page OCR, with the second exponent still refused.
The separate authority replay distinguishes a repaired real retained edge from
three unresolved title-conflicted observations. This supersedes the older
pilots' missing real-embedding/German-caption checks, without turning selected
pages into a full-paper, full-corpus or fresh citation-graph claim.

A CPU extraction/chunking pilot on 2026-09-18 used copies of Tung2003 physical
pages 30 and 84 (mapped to pilot pages 1 and 2), preserving their source
geometry. It passed source assertions for the Chinese heading and figure-1
caption, `Nanomia bijuga`, `5301±8525`, and the exponent in `ind./100m³` in
saved text/chunks. Encoding receipts survived into chunk metadata without
being serialized into source prose. This is a two-page production-path pilot,
not a full-paper build, embedding evaluation, or corpus acceptance.
Bounded chunk discovery/fetch and figure-record/image calls now replay those
saved artifacts: the returned caption matches the source and returned image
bytes match the saved figure, with the complete pilot inventory unchanged.
A source-verified, already-legible Liu2012 title survives encoding recovery,
serialization, chunking and serving without edits. That control reconstructs
one text item from captured source geometry; it makes no fresh conversion or
whole-paper fidelity claim. See `tests/fixtures/text_encoding/README.md`.

The current default classifier independently identifies Boysen-Ennen1987 and
Mapstone2009 as scans and selects force-OCR. One-page pilots preserving those
full-document decisions remove the old German mojibake but still misread
`Fänge` as `Finge` and `Alvariño` as `Alvarifio`. Do not count the routing change
alone as repair.

The integrated #312 policy retains candidate regions from the original native
layer before preparation. It requires agreeing 300/600-dpi regional readings
and aligned prepared-word geometry. Replaying Boysen-Ennen physical page 12
through atomic figure/text materialization, real HybridChunker and `get_chunks`
returns `RMT-8-Fänge` with its original evidence retained. Eight other suspect
tokens remain explicitly unresolved; this is not a whole-paper German accuracy
claim. The replay also verifies receipt transfer into temporary extraction
outputs. Resume tests cover producer/model changes, retirement of obsolete
receipts, clean-build equivalence and an unchanged rerun.

Saved source cases also pass through real HybridChunker, the production
embedding callback, temporary LanceDB storage and bounded `get_chunks` calls.
The exact callback inputs and stored row texts retain the reviewed signs,
units, exponents, caption equations and German repair. An inexpensive recording
backend checks the input contract; it does not measure semantic model quality.
Unchanged notation and date-range controls survive, changed source text makes
the old embedding receipt stale, and an unchanged rerun skips the backend.
The German pilot has no figures, so it establishes no German caption result.
See `tests/fixtures/text_integrity/materialized/README.md` for captured-artifact
origins and limits; broader scientific precision/recall remains pending.

The #315 surname policy uses curated author/year candidates and requires two
unhinted regional OCR readings to agree. Its selected source sample covers
14 proposals across six papers: ten repairs agree with visual source review,
and four remain unresolved. Two already-correct source cases and four separate
synthetic controls are retained. This selected sample does not estimate
corpus-wide precision or recall. The compact crops, full source hashes and
reading receipts live in `tests/fixtures/surname_recovery`; see
[SURNAME_RECOVERY.md](SURNAME_RECOVERY.md) for the method and remaining limits.
Atomic extraction using the complete current curated catalog, real chunking
and `get_chunks` recovers Mapstone physical page 47's two tested `Alvariño`
citations (1971 and 1991). A separate retained-reference test demonstrates
supported edge rematerialization without changing the original observations;
that is not a fresh full-corpus Grobid/graph validation.

Those atomic surname pilot receipts predate policy v3. The strict current crop
fixtures retain their source-reviewed decisions, and the v3 producer invalidates
extraction and its consumers. Candidate build acceptance must use the updated
policy; a previous successful pilot does not establish that later rebuild.

The four named #314 observations have separate retained v1.2.1 rows and fresh
Grobid 0.8.1 source-page captures with raw citations enabled. A compact real
authority refresh repairs the line-wrap and article-omission cases, preserves
the already-correct fresh case, and leaves historical observations unchanged.
All four compact citing-document cases resolve to the source-supported publication,
with clean/incremental edge and missing-reference agreement. The subsequent
[full Tung replay](../tests/fixtures/bibliographic_integrity/publication_year_full_tung_2026_09_18/README.md)
uses actual normal preparation and full Grobid parsing: all 121 PDF pages are
preserved byte-for-byte, all 68 references are materialized, and the named
observation maps to Pugh1974 while its original parsed 1965 remains inspectable.
The [subsequent scanned-document replay](examples/siphonophore_citation_year_2026_09_21.json)
uses normal preparation and full Grobid parsing for Pugh1990 (76 pages/28
references), Pugh1992a (42/40), and Alvarino1990 (440/110). Every page is retained;
all three OCR subprocesses exit zero, without blanked, timed-out or textless
pages. All 178 reference observations survive real authority materialization.
The named observations map to Pugh1974 and are absent from the compact missing
list; Alvarino's fresh XML identifier is b83, distinct from historical b80.
Unchanged refresh is a no-op. The original report SQL used the wrong column
name; saved-artifact replay fixes reporting without repeating OCR/Grobid.
A broad report selector also included unrelated Kinzer raw1977/parsed1972;
that conflict remains unresolved and separate from the named Pugh acceptance.
The receipt preserves both failed reporting statuses and the final scoped pass.
Full-corpus overlapping citation sets and missing-work ranking remain unproved.
See the source crops, TEI fragments and receipts in
`tests/fixtures/bibliographic_integrity/publication_year_sources`.

Read-only legacy-bundle transport acceptance passes all registered MCP tools,
real query embedding, strict image crops, signed HTTP whole/panel delivery and
authentication checks, with the artifact inventory unchanged. It uses the
retained v1.2.1 snapshot identified below and proves compatibility/immutability,
not correctness of freshly rebuilt evidence.

Preparation for #305 copies three named source rasters (Siebert,
Deevey–Brooks and Hays), verifies original PDF hashes and renders their physical
source pages. `tools/qc/vision_source_pilot.py` runs these through the production
local backend while capturing processor patch tensors, reconstructed input
frames, prompts, raw responses, coordinate provenance and proposed crops.
Its numeric checks leave scientific source review explicitly pending. No fresh
model inference is established by preparation or the processor inversion tests.
The available v1.2.1 Hays receipt and panel population differ from the audited
snapshot; the original six-case seed/control manifest is still needed to claim
reproduction of that earlier pilot.

An [adapted operator tour](examples/siphonophore_operator_acceptance_2026_09_18.json)
passes 39 public CLI calls in the installed environment against a separate copy
of completed Hosia artifacts: scaffolding, all verb help, shell completions,
citation output, preflight, dry-run variants, status reports and BibTeX
export/import. Bibliographic fields and ordered authors survive unchanged
imports, and a controlled title edit exports correctly and can be restored.
Import deliberately refreshes provenance under #100, so the database is not
byte-identical after an explicit unchanged import; dry-run leaves it untouched.
This is partial operator acceptance, not a fresh environment, full build/resume,
live-server walkthrough or final-candidate validation.

## Independent review and fresh gold findings

The [intermediate gold baseline](examples/siphonophore_gold_baseline_2026_09_18/README.md)
completed extraction at `6fbf4e0`: all 35 document summaries and required stage
records report success, with no stage failures. The application logged successful
completion; the retained session wrapper returned 143, whose cause is unknown.
The receipt records both rather than claiming a clean shell exit. Three read-only
scorers then exited zero. Curator page selection excludes 86 of the 761 gold pages;
all 675 included pages are scored. Median prose-token coverage is 0.9448, physical
figure-count recall/precision is 0.8803/0.8642, and caption-identity binding is
0.5757/0.9817. These are intermediate CPU/OCR measurements, not final-producer,
Qwen or full-corpus acceptance. Token coverage cannot validate scientific signs,
and figure-count agreement cannot validate object identity.

The [independent retrieval sample and source decisions](examples/siphonophore_retrieval_review_2026_09_18/README.md), frozen against this completed
build before source review, found only five eligible units from two papers.
That cannot meet the previously defined minimum of ten queries from five papers.
Source inspection finds three usable units, all from Totton; two carry incorrect
stored treatment names. The other two selections are a bibliography fragment
and a nematocyst measurement table falsely classified as identification keys.
Keep the insufficient population and false-positive findings explicit; do not
replace them with convenient cases or claim a passing independent evaluation.

The subsequent table-classification correction requires destinations that resolve
to observed couplet numbers, supported by dotted leaders or paired alternatives
and separate destination cells. Both false-key source tables now remain ordinary
tables, preserving their cells and row metadata. The two actual Totton keys,
two Pugh–Haddock source keys and synthetic short-key controls remain keys through
serialization/chunk checks. This corrects the demonstrated classification defect;
it does not repair the two inherited treatment names, establish general key
recall or fill the missing independent evaluation coverage.

The subsequent treatment-context v3 correction addresses those two names.
The Prayinae key is a structural `TableItem` labelled `document_index`; it now
receives the same unknown treatment and separate merge role as other tables.
The printed `Sub-family` and key headings also end the preceding species scope.
For the other case, the entire plain-text item `Genus: PRAYOIDES Leloup 1934`
provides an explicit genus and dated authority. Its diagnosis is served under
`Prayoides` with `rank="genus"`, until the following explicit species heading
establishes `Prayoides intermedia`. No species is inferred from monotypic prose.
Both original saved-document contexts and compact save/chunk/serve replays pass;
the v3 resume test refreshes chunks and preserves extraction. This does not
establish complete treatment-boundary recall or replace integrated evaluation.

A source-first notation review froze 24 selections from 14 gold papers before
inspecting their candidate outputs. It contains 23 printed expressions and one
unscorable transcription-note range; the latter was retained rather than
replaced after review. At the recorded interim build snapshot, ten selections
were on included pages of completed documents: five were faithful in Docling,
three were corrupted, and two were omitted (a figure axis and a publisher
footer). Four other source pages were deliberately excluded by `keeppages`;
nine selections awaited document completion. Flattened text adds uncertainty
about one footnote's role. This is occurrence-level evidence, not a complete
corpus score or a repair-precision estimate.

The review exposed corrupted Chen quantities and a Russian cubic unit despite
passing named notation regressions. Two selected Mańko exponent repairs agree
with their source, but both come from the same paper. The sample contains no
selected ±, µm or subtraction-equation case. Labels remain frozen; subsequent
fixes informed by these findings make it a regression sample, not a new unseen
evaluation of those fixes.

The subsequent [Chen regional repair](../tests/fixtures/pdf_cmap/README.md)
requires exact same-region re-encoding through the original Type1 font map.
Three `mg/m³` occurrences additionally have source geometry and agreeing unhinted
digit OCR; reported statistical values retain their actual source typography.
Original-PDF replay and the saved-document → chunk → production embedding-input
→ vector-row → bounded-response path pass. Six hyphen candidates remain raw and
unresolved where the raster evidence is insufficient. Normally decoded text and
malformed source maps are nonmutating controls. This repair does not establish
full-paragraph accuracy, fix the Russian exponent or complete broader #303
source acceptance. The frozen [notation review](../tests/fixtures/text_integrity/source_review/README.md)
and its historical denominators are retained unchanged.

A separate [selected spacing review](../tests/fixtures/table_structure/spacing_review_2026_09_18/README.md)
retains three available produced-repair
records: two source-correct and one unresolved at a tight italic-to-roman
authority boundary. That is 16 correct and one unresolved boundary; dropping
the uncertainty would overstate precision. Five unchanged/rejected controls
and a separately labelled historical reconstruction do not enter that repair
denominator. Three additional geometric proposals have source-correct gaps,
but missing historical per-proposal OCR receipts prevent classifying them as
accepted repairs. The newly examined mixed-font DuClos run remains unresolved
because actual OCR modes disagree at `at | Friday`.

The [fresh Hosia acceptance record](../tests/fixtures/figure_integrity/hosia_gold_2026_09_18/README.md)
uses completed artifacts from the frozen `6fbf4e0` gold build. These now pass
real authority/taxon materialization, packaging and live serving under `9d1eac1`;
the relevant figure/rights/query code is identical at those revisions.
Twenty-seven direct/boundary cases and 25 live MCP calls plus three HTTP
downloads verify caption exclusions, supported abbreviation links, panel
fallbacks and matching stored/served record counts. Strict profiles refuse
the excluded Figure 1; the report profile permits it. The bundle and source
artifact inventories remain unchanged. All count routes agree on 15 records
for this paper, which is not a claim of 15 scientifically correct objects:
Figure 2 is missing, a Table 1 image is classified as a figure, and the missing
figure queue includes an apparent cross-publication number. This source replay
does not replace later-producer or full-corpus validation.

## Acceptance still required

The [Erenna Figure 51 delivery replay](../tests/fixtures/figure_integrity/erenna_fig51_delivery_2026_09_21/README.md)
completes #329's named crop check using saved real extraction geometry and the
original PDF. Production rendering, bundling and in-process MCP image conversion
preserve identical PNG bytes with the complete right/bottom species label, five
source panels and five scale bars, excluding neighboring prose. Fourteen
raster-edge controls pass, including rotation, cropbox margins and a prose veto.
This does not replace fresh layout extraction or corpus-wide figure recall.

Rebuild the existing gold set and score text/figure/caption fidelity, including
issue-specific assertions; do not substitute fixture presence or token coverage
for scientific correctness. Verify clean/incremental equivalence and downstream
invalidation. Inspect fresh Qwen source-pilot crops: coordinate arithmetic alone
does not prove the model selected the right panel. Finally rebuild the full
reference corpus, compare counts/mappings and replay the audited MCP calls.

The retained local served snapshot identifies itself as v1.2.1. It is useful
for historical reproduction and query controls, but it is not the audited
v1.4.0.dev0 bundle and cannot establish that deployment's reconciliation history
or the candidate release's correctness.
