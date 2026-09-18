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

A CPU extraction/chunking pilot on 2026-09-18 used copies of Tung2003 physical
pages 30 and 84 (mapped to pilot pages 1 and 2), preserving their source
geometry. It passed source assertions for the Chinese heading and figure-1
caption, `Nanomia bijuga`, `5301±8525`, and the exponent in `ind./100m³` in
saved text/chunks. Encoding receipts survived into chunk metadata without
being serialized into source prose. This is a two-page production-path pilot,
not a full-paper build, embedding evaluation, or corpus acceptance.

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

Read-only legacy-bundle transport acceptance passes all registered MCP tools,
real query embedding, strict image crops, signed HTTP whole/panel delivery and
authentication checks, with the artifact inventory unchanged. It uses the
retained v1.2.1 snapshot identified below and proves compatibility/immutability,
not correctness of freshly rebuilt evidence.

## Acceptance still required

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
