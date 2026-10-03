# Full Tung preparation and reference replay, 2026-09-18

This dated acceptance record advances #314 for one of its four named citing
sources. The complete 121-page Tung2003 PDF passed actual normal page preparation,
full Grobid parsing and authority materialization. **#314 remains open** for
three complete scanned documents and the full-corpus degree/overlap/ranking gates.
This is not a Docling/full build, a deployed-corpus result or historical v1.4 proof.

Code was pinned to `a3e638dc2003ffd7ebc0a32cefbc359d97c8461a`.
The original source hash is
`db7338ea18673a705c43f71c1833e980e3c7ee8de14879c81fc9a49a3955cf1f`.
Actual library BibTeX supplied the header and had no `keeppages`/OCR override.
Source, Bib, config, module and artifact identities are in [acceptance.json](acceptance.json).

## Measured path

1. Actual `detect_scan_type` selected `born_digital`, `clean_text_layer`,
   `needs_ocr=false`, with the installed English pack available. Actual
   `prepare_pdf` found no native recovery candidates and copied all 121 pages
   byte-for-byte. No synthetic detection result or raw-PDF substitution was used.
   Detection took 0.28 seconds and preparation 0.61 seconds.
2. Actual `extract_metadata` sent that complete prepared PDF to Grobid 0.8.1,
   with header consolidation 1, citation consolidation 0, raw citations enabled
   and ref coordinates. A fresh request took 4.06 seconds and produced **68
   references**, full TEI and the normal hash/input provenance receipt. The
   running container/image identity is retained. No OCR or model build ran.
3. Production `phase1_corpus_papers` and `phase2_references` materialized all 68
   observations in a compact authority containing Tung and the actual Pugh1974
   BibTeX target header. No references were filtered out before materialization,
   and no derived observation/edge row was manually inserted.
4. Full-document XML `b40` retains parsed year 1965 and the original raw 1974
   citation. Its raw string exactly matches the earlier source-page capture
   grounded in original physical page 60. The authority maps it to
   `10.1017/s0025315400022086` via `raw_publication_year_title_authors`, with
   publication-year decision evidence and the original observation intact.
5. The actual graph contains one Tung→Pugh1974 edge. Direct bibliography calls
   retrieve all 68 unique XML IDs in pages of 50, 18 and 0. Bibliography retains
   raw year 1965 and exposes the DOI/decision provenance; the canonical formatter
   reports 1974. The compact missing list has 60 other leads and no false Pugh1965
   title. This is explicitly not the full-corpus missing-reference ranking.
6. Unchanged authority refresh is a no-op. Serving leaves authority bytes
   unchanged and the source PDF is unchanged. These are production tool-function
   calls, not an MCP transport or bundle replay.

The first harness run omitted the conda executable PATH and marked Tesseract
availability unknown. That run was retained separately and superseded by this
fresh run with the correct PATH and available pack. Its TEI was not reused.

## Files and source identity

- `acceptance.json` records source/producer identities, normal routing,
  timings, all generated artifact hashes, raw `b40`, source page/bbox anchor,
  source-module hashes, tool checks and limitations.
- `complete_reference_observations.json` preserves all **68 actual full-document
  observations** and resolved mappings, including raw strings and parsed years.
  Its SHA-256 is
  `c0af962c209a5319a6191764a5fa80fb7d53798abefb7667b5810213a3b45e4f`.
  `complete_quality_decisions.json` and `target_authority_mapping.json` retain
  the corresponding quality decisions and named edge.
- `execution_config.json` pins the resolved OCR/Grobid/timeout configuration.
  Its service URL identifies local execution, not a public service dependency.
- `queue.json` pins each original PDF, exact actual Bib entry and named source
  page, along with the preparation stages still required for the three scans.

PDFs, full TEI, source images and databases are deliberately not distributed
here. Their original hashes identify the measured artifacts; hashes cannot
reconstruct them. The source anchor and existing source crops are retained in
[publication_year_sources](../publication_year_sources/). The observations above
come from this full-document run, not copies of those earlier one-page results.

## Remaining source preparation and reproduction

| Source | Full pages | Named reference physical page | Pinned preparation |
| --- | ---: | ---: | --- |
| Pugh1990.pdf | 76 | 75 | eng, no keeppages |
| Pugh1992a_Nectopyramidinae.pdf | 42 | 42 | eng, no keeppages |
| Alvarinoetal1990.pdf | 440 | 439 | eng, no keeppages |

All three have first-five-page raster fraction 1.0. Current normal policy
re-OCRs their existing scanned text layers, so **558 pages** require actual
preparation. These are read-only geometry findings, not claims that their
normal detection probes/OCR ran. No provenance-backed complete prepared PDFs
were available in the inspected retained artifacts. Do not substitute their
raw embedded text or one-page captures for complete prepared-document evidence.
The scans remain queued for coordinated CPU/OCR scheduling after gold; no scan
was launched by this replay.

To reproduce Tung, use its exact original PDF and the pinned Bib entry/config
with the recorded checkout. Run production page selection (none is selected out
here), `detect_scan_type` and `prepare_pdf`; retain actual routing and prepared
hash. Call `extract_metadata` on the complete prepared PDF with the real Bib
entry and live Grobid version. Retain the entire returned reference list. Create
an authority using production schema and phases, including the actual Pugh1974
Bib header, then inspect source-supported `b40`, all observation rows, refresh
behavior, graph, paginated bibliography, canonical formatter and missing list.
New TEI timestamps and DB observation timestamps can change output hashes;
record a new dated receipt rather than rewriting this one.

For the queued scans, follow `queue.json` and run actual full OCR preparation
before Grobid. Validate each named source observation independently. Final
acceptance still needs the intended complete corpus rebuilt, true/false Pugh
citing-paper sets and their overlap measured, and acquisition/ranking surfaces
replayed. Historical 59 and 23 cannot simply be added. Existing edition,
legitimate title-date and genuinely missing-reference tests remain separate
controls. No issue was closed and no production code changed in this replay.
