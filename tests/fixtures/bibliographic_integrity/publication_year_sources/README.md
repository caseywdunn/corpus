# Named Pugh publication-year observations (#314)

These compact reference regions supplement the existing siphonophore gold
subsample; they are not another corpuscle or a full-paper fixture. The five
PNGs together occupy about 170 KB. Capture and source review: 2026-09-18.
Library revision: `5ad0164e6840ea16adb6e54a3a5712e20b234feb`.

`cases.json` keeps three distinct kinds of evidence:

1. Exact parsed reference rows retained in the local August 31 corpuscle
   (`_serve`, v1.2.1), including their empty `raw` fields. This is **not** the
   September 9 v1.4.0.dev0 deployment cited by the issue.
2. Original PDF hashes, physical page numbers, bounded reference crops and
   uncorrected native text lines/coordinates. The crops were visually checked.
   `Pugh1974-title.png` independently shows the 1974 publication, title ending
   in 1965, author, volume 54 and pages 25–90.
3. Fresh Grobid 0.8.1 outputs from each complete source page, with
   `consolidateHeader=0`, `consolidateCitations=0`, `includeRawCitations=1` and
   `teiCoordinates=ref`, through the production `GrobidClient`. No OCR, Docling,
   external consolidation or full-document build ran for this capture. The
   saved TEI fixtures contain the unchanged selected `biblStruct` inside a
   minimal generated TEI wrapper. Full-response and fragment hashes are
   recorded separately. Page-only XML IDs are new, not replacements for the
   historical IDs.

| Historical observation | Source PDF / physical page | Fresh ID / parsed year | Source and mapping evidence |
| --- | --- | --- | --- |
| `db7338ea1867/b40` | `T/Tung2003.pdf`, 60 | b2 / 1965 | Raw citation explicitly begins Pugh 1974; exact title mapping. |
| `d0894c24715f/b26` | `P/Pugh1992a_Nectopyramidinae.pdf`, 42 | b5 / 1965 | Printed line wraps `siphono-` / `phores`; raw Grobid joins lines with a space. Comparison-only joining must reproduce the curated title. |
| `b8ae22a47c6a/b80` | `A/Alvarinoetal1990.pdf`, 439 | b18 / 1974 | Source says Cruise. Historical parsed Crise is retained unchanged. Fresh parsing already yields 1974, with a journal fragment in its title; the existing title/year/author resolver maps it. No guessed spelling correction was added. |
| `a16337443af7/b15` | `P/Pugh1990.pdf`, 75 | b13 / 1965 | The printed citation omits “the” before siphonophores. Mapping tolerates that one article omission only with the independent exact raw volume/page pair, authors and publication year. The source citation is not rewritten. |

Production obtains new raw evidence by rerunning metadata extraction with
`includeRawCitations=1`; verified TEI cache receipts must match those options and
the prepared PDF. A bibliography-only rebuild cannot supply raw text that is
absent from an old observation. The acceptance test first materializes all four
retained rows, proves a producer-only refresh cannot repair them, then supplies
the four **actual new parser outputs** and preserves the older rows as history.
No source transcript is silently inserted into an old raw field.

The actual compact authority replay (five curated source works and these four
references) gave:

| Materialization | Citing documents at Pugh 1974 | At false 1965 work | Overlap | Missing rows |
| --- | ---: | ---: | ---: | ---: |
| Retained empty-raw observations | 0 | 4 | 0 | 1 |
| Fresh page parses, before this fix | 2 | 2 | 0 | 1 |
| Fresh page parses, after this fix | 4 | 0 | 0 | 0 |

Clean and incremental builds agree on citation edges and missing rows. The
incremental authority retains eight observation rows (four historical plus four
fresh); an unchanged rerun makes no writes. Separate existing tests retain a
genuinely missing synthetic work and prove deduplication when a citing document
has both forms. The compact counts above are actual named cases; they do not
estimate full-corpus precision or the overlap behind the issue's 59/23 counts.

Reproduce the compact source acceptance without a running service:

```sh
python -m pytest tests/test_reference_year_sources.py tests/test_reference_year.py -q
```

To recapture with the original library and local service, use PyMuPDF
`insert_pdf(source, from_page=physical_page-1, to_page=physical_page-1)` for each
manifest case, verify `source_sha256` first, then call:

```python
tei = client.process_fulltext(page_pdf, consolidate_header=0,
                             consolidate_citations=0, include_raw_citations=True)
references = parse_tei_references(tei)
```

Save the full raw TEI and actual service version before selecting the one
reference containing both “vertical distribution” and “siphono”. The original
capture's full responses, page PDFs, commands and regenerated authorities are
under `/tmp/corpus-issue314-source/`; the portable regression needs only the
committed fragments and curated values.

Remaining release gates: fresh full prepared-document metadata; full authority
and citation rebuild; actual citing-document sets and overlap for the reported
1965/1974 nodes; the complete missing-reference ranking; and those counts and
provenance through the candidate served bundle. These page pilots do not show
that the deployed historical bundle has changed, or recover its missing TEI.
