# Existing gold set: completed CPU build and source acceptance, 2026-09-22

The three independent gold scorers completed successfully after the normal
35-document CPU extraction refresh finished. Extraction used `ab6353f`; scoring
used integration revision `bb84b4b55f9f1611afd695c22cf7aaf4c6a5cd25`. This packet
covers source scoring plus the subsequent normal embedding, post-processing,
bundling, live serving and unchanged-resume checks. Those phases have separate
receipts so successful execution is not confused with source fidelity.

The comparison preserves the same **35 documents and 675 scored gold pages**.
The gold contains 761 pages; configured `keeppages` excludes the same 86 pages
across 20 documents. There are no unmatched documents. Empty extraction remains
in the denominator: **48 extraction-empty pages and one script-missing page**
persist. These results do not establish that every source page is faithful.

| Text measure | Earlier `6fbf4e0` build | Current `ab6353f` build |
| --- | ---: | ---: |
| Median prose coverage | 0.9448 | 0.9482 |
| Median combined coverage | 0.9008 | 0.9015 |
| Mean taxon coverage | 0.8902 | 0.9035 |
| 10th-percentile taxon coverage | 0.6667 | 0.7436 |
| Median recall | 0.9329 | 0.9376 |
| Pages below 0.5 coverage | 87 | 86 |

Fourteen pages improve on at least one reported coverage/recall/similarity
measure: all 12 Chen pages, Stepanjants2014 page 8 and Tilesius1814 page 10. The
other 661 page metric records are unchanged. None of the six compared measures
(combined/prose/figure/taxon coverage, recall, similarity) decreases. Only Chen
and Tilesius document aggregates change; the small Stepanjants page improvement
does not move its aggregate. Chen median prose coverage rises from 0.8118 to 0.9637
and combined coverage from 0.6673 to 0.8135, while its page 6 combined coverage
remains low at 0.231. The document is not treated as wholly repaired.

Figure-detection reports are identical. Per-page count matching finds 331 of 376
gold figure/plate blocks among 383 physical detections (recall 0.8803,
precision 0.8642). The default served-type count surface matches 310/376 with no
count surplus (recall 0.8245, precision 1.0). This is a **count approximation**:
it cannot show that a crop is the correct object merely because counts agree.

Caption identity and panel totals are also unchanged: 483/839 gold identities
matched, 492 reported identities, recall 0.5757 and precision 0.9817. Of 98 gold
panelled figures, 84 have reported declarations and 66 have exact label sets.
The baseline command forced `--figure-panels ocr`; the current configuration
uses `vision-local` with unavailable vision deliberately deferred. These are
**not a clean panel-model comparison or fresh Qwen acceptance**.

Five caption previews change in Chen after text recovery, while their records
remain unnumbered and therefore do not change identity metrics. One inherited
association, page 6 `docling_6`, remains a URL marked as a bound caption; before
recovery its same source string was damaged. This is an outstanding source
association limitation, not a newly introduced regression or a reason to infer
that all captions are correct from the aggregate precision.

The gold's `CROSSCHECK_REPORT.md` was read before interpreting these measures.
Here the independent transcription is the reference and the extraction is on
trial. The primary coverage measure asks how much source text extraction
recovers; missing/garbled extraction must not be dropped as unscorable. The
scorer normalizes case, diacritics and punctuation, so these token scores cannot
establish exact scientific signs, exponent roles or freedom from invented text.
The separate source-region notation review remains necessary.

## Complete pipeline, serving and unchanged resume

`pipeline-acceptance.json` records the normal CLI phases and input identities.
Extraction finishes with a clean exit after 5,203.7 seconds: all 35 summaries
and required stages are complete, current producer checks pass, and there are
no stage failures. Existing source-integrity, low-text and missing-reference
warnings remain visible. Real CPU BGE-M3 embedding completes all 35 documents
with zero failures in 1,319.4 seconds. All four post-processing steps and bundling
pass. The bundle holds 35 papers, 3,336 chunks and 605 figure records; path
scrubbing reports a clean audit. Production modules remain identical to
`ab6353f` throughout the later receipt-only commits.

`serving-summary.json` records an actual stdio MCP replay at `7b4343f`. All
3,336 unique 1,024-dimensional vectors are finite and nonzero, match exact
materialized payloads, and have valid current producer/generation markers.
The same 3,336 chunk texts/headings survive 152 bounded `get_chunks` calls;
the largest response is 257,019 bytes against a 524,288-byte verifier limit.
All fourteen clearly faithful source-sample expressions survive at their exact
reviewed occurrences. The ten other source selections remain in the separate
[notation review](../../../tests/fixtures/text_integrity/source_review/current_ab6353f/README.md),
not silently excluded from its fidelity denominator. Six report-profile image
deliveries preserve bytes. Manuscript mode permits one licensed-open figure
and issues five structured refusals for absent clearance records. The full
bundle file inventory remains unchanged during verification.

Unchanged extraction skips all 35 papers; unchanged embedding reports
`embedded=0, skipped=35, failed=0`. Document artifacts and embedding markers
retain their hashes. The first, deliberately strict physical-file comparison
returns false: normal orphan pruning creates a no-row deletion transaction,
adding a transaction and manifest and updating the version hint. Both strict
receipts are preserved. Direct comparison of every row in the pre-resume bundle
(table version 36) and resumed build (version 37) proves exact equality of
keys, text, full metadata, vector values and generations. This meets semantic
unchanged-resume acceptance without claiming byte-identical database files or
a complete clean/incremental change-class matrix.

The initial live verifier attempt also remains recorded: it stopped after 135
successful chunk calls because its decoder unwrapped a singleton list. Using
the SDK's structured list envelope fixes that scratch harness; the complete
second attempt passes. No product change or model rerun was needed.

These results complete CPU gold acceptance. Fresh Qwen inference, a complete
reference-corpus build and matching-population retrieval evaluation remain
separate release work.

## Evidence and reproduction

- `scoring-receipt.json` pins commands, extraction/scoring revisions, raw report
  hashes, source/config identities, completion-receipt pointers, scorer resources
  and validation. The three scorer files and their directly imported production
  declarations are unchanged from the baseline implementation.
- `scoring-summary.json` preserves both corpus summaries, all changed page
  metrics, figure/caption counts and the five caption-preview changes.
- `scoring-input-hashes.json` pins all 245 tracked build artifacts and the digest
  of the full 761-page source-transcription hash inventory. All stayed unchanged
  during scoring. The complete per-page inventory remains in the hashed scratch
  launch receipt; no PDF or transcription corpus is vendored here.
- `SHA256SUMS` protects this scoring packet. Raw full reports and logs remain in
  `/tmp/corpus-v15-acceptance/gold/scores_ab6353f_20260922/`; prior reports remain
  in `scores_6fb_2026_09_18/`. Absolute paths are historical provenance and must
  be adapted on another machine. The earlier report's recorded `output` path
  predates its preservation as `baseline-before-sep22-refresh/`.
- The coordinator's `refresh-20260922-completion.json` and
  `refresh-20260922-artifact-comparison.json` are cross-linked by exact hash in
  the scoring receipt. Their extraction, embedding and serving records remain
  distinct from this scoring result.

At the pinned scoring revision, substitute paths to the same source gold,
complete build, and configured bibliography, then run these independent
read-only CPU jobs. They load no extraction or query model:

```bash
python tools/qc/fidelity.py \
  --gold /path/to/transcriptions --corpuscle /path/to/output \
  --bib /path/to/inputs/siphonophores.bib --out /tmp/fidelity.json
python tools/qc/figure_detection.py \
  --gold /path/to/transcriptions --corpuscle /path/to/output \
  --out /tmp/figure_detection.json
python tools/qc/caption_binding.py \
  --gold /path/to/transcriptions --corpuscle /path/to/output \
  --out /tmp/caption_binding.json
```

The completed jobs took 11.41s, 0.69s and 1.02s respectively, with peak RSS 334,320,
26,720 and 32,276 KiB. All exited 0. They made no source or build changes, launched
no OCR/model process, and did not inspect sealed retrieval-evaluation outcomes.
