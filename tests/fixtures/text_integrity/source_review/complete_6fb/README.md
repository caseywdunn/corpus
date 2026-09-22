# Completed source-sample comparison for #303

The original 24-expression/14-paper source labels remain byte-for-byte unchanged.
This report finishes grading the **preserved complete `6fbf4e0` CPU baseline**,
including the nine rows that were pending in `../interim_6fb`. It does not present
the baseline as the current September 22 refresh. Artifact SHA-256s, completed
stage receipts, source pages, item references and exact local contexts identify
the evidence. All four graded artifact files for previously graded documents
match their interim captures exactly. No OCR, extraction or models were run.

| Surface / measure | Complete baseline result |
| --- | --- |
| Docling semantic fidelity among emitted expressions | 14/17 |
| Docling source-expression recall, consumed pages | 14/19 |
| Docling source-expression recall, all verified source selections | 14/23 |
| Markdown/chunks semantic fidelity among adjudicated emitted expressions | 12/15; one additional footnote role is indeterminate (bounds 12/16–13/16) |
| Markdown/chunks source-expression recall, consumed pages | 12/19–13/19 |
| Markdown/chunks source-expression recall, all verified selections | 12/23–13/23 |
| Admitted scientific-repair precision in selected occurrences | 2/2, both inverse-exponent repairs in one Mańko paper |

Four source selections lie outside the library's consumed `keeppages`; they are
coverage exclusions, not demonstrated notation failures. One transcriber-note
lead is not printed and remains unscorable. No pending baseline rows remain.
`summary.json` includes all omissions, corruption, exclusions and uncertainty;
the favorable emitted-expression denominator is never the only reported one.

The nine newly graded rows are 5, 6, 7, 9, 11, 18, 19, 22 and 24. The three
Totton magnifications preserve their clear caption meaning using `X`/`x` for
`×` and a decimal point for the middle dot. Those literal substitutions are
recorded, not generalized to scientific variables or arbitrary operators.
Eight rows survive markdown/chunks. Row 22's station-1772 `750–500 m` table
cell survives in Docling text node 121, but its parent is a picture and the text
is absent from markdown/chunks. That is a baseline surface omission; it must
be checked on the completed current refresh before describing current code.
Yamamori's selected `31.5 m` survives, while substantial neighboring Japanese
prose corruption is explicitly outside this narrow expression grade.

The two repairs are source-correct, but 2/2 does not establish a corpus-wide
precision estimate. Five aligned selected expressions show repair-relevant
damage in their observations: the two Mańko exponents were repaired, whereas
Chen's unit/relation and the Russian cubic unit were still corrupt in this
baseline. Missing figure-axis/footer text and the table-to-picture loss are
reported separately as missing-text recall, not false repair operations.
All spurious notation outputs or rejected proposals have not been independently
enumerated, so output/proposal precision remains unmeasured. Source-expression
fidelity is not renamed repair precision.

## Separately identified later evidence

* `chen-current-helper-overlay.json`: current `ab6353f` pure-helper replay of
  the two pinned real Chen regions, actual font maps/native glyphs and captured
  source-raster decisions restores `mg/m³` and `R = 0.596`. It is captured-region
  regression evidence, not a fresh normal full-document extraction. Other
  unsupported hyphens remain unresolved. The existing fixture is
  `tests/fixtures/pdf_cmap/chen_regions.json`.
* `russian-fresh-overlay.json`: the fresh normal page-8 build, real BGE-M3 vector
  text and actual MCP `get_chunks` at bundle revision `bfbf90c` all retain the
  selected **`0.3–10 мм³` as `0.3-10 мм'`**. It is still a fidelity failure,
  now explicitly marked `digit_ocr_disagreement`. The adjacent `2000 мм³` is
  repaired and served correctly, but is a different source occurrence and
  earns no credit for frozen row 16. Full original → selected page → prepared
  page identities are retained. No guessed exponent was introduced.

These overlays are not pooled into a fictitious single-build score. The
September 22 full-gold refresh was still running when this report was prepared;
completed source/artifact identity must be checked before replacing any baseline
row. No full-sample embedding/serving score is inferred from CPU chunks.

## Scope and reproduction

The frozen selection was source-first and separate from the issue's named
regressions. Chen rows 4/8 and Russian row 16 subsequently informed fixes, so
later comparisons on these unchanged labels are regression scores, not an
unseen holdout. This sample has no ±/∓, µm or subtraction equation; all three
inverse-unit selections come from one paper. The named regression fixtures
cover additional classes but cannot substitute for independent broad coverage.
These limitations constrain conclusions; the measurement requirement does not
imply zero errors across every kind of source text.

An independent agent re-viewed and hash-checked source crops 5, 6, 7 and 22,
checked their candidate regions, and verified the whole-sample arithmetic.
`independent-review-fix-figures.json` records that bounded review; it does not
claim a second source grade for all 24 selections.

The frozen source render/selection reproduction remains `../regenerate.py`.
The complete read-only local scoring packet is
`/tmp/corpus-v15-notation-comparison-complete-20260922/`; its `capture.py`,
`align.py`, `grade.py` and `current_overlays.py` are hashed in
`comparison-receipt.json`. The large raw document copies are kept there, not
committed. To audit each row, resolve the original source hash/page and frozen
crop, inspect the named candidate item and bbox, then inspect **all chunks
linked to that item** for the selected local context. A repeated matching token
elsewhere does not count. `grade.py` asserts the exact nine new contexts,
their markdown/chunk presence or absence, no admitted scientific repair for
those occurrences, and identity of all reused historical artifacts.
