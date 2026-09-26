# Current source-sample comparison for #303

This receipt grades the completed September 22 CPU gold refresh at
`ab6353f8d4a8ef59b1be314458d32437913635ec`, using the unchanged 24-expression,
14-paper source sample. All 35 gold documents passed the pinned completed-stage
fingerprint checks before the final capture. The scorer only read saved artifacts;
it ran no OCR, model or extraction pipeline. This is separate from the preserved
`6fbf4e0` baseline and earlier captured-region/source-page overlays.

| Surface / measure | Current result |
| --- | --- |
| Docling semantic fidelity among emitted expressions | 16/17 |
| Docling source-expression recall, consumed pages | 16/19 |
| Docling source-expression recall, all verified source selections | 16/23 |
| Markdown/chunks semantic fidelity among adjudicated emitted expressions | 14/15; one additional footnote role is indeterminate (bounds 14/16–15/16) |
| Markdown/chunks source-expression recall, consumed pages | 14/19–15/19 |
| Markdown/chunks source-expression recall, all verified selections | 14/23–15/23 |
| Selected current atomic-repair precision | Mańko signed exponents 2/2; Chen first-unit raised-digit formatting 1/1 |

The source sample contains 23 visually verified expressions. Four fall outside
the recorded library `keeppages`; they remain coverage exclusions. One retained
transcriber-note lead is not actually printed and remains unscorable. Neither
category is silently replaced. The consumed-page denominator is 19, and all 23
verified source selections are also reported. No incomplete document is scored
as a failure or replaced with an older artifact.

Fresh normal Chen extraction now retains the selected first/surface `mg/m³`
unit and `R=0.596` relation. Their success grades only those exact source
occurrences. It does not establish the accuracy of whole Chinese/English
paragraphs, all font mappings, or the surrounding text. The selected Russian
`0.3–10 мм³` still appears as `0.3-10 мм'` and is explicitly marked
`digit_ocr_disagreement`; a neighboring repaired `2000 мм³` earns no credit for
that frozen row. Both kinds of result remain in the source-recall denominator.

Current Totton1965b row 22 remains a table-to-picture omission: source-aligned
Docling text 121 retains station 1772's `750–500 m` depth under picture 3, but
markdown and chunks omit it. The selected Mańko figure-axis expression and Chen
footer are also omitted. These are missing-text recall outcomes, not wrong
admitted scientific repairs. The Beklemishev footnote marker's role remains
indeterminate in markdown/chunks and is represented by bounds. The three Totton
caption magnifications preserve meaning via `X`/`x` for `×` and decimal points
for middle dots; those literal differences are recorded and are not generalized
to scientific variables or other operators. Yamamori's selected `31.5 m` success
does not validate the corrupted Japanese prose around it.

## Precision, recall and scope

`selected-current-operations.json` records actual current operations, rather
than inheriting the historical 2/2: two selected Mańko inverse-exponent repairs
are source-correct (2/2), and the first/surface Chen unit's isolated raised-digit
formatting entry is source-correct (1/1). These are three reviewed atomic
exponent operations in two producer strata and two papers. The Chen paragraph's
other scalar mappings and its middle/bottom unit formatting entries are outside
that selected-operation denominator. `independent-selected-operations.json`
records a second agent's source-crop, charspan/bbox and item/chunk review for
these three operations. It is not a human domain review or a review of every
character in the source regions.

The two source-selected Chen expressions count separately as expression
recoveries; they are not two additional fully verified paragraph operations.
Four of the five selected expressions with repair-relevant damaged observations
now preserve their source meaning. That 4/5 is conditional expression recovery,
not independently enumerated detector/proposal recall. Semantic fidelity among
emitted expressions is also kept separate from admitted-repair precision.
Spurious outputs across the corpus and all rejected proposals were not fully
enumerated; output/proposal precision remains unmeasured. None of these bounded
counts is a corpus-wide error rate.

The frozen selection was independent of candidate inspection and of the issue's
named regression fixtures. Chen rows 4/8 and Russian row 16 subsequently informed
fixes, so the current comparison is a regression-sample result, not an unseen
holdout score. `graded-summary.json` keeps these three rows separate from the
other 21 frozen selections. This broader sample has no ±/∓, µm or subtraction
equation, and all inverse-unit examples come from one paper. The named source
fixtures separately cover those additional classes through extraction, real
HybridChunker/embedding-input code and a bounded served-text replay; the numerical
embedding double in that focused test is documented rather than described as a
real BGE build. The current CPU receipt itself makes no vector or live-server
claim; those surfaces have a separate verifier and receipt.

## Evidence and reproduction

`results.json` records every source identity, physical/prepared page, grade,
selected item/chunk reference and current artifact key. Exact unchanged selected
item/chunk evidence plus identical whole-document markdown can reuse an earlier
source grade; altered surfaces require an explicit new judgment, keyed to all
four current artifact hashes. Source regions 5/6/7, 14, 16, 22 and Chen 4/8/13 were
explicitly re-reviewed where changed. `capture-receipt.json` pins every captured
file. A second agent checked the final 24-row arithmetic and the current row-22
parent/omission evidence; `independent-final-denominator.json` records that
bounded review, not a second visual grade for every row.
`comparison-receipt.json` pins the scripts, input/config/producer identities
and source-label immutability checks; full saved snapshots remain in scratch.

Recreate the original source selection/crops using `../regenerate.py` and its
pinned source hashes. From the repository root, with the exact completed artifacts
retained, replay using:

```sh
/home/claude/miniforge3/envs/corpus/bin/python \
  tests/fixtures/text_integrity/source_review/current_ab6353f/replay.py \
  --gold /tmp/corpus-v15-acceptance/gold \
  --out /tmp/corpus-v15-notation-current-final \
  --require-complete
/home/claude/miniforge3/envs/corpus/bin/python \
  tests/fixtures/text_integrity/source_review/current_ab6353f/summarize_final.py \
  /tmp/corpus-v15-notation-current-final
```

Use a new output directory for any later repetition; snapshots are never
overwritten. Exit 2 rejects incomplete or prior-producer stages, and exit 3 means
changed regions still require source adjudication. Producer requirements and
adjudications are preserved with this receipt. The exact local scripts are
hashed, and their compact read-only forms are retained here for reproducibility;
paths identify this worked example, not requirements for other corpora. Source
labels/crops, gold artifacts, library files and configuration were never edited.
