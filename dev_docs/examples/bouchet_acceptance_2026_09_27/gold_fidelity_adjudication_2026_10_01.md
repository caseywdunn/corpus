# v1.5 gold fidelity adjudication — 2026-10-01

This checks the earlier full candidate's extraction against the v1.4 reference
on the same 35 documents and 675 selected gold pages. It does **not** accept the
latest-code full-library build running from `792f478`; that build needs its own
source and served-bundle checks. The [September 29 comparison](gold_comparison_2026_09_29.json)
remains the historical raw receipt.

The September 29 fidelity scorer assigned every cross-page Docling text item to
its first provenance page. One Beklemishev item starts with seven stray
characters from a plate and continues with 186 words on the next page. The
scorer reported nearly all those words on the plate. `tools/qc/fidelity.py` now
uses each provenance `charspan` to score the text on the page where it was
extracted. The regression test in `tests/test_fidelity_harness.py` covers this
case. All 40 fidelity harness checks pass, and Ruff passes. The corrected raw
reports are retained under
`scratch/v15-bouchet-20260925/gold-acceptance-20260929/build-scores/` as
`reference-fidelity-corrected.json` and `candidate-fidelity-corrected.json`.
The original input-hash audit still matches for every `docling_doc.json` used
in this comparison; only two `figures.json` files have since changed while the
latest build runs, and fidelity scoring does not read those files.

| Measure | v1.4 reference | Earlier v1.5 candidate |
| --- | ---: | ---: |
| Median total coverage | 0.9388 | 0.9398 |
| Median prose coverage | 0.9789 | 0.9806 |
| Mean taxon coverage | 0.8847 | 0.9021 |
| Median figure-text coverage | 0.7644 | 0.7638 |
| Pages below 0.5 total coverage | 73 | 74 |

The correction removes most supposed Beklemishev prose losses, but it does not
make the comparison a blanket pass. On two landscape plates in that book,
re-OCR reads the rotated captions upside down: source page 17 coverage drops
0.9236 → 0.0417 and page 33 drops 0.9348 → 0.0217. The original text layer is
readable on those pages. This requires a source-backed OCR or page-preservation
fix, or an explicit release deferral with a known limitation. Issue #346 tracks
it.

Yamamori 2014 is a real tradeoff: median prose coverage drops 0.8819 → 0.7595
on Japanese pages, while taxon coverage rises 0.0 → 0.9444 and figure-text
coverage rises 0.6951 → 0.8333. The original layer has spaced/broken Japanese
and Latin characters; the candidate's re-OCR gives different, often cleaner
words but omits some gold prose. Page-level source review is still needed.
Ahuja page 1 drops only 0.0229 in prose coverage; its document median is
unchanged after the scorer correction. The Quoy–Gaimard plate aggregate is
unchanged after correction.

Physical figure scoring and caption-binding results in the September 29
comparison were not changed by this fidelity-scorer correction. They still
need source review, as does a fresh gold score against the latest-code build.
