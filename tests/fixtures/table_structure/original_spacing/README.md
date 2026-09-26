# Original spaces lost during OCR (#334)

The three examples were found in a fresh normal pipeline build of selected
Mapstone2009 source pages on 2026-09-22, after the earlier native-text spacing
tests had passed. They are regression cases selected after observing failures,
not an independent precision sample.

`mapstone.json` contains original/prepared character observations around each
affected line, an actual prepared Docling text item or table cell, and the
bounded source-crop readings. It pins the original full PDF, kept-page PDF and
prepared PDF separately. Kept pages 1/2/3 correspond to physical source pages
47/68/200. No full PDF or corpuscle is vendored. The PNGs are small source crops
recaptured from a separate PDF handle; their hashes are recorded separately from
the original OCR-input crop hashes. All three recaptured images were visually
checked by the implementing agent, not a human domain reviewer.

| Physical page | Prepared text | Original wording, confirmed by both OCR modes |
| --- | --- | --- |
| 47 | `nectophoreslisted` | `nectophores listed` |
| 68 | `alonedeveloped` | `alone developed` |
| 68 | `nectophoreonly` | `nectophore only` |

Only missing original spaces are inserted. Exact letters and case must agree,
two nearby words must independently align the page coordinates, and both
unhinted Tesseract modes 6/7 must confirm all original boundaries. OCR-only
splits are ignored. Page identity, page count, rotation/crop geometry, finite
candidate/time budgets and a unique structured owner constrain application.
The source/prepared hashes, page mapping, original token, cell index, boxes,
anchors, crop readings and producer survive into chunk provenance. Applying
the same receipt again changes nothing.

This route does not lower the existing geometric fallback's 20-letter or
three-gap thresholds. Existing source-printed closed runs, legitimate compounds,
ambiguous owners, altered letters and uncorroborated OCR remain unchanged.
The existing Mapstone2009 physical-page-200 fixture remains a negative control.
The direct replay repaired all three observations in the actual prepared
document. A fresh normal-pipeline refresh and unchanged-resume check are
recorded separately when completed; the captured replay is not that refresh.
