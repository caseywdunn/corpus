# Contextual exponent experiment, 2026-09-22

**Result: no additional exponent verified. Keep the second expression unresolved;
leave the production verifier and producer identities unchanged.**

This bounded experiment tested whether including the printed unit letters could
independently corroborate the second exponent in Stepanjants2014, physical page 8.
It also applied the same crop rule to the first exponent as a positive control.
It did not tune a repair, update frozen labels, or run full OCR/extraction/models.

Source: `S/Stepanjants2014.pdf`, SHA-256
`d183f5012c5438cd99f3016b1ade97241d9207bdcd7a0d92eb844ba5ff151b4a`.

Settings declared before observing the OCR outputs:

- Original-source box containing the two adjacent unit letters and encoded
  digit, including an immediately following closing parenthesis if present;
  pad the union by 1.5 points and clip to the page.
- Render directly in grayscale at 300 and 600 DPI; no image preprocessing.
- Tesseract 5.5.2, English model, PSM 7, no candidate text/whitelist/dictionary.
- English traineddata SHA-256:
  `8280aed0782fe27257a68ea10fe7ef324ca0f8d85bd2fd145d1c2b560bcb66ba`.
- Four total OCR calls, each bounded to three seconds and 40,000 crop pixels;
  `OMP_THREAD_LIMIT=1`. No retries or follow-up crop/settings tuning.

| Case | 300 DPI raw stdout | 600 DPI raw stdout |
|---|---|---|
| First exponent, positive control | `MM?\n` | `MMS\n` |
| Second exponent, unresolved target | `VIVE )\n` | `MMS )\n` |

All four calls exited 0 with empty stderr. Exact original-page crop coordinates
in PDF points (left, top, right, bottom):

- First: `[174.2631072998047, 626.840576171875, 193.4008026123047, 639.897216796875]`.
- Second: `[86.86289978027344, 648.7823486328125, 109.41612243652344, 663.6718139648438]`.

The 600-DPI source crops were visually checked by the same agent conducting this
experiment. They retain the printed unit context and raised glyph. This is an
agent visual check, not independent human/domain review. It does not replace the
required unhinted OCR corroboration. In particular, `S` is not accepted as `3`.

Neither pair supplies two readings agreeing with the encoded digit. The first
positive control also loses its successful digit-only corroboration when this
context is included. These two observations therefore do not justify a general
contextual fallback. They do not establish that every possible wider crop or OCR
configuration must fail; no further configuration was tried in this experiment.

`report.json` is the immutable raw receipt, including source and baseline fixture
identities, exact producer, commands, timings, crops, hashes and outputs.
Its SHA-256 is
`fed1b5081d76b124d1bbfb693fd90b4e15954bc67b51679464cb6d75b7806f32`.
The receipt records the exploratory driver hash; that path-specific script is
not part of this portable fixture. The exact crop boxes, rendering parameters
and OCR commands are retained, and the four PNGs preserve every OCR input.
No production code, gold build, library source or existing capture was modified.

The raw receipt and all four PNGs are byte-identical to the experiment outputs.
`manifest.json` inventories their hashes and the documentation hash. The parent
`../capture.json` is the original digit-only evidence; its identity at experiment
time is retained in `report.json["baseline_capture_sha256"]`. This failed-context
experiment neither replaces that evidence nor changes its source decisions.
