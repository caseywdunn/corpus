# Frozen source review for #303

These are the exact policy, selection, labels and freeze receipts from the
2026-09-18 source-first review. Labels were frozen at 20:02:12 UTC before any
candidate extraction, native text or OCR was read. SHA-256 identities and actual
physical pages/region coordinates are retained. No PDFs, full-page rasters or
2 MB proposal population are committed.

The deterministic, stratified selection took 24 expressions from 14 papers out
of 2,177 regex occurrences across 688 eligible transcribed pages. One reviewer
checked every label against the original page and region raster: 23 printed
expressions were verified; one transcriber-note-only range remains unscorable.
It is a capped source sample, not a representative corpus sample or a measured
inter-reviewer agreement exercise. It does not cover ±/∓, µm or a subtraction
equation; all three inverse powers come from one paper. Transcriptions proposed
regions; they were not treated as source truth.

The subsequent interim comparison against actual build `6fbf4e0` found Chen
CMap decoding and a Russian prepared-OCR exponent failure. Those observations
informed fixes, so this sample is now a **regression sample for subsequent
candidates**, not an unseen holdout. Labels remain unmodified. No final
precision/recall number is implied by these source labels or the targeted Chen
regressions. The historical comparison remains at
`/tmp/corpus-v15-notation-comparison-6fb-interim/` with a before-inspection
counting policy, artifact hashes, pending/excluded/unscorable cases and distinct
expression-fidelity versus admitted-repair denominators.

To reproduce selection from the pinned library/transcriptions and render the
24 review regions into a new directory:

```sh
python tests/fixtures/text_integrity/source_review/regenerate.py LIBRARY_REPO NEW_OUTPUT
```

This verifies the source-manifest and every transcription hash and requires an
identical `selected.json`. Recreated policy timestamps/absolute paths naturally
differ; the historical receipts remain unchanged. PyMuPDF raster bytes can vary
with version; source PDFs, physical pages, region geometry and independently
reviewed meaning are the stable evidence. `label-freeze-receipt.json` also lists
hashes of the original local full-page/crop packets and population, intentionally
not all included here. Original packet: `/tmp/corpus-v15-notation-holdout/`.

`interim_6fb/` preserves the compact historical counting policy, per-expression
adjudications, denominators, uncertainty clarification and artifact-hash receipt.
The large build-artifact copies remain external. Those historical counts are
not updated to credit the new Chen regression: candidate comparisons must have
separate receipts and their actual build identity.
