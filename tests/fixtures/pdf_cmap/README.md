# Chen regional font evidence (#303)

`chen_regions.json` retains two actual Docling text items from the frozen
`6fbf4e0` gold extraction and native character geometry from the original
`library/C/Chenetal2015.pdf`. Source and observed-artifact SHA-256 identities,
physical pages, font IDs, exact ToUnicode scalar maps and CMap/Encoding stream
hashes are recorded. Text/orig/provenance/labels are actual fields; Pydantic has
normalized the item serialization. The original and prepared PDF hashes agree.
This is a regional capture, not a fresh Docling extraction or complete rebuilt
corpuscle.

The source font is Type1; Docling exposes encoded bytes (`HDYH !`) where the
source CMap maps them to `mg/m3`. The helper requires exact same-region,
one-to-one re-encoding of every nonspace scalar; only the installed Docling
curly-single-quote sanitization and scalar CJK compatibility equivalents are
accounted for. It does not use fuzzy text replacement or language guessing.

Three distinct `mg/m³` occurrences in the actual English-abstract item have a
small raised native `3` adjacent to the metre glyph. The retained independent
600-dpi raster crops show that digit; unhinted Tesseract eng PSM7 and PSM10 both
read `3`. The actual producer/model receipt and raw outputs are pinned. The
second item yields the source `R=0.596` and `P=0.001` with original fullwidth
letter/decimal typography retained; the equals crops pass a two-bar check.

Six unrelated hyphen candidates across the two paragraphs remain raw and
explicitly unresolved because their glyph boxes contain neighboring strokes.
The fixture keeps both the damaged observation and the candidate result. The
candidate is not complete paragraph proofreading: for example `ind /m 3` is
outside the narrow recovered unit-formatting rule and remains baseline.
Russian original-to-prepared-layer recovery is a separate unresolved problem.

`test_pdf_cmap_recovery.py` replays source rows and independent captured glyph
evidence, plus synthetic adversarial controls. Its driver double has an
explicitly synthetic PDF identity. The optional original-PDF replay uses the
actual parser, rendering and regional OCR:

```sh
CORPUS_LIBRARY_ROOT=LIBRARY_REPO python -m pytest tests/test_pdf_cmap_recovery.py -k original_pdf
```

Put the installed Tesseract executable on PATH. The materialization test runs
actual saved Docling serialization, HybridChunker, production embedding-input
selection, LanceDB writes and bounded chunk fetches. Only source-row/OCR replay,
token counting and the numeric embedding backend are instrumented. No embedding
quality, global proposal precision or full-corpus recall is claimed.
