These are small source-derived excerpts for issues #303, #304 and #319. They
contain selected text items, bounding boxes, and glyph evidence from the existing
siphonophore source cases, rather than PDFs or a built corpuscle. The JSON records
source paths/hashes, physical pages, and the Docling capture settings. `[…]`
marks deliberate shortening of prose used only for layout tests; notation and
diagnosis assertions use unshortened fragments.

`glyph_crops.json` records the source coordinates of four tiny raster crops.
The Kidwai date-range glyph is printed as a dash despite being encoded as ±.
That negative control prevents treating native PDF Unicode as ground truth.
Superscripts require raised source spans. Symbol replacements additionally
require matching local text and independently checked source glyph geometry.
The Haddock unit uses a source descender, a different font, and rendered OCR;
missing OCR leaves an explicit unresolved observation.

The Pakhomov current TableFormer cells correctly preserve ±. The external-source
replay applies the historical served `+` corruption as a controlled input, then
checks correction against the rotated source page; the already-correct cells
must remain unchanged.

Set `CORPUS_LIBRARY_DIR` to the read-only source library to replay glyph repairs
against the hashed PDFs. The local tests use only these fragments. They establish
the reported failure mechanisms and negative controls, not full-corpus recall or
the complete build/embedding/serve acceptance gate.

`test_source_layout_retrieval.py` extends the Church pp. 7–8 and Haddock p. 1
fragments through the production caption-role/order functions, Docling JSON
save/reload, real HybridChunker splitting/merging, and metadata-first bounded
`get_chunks_by_section` → `get_chunks` retrieval. A deterministic local word
counter avoids downloading a tokenizer; it does not replace the serializer or
chunking algorithm. Church is tested with separate, merged and split
continuation chunks: all keep *Physalia minuta*, precede the adjacent
*P. physalis* treatment, and exclude captions. Removing the preceding heading
leaves the continuation unknown despite its literal `P. minuta` mention.
Haddock's `inside two of / our specimens` and `yellow to red / (583, 620, …)`
boundaries remain contiguous, with its caption in a separate chunk. The served
projection explicitly reports a 6-of-9 source-record limit for its merged
prose; all nine records remain in the unchanged artifact.

The same save/chunk/serve path with reading-order repair disabled reproduces
both named failures. This establishes the captured-fragment acceptance path,
including provenance and adjacent-species controls. The deliberate `[…]`
shortening remains in these fixtures, so these tests do not claim a fresh
complete-document parse, default-model token boundaries, embedding retrieval
quality, corpus-wide reading-order recall or a deployed bundle repair.


`totton_treatment_boundaries.json` captures exact selected item text, original
layout labels, table cells and page provenance from the frozen `6fbf4e0` Totton
extraction. Minimal document wrappers omit intervening prose; original refs and
the complete extraction/PDF hashes remain recorded. The source-reviewed cases
are the Prayinae key (physical p120) and Prayoides genus/diagnosis transition
(pp129–130), retaining prior and next species headings as boundary controls.
The table's `document_index` label is intentionally preserved. Tests use actual
save/reload, production chunking and bounded serving; only the tokenizer is
instrumented. No PDF, gold artifact or frozen source label was modified.
