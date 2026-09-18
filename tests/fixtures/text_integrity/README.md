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
