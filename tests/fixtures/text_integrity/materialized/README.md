# Materialized source text through embedding inputs (#303, #312)

These seven compact saved Docling artifacts add the missing downstream check to
existing source/glyph evidence. Their combined size is about 136 KiB, plus the
manifest; they contain no PDF or vector database. Redundant computed table grids
are omitted where present, while the actual cells and document structure remain.

Five scientific-text artifacts start from the actual saved Docling source-page
extractions in the original notation investigation. `prepare_scientific_text`
was run against the original PDFs again at `2a4abc5`, with each full PDF SHA-256
checked against the existing source manifest. The operation used source glyph
geometry and rendered glyph checks, including the bounded English OCR check for
the Haddock micro-unit. It did not rerun Docling or load an extraction model.
`manifest.json` records source and artifact hashes, source pages, the original
Docling artifact hashes, policy, OCR runtime/model identity and actual repair or
unresolved receipts. Its repair counts include non-scientific superscripts, so
those counts are not scientific-symbol accuracy measures.

The Pakhomov artifact contains only the two previously captured TableFormer
cells. Both already contain the correct ± signs. This is the unchanged real
source control, rather than the separate controlled historical `+` corruption.
The Kidwai source artifact retains the printed `1857-1859` date range and its
unresolved misleading-native-glyph observations.

The German artifact is the byte-exact `docling_doc.json` from the earlier actual
Boysen-Ennen physical-page-12 pilot: `rebuild_figure_base` →
`extract_docling_content` after full-page German/English OCR, with original-layer
regional evidence passed through preparation. Physical page 12 is prepared page
1. Its provenance retains the original `RMT-8-FÃ¤ng`, the prepared
`RMT-8-Finge`, the confirmed `RMT-8-Fänge`, and other unresolved observations.
The manifest distinguishes the full source PDF, original page slice and prepared
PDF hashes. This saved pilot contains no figure/caption objects; it supplies no
German caption-propagation evidence.

`tests/test_source_text_embedding.py` loads these saved artifacts through the
production `DoclingDocument` reader, serializes them, and runs the real
`HybridChunker` through `chunk_text`. Only its token counter is replaced with a
small local word counter. The test then calls production `embed_document` with a
recording two-dimensional backend and real temporary LanceDB. The recorder
captures the exact strings received by the embedding backend; it does not
reconstruct an alternative embedding input. The committed vector-row text and
one-chunk-at-a-time `get_chunks` responses must be byte-for-byte equal to those
strings. Source integrity metadata remains separate from embedded prose.

The assertions cover uncertainty signs, superscript units/exponents, both copies
of the corrected Sutherland caption equation, the Haddock `200 µm` caption,
Pakhomov's existing correct values, Kidwai's date range, and the German corrected
word plus unchanged unresolved prose. Earlier-generation controls reconstruct
the exact persisted `orig` text for Sutherland and the German page. Repaired
chunk text must invalidate the old embedding marker, replace its rows, and then
pass unchanged `--resume` without another backend call. Policy fingerprint
coverage separately checks scientific extraction consumers; the existing native
OCR producer test covers preparation and its descendants.

This is a text-preservation test. The numerical vectors have no semantic meaning,
and the local token counter does not establish production model token boundaries
or retrieval quality. It is not a fresh full-document build, bundle/SSE replay,
German figure-route check, or the independently graded broader precision/recall
sample required by #303. The existing source-reviewed sample limitations remain.
