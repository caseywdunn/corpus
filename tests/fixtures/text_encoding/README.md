# Source encoding and accent evidence

`source_regions.json` contains nine small source regions captured on 2026-09-18
with Docling 2.94.0, CPU, OCR disabled. Source SHA-256, library revision, physical
page, item bounds and original extracted strings are included. Full PDFs remain
in the external library. This adds no fixture corpus.

The cases cover Tung2003 physical page 30 (heading and scientific paragraph),
page 84 (figure 1 caption), Hung2002 page 3 and Yu2006 page 2 (short headings and
prose), and Fernandezetal2007 page 1 (author and affiliation lines). The source
Chinese heading and Fernández author line were also inspected as rendered
images. Accent cases retain individual native character boxes; Chinese cases
retain the independent Unicode text from the same region. Tests can replay the
original geometry with `CORPUS_LIBRARY_DIR=/path/to/siphonophores/library`.

Big5-compatible corruption is freshly reproduced in all three reported Chinese
sources. Automatic recovery requires every observed non-whitespace byte to
agree, in order, with the same region's native Unicode text re-encoded as Big5
or CP950. The native region must have Chinese characters and at least 65% of its
encoded bytes must be observed. Source-only insertions account for independently
verified line tails lost during extraction; the replacement and coverage are
recorded. Any conflicting byte—including a year or quantity—prevents automatic
recovery. Correct Unicode/Latin-1, literal spacing symbols and legitimate `fi`
sequences are negative controls. A successful second pass does not repeat edits.

Accent recovery uses source overlap and alignment, not a name dictionary:
`Ferna ´ndez`, `Espan ˜ol` and dotless-i plus acute become the printed spellings.
The exact pre-repair substrings/offsets, base/accent boxes and method survive in
metadata; Docling `orig` remains unchanged. This does not repair independent
Grobid author parsing. Orphan surname accents instead produce review warnings
and observation-specific unresolved identities unless an independent DOI is
available. They cannot acquire guessed surname aliases or attract clean
references through reverse-order fuzzy matching/reconciliation.

`already_legible_chinese.json` adds a separate genuine negative control from
Liu et al. 2012, physical page 2: the Chinese title
`哥斯达黎加外海夏季表层浮游动物种类组成及分布`. The retained v1.2.1 served title matches
the native glyphs and the separately rendered, visually checked source crop
`liu2012_legible_title.png`. Its JSON records the source PDF hash, exact title
bounds, native character boxes, retained artifact hash and heading offsets.
The one-item Docling wrapper used in the test is an explicit reconstruction;
its `text` label is a test choice, not a claimed new layout-model prediction.
No fresh Docling conversion was run for this control. With
`CORPUS_LIBRARY_DIR` set, the test replays the actual source region through
encoding recovery, serialization, real HybridChunker with a local word counter,
and bounded chunk retrieval. It requires no repair, no added repair receipt and
the same readable title throughout. This establishes only this title region;
it does not establish Liu's entire paper or the wider Chinese cluster.

These tests do not establish complete-document or full-corpus fidelity. Fresh
current-policy OCR pilots of Boysen-Ennen1987 page 12 and Mapstone2009 page 47
still produce `Finge` and `Alvarifio`; #312/#315 remain under investigation.
Regional OCR with appropriate packs can read the printed forms, but that pilot
alone does not establish a safe automatic recovery policy. Source uncertainty
is surfaced through `text.json` receipts and the `source_text_integrity` quality
gate. The independent gold rebuild, scientific notation checks, authority-edge
acceptance and full-corpus audit replay remain release gates.
