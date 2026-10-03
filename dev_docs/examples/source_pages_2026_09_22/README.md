# Selected source-page acceptance, 2026-09-22

These receipts follow normal `corpus run` preparation, extraction, annotation,
embedding, post-processing and bundling. Only copied BibTeX `keeppages` values
and input symlinks select the pages; original library files remain unchanged.
The complete copied author catalog supports surname verification. Grobid is
disabled, so this build makes no fresh reference-edge or bibliography-completeness
claim. Page selection precedes scan classification; it does not test full-paper
classification. This is a bounded source pilot, not a second fixture corpus.

`baseline-serving.json` records eight documents/12 physical pages at `4d0ecfd`:
Tung 30/84, Hung 3, Yu 2, Liu 2, Boysen-Ennen 12/16, Mapstone 47/68/200,
Hissmann 7 and DuClos 6. It verifies 62 real BGE-M3 vector rows, Chinese/German
chunks and captions, surnames, two byte-identical MCP image deliveries, and an
unchanged bundle inventory. The Liu page is an already-legible Chinese control.
Separate actual-source tests passed 153 checks with no skips.

`integrated-serving.json` repeats the normal build with the spacing/exponent
producers at `b260406` and adds Stepanjants physical page 8. Its nine documents
produce 75 nonzero 1024-dimensional BGE-M3 rows; every stored row's text and
heading metadata agrees with its chunk. Six named documents pass actual stdio
MCP chunk checks; Chinese and German figure delivery preserves exact bundled
image bytes. The unit-fixture-only follow-up `bfbf90c` does not alter production
code. All 99 compared materialized files, including pipeline state, remain
identical on unchanged extraction; embedding refresh changes two documents,
skips seven, and then skips all nine on the unchanged follow-up.

For #306/#312, the Chinese heading, caption, taxa and signed quantity are
preserved, as are German `RMT-8-Fänge` and the page-16 figure caption. Recovery
and refusal metadata are retained. Existing literal-Latin-1 and literal-`Ã`
controls remain tested. This does not assert complete character fidelity for
every page of those papers.

For #334, the fresh force-OCR route now preserves `Anterior nectophore alone
developed`; its `nectophore` annotation has the exact whole-word span in chunk7.
The Hissmann holotype sentence and DuClos siphosome relationship also survive
chunking, real vectors and serving. Source-printed closed typography on Mapstone
page200 remains unchanged. The third source-confirmed proposed repair,
`nectophoreonly`, applies in the earlier saved-Docling replay but is refused in
the fresh build because its structured owner cannot be uniquely aligned. The
additional DuClos G1 run remains unresolved. Keep the separate frozen spacing
review's **2 correct, 0 incorrect, 1 unresolved** accepted decisions; neither
this targeted pilot nor its three selected proposals establishes corpus-wide
precision or recall.

For #303, original scan evidence survives fresh full-page OCR and changes only
the verified `2000 мм?` to `2000 мм³`. The second Russian exponent remains an
explicit OCR disagreement. Broader scientific-notation grading remains open.

`authority.json` separates #315's source text from citation identity. Fresh v3
Mapstone source evidence, applied to the real retained-v1.2.1 reference b166,
changes its actual authority edge to curated Alvariño1971a; raw observations
remain unchanged and refresh is idempotent. Three references from a separate
prior complete source TEI remain unresolved because their titles do not satisfy
the conservative identity gate. They are not claimed as repaired edges. The
selected surname sample and its limits remain in `dev_docs/SURNAME_RECOVERY.md`.

Full PDFs, model weights, vector values and duplicated extracted prose are not
vendored. Hashes identify the source/artifact producers and retained scratch
captures under `/tmp/corpus-original-text-acceptance-20260922` and
`/tmp/corpus-original-biblio-evidence/source-pilot`. Normal CLI logs there retain
launch failures as well as successful runs. Final gold/full-corpus, Qwen and
retrieval evaluation remain release gates.
