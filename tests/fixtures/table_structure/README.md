# Source table, identification-key and word-boundary regressions

These are small **source-derived item/cell fragments**, not whole corpuscles or
PDF fixtures. They were captured from local library PDFs with the installed
Docling native-text CPU route (`do_ocr=False`, TableFormer enabled) on
2026-09-18. Each JSON records the original PDF SHA-256, physical page numbers,
unchanged selected text/cell geometry/spans, and independent PyMuPDF source-line
observations. Reference indexes were renumbered after selecting the relevant
items; derived grid copies were omitted. Full unmodified page captures were
retained locally under `/tmp/table-integrity-source/` during investigation.

| Source / physical page | Hand-checked expectation and scope |
| --- | --- |
| Daniel 1974 / 20 | The `Genus ?Epibulia Eschscholtz, 1829` heading is one merged cell spanning four columns. Upstream Markdown repeats it four times. Logical serialization emits it once; distinct equal-valued cells are still separate observations. The surrounding OCR spelling and TableFormer row assignments are not claimed corrected. |
| Pugh & Haddock 2016 / 42 | All eight nectophore-key leads retain their destinations, including `Erenna cornuta Pugh, 2001` as one endpoint. |
| Pugh & Haddock 2016 / 44 | The bract key maps largest bracts >35 mm to 2 and <35 mm to 3; the complete transverse-ridge branch maps to `Erenna cornuta Pugh, 2001`. All eight leads retain endpoints within the tested 80-word chunk budget. Larger individual rows receive explicit partial-row metadata if they must split. |
| Daniel 1985 / 271–272 | Geometry binds the sharply conical hydroecium lead to `kochi`, shallower/nearly horizontal lead to `delsmani`, beyond-apex/deep-hydroecium lead to `atlantica`, and rounded-apex/sausage-shaped lead to the final destination. **Spelling remains unresolved:** the native PDF and Docling say `hargmannae`, but the included p272 crop visibly prints `bargmannae`. The association records `geometry_verified_spelling_unverified`; no global b/h substitution or bibliography-based renaming is made. |
| Hosia et al. 2024 / 5 | Positive control: the usable comparative Nanomia character table retains every cell value, column position, and source header. It is not classified as an identification key. |
| Mapstone 2009 / 68 | Source lines print `Anterior nectophore alone developed`; current Docling preserves these spaces. Replaying the reported old collapsed value against this same source evidence regenerates the same source text/chunks as a clean extraction. |
| Hissmann 2005 / 7 | Fresh Docling still emits `Theholotypepossessedninenectophoresandaboutnine`; exact native source letters plus spaces recover the sentence. An independent whole-word annotation check recovers the missing nectophore mention. |
| DuClos et al. 2022 / 6 | The PDF visibly separates the primary phrase's words, but native text omits spaces. Source boxes have uniform ~1.085 pt interword jumps and zero intraword gaps. Both recorded crop-OCR modes support these boundaries. OCR's extra `cou nt` split has **no** geometric support and is discarded. |
| Mapstone 2009 / 200 | **Adjudicated source typography, not an extraction error:** the included rendered crop actually prints `Anteriornectophorewithsevencompletelongitudinal`. Source character gaps are uniform. This remains unchanged and is a negative control against invented segmentation. |

The DuClos observation includes its crop hash, both OCR outputs, Tesseract /
Leptonica identity, English traineddata SHA-256, PyMuPDF version and fixed crop
policy. Tests replay these captured observations without OCR executables or
model downloads. Actual source-page replay on the same producer repaired the
primary phrase and two additional gap-supported runs; a fourth proposed run
was retained because both OCR modes omitted a geometric boundary. That is a
reviewable disagreement, not a successful repair or a corpus-wide error count.
The extra three runs' individual OCR receipts were not retained; their earlier
aggregate outcome is not a source-graded accepted-repair precision denominator.
The [frozen 13-decision spacing review](spacing_review_2026_09_18/README.md)
preserves the separate source grades, denominators, uncertainty and crop hashes,
with instructions to re-render every region from hash-verified originals.

The additional `DuClos_etal2022-6-mixed-font.json` retains an unchanged item
from that actual page capture and fresh regional OCR observations. Its native
PDF joins italic `N.bijuga` to upright
`werecollectedatFridayHarborLaboratoriesbetween`. The v2 spacing policy proposes
only internal gaps of existing uniform-font portions, with the same numeric
thresholds and both-mode exact-letter agreement. It never proposes a space at
the font transition. Multiple portions still share one crop and two OCR calls
within the original eight-line page budget; the actual crop bounds and hash
are recorded alongside full native-line offsets.

**This named phrase remains unresolved:** fresh full-line and narrowed-region
mode 7 OCR both omit the `at | Friday` boundary. The original Docling text and
`N.bijuga` spacing stay unchanged, and `ocr_disagreement` is retained. The
controlled agreement tests exercise safe application and preservation of all
letters/outer boundaries; they are explicitly not passing source OCR evidence.
These captures used the existing stored Docling item, not a fresh extraction.
Unresolved long-run receipts now also produce the existing review warning;
they are not automatically counted as spacing errors or forcibly segmented.

`pipeline.table_structure` runs before Markdown and chunk production. Logical
cells retain row/column indexes and spans; extraction Markdown uses actual HTML
spans, while chunk prose uses one logical row per line. Dotted key leaders are
layout-normalized so they do not consume the token budget. Headers, row spans,
couplet context and partial-row flags are chunk metadata, not repeated source
mentions. Custom `corpus__*` repair metadata is blocked from Docling's prose
serializer; raw observations remain in the structured build artifact.

Missing spaces use exact same-letter source lines first. The optional OCR route
requires a strongly separated native character-gap pattern plus confirmation
from both OCR segmentation modes. It accepts only existing geometric gaps,
never OCR-only word splits or new letters. At most eight small line crops per
page are considered, with a pixel ceiling and subprocess deadlines. Missing
OCR, timeouts and disagreement keep source text unchanged and retain quality
observations. Ordinary letter spacing, intact compounds, multilingual words
and morphological variants are explicit negative controls.

Deployment requires extraction/chunk regeneration, annotation rebuilding and
new embeddings/bundling. The producer policy and OCR/model identity belong in
the extraction input fingerprint, inherited by downstream stages; merely
rechunking old `text.json` cannot recover source boundaries or cell spans.
The tests include save/reload stability and clean-versus-old-artifact repair
agreement. These fixtures cover the named source cases, not every table, key,
OCR spelling error or collapsed run in the full library.

`tests/test_key_branch_retrieval.py` exercises the Daniel p271–272 key through
the actual geometry-association helper, JSON save/reload, production
`chunk_text`/HybridChunker, `CorpusIndex` and bounded `get_chunks` calls.
Local word-count budgets of 12, 25 and 2000 force split and merged variants
without downloading a tokenizer. Every captured branch retains its destination,
source geometry and unverified-spelling status; partial branches additionally
expose completeness, endpoint presence and adjacent fragment IDs. Removing
geometry association leaves names as unassociated text. A missing-scope legacy
artifact rechunks to the same result as fresh materialization, and policy
changes invalidate chunk consumers without rerunning extraction. The test
retains per-budget `acceptance.json` receipts in pytest's temporary directory.
`tests/test_key_branch_resume.py` additionally exercises both production Stage 1
resume gates after a policy change and compares incremental and clean chunks.
The key metadata leaves the exact embedding payload and its fingerprint
unchanged; a chunk text change still changes that fingerprint. Stage 2 does not
store this scope in vector rows, so metadata-only rechunking requires no model
rerun.

This path begins at the captured Docling items; it does not rerun full-PDF
extraction, OCR, the production embedding tokenizer or a live MCP transport.
The served `hargmannae` spelling remains explicitly unverified against the
rendered `bargmannae` crop; this acceptance proves the geometric association and
fragment route, not a spelling correction.


### Source-reviewed key classification (#320)

`key_classification_source.json` captures exact table cells and page provenance
from the completed `6fbf4e0` gold build in minimal Docling wrappers. Original
PDF identities, source table references and extracted-document hashes are
retained; the original observations and frozen retrieval selection are unchanged.
The independent source review is recorded in
`dev_docs/examples/siphonophore_retrieval_review_2026_09_18/`.

Ahuja physical p25 is numbered bibliography, and Totton p153 is a nematocyst
measurement table. Both must remain ordinary tables. Totton pp57/120 are
genuine keys and remain classified as keys, alongside the Pugh–Haddock
fixtures. Source tables pass actual save/reload, chunk serialization and bounded
row metadata checks. Synthetic controls cover quantities, years, page numbers,
missing/self links, a short two-couplet key and a separate destination column
without leaders. Neither key classification nor retained cell text establishes
that the extractor recovered every source relationship. This change does not
repair wrong treatment names or satisfy the independent retrieval sample gate.
