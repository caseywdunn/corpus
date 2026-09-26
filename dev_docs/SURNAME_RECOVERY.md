# Source-supported citation surname recovery

Citation surname recovery (#315) runs during extraction. It uses the reviewed
BibTeX author/year catalog to propose near-name readings, then checks a crop of
the same source citation before changing text. It does not replace character
sequences globally. A surname containing `fi` is not itself evidence of damage.

The current bounded policy covers single-token surnames of at least five
letters, a nearby publication year, an exact first-three-letter match, and one
or two Unicode edit operations against a catalog surname containing a non-ASCII
character. Already known surname/year pairs are left alone. Multiple possible
names, quoted or `[sic]` text, missing or ambiguous native word/year geometry,
unavailable OCR models, and OCR disagreement remain recorded review leads.
This deliberately limited detector does not establish the absence of other OCR
errors, identify every compound surname, or reinterpret source-printed names.

Each proposal needs a unique source word/year anchor inside the corresponding
Docling item's page and bounding box. Native PDF word coordinates are transformed
to the displayed page; the crop is rendered upright even for rotated pages.
A small 600 DPI crop is read in two
Tesseract segmentation modes. Both readings must independently contain exactly
the complete proposed surname next to a number token. At least one reading must
contain the complete expected four-digit year, optionally followed by a single
lowercase citation suffix. Any other valid year in either reading blocks the
repair. A malformed token such as `19814` is retained as invalid evidence; it
never supplies `1981` by taking a prefix. This check changes only the surname,
not the year or its suffix. No candidate is supplied as a user
word or character whitelist, and system/frequency dictionary hints are disabled.
The two modes share a model; their agreement is corroboration, not an independent
human judgment. A printed alternate spelling that OCR reproduces stays unchanged.

The language pack is chosen before OCR, from declared non-English OCR languages
of the curated author's same-year publications, then the author's other curated
publications, then English. Frequencies in that reviewed catalog break choices,
with a deterministic lexical tie break. The policy never tries successive packs
until one produces a desired name. This can leave a correct proposed diacritic
unresolved when the chosen pack lacks suitable character coverage.

The producer receipt includes the policy, exact catalog digest, PyMuPDF version,
Tesseract version and selected traineddata hashes, crop settings, OCR modes,
timeouts and candidate budgets. Processing is limited to 12 candidates per page,
64 per document, two OCR calls per unique candidate, 15 seconds per call and two
million pixels per crop. Original Docling text remains in `orig`; corrected prose
and structured `corpus__surname_recovery` decisions are separate. The extraction
report belongs in `text.json` under `source_text_integrity.surnames`, and structured
chunk integrity carries the source decisions without placing them in chunk prose.
Changing the source/catalog/model inputs must invalidate extraction and its
consumers; the server only reads these materialized artifacts. The current
receipt retains only surname/year exclusions for ASCII names and names too short
to be candidates. Title or BibTeX-key changes to those records do not invalidate
extraction. Eligible non-ASCII candidate records still include all title and
language evidence consumed by recovery and reference corroboration, so edits to
those fields invalidate extraction across the corpus. Narrowing that remaining
cost safely requires an independent current curated-title lookup in reference
materialization, rather than ignoring a field that the stored source report
still consumes.

The report retains every decision and also materializes an `unresolved` subset
for the extraction quality gate. Confirmed repairs, source-confirmed alternate
spellings, and quoted/`[sic]` contexts do not enter that subset. Ambiguous anchors,
unavailable OCR, disagreements and other unverified candidates produce a source
integrity warning with the review count.

## Bibliographic identity

A corrected citation in prose does not by itself identify a reference. The
bibliography materializer also requires an exact normalized full curated title
(at least 25 normalized characters), the same year, and one unambiguous verified
surname replacement. Missing or damaged titles remain unresolved with a review
reason. Raw `reference_observations` never change. The derived author and edge
record the observed author, curated keys, source PDF/crop hashes, page and producer
policy. The exact consumed extraction reports participate in materialization
identity, so changed or removed source evidence rederives unchanged observations.

## Small source replay sample

`tests/fixtures/surname_recovery/sample.json` records source filenames, SHA-256
identities, physical pages, crop boxes, rendered image hashes, visual labels,
recorded OCR outputs and producer identity. The compact crops come from a fixed
library revision. Fourteen candidate occurrences were selected from six papers:
the two reported Mapstone page 47 cases and the first up to four uniquely anchored
near-name citations per paper in a fixed six-file search list. A paper in that
search list supplied no eligible candidate. Two correctly accented source
citations from two additional papers are negative controls.

On this selected source sample, all 14 proposed names are visibly damaged in the
extracted reading; 10 are repaired correctly and four Forskål readings remain
unresolved because the selected Latin OCR output disagrees with the source's
printed diacritic. Both correctly accented source controls remain untouched.
These are 14/14 correct proposals and 10/10 correct admitted repairs **within the
selected cases**, with four missed repairs. Selection was not random and these
counts cannot estimate corpus-wide precision, recall or the number of damaged
citation edges.

Four separately labelled synthetic crops exercise a source-printed misspelling,
a legitimate `fi` surname, an unrelated name and an already correct accented
name. The printed misspelling produces a review proposal but no repair; the other
three produce no proposal. Synthetic controls are not counted in source rates.

The reference-edge replay uses the actual retained Mapstone `b166` parsed record
from a v1.2.1 bundle and a freshly verified source crop. It demonstrates that the
wrong derived ghost edge changes to the curated Alvariño 1971 work, that raw
observations remain byte-for-byte unchanged, and that removing verification
rederives the unresolved edge. It is not a fresh Grobid result, a reconstruction
of the later audited v1.4 bundle, or a corpus-wide acceptance run. An optional
original-PDF replay (`CORPUS_LIBRARY_DIR`) exercises source recovery through the
real HybridChunker and checks that provenance stays out of the prose. Extraction,
atomic figure resets and the separately scheduled vision phase share the same
catalog and producer receipt. Both resume gates and configuration/source audits
track these inputs; tests cover catalog changes/removal, unchanged reruns and
clean-build equivalence. A direct module invocation also resolves the configured
BibTeX relative to its config file, with an explicit CLI argument taking
precedence. Full rebuilt-corpus acceptance remains required.

## Hosted OCR year variation

The v3 policy was prompted by an actual hosted Tesseract 5.5.3 result on the
Alvarinoetal1990 page 14 source crop: mode 6 read `[Alvariño, 19814].` and mode 13
read `[Alvariño, 1981a].`. Both modes recovered the complete surname; only the
second supplied a valid matching year. The previous policy discarded the first
name along with its damaged date and left the surname unresolved. The hosted
outputs and selected model identity are preserved in
`tests/fixtures/surname_recovery/hosted_ocr_553.json`; the regression replays those
outputs, rather than claiming to run that hosted binary locally.

The source anchor still supplies the expected year, and OCR must independently
corroborate it. Swapping the modes does not change the decision. Conflicting valid
years, two malformed dates, different or incomplete names, multiple name/date
readings, and a single available segmentation mode cannot authorize a repair.
Raw OCR outputs and the matching, conflicting and invalid year evidence remain
in the receipt. The producer policy changes from v2 to v3, invalidating extraction
and its existing chunk, annotation and figure consumers; builds made with v2 do
not acquire v3 acceptance evidence without re-extraction.
