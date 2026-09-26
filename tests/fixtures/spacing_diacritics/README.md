# Niño/Niña source acceptance (#316)

`nino_source.json` and its small rendered crop capture the first-page abstract
of Keister & Peterson (2003), *Zonal and seasonal variations in zooplankton
community structure off the central Oregon coast, 1998–2000*, physical page 1 /
printed page 341. The source PDF SHA-256, library revision, exact native glyph
boxes and crop settings are recorded. The rendered page visibly prints
**El Niño/La Niña**; its native layer contains separate overlapping tilde glyphs.

The exact damaged substring **El Nin ˜o/La Nin ˜a** occurs in `chunk_2` of the
retained **v1.2.1** bundle for `df1447c3d475`. Its byte hash and substring offsets
are retained. Issue #316 names this string family without assigning it a paper;
this is an actual matching source case, not a claim to have identified the
specific unnamed chunk in the later v1.4 audit. The regression checks current
geometry repair directly, and optionally against the original PDF with
`CORPUS_LIBRARY_DIR`. It does not rerun a full extraction model or assert a
corpus-wide error rate. Shifted nonoverlapping source accents are a negative
control: a familiar phrase alone must not authorize correction.

The adjoining authority regression uses the pinned library's real Kolliker1853
bibliographic fields, but its truncated and explicitly corrected parsed-author
states are **controlled input states**. It verifies that rerunning materialization
after an upstream adjudicated correction replaces the active unresolved edge,
retains both immutable observations and their identical raw citation, preserves
the complete curated author, and is a no-op on an unchanged refresh. It does not
guess that every `¨lliker` is Kölliker or claim to have regenerated any particular
production Grobid record. Actual corpus adjudications still require source review.
